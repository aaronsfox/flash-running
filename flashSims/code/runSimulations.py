# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This script loops through the various sprint speeds using the Flash simulated model.
    
    TODO:
        > Add detailed notes...

"""

# %% Import packages

import opensim as osim
import osimHelper as helper
import numpy as np
import os
import shutil
from scipy.signal import butter, filtfilt
# import matplotlib.pyplot as plt

# %% Set-up

#Add OpenSim geometry path (weird issues with this on new laptop)
osim.ModelVisualizer.addDirToGeometrySearchPaths('C:\\OpenSim 4.3\\Geometry')

#Set the desired speed list
speedList = [10,20,40]

#Set a list of joint coordinates we want to track in later simulations
trackingCoords = ['pelvis_tx', 'pelvis_ty', 'pelvis_tz',
                  'pelvis_tilt', 'pelvis_list', 'pelvis_rotation', 
                  'hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r',
                  'knee_angle_r', 'ankle_angle_r', 'mtp_angle_r',
                  'hip_flexion_l', 'hip_adduction_l', 'hip_rotation_l', 
                  'knee_angle_l', 'ankle_angle_l', 'mtp_angle_l', 
                  'lumbar_extension', 'lumbar_bending', 'lumbar_rotation',
                  'arm_flex_r', 'arm_add_r', 'arm_rot_r',
                  'elbow_flex_r', 'pro_sup_r',
                  'arm_flex_l', 'arm_add_l', 'arm_rot_l',
                  'elbow_flex_l', 'pro_sup_l']

# %% Create updated coordinates and GRF profiles from experimental data

"""

This section extracts the GRF data from the experimental section of interest and
scaled up the magnitude of forces based on the increase in speed. We fit a logarithmic
function to the increase in vGRF with speed based on data presented in a PhD thesis
by Hierholzer (https://dc.etsu.edu/cgi/viewcontent.cgi?article=5114&context=etd) and
apply this scaled change to the experimental data.

"""

# #Set the experimental speeds from Hierholzer (between 10.1 km/h to 21.7)
# #The speeds where a larger sample size was tested have been used
# testSpeedsKMH = np.array((10.1,11.4,12.7,14.0,15.2,16.5,17.8,19.1, 20.4, 21.7))

# #Convert speeds to metres per second for consistency
# testSpeedsMS = testSpeedsKMH / 3.6

# #Set vGRF in BW
# vGRF = np.array((2.38, 2.49, 2.56, 2.62, 2.66, 2.69, 2.76, 2.83, 2.96, 2.92))

# #Take log values of the speed data
# testSpeedsMS_log = np.log(testSpeedsMS)

# #Get the logarithmic curve fit to the log_x data and vGRF
# #This gives the a & b coefficients to create the linear equation
# logCurve = np.polyfit(testSpeedsMS_log, vGRF, 1)
# a = logCurve[0]
# b = logCurve[1]

# #Create an x and y fit for the logarithmic curve
# #This encompasses speed values up to the max that we'll be testing
# x_fitted = np.linspace(np.min(testSpeedsMS), np.max(speedList), 100)
# y_fitted = a * np.log(x_fitted) + b

# #Visualise fit
# fig, ax = plt.subplots(nrows = 1, ncols = 1)
# ax.scatter(testSpeedsMS, vGRF, c = 'k')
# plt.plot(x_fitted, y_fitted, c = 'r')

#Get the gait timings from the experimental data to extract from the GRF file
startTime, endTime = helper.getGaitTimings(grfFile = '..\\..\\expSims\\data\\sprint_grf.mot',
                                           extLoads = '..\\..\\expSims\\data\\sprint_grf.xml',
                                           startForceName = 'RightGRF1',
                                           stopForceName = 'LeftGRF1',
                                           forceThreshold = 20)

#Get the experimental kinematics
expTrackingTrajectory = osim.MocoTrajectory('..\\..\\expSims\\fullTrackingSim\\fullTrackingSimSolution.sto')

#Determine the original sprinting speed from this time frame of the experimental data
#This is the entire period examined in the experimental marker tracking solution
expSpeedMS = expTrackingTrajectory.getState('/jointset/ground_pelvis/pelvis_tx/speed').to_numpy().mean()

# #Extract the vGRF in BW from the fitted logarithmic curve
# expGrfBW = a * np.log(expSpeedMS) + b

#Convert tracking trajectory to states table
expTrackingStates = expTrackingTrajectory.exportToStatesTable()

#Read in the experimental GRF data to scale up
grfTable = osim.TimeSeriesTable('..\\..\\expSims\\data\\sprint_grf.mot')

#Get the time data from the grf table
origForceTime = np.array(grfTable.getIndependentColumn())

#Find the indices closest to the start and end times
#Extract the experimental GRF time from here
grfStartInd = np.absolute(origForceTime - startTime).argmin()
grfEndInd = np.absolute(origForceTime - endTime).argmin()
expGrfTime = origForceTime[grfStartInd:grfEndInd+1]

#Extract the right lim ground reaction force over this time
#These data were on force plate 2 and related to a right foot strike
grfX = grfTable.getDependentColumn('ground_force_2_vx').to_numpy()[grfStartInd:grfEndInd+1]
grfY = grfTable.getDependentColumn('ground_force_2_vy').to_numpy()[grfStartInd:grfEndInd+1]
grfZ = grfTable.getDependentColumn('ground_force_2_vz').to_numpy()[grfStartInd:grfEndInd+1]

#Create a simplified set of table labels for right and left force values
tableLabels = osim.StdVectorString()
for currLabel in ['ground_force_r_vx', 'ground_force_r_vy', 'ground_force_r_vz',
                  'ground_force_l_vx', 'ground_force_l_vy', 'ground_force_l_vz']:
    tableLabels.append(currLabel)

#Create upscaled ground reaction force data for the desired sim speeds

#Loop through speed list
for sprintSpeed in speedList:
    
    #Create folder to store data
    try:
        os.mkdir(f'..\\flashSimSpeed{str(sprintSpeed)}')
    except:
        print('Flash sim folder for current speed already detected...')
        
    ########## CREATE MODELS FOR EACH SPEED ##########
    
    #Create the model to use in this simulation
    osimModel = helper.createFlashModel(inputModelFile = '..\\..\\expSims\\data\\JA1_SCALED_Osim40_Muscles.osim',
                                        outputModelFile = f'..\\flashSimSpeed{str(sprintSpeed)}\\osimModel.osim',
                                        unilateralMuscles = False,
                                        jointsToWeld = ['subtalar_r', 'subtalar_l', 'radius_hand_r', 'radius_hand_l'],
                                        addMetabolicsModel = True)
        
    ########## CREATE UPDATED COORDINATES FOR SPEED ##########
    
    #Create a copy of the tracking states for export
    if speedList.index(sprintSpeed) == 0:
        
        #Grab the experimental tracking data
        speedTrackingStates = expTrackingTrajectory.exportToStatesTable()
        
    else:
        
        #Get the previous solution
        previousSolution = osim.MocoTrajectory(f'..\\flashSimSpeed{str(speedList[speedList.index(sprintSpeed)-1])}\\flashSimSpeed{str(speedList[speedList.index(sprintSpeed)-1])}Solution.sto')
        speedTrackingStates = previousSolution.exportToStatesTable()
    
    #Determine the new time duration based on change in speed
    #Time
    initialTime = speedTrackingStates.getIndependentColumn()[0]
    finalTime = speedTrackingStates.getIndependentColumn()[-1]
    timeDur = finalTime - initialTime
    #Speed
    speedMS = speedTrackingStates.getDependentColumn('/jointset/ground_pelvis/pelvis_tx/speed').to_numpy().mean()
    #New time duration
    newTimeDur = timeDur * (speedMS / sprintSpeed)
    
    #Create a new time array to apply to coordinates
    newCoordinatesTime = np.linspace(initialTime, initialTime + newTimeDur, len(speedTrackingStates.getIndependentColumn()))
    
    #Set in the new states table
    for ii in range(len(newCoordinatesTime)):
        speedTrackingStates.setIndependentValueAtIndex(ii, newCoordinatesTime[ii])
        
    #Write to file
    osim.STOFileAdapter().write(speedTrackingStates, f'..\\flashSimSpeed{str(sprintSpeed)}\\coordinatesToTrack.sto')
    
    #Convert to a kinematics file for inverse dynamics
    for colLabel in speedTrackingStates.getColumnLabels():
        if colLabel.endswith('/value'):
            #Rename to just joint coordinate
            speedTrackingStates.setColumnLabel(list(speedTrackingStates.getColumnLabels()).index(colLabel),
                                               colLabel.split('/')[3])
        else:
            #Remove the column
            speedTrackingStates.removeColumn(colLabel)
            
    #Export pure kinematics file
    osim.STOFileAdapter().write(speedTrackingStates, f'..\\flashSimSpeed{str(sprintSpeed)}\\kinematicsToTrack.sto')
    
    ########## CREATE UPDATED GRF FOR SPEED ##########
    
    #Read in the inverse dynamics template
    idTool = osim.InverseDynamicsTool('..\\resources\\inverseDynamicsTemplate.xml')
    
    #Set options in ID tool
    idModel = osim.Model('..\\..\\expSims\\data\\JA1_SCALED_Osim40_Muscles.osim')
    idTool.setModel(idModel)
    idTool.setCoordinatesFileName(f'..\\flashSimSpeed{str(sprintSpeed)}\\kinematicsToTrack.sto')
    idTool.setStartTime(newCoordinatesTime[0])
    idTool.setEndTime(newCoordinatesTime[-1])
    idTool.setResultsDir(f'..\\flashSimSpeed{str(sprintSpeed)}\\')
    
    #Run ID tool
    idTool.run()
    
    #Read in body forces from inverse dynamics
    bodyForces = osim.TimeSeriesTable(f'..\\flashSimSpeed{str(sprintSpeed)}\\body_forces_at_joints.sto')
    
    #Identify force offset by inspecting average data 60-90% through phase, keeping
    #in mind we had a 50 node solution. This period equates to the flight phase
    forceOffsetX = bodyForces.getDependentColumn('ground_pelvis_pelvis_offset_FX').to_numpy()[30:45].mean()
    forceOffsetY = bodyForces.getDependentColumn('ground_pelvis_pelvis_offset_FY').to_numpy()[30:45].mean()
    forceOffsetZ = bodyForces.getDependentColumn('ground_pelvis_pelvis_offset_FZ').to_numpy()[30:45].mean()
    
    #Get time array from body forces file
    bodyForcesTime = np.array(bodyForces.getIndependentColumn())
    
    #Create new table to fill
    newGrfTable = osim.TimeSeriesTable()
    
    #Set the column labels
    newGrfTable.setColumnLabels(tableLabels)
    
    #Create an array that stores the column data
    newForceVals = np.zeros((len(bodyForcesTime),len(newGrfTable.getColumnLabels())))
    
    #Set the X, Y and Z force components for right limb in the new GRF data
    #We leave the left limb force values at zero
    newForceVals[:,0] = bodyForces.getDependentColumn('ground_pelvis_pelvis_offset_FX').to_numpy() - forceOffsetX
    newForceVals[:,1] = bodyForces.getDependentColumn('ground_pelvis_pelvis_offset_FY').to_numpy() - forceOffsetY
    newForceVals[:,2] = bodyForces.getDependentColumn('ground_pelvis_pelvis_offset_FZ').to_numpy() - forceOffsetZ
    
    #Define low-pass Butterworth filter
    fs = 1 / np.diff(bodyForcesTime).mean()
    filtFreq = 100
    nyq = 0.5 * fs
    normCutoff = filtFreq / nyq
    b, a = butter(4, normCutoff, btype = 'low', analog = False)
    
    #Filter data
    newForceVals[:,0] = filtfilt(b, a, newForceVals[:,0])
    newForceVals[:,1] = filtfilt(b, a, newForceVals[:,1])
    newForceVals[:,2] = filtfilt(b, a, newForceVals[:,2])
    
    # #Set any force values past 60% back to zero
    # newForceVals[30::] = 0
    
    #Find first point where vertical force becomes negative > 30% in (i.e. 15 samples)
    zeroInd = [ii for ii, xx in enumerate(newForceVals[15:,1] < 0) if xx][0] + 15
    
    #Set any force values past this point to zero
    newForceVals[zeroInd::] = 0
    
    #Fill the table with rows by looping through the rows and columns
    for iRow in range(len(bodyForcesTime)):
        #Create the current row
        osimRow = osim.RowVector.createFromMat(newForceVals[iRow,:])
        #Append to table
        newGrfTable.appendRow(iRow, osimRow)
        #Set time for current row
        newGrfTable.setIndependentValueAtIndex(iRow, bodyForcesTime[iRow])
    
    # #Extract the predicted vGRF in BW for the current sprint speed
    # simGrfBW = a * np.log(sprintSpeed) + b

    # #Determine force scaling factor by increase in speed
    # upscaleFac = simGrfBW / expGrfBW        
    
    # #Create new table to fill
    # newGrfTable = osim.TimeSeriesTable()
    
    # #Set the column labels
    # newGrfTable.setColumnLabels(tableLabels)
    
    # #Create an array that stores the column data
    # newForceVals = np.zeros((len(expGrfTime),len(newGrfTable.getColumnLabels())))
    
    # #Set the X, Y and Z force components for right limb in the new GRF data
    # #This corresponds to the first, second and third columns
    # #We leave the left limb force values at zero
    # newForceVals[:,0] = grfX * upscaleFac
    # newForceVals[:,1] = grfY * upscaleFac
    # newForceVals[:,2] = grfZ * upscaleFac

    # #Fill the table with rows by looping through the rows and columns
    # for iRow in range(len(expGrfTime)):
    #     #Create the current row
    #     osimRow = osim.RowVector.createFromMat(newForceVals[iRow,:])
    #     #Append to table
    #     newGrfTable.appendRow(iRow, osimRow)
    #     #Set time for current row
    #     newGrfTable.setIndependentValueAtIndex(iRow, expGrfTime[iRow])
        
    #Write to mot file format
    osim.STOFileAdapter().write(newGrfTable, f'..\\flashSimSpeed{str(sprintSpeed)}\\grfToTrack.mot')
    
    #Create the external loads .xml file
    forceXML = osim.ExternalLoads()
    
    #Create and append the right GRF external forces
    rightGRF = osim.ExternalForce()
    rightGRF.setName('RightGRF')
    rightGRF.setAppliedToBodyName('calcn_r')
    rightGRF.setForceExpressedInBodyName('ground')
    rightGRF.setPointExpressedInBodyName('ground')
    rightGRF.setForceIdentifier('ground_force_r_v')
    forceXML.cloneAndAppend(rightGRF)
    
    #Create and append the left GRF external forces
    leftGRF = osim.ExternalForce()
    leftGRF.setName('LeftGRF')
    leftGRF.setAppliedToBodyName('calcn_l')
    leftGRF.setForceExpressedInBodyName('ground')
    leftGRF.setPointExpressedInBodyName('ground')
    leftGRF.setForceIdentifier('ground_force_l_v')
    forceXML.cloneAndAppend(leftGRF)
    
    #Set GRF datafile
    forceXML.setDataFileName('grfToTrack.mot')
    
    #Write to file
    forceXML.printToXML(f'..\\flashSimSpeed{str(sprintSpeed)}\\grfToTrack.xml')

# %% Predictive Flash simulations

"""

TODO:
    > Add notes...
    > Does logarithmic GRF scaling work well? - no...
    > Should we always go back and use kinematics from original tracking sim?
        >> The initial guess shouldn't always be this though...
        >> If this is done then the tracking states don't need to be called every time
    > Muscle force (and any other analyses) calculation after solution
    > Possibly need higher weights on certain coordinates (e.g. arms)?
    > Can muscles produce speedy movements with no changes to activation dynamics?
    > Going from coarse to fine problem works very differently than just starting with fine mesh...
        >> Seems like the coarse solution isn't actually helping much?
        
    > If tracking kinematics can't be achieved as speed increases - is targeting these feasible?
        >> Whatever the bounds, the tracking of prior solutions seems more feasible...
        >> This goes OK with the standard states tracking and 1e-3 GRF weights (up to 20 m/s)
            >> Maybe need to increase GRF tracking weights though?

"""



#Loop through speeds
#### TODO: add loop
sprintSpeed = speedList[0]
    
# #Create the model to use in this simulation
# osimModel = helper.createFlashModel(inputModelFile = '..\\..\\expSims\\data\\JA1_SCALED_Osim40_Muscles.osim',
#                                     outputModelFile = f'..\\flashSimSpeed{str(sprintSpeed)}\\osimModel.osim',
#                                     unilateralMuscles = False,
#                                     jointsToWeld = ['subtalar_r', 'subtalar_l', 'radius_hand_r', 'radius_hand_l'],
#                                     addMetabolicsModel = True)

#Read in model for this simulation
osimModel = osim.Model(f'..\\flashSimSpeed{str(sprintSpeed)}\\osimModel.osim')

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName(f'flashSimSpeed{str(sprintSpeed)}')

#Set model in tracking tool
trackModelProcessor = osim.ModelProcessor(osimModel)
track.setModel(trackModelProcessor)

#Load tracking states for later reference
trackingStates = osim.TimeSeriesTable(f'..\\flashSimSpeed{str(sprintSpeed)}\\coordinatesToTrack.sto')

# #Set the coordinates reference in tracking tool from previous simulations
# #Check for initial sprint speed or subsequent, and choose the appropriate tracking states
# if speedList.index(sprintSpeed) == 0:
#     trackingStates = osim.MocoTrajectory('..\\..\\expSims\\fullTrackingSim\\fullTrackingSimSolution.sto').exportToStatesTable()
# else:
#     trackingStates = osim.MocoTrajectory(f'..\\flashSimSpeed{str(speedList[speedList.index(sprintSpeed)-1])}\\flashSimSpeed{str(speedList[speedList.index(sprintSpeed)-1])}Solution.sto').exportToStatesTable()
    
# #Export to a coordinates file for tracking
# osim.STOFileAdapter().write(trackingStates, f'..\\flashSimSpeed{str(sprintSpeed)}\\coordinatesToTrack.sto')

#Set coordinates file as reference in tool
track.setStatesReference(osim.TableProcessor(f'..\\flashSimSpeed{str(sprintSpeed)}\\coordinatesToTrack.sto'))
track.set_states_global_tracking_weight(1)
# track.set_states_global_tracking_weight(0.1)
 
#Set allow unused references in case there are extra markers
track.set_allow_unused_references(True)

#Set speeds a derivatives of coordinate references
track.set_track_reference_position_derivatives(True)

#Set tracked states to guess
track.set_apply_tracked_states_to_guess(True)

#Create weight set for state tracking
stateWeights = osim.MocoWeightSet()

#Create a dictionary that provides the kinematic task weights for function
# taskWeights = {'pelvis_tx': 1, 'pelvis_ty': 1, 'pelvis_tz': 1,
#                 'pelvis_tilt': 1, 'pelvis_list': 1, 'pelvis_rotation': 1, 
#                 'hip_flexion_r': 1, 'hip_adduction_r': 1, 'hip_rotation_r': 1, 
#                 'knee_angle_r': 1, 'ankle_angle_r': 1, 'mtp_angle_r': 1,
#                 'hip_flexion_l': 1, 'hip_adduction_l': 1, 'hip_rotation_l': 1, 
#                 'knee_angle_l': 1, 'ankle_angle_l': 1, 'mtp_angle_l': 1,
#                 'lumbar_extension': 1, 'lumbar_bending': 1, 'lumbar_rotation': 1,
#                 'arm_flex_r': 1, 'arm_add_r': 1, 'arm_rot_r': 1,
#                 'elbow_flex_r': 1, 'pro_sup_r': 1,
#                 'arm_flex_l': 1, 'arm_add_l': 1, 'arm_rot_l': 1,
#                 'elbow_flex_l': 1, 'pro_sup_l': 1}
taskWeights = {'pelvis_tx': 1, 'pelvis_ty': 0.5, 'pelvis_tz': 1,
                'pelvis_tilt': 10, 'pelvis_list': 10, 'pelvis_rotation': 10, 
                'hip_flexion_r': 5, 'hip_adduction_r': 5, 'hip_rotation_r': 2.5, 
                'knee_angle_r': 5, 'ankle_angle_r': 2.5,
                # 'subtalar_angle_r': 0,
                'mtp_angle_r': 2.5,
                'hip_flexion_l': 5, 'hip_adduction_l': 5, 'hip_rotation_l': 2.5, 
                'knee_angle_l': 5, 'ankle_angle_l': 2.5,
                # 'subtalar_angle_l': 0,
                'mtp_angle_l': 2.5,
                'lumbar_extension': 10, 'lumbar_bending': 5, 'lumbar_rotation': 5,
                'arm_flex_r': 5, 'arm_add_r': 1, 'arm_rot_r': 1,
                'elbow_flex_r': 5, 'pro_sup_r': 1,
                'arm_flex_l': 5, 'arm_add_l': 1, 'arm_rot_l': 1,
                'elbow_flex_l': 5, 'pro_sup_l': 1}

#Loop through coordinates to apply weights
#Only apply to joint coordinate values, as the speeds will change
for coordInd in range(osimModel.updCoordinateSet().getSize()):
    
    #Get name and absolute path to coordinate
    coordName = osimModel.updCoordinateSet().get(coordInd).getName()
    coordPath = osimModel.updCoordinateSet().get(coordInd).getAbsolutePathString()

    #If a task weight is provided, add it in
    if coordName in list(taskWeights.keys()):
        #Append state into weight set
        #Coordinate value
        stateWeights.cloneAndAppend(osim.MocoWeight(f'{coordPath}/value',
                                                    taskWeights[coordName]))
        #Coordinate speed
        stateWeights.cloneAndAppend(osim.MocoWeight(f'{coordPath}/speed', 0))
        
#Add to tracking problem
track.set_states_weight_set(stateWeights)

#Initialise the Moco study and problem
study = track.initialize()
problem = study.updProblem()
        
#Update the weight and exponent on the default control effort goal
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
effort.setWeight(0.001)
effort.setExponent(3)

#Update individual weights in control effort goal to be relative to
#actual muscle and reserve actuator names
#Set appropriate patterns in the weight set
#Muscles
effort.setWeightForControlPattern('/forceset/.*/activation', 0.01) #lower given higher muscle activation scope
#Idealised actuators
effort.setWeightForControlPattern('/forceset/.*_actuator', 0.1)
#Reserves
effort.setWeightForControlPattern('/forceset/.*_reserve', 10)
    
#Add periodicity constraint
periodicityGoal = osim.MocoPeriodicityGoal('symmetryGoal')

#Opposite coordinates are periodic except pelvis anterior-posterior translation
for coordName in taskWeights.keys():
    
    if coordName != 'pelvis_tx':
        
        #Get full path to state
        stateName = osimModel.updCoordinateSet().get(coordName).getAbsolutePathString()+'/value'
        
        #Check for singular periodicity conditions
        if coordName in ['pelvis_ty', 'pelvis_tz', 'pelvis_tilt', 'lumbar_extension']:
            
            #Add state to goal
            periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair(stateName))

        #Check for singular negated conditions
        elif coordName in ['lumbar_bending', 'lumbar_rotation', 'pelvis_list', 'pelvis_rotation']:
            
            #Create goal and negate
            pg = osim.MocoPeriodicityGoalPair(stateName)
            pg.set_negate(True)
            periodicityGoal.addStatePair(pg)
            
        #Set remaining periodicity goals to oppose one another left-right
        elif coordName.endswith('_r'):
            
            #Add state pair to goal for right side
            periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair(stateName, stateName.replace('_r/','_l/')))
            
        elif coordName.endswith('_l'):
            
            #Add state pair to goal for right side
            periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair(stateName, stateName.replace('_l/','_r/')))
            
#Opposite muscle activations are periodic
for stateName in trackingStates.getColumnLabels():
    
    #Right side muscles
    if stateName.endswith('_r/activation'):
        
        #Add state pair to goal for right side
        periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair(stateName, stateName.replace('_r/','_l/')))
        
    #Left side muscles
    elif stateName.endswith('_l/activation'):
        
        #Add state pair to goal for right side
        periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair(stateName, stateName.replace('_l/','_r/')))
        
#Add lumbar and arm controls to periodicity constraint

#Set control list
periodicControls = ['/forceset/lumbar_extension_actuator',
                    '/forceset/lumbar_bending_actuator',
                    '/forceset/lumbar_rotation_actuator',
                    '/forceset/arm_flex_r_actuator',
                    '/forceset/arm_add_r_actuator',
                    '/forceset/arm_rot_r_actuator',
                    '/forceset/elbow_flex_r_actuator',
                    '/forceset/pro_sup_r_actuator',
                    '/forceset/arm_flex_l_actuator',
                    '/forceset/arm_add_l_actuator',
                    '/forceset/arm_rot_l_actuator',
                    '/forceset/elbow_flex_l_actuator',
                    '/forceset/pro_sup_l_actuator']

#Add periodic control parameters
for controlName in periodicControls:
    #Right side
    if '_r_' in controlName:
        periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair(controlName, controlName.replace('_r/','_l/')))
    #Left side
    elif '_l_' in controlName:
        periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair(controlName, controlName.replace('_l/','_r/')))

#Add to problem
problem.addGoal(periodicityGoal)

#Create a speed goal to match the desired speed
speedGoal = osim.MocoAverageSpeedGoal('speed')
speedGoal.setWeight(1)
speedGoal.set_desired_average_speed(sprintSpeed)

#Add to the problem
problem.addGoal(speedGoal)

#Set the time bounds in the problem
#Here we constrain the end time to be relative to the increase in speed

# #Get average speed from pelvis_tx speed
# previousSpeed = trackingStates.getDependentColumn('/jointset/ground_pelvis/pelvis_tx/speed').to_numpy().mean()
# speedFactor = previousSpeed / sprintSpeed

# #Calculate desired new time duration
# initialTime = trackingStates.getIndependentColumn()[0]
# finalTime = trackingStates.getIndependentColumn()[-1]
# timeDur = finalTime - initialTime
# newTimeDur = timeDur * speedFactor

#Set time bound
# problem.setTimeBounds(initialTime, initialTime + newTimeDur)
problem.setTimeBounds(trackingStates.getIndependentColumn()[0],
                      trackingStates.getIndependentColumn()[-1])

#Constrain joint coordinates to bounds relative to tracking targets
#Constrain joint speeds to infinite values

#### TODO: if using previous solution as tracking - do we need to use wider total bounds?

#Loop through coordinates
for coordInd in range(osimModel.updCoordinateSet().getSize()):
    
    #Get coordinate name and absolute path
    coordName = osimModel.updCoordinateSet().get(coordInd).getName()
    coordPath = osimModel.updCoordinateSet().get(coordName).getAbsolutePathString()
    
    #Set joint speeds to infinite total bounds
    problem.setStateInfo(coordPath+'/speed', [np.inf*-1, np.inf])
    
    #Check for pelvis_tx
    if coordName != 'pelvis_tx':
    
        #Get minimum and maximum value from the tracking states and calculate range
        minVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy().min()
        maxVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy().max()
        rangeVal = np.diff((minVal,maxVal))[0]
        
        #Get the starting value from the tracked states
        initialVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy()[0]
        
        #Calculate total and initial bound ranges
        #Calculate and set the 25% range either side of the min and maximum
        #Set the initial value to be within 10% of the starting value
        totalBounds = [minVal - (rangeVal * 0.25), maxVal + (rangeVal * 0.25)]
        initialBounds = [initialVal - (np.abs(initialVal) * 0.1), initialVal + (np.abs(initialVal) * 0.1)]
        
        #Put check in place to correct initial bounds if outside of total ranges
        correctedInitialBounds = []
        #Lower bound
        if initialBounds[0] > totalBounds[0]:
            correctedInitialBounds.append(initialBounds[0])
        else:
            correctedInitialBounds.append(totalBounds[0])
        #Upper bound
        if initialBounds[1] < totalBounds[1]:
            correctedInitialBounds.append(initialBounds[1])
        else:
            correctedInitialBounds.append(totalBounds[1])
    
        #Set in problem
        problem.setStateInfo(coordPath+'/value',
                              #Total bounds range
                              totalBounds,
                              #Initial bounds range
                              correctedInitialBounds
                              )

#Add contact tracking goal using upscaled GRF data
#Uses three separate vectors to apply different contact tracking weights

#Set right and left contact sphere groups
rightFootContacts = []
leftFootContacts = []
for forceInd in range(osimModel.updForceSet().getSize()):
    if 'contact_' in osimModel.updForceSet().get(forceInd).getAbsolutePathString():
        if osimModel.updForceSet().get(forceInd).getAbsolutePathString().endswith('_r'):
            rightFootContacts.append(osimModel.updForceSet().get(forceInd).getAbsolutePathString())
        elif osimModel.updForceSet().get(forceInd).getAbsolutePathString().endswith('_l'):
            leftFootContacts.append(osimModel.updForceSet().get(forceInd).getAbsolutePathString())

#Set left foot tracking parameters for cut plant foot
forcesLeftFoot = osim.StdVectorString()
for contactLabel in leftFootContacts:
    forcesLeftFoot.append(contactLabel)
trackLeftGRF = osim.MocoContactTrackingGoalGroup(forcesLeftFoot, 'LeftGRF')
trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l')        

#Set right foot tracking parameters for no GRFs
forcesRightFoot = osim.StdVectorString()
for contactLabel in rightFootContacts:
    forcesRightFoot.append(contactLabel)
trackRightGRF = osim.MocoContactTrackingGoalGroup(forcesRightFoot, 'RightGRF')
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r')

#Set list to create tracking goals in to avoid overwriting
contactTracking = []

#Create contact settings dictionary
contactTrackingSettings = {'vector': [osim.Vec3(1,0,0),
                                      osim.Vec3(0,1,0),
                                      osim.Vec3(0,0,1)],
                            # 'weight': [1e-3, 1e-3, 1e-3]
                            'weight': [1e-2, 1e-2, 1e-2]
                           }

#Loop through contact tracking settings to create separate weighted goals
for ii in range(len(contactTrackingSettings['vector'])):
    
    #Create tracking goal
    contactTracking.append(osim.MocoContactTrackingGoal(f'contact{ii}',
                                                        contactTrackingSettings['weight'][ii]))
    
    #Set external loads
    contactTracking[ii].setExternalLoadsFile(f'..\\flashSimSpeed{str(sprintSpeed)}\\grfToTrack.xml')
    
    #Add the left and right tracking groups to the goal
    contactTracking[ii].addContactGroup(trackLeftGRF)
    contactTracking[ii].addContactGroup(trackRightGRF)
    
    #Set projection vector in problem
    contactTracking[ii].setProjection('vector')
    contactTracking[ii].setProjectionVector(contactTrackingSettings['vector'][ii])

    #Add to problem
    problem.addGoal(contactTracking[ii])

#Define the solver for the coarse first run through
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_optim_constraint_tolerance(1) #### relaxed this a little bit given scenario
solver.set_optim_convergence_tolerance(0.01)
solver.set_multibody_dynamics_mode('explicit')
solver.set_num_mesh_intervals(5)
solver.resetProblem(problem)

#Set guess from previous tracking simulation

##### TODO: get guess and fill activation states - this would reset the tracked
##### states being applied just by setting a guess file here
##### Get the guess, edit and then save to file

#Check for first sprint speed
if speedList.index(sprintSpeed) == 0:
    
    #Grab and save the full tracking simulation as the guess
    initialGuess = osim.MocoTrajectory('..\\..\\expSims\\fullTrackingSim\\fullTrackingSimSolution.sto')
    
    
else:
    
    #Grab the current guess and fill the activations and controls states
    initialGuess = solver.getGuess()
    
    #Get the previous solution
    previousSolution = osim.MocoTrajectory(f'..\\flashSimSpeed{str(speedList[speedList.index(sprintSpeed)-1])}\\flashSimSpeed{str(speedList[speedList.index(sprintSpeed)-1])}Solution.sto')
    
    #Fill the additional states from the previous solution
    #States
    for stateName in initialGuess.getStateNames():
        if stateName.endswith('/activation'):
            initialGuess.setState(stateName, previousSolution.getState(stateName).to_numpy())
    #Controls
    for controlName in initialGuess.getControlNames():
        initialGuess.setControl(controlName, previousSolution.getControl(controlName).to_numpy())

#Write the initial guess to file and set in problem
initialGuess.write(f'..\\flashSimSpeed{str(sprintSpeed)}\\initialGuess.sto')
solver.setGuessFile(f'..\\flashSimSpeed{str(sprintSpeed)}\\initialGuess.sto')
   
#Final reset of problem after setting guess
solver.resetProblem(problem)

#Solve the tracking problem
solution = study.solve()
# study.visualize(solution)

#Save solution to file
solution.write(f'..\\flashSimSpeed{str(sprintSpeed)}\\flashSimSpeed{str(sprintSpeed)}Solution_coarse.sto')

#Extract predicted GRFs

#Set right and left contact sphere groups
forcesRightFoot = osim.StdVectorString()
forcesLeftFoot = osim.StdVectorString()
for forceInd in range(osimModel.updForceSet().getSize()):
    if 'contact_' in osimModel.updForceSet().get(forceInd).getAbsolutePathString():
        if osimModel.updForceSet().get(forceInd).getAbsolutePathString().endswith('_r'):
            forcesRightFoot.append(osimModel.updForceSet().get(forceInd).getAbsolutePathString())
        elif osimModel.updForceSet().get(forceInd).getAbsolutePathString().endswith('_l'):
            forcesLeftFoot.append(osimModel.updForceSet().get(forceInd).getAbsolutePathString())

#Extract forces from contact spheres in solution
externalForcesTableFlat = osim.createExternalLoadsTableForGait(osimModel, solution,
                                                               forcesRightFoot, forcesLeftFoot)

#Write table to file
osim.STOFileAdapter().write(externalForcesTableFlat, f'..\\flashSimSpeed{str(sprintSpeed)}\\flashSimSpeed{str(sprintSpeed)}Solution_coarse_grf.sto')

#Reset the solver to progress to a finer mesh

#Reset solver intervals
solver.set_num_mesh_intervals(25)

#Set initial guess to coarse solution
solver.setGuessFile(f'..\\flashSimSpeed{str(sprintSpeed)}\\flashSimSpeed{str(sprintSpeed)}Solution_coarse.sto')
solver.resetProblem(problem)

#Re-solve on finer grid
solutionFiner = study.solve()
# study.visualize(solutionFiner)

#Save solution to file
solutionFiner.write(f'..\\flashSimSpeed{str(sprintSpeed)}\\flashSimSpeed{str(sprintSpeed)}Solution.sto')

#Extract forces from contact spheres in solution
externalForcesTableFlat = osim.createExternalLoadsTableForGait(osimModel, solutionFiner,
                                                               forcesRightFoot, forcesLeftFoot)

#Write table to file
osim.STOFileAdapter().write(externalForcesTableFlat, f'..\\flashSimSpeed{str(sprintSpeed)}\\flashSimSpeed{str(sprintSpeed)}Solution_grf.sto')

#Copy tracked states file over to main directory
shutil.move(f'flashSimSpeed{str(sprintSpeed)}_tracked_states.sto',
            f'..\\flashSimSpeed{str(sprintSpeed)}\\flashSimSpeed{str(sprintSpeed)}_tracked_states.sto')

#Extract muscle forces from solution

#Create output paths for muscle forces
outputPaths = osim.StdVectorString()
outputPaths.append('.*tendon_force')

#Analyze study to get muscle forces through tendon
outputTable = study.analyze(solutionFiner, outputPaths)

#Write muscle forces to STO file
osim.STOFileAdapter().write(outputTable, f'..\\flashSimSpeed{str(sprintSpeed)}\\flashSimSpeed{str(sprintSpeed)}Solution_muscleForces.sto')

# %% ----- End of runSimulations.py ----- %% #