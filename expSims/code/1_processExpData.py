# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 20:25:31 2020

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This script:
        1) Converts the data provided by Dorn et al. (2012) into 
           OpenSim format for futher processing
           
        2) Processes the 3D marker data from the sprint trial to get appropriate
           sprinting kinematics
    
"""

# %% Import packages

import opensim as osim
import numpy as np
from scipy.signal import butter, filtfilt
import os
import osimfunctions as helper

# %% Modify any details in here

"""

NOTE: this section of code can be modified to suit another computer set-up

"""

#Set opensim install path directory (e.g. C:\\OpenSim 4.4)
opensimPath = os.path.join('C:', os.sep, 'OpenSim 4.4')

# %% Set-up

#Add OpenSim geometry path (weird issues with this on new laptop)
osim.ModelVisualizer.addDirToGeometrySearchPaths(os.path.join(opensimPath, 'Geometry'))

#Set-up initial logger for converting experimental data
osim.Logger.removeFileSink()
osim.Logger.addFileSink(os.path.join('..', 'data', 'processDataLog.log'))

# %% 1) Convert c3d data to OpenSim formats

# %% Static trial

#Set static file
staticFile = os.path.join('..','data','JA1Static05.c3d')

#Construct opensim 3d object
c3dFile = osim.C3DFileAdapter()
c3dFile.setLocationForForceExpression(osim.C3DFileAdapter.ForceLocation_CenterOfPressure)

#Read in the static trial
staticC3D = c3dFile.read(staticFile)

#Get markers table
staticMarkers = c3dFile.getMarkersTable(staticC3D)

#Rotate marker data
#Create the two rotations needed
markerRot1 = osim.Rotation(np.deg2rad(-90), osim.Vec3(0,0,1))
markerRot2 = osim.Rotation(np.deg2rad(-90), osim.Vec3(1,0,0))
#Rotate the data
for iRow in range(staticMarkers.getNumRows()):
    #Apply the two rotations
    staticMarkers.setRowAtIndex(iRow, markerRot1.multiply(staticMarkers.getRowAtIndex(iRow)))
    staticMarkers.setRowAtIndex(iRow, markerRot2.multiply(staticMarkers.getRowAtIndex(iRow)))
    
#Remove black markers with '*' notation
#Loop through marker labels
markerLabels = staticMarkers.getColumnLabels()
for mm in range(len(markerLabels)):
    #Check if current marker starts with *
    if markerLabels[mm].startswith('*'):
        #Remove it
        staticMarkers.removeColumn(markerLabels[mm])
    
#Write static markers to TRC file
osim.TRCFileAdapter().write(staticMarkers,
                            os.path.join('..','data','static.trc'))

# %% Dynamic trial

#Set the dynamic file name
dynamicFile = os.path.join('..','data','JA1Gait35_9ms.c3d')
    
#Construct opensim 3d object
c3dFile = osim.C3DFileAdapter()
c3dFile.setLocationForForceExpression(osim.C3DFileAdapter.ForceLocation_CenterOfPressure)

#Read in the c3d file
dynamicC3D = c3dFile.read(dynamicFile)

#Get markers table
dynamicMarkers = c3dFile.getMarkersTable(dynamicC3D)

#Rotate the data
#Use the same two rotations as earlier
for iRow in range(dynamicMarkers.getNumRows()):
    #Apply the two rotations
    dynamicMarkers.setRowAtIndex(iRow, markerRot1.multiply(dynamicMarkers.getRowAtIndex(iRow)))
    dynamicMarkers.setRowAtIndex(iRow, markerRot2.multiply(dynamicMarkers.getRowAtIndex(iRow)))
    
#Remove black markers with '*' notation
#Loop through marker labels
markerLabels = dynamicMarkers.getColumnLabels()
for mm in range(len(markerLabels)):
    #Check if current marker starts with *
    if markerLabels[mm].startswith('*'):
        #Remove it
        dynamicMarkers.removeColumn(markerLabels[mm])
    
#Write markers to TRC file
osim.TRCFileAdapter().write(dynamicMarkers,
                            os.path.join('..','data','sprint.trc'))

#Extract the GRF data

#Get forces table
dynamicForces = c3dFile.getForcesTable(dynamicC3D)

#Rotate forces data
#Use the same rotations as earlier
for iRow in range(dynamicForces.getNumRows()):
    dynamicForces.setRowAtIndex(iRow, markerRot1.multiply(dynamicForces.getRowAtIndex(iRow)))
    dynamicForces.setRowAtIndex(iRow, markerRot2.multiply(dynamicForces.getRowAtIndex(iRow)))
    
#Flatten forces data
forcesFlat = dynamicForces.flatten()

#Convert to numpy array
#Pre-allocate numpy array based on data size
dataArray = np.zeros((forcesFlat.getNumRows(),
                      forcesFlat.getNumColumns()))
#Extract data
for forceInd in range(forcesFlat.getNumColumns()):
    dataArray[:,forceInd] = forcesFlat.getDependentColumn(forcesFlat.getColumnLabels()[forceInd]).to_numpy()
    
#Replace nan's for COP and moment data with zeros
np.nan_to_num(dataArray, copy = False, nan = 0.0)

#Convert force point data from mm to m
for forceName in list(forcesFlat.getColumnLabels()):
    if forceName.startswith('p') or forceName.startswith('m'):
        #Get force index
        forceInd = list(forcesFlat.getColumnLabels()).index(forceName)
        #Convert to m units in data array
        dataArray[:,forceInd] = dataArray[:,forceInd] / 1000

#Filter force data
    
#Get the sampling rate
fs = float(dynamicForces.getTableMetaDataAsString('DataRate'))

#Define low-pass Butterworth filter
filtFreq = 50
nyq = 0.5 * fs
normCutoff = filtFreq / nyq
b, a = butter(4, normCutoff, btype = 'low', analog = False)

#Apply lowpass filter to data
for forceName in list(forcesFlat.getColumnLabels()):
    #Get force index
    forceInd = list(forcesFlat.getColumnLabels()).index(forceName)
    #Apply filter
    dataArray[:,forceInd] = filtfilt(b, a, dataArray[:,forceInd])
    
#Build the new time series table
forcesStorage = osim.Storage()

#Get the time data
time = forcesFlat.getIndependentColumn()

#Create maps to replace text from force labels with
#Force plate and type identifiers
forceType = {}
for ii in range(1,9):
    forceType[f'f{ii}'] = f'ground_force_{ii}_v'
    forceType[f'p{ii}'] = f'ground_force_{ii}_p'
    forceType[f'm{ii}'] = f'ground_force_{ii}_m'
#Axis identifiers
forceAxis = {'1': 'x',
             '2': 'y',
             '3': 'z'}

#Set labels in table
newLabels = osim.ArrayStr()
newLabels.append('time')
for forceLabel in forcesFlat.getColumnLabels():
    #Split the label to get parts
    labelSplit = forceLabel.split('_')
    #Create new label
    forceLabel = f'{forceType[labelSplit[0]]}{forceAxis[labelSplit[1]]}'
    #Append to labels vector
    newLabels.append(forceLabel)
forcesStorage.setColumnLabels(newLabels)

#Add data
for iRow in range(dataArray.shape[0]):
    row = osim.ArrayDouble()
    for iCol in range(dataArray.shape[1]):
        row.append(dataArray[iRow,iCol])
    #Add data to storage
    forcesStorage.append(time[iRow], row)

#Set name for storage object
forcesStorage.setName('sprint_grf')

#Write to file
forcesStorage.printResult(forcesStorage, 'sprint_grf', os.path.join('..','data'), 0.001, '.mot')

#Create the external loads .xml file
#Set the different force plates to the varying foot contacts
#have which contacts
#Right foot = fp2,fp7
#Left foot = fp4,fp8
forceXML = osim.ExternalLoads()

#Create and append the right GRF external forces
#FP2
rightGRF1 = osim.ExternalForce()
rightGRF1.setName('RightGRF1')
rightGRF1.setAppliedToBodyName('calcn_r')
rightGRF1.setForceExpressedInBodyName('ground')
rightGRF1.setPointExpressedInBodyName('ground')
rightGRF1.setForceIdentifier('ground_force_2_v')
rightGRF1.setPointIdentifier('ground_force_2_p')
rightGRF1.setTorqueIdentifier('ground_force_2_m')
forceXML.cloneAndAppend(rightGRF1)
#FP7
rightGRF2 = osim.ExternalForce()
rightGRF2.setName('RightGRF2')
rightGRF2.setAppliedToBodyName('calcn_r')
rightGRF2.setForceExpressedInBodyName('ground')
rightGRF2.setPointExpressedInBodyName('ground')
rightGRF2.setForceIdentifier('ground_force_7_v')
rightGRF2.setPointIdentifier('ground_force_7_p')
rightGRF2.setTorqueIdentifier('ground_force_7_m')
forceXML.cloneAndAppend(rightGRF2)

#Create and append the left GRF external forces
#FP4
leftGRF1 = osim.ExternalForce()
leftGRF1.setName('LeftGRF1')
leftGRF1.setAppliedToBodyName('calcn_l')
leftGRF1.setForceExpressedInBodyName('ground')
leftGRF1.setPointExpressedInBodyName('ground')
leftGRF1.setForceIdentifier('ground_force_4_v')
leftGRF1.setPointIdentifier('ground_force_4_p')
leftGRF1.setTorqueIdentifier('ground_force_4_m')
forceXML.cloneAndAppend(leftGRF1)
#FP8
leftGRF2 = osim.ExternalForce()
leftGRF2.setName('LeftGRF2')
leftGRF2.setAppliedToBodyName('calcn_l')
leftGRF2.setForceExpressedInBodyName('ground')
leftGRF2.setPointExpressedInBodyName('ground')
leftGRF2.setForceIdentifier('ground_force_8_v')
leftGRF2.setPointIdentifier('ground_force_8_p')
leftGRF2.setTorqueIdentifier('ground_force_8_m')
forceXML.cloneAndAppend(leftGRF2)

#Set GRF datafile
forceXML.setDataFileName('sprint_grf.mot')

#Write to file
forceXML.printToXML(os.path.join('..','data','sprint_grf.xml'))

# %% 2) Process the 3D sprinting data

"""

Generate a torque-driven marker tracking problem that generates dynamically consistent
sprint kinematics considering the marker and GRF data.

"""

##### SETTINGS ####

#Create a dictionary that provides the task weights for tracking markers
markerWeightVals = {'C7': 2.0, 'RSH': 2.5, 'LSH': 2.5, 'MAN': 2.0, 'T7': 2.0, 'LARM': 1.0,
                    'LELB': 2.0, 'LFOREARM': 1.0, 'LWR': 2.0, 'RARM': 1.0, 'RELB': 2.0, 
                    'RFOREARM': 1.0, 'RWR': 2.0, 'RASI': 2.0, 'LASI': 2.0, 'SACR': 2.5,
                    'LTHLP': 1.0, 'LTHLD': 1.0, 'LTHAP': 1.0, 'LTHAD': 1.0, 'LLEPI': 2.0,
                    'LTIAP': 1.0, 'LTIAD': 1.0, 'LTILAT': 1.0, 'LLMAL': 2.5,
                    'LHEEL': 1.0, 'LMFS': 1.0, 'LMFL': 1.0, 'LP1MT': 2.0, 'LTOE': 2.0,
                    'LP5MT': 2.0, 'RTHLP': 1.0, 'RTHLD': 1.0, 'RTHAP': 1.0, 'RTHAD': 1.0,
                    'RLEPI': 2.0, 'RTIAP': 1.0, 'RTIAD': 1.0, 'RTILAT': 1.0,
                    'RLMAL': 2.5, 'RHEEL': 1.0, 'RMFS': 1.0, 'RMFL': 1.0, 'RP1MT': 2.0,
                    'RTOE': 2.0, 'RP5MT': 2.0,
                    }

#Create dictionary that sets optimal forces for torque model actuators
optForces = {
            #upper body
            'lumbar_extension': 300, 'lumbar_bending': 300, 'lumbar_rotation': 300,
            'arm_flex_r': 300, 'arm_add_r': 300, 'arm_rot_r': 300,
            'elbow_flex_r': 100, 'pro_sup_r': 100,
            'arm_flex_l': 300, 'arm_add_l': 300, 'arm_rot_l': 300,
            'elbow_flex_l': 100, 'pro_sup_l': 100,
            #left limb
            'hip_flexion_l': 300, 'hip_adduction_l': 200, 'hip_rotation_l': 100,
            'knee_angle_l': 300, 'ankle_angle_l': 300, 'subtalar_angle_l': 100, 'mtp_angle_l': 50,
            #right limb
            'hip_flexion_r': 300, 'hip_adduction_r': 200, 'hip_rotation_r': 100,
            'knee_angle_r': 300, 'ankle_angle_r': 300, 'subtalar_angle_r': 100, 'mtp_angle_r': 50,
            #pelvis
            'pelvis_tx': 1, 'pelvis_ty': 1, 'pelvis_tz': 1,
            'pelvis_tilt': 1, 'pelvis_list': 1, 'pelvis_rotation': 1
            }

#Set joints to weld for simulation
jointsToWeld = ['radius_hand_r', 'radius_hand_l']
#Create vector string object
weldVectorStr = osim.StdVectorString()
[weldVectorStr.append(joint) for joint in jointsToWeld]

#Set start and end times in problem
#Get the gait timings from helper function for a full gait cycle
startTime, endTime = helper.getGaitTimings(grfFile = os.path.join('..', 'data', 'sprint_grf.mot'),
                                           extLoads = os.path.join('..', 'data', 'sprint_grf.xml'),
                                           startForceName = 'RightGRF1',
                                           stopForceName = 'RightGRF2',
                                           forceThreshold = 20)

#Calculate mesh interval required for an approximate 0.01 sample step
numMeshInterval = int((np.round((endTime - startTime) / 0.01) - 1) / 2)

##### ADJUST MODEL FOR PROBLEM #####

#Read in model for processing
osimModelProc = osim.ModelProcessor(osim.Model(os.path.join('..', 'model', 'JA1_SCALED_Osim4_Muscles.osim')))

#Append model processors
#Remove muscles
osimModelProc.append(osim.ModOpRemoveMuscles())
#Add external loads
osimModelProc.append(osim.ModOpAddExternalLoads(os.path.join('..', 'data', 'sprint_grf.xml')))
#Weld joints
osimModelProc.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))

#Process model
osimModel = osimModelProc.process()

#Add torque actuators to coordinates
osimModel = helper.addTorqueActuators(osimModel, optForces, np.inf*-1, np.inf)
    
#Finalise model connections
osimModel.finalizeConnections()

##### CREATE TRACKING PROBLEM #####
    
#Define the motion tracking problem
track = osim.MocoTrack()
track.setName('markerTrackingSim')

#Set model
modelProcessor = osim.ModelProcessor(osimModel)
track.setModel(modelProcessor)

#Set the markers reference from TRC file
#Note that this needs to be done using a flattened table due to a bug
#See: https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=1815&t=14612&p=43550&start=0&view=
markerData = osim.TimeSeriesTableVec3(os.path.join('..', 'data', 'sprint.trc'))
markerTableProcessor = osim.TableProcessor(markerData.flatten())
markerTableProcessor.append(osim.TabOpLowPassFilter(12))
track.setMarkersReference(markerTableProcessor)

#Set allow unused references in case there are extra markers
track.set_allow_unused_references(True)

#Set the global markers tracking weight
track.set_markers_global_tracking_weight(10.0)

#Set initial and final times
track.set_initial_time(startTime)
track.set_final_time(endTime)

#Set the tracking weights for markers in a weights object
markerWeights = osim.MocoWeightSet()
for marker in markerWeightVals.keys():
    markerWeights.cloneAndAppend(osim.MocoWeight(marker, markerWeightVals[marker]))

#Add to tracking problem
track.set_markers_weight_set(markerWeights)

#Convert the tracking tool to a Moco study and problem
study = track.initialize()
problem = study.updProblem()

##### GOALS #####

#Set the weight on the default control effort goal
osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort')).setWeight(0.001)

##### CONSTRAINTS #####

#Add constraints to ensure that ground markers do not penetrate floor

#Create a dictionary to store constraints in to avoid over-writing
contactConstraints = {}

#Loop through markers and identify those to add constraints for
for markerInd in range(osimModel.updMarkerSet().getSize()):
    if '_fp_' in osimModel.updMarkerSet().get(markerInd).getName():
                
        #Set marker name to variable
        markerName = osimModel.updMarkerSet().get(markerInd).getName()
        
        #Add the constraint to dictionary
        contactConstraints[markerName] = osim.MocoOutputConstraint()
        
        #Set the constraint name
        contactConstraints[markerName].setName('groundConstraint_'+markerName)
        
        #Set output path to marker location
        contactConstraints[markerName].setOutputPath('/markerset/'+markerName+'|location')
        
        #Set axis index of Vec 3
        #Index 1 relates to the Y-component of the marker location
        contactConstraints[markerName].setOutputIndex(1)
        
        #Set the bounds on the output constraint
        #Create the bounds
        verticalMarkerBounds = osim.StdVectorMocoBounds()
        verticalMarkerBounds.append(osim.MocoBounds(0.0,1.0))
        #Set in constraint
        contactConstraints[markerName].updConstraintInfo().setBounds(verticalMarkerBounds)
        
        #Add to problem
        problem.addPathConstraint(contactConstraints[markerName])

##### SOLVER #####

#Define and reset the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.resetProblem(problem)

#Set solver options
solver.set_optim_max_iterations(1000)
solver.set_num_mesh_intervals(numMeshInterval) #would be too small for contact tracking
# solver.set_minimize_implicit_multibody_accelerations(True) #would this help with contact tracking?
# solver.set_scale_variables_using_bounds(True) #https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=1815&t=17154&p=0&start=0&view=&sid=d92da2432f6fc11cc0ce2beb83ad7672
solver.set_optim_constraint_tolerance(1e-4)
solver.set_optim_convergence_tolerance(1e-2)

##### RUN SOLVER #####
    
#Solve!
solution = study.solve()

#Check to unseal
if solution.isSealed():
    solution.unseal()

# #Review visualisation
# study.visualize(solution)

#Remove the original tracked states file
os.remove('markerTrackingSim_tracked_markers.sto')

#Save tracking sim to file
os.makedirs(os.path.join('..', 'results', 'markerTracking'), exist_ok = True)
solution.write(os.path.join('..', 'results', 'markerTracking', 'markerTracking_fullGaitCycle_solution.sto'))

# %% ----- End of 1_processExpData.py ----- %% #