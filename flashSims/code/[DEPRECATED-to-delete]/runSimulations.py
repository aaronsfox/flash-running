# -*- coding: utf-8 -*-
"""

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    TODO: add notes


"""

# %% Import packages

import opensim as osim
import osimHelper as helper
import numpy as np
import os
import shutil

# %% Set-up

#Add OpenSim geometry path (weird issues with this on new laptop)
osim.ModelVisualizer.addDirToGeometrySearchPaths('C:\\OpenSim 4.3\\Geometry')

# %% Tracking simulation with Flash muscle parameters

"""

TODO: 
    > add notes
    > how to control activation states if not included?

"""

#Create folder to store data
try:
    os.mkdir('..\\flashTrackingSim')
except:
    print('Flash tracking sim folder already detected...')
    
#Create the model to use in this simulation
osimModel = helper.createFlashModel(inputModelFile = '..\\..\\expSims\\data\\JA1_SCALED_Osim40_Muscles.osim',
                                    outputModelFile = '..\\flashTrackingSim\\osimModel.osim',
                                    unilateralMuscles = False,
                                    jointsToWeld = ['subtalar_r', 'subtalar_l', 'radius_hand_r', 'radius_hand_l'],
                                    addMetabolicsModel = True)

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName('flashTrackingSim')

#Set model in tracking tool
trackModelProcessor = osim.ModelProcessor(osimModel)
track.setModel(trackModelProcessor)

#Set the coordinates reference in tracking tool

#Convert the solution data from the tracking sim to coordinates
trackingStates = osim.MocoTrajectory('..\\..\\expSims\\fullTrackingSim\\fullTrackingSimSolution.sto').exportToStatesTable()
osim.STOFileAdapter().write(trackingStates, '..\\flashTrackingSim\\coordinatesToTrack.sto')

#Set coordinates file as reference in tool
track.setStatesReference(osim.TableProcessor('..\\flashTrackingSim\\coordinatesToTrack.sto'))
track.set_states_global_tracking_weight(1)
 
#Set allow unused references in case there are extra markers
track.set_allow_unused_references(True)

#Set speeds a derivatives of coordinate references
track.set_track_reference_position_derivatives(True)

#Set tracked states to guess
track.set_apply_tracked_states_to_guess(True)

#Create weight set for state tracking
stateWeights = osim.MocoWeightSet()

#Create a dictionary that provides the kinematic task weights for function
taskWeights = {'pelvis_tx': 5, 'pelvis_ty': 0.1, 'pelvis_tz': 2.5,
                'pelvis_tilt': 10, 'pelvis_list': 10, 'pelvis_rotation': 10, 
                'hip_flexion_r': 10, 'hip_adduction_r': 5, 'hip_rotation_r': 2.5, 
                'knee_angle_r': 10, 'ankle_angle_r': 2.5,
                # 'subtalar_angle_r': 0,
                'mtp_angle_r': 2.5,
                'hip_flexion_l': 10, 'hip_adduction_l': 5, 'hip_rotation_l': 2.5, 
                'knee_angle_l': 10, 'ankle_angle_l': 2.5,
                # 'subtalar_angle_l': 0,
                'mtp_angle_l': 2.5,
                'lumbar_extension': 15, 'lumbar_bending': 10, 'lumbar_rotation': 5,
                'arm_flex_r': 2.5, 'arm_add_r': 1, 'arm_rot_r': 1,
                'elbow_flex_r': 5, 'pro_sup_r': 1,
                'arm_flex_l': 2.5, 'arm_add_l': 1, 'arm_rot_l': 1,
                'elbow_flex_l': 5, 'pro_sup_l': 1}

#Set constant weight to scale tracking error speeds by
w = 0.001

#Loop through coordinates to apply weights
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
        stateWeights.cloneAndAppend(osim.MocoWeight(f'{coordPath}/speed',
                                                    w * taskWeights[coordName])) 

#Add to tracking problem
track.set_states_weight_set(stateWeights)

#Set start and end times in problem
track.set_initial_time(trackingStates.getIndependentColumn()[0])
track.set_final_time(trackingStates.getIndependentColumn()[-1])

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
effort.setWeightForControlPattern('/forceset/.*/activation', 0.1)
#Idealised actuators
effort.setWeightForControlPattern('/forceset/.*_actuator', 0.1)
#Reserves
effort.setWeightForControlPattern('/forceset/.*_reserve', 10)
#Residuals
effort.setWeightForControlPattern('/forceset/.*_residual', 10)

#Add contact tracking goal
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
trackLeftGRF = osim.MocoContactTrackingGoalGroup(forcesLeftFoot, 'LeftGRF1')
trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l')        

#Set right foot tracking parameters for no GRFs
forcesRightFoot = osim.StdVectorString()
for contactLabel in rightFootContacts:
    forcesRightFoot.append(contactLabel)
trackRightGRF = osim.MocoContactTrackingGoalGroup(forcesRightFoot, 'RightGRF1')
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r')

#Set list to create tracking goals in to avoid overwriting
contactTracking = []

#Create contact settings dictionary
contactTrackingSettings = {'vector': [osim.Vec3(1,0,0),
                                      osim.Vec3(0,1,0),
                                      osim.Vec3(0,0,1)],
                           'weight': [1e-2, 1e-2, 1e-2]
                           }

#Loop through contact tracking settings to create separate weighted goals
for ii in range(len(contactTrackingSettings['vector'])):
    
    #Create tracking goal
    contactTracking.append(osim.MocoContactTrackingGoal(f'contact{ii}',
                                                        contactTrackingSettings['weight'][ii]))
    
    #Set external loads
    contactTracking[ii].setExternalLoadsFile('..\\\..\\expSims\\data\\sprint_grf.xml')
    
    #Add the left and right tracking groups to the goal
    contactTracking[ii].addContactGroup(trackLeftGRF)
    contactTracking[ii].addContactGroup(trackRightGRF)
    
    #Set projection vector in problem
    contactTracking[ii].setProjection('vector')
    contactTracking[ii].setProjectionVector(contactTrackingSettings['vector'][ii])

    #Add to problem
    problem.addGoal(contactTracking[ii])
    
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
            
##### TODO: consider periodic excitations instead if activation dynamics ignored...            

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
    # #Remaining controls
    # else:
    #     periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair(controlName))

#Add to problem
problem.addGoal(periodicityGoal)

#Set kinematic bounds in problem

#Set the joint coordinates to place bounds on
boundedCoords = ['pelvis_tx', 'pelvis_ty', 'pelvis_tz', 'pelvis_tilt', 'pelvis_list', 'pelvis_rotation', 
                  'hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r', 
                  'knee_angle_r', 'ankle_angle_r', 
                  'mtp_angle_r', # 'subtalar_angle_r',
                  'hip_flexion_l', 'hip_adduction_l', 'hip_rotation_l', 
                  'knee_angle_l', 'ankle_angle_l',
                  'mtp_angle_l', # 'subtalar_angle_l',
                  'lumbar_extension', 'lumbar_bending', 'lumbar_rotation',
                  'arm_flex_r', 'arm_add_r', 'arm_rot_r',
                  'elbow_flex_r', 'pro_sup_r',
                  'arm_flex_l', 'arm_add_l', 'arm_rot_l',
                  'elbow_flex_l', 'pro_sup_l']

#Also create a dictionary which limits initial bounds to be relative to starting point
#from the experimental data. We don't set total bounds here as it's possible that
#to meet the constraints the joint coordinates will need to exceed those.

#Loop through coordinates
for coordName in boundedCoords:
    
    #Get the full path to coordinate
    coordPath = osimModel.updCoordinateSet().get(coordName).getAbsolutePathString()
    
    #Get the starting value from the tracked states
    initialVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy()[0]
    
    #Calculate total and initial bound ranges
    #Set the initial value to be within 10% of the starting value
    initialBounds = [initialVal - (np.abs(initialVal) * 0.1), initialVal + (np.abs(initialVal) * 0.1)]
    
    #Set in problem
    problem.setStateInfo(coordPath+'/value',
                          #Total bounds range
                          [],
                          #Initial bounds range
                          initialBounds
                          )
    
#Define the solver for the coarse first run through
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_optim_constraint_tolerance(0.01)
solver.set_optim_convergence_tolerance(0.01)
solver.set_multibody_dynamics_mode('explicit')
solver.set_num_mesh_intervals(5)
solver.resetProblem(problem)

# #Recreate guess file, which requires removing activation states given they no
# #longer exist in the problem with activation dynamics ignored
# newGuess = osim.TimeSeriesTable('..\\..\\expSims\\fullTrackingSim\\fullTrackingSimSolution.sto')
# for colLabel in newGuess.getColumnLabels():
#     if colLabel.endswith('/activation'):
#         newGuess.removeColumn(colLabel)
        
# #Write new guess to file
# osim.STOFileAdapter.write(newGuess, '..\\flashTrackingSim\\initialGuess.sto')

# #Set guess from marker tracking simulation
# solver.setGuessFile('..\\flashTrackingSim\\initialGuess.sto')
# solver.resetProblem(problem)

#Solve the tracking problem
solution = study.solve()
# study.visualize(solution)

#Save solution to file
solution.write('..\\flashTrackingSim\\flashTrackingSimSolution_coarse.sto')

#Extract predicted GRFs

#Extract forces from contact spheres in solution
externalForcesTableFlat = osim.createExternalLoadsTableForGait(osimModel, solution,
                                                               forcesRightFoot, forcesLeftFoot)

#Write table to file
osim.STOFileAdapter().write(externalForcesTableFlat, '..\\flashTrackingSim\\flashTrackingSimSolution_coarse_grf.sto')

#Reset the solver to progress to a finer mesh

#Reset solver intervals
solver.set_num_mesh_intervals(25)

#Set initial guess to coarse solution
solver.setGuessFile('..\\flashTrackingSim\\flashTrackingSimSolution_coarse.sto')
solver.resetProblem(problem)

#Re-solve on finer grid
solutionFiner = study.solve()
# study.visualize(solutionFiner)

#Copy tracked states file over to main directory
shutil.move('flashTrackingSim_tracked_states.sto',
            '..\\flashTrackingSim\\flashTrackingSim_tracked_states.sto')

# %% Run the Flash simulations - TEST TRACKING APPROACH...

"""

Notes?
    > Seems that tracking all states loosely-ish might be the best way - even pelvis_tx
        >> This just basically means the similar movement will be done - but just quicker
    > Approach didn't really work for 10 up to 20 m/s - maybe because an initial guess was used?
        >> What happens without a guess?
            >> Mesh interval of 5 solved in 18 minutes
            >> Solved in a similar way
        >> What if we constraint pelvis_tx to start at the same point but progress to the point based on speed?
            >> This seems to work OK for 20m/s with a coarse mesh (20 minute solution)
                >>> Progress to finer grid?
                    >>>> Seems to solve similarly (30 minutes)
            >> Possibly the most workable solution - but does somewhat constrain simulations to longer, rather than faster strides
        >> The other option is to constraint pelvis_tx motion to the same distance travelled, but reduce the time duration
            >> This subsequently enforces stride rate increases...
            >> Solves probably easier, but kinematics can go off
                >>> Maybe higher state weights as per tracking?
                >>> And probably take out speeds too...
                    >>>> Seemed to solve easier, possibly due to not requiring coordinate speed tracking
                    >>>> Looked more like original too...
                >>> coarse mesh solution as initial guess seemed to work better for finer mesh
            

TODO:
    > Evaluate outcomes of this approach...
    > Does this continually need a coarse mesh - finer meashes seem to struggle?
    > If this works, do we just start with this at matches speed as the "tracking sim"?

"""

#Set a desired speed
sprintSpeed = 20

#Create folder to store data
try:
    os.mkdir(f'..\\flashSimSpeed{str(sprintSpeed)}')
except:
    print('Flash sim folder for current speed already detected...')
    
#Create the model to use in this simulation
osimModel = helper.createFlashModel(inputModelFile = '..\\..\\expSims\\data\\JA1_SCALED_Osim40_Muscles.osim',
                                    outputModelFile = f'..\\flashSimSpeed{str(sprintSpeed)}\\osimModel.osim',
                                    unilateralMuscles = False,
                                    jointsToWeld = ['subtalar_r', 'subtalar_l', 'radius_hand_r', 'radius_hand_l'],
                                    addMetabolicsModel = True)

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName(f'flashSimSpeed{str(sprintSpeed)}')

#Set model in tracking tool
trackModelProcessor = osim.ModelProcessor(osimModel)
track.setModel(trackModelProcessor)

#Set the coordinates reference in tracking tool

#Convert the solution data from the tracking sim to coordinates
trackingStates = osim.MocoTrajectory('..\\flashTrackingSim\\flashTrackingSimSolution_coarse.sto').exportToStatesTable()
osim.STOFileAdapter().write(trackingStates, f'..\\flashSimSpeed{str(sprintSpeed)}\\coordinatesToTrack.sto')

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
taskWeights = {'pelvis_tx': 1, 'pelvis_ty': 1, 'pelvis_tz': 1,
                'pelvis_tilt': 1, 'pelvis_list': 1, 'pelvis_rotation': 1, 
                'hip_flexion_r': 1, 'hip_adduction_r': 1, 'hip_rotation_r': 1, 
                'knee_angle_r': 1, 'ankle_angle_r': 1,
                # 'subtalar_angle_r': 0,
                'mtp_angle_r': 1,
                'hip_flexion_l': 1, 'hip_adduction_l': 1, 'hip_rotation_l': 1, 
                'knee_angle_l': 1, 'ankle_angle_l': 1,
                # 'subtalar_angle_l': 0,
                'mtp_angle_l': 1,
                'lumbar_extension': 1, 'lumbar_bending': 1, 'lumbar_rotation': 1,
                'arm_flex_r': 1, 'arm_add_r': 1, 'arm_rot_r': 1,
                'elbow_flex_r': 1, 'pro_sup_r': 1,
                'arm_flex_l': 1, 'arm_add_l': 1, 'arm_rot_l': 1,
                'elbow_flex_l': 1, 'pro_sup_l': 1}

# #Set constant weight to scale tracking error speeds by
# w = 0.001

#Loop through coordinates to apply weights
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
        # stateWeights.cloneAndAppend(osim.MocoWeight(f'{coordPath}/speed',
        #                                             w * taskWeights[coordName])) 

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
effort.setWeightForControlPattern('/forceset/.*/activation', 0.1)
#Idealised actuators
effort.setWeightForControlPattern('/forceset/.*_actuator', 0.1)
#Reserves
effort.setWeightForControlPattern('/forceset/.*_reserve', 10)
#Residuals
effort.setWeightForControlPattern('/forceset/.*_residual', 10)
    
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
    # #Remaining controls
    # else:
    #     periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair(controlName))

#Add to problem
problem.addGoal(periodicityGoal)

#Create a speed goal to match the desired speed
speedGoal = osim.MocoAverageSpeedGoal('speed')
speedGoal.setWeight(1)
speedGoal.set_desired_average_speed(sprintSpeed)

#Add to the problem
problem.addGoal(speedGoal)

#Set the time bounds in the problem
#The initial time is locked to the same as the tracking simulation
#No final time bounds are set as this would limit the stride style (i.e. if running
#a faster speed, locking the time does not permit shorter strides)
# problem.setTimeBounds(trackingStates.getIndependentColumn()[0], [])

####################### CONSTRIANED TIME BOUNDS APPROACH ######################

#Set previous speed and calculate factor
previousSpeed = 10
speedFactor = previousSpeed / sprintSpeed

#Calculate desired new time duration
initialTime = trackingStates.getIndependentColumn()[0]
finalTime = trackingStates.getIndependentColumn()[-1]
timeDur = finalTime - initialTime
newTimeDur = timeDur * speedFactor

#Set time bounds
problem.setTimeBounds(initialTime, initialTime + newTimeDur)

####################### CONSTRIANED TIME BOUNDS APPROACH ######################

#Set kinematic bounds in problem
#Here we allow a 25% window of the total range either side of the max and min value
#from the tracking simulation. We also set the initial bounds to be within 10% of
#the initial values from the tracking simulation. This needs to ignore pelis_tx,
#or at least provide wider bounds when speed is increased.
trackingData = osim.MocoTrajectory('..\\flashTrackingSim\\flashTrackingSimSolution_coarse.sto')

#Loop through coordinates
for coordInd in range(osimModel.updCoordinateSet().getSize()):
    
    #Get coordinate name and absolute path
    coordName = osimModel.updCoordinateSet().get(coordInd).getName()
    coordPath = osimModel.updCoordinateSet().get(coordName).getAbsolutePathString()
    
    #Check for pelvis_tx
    if coordName != 'pelvis_tx':
    
        #Get minimum and maximum value from the tracking states and calculate range
        minVal = trackingData.getState(coordPath+'/value').to_numpy().min()
        maxVal = trackingData.getState(coordPath+'/value').to_numpy().max()
        rangeVal = np.diff((minVal,maxVal))[0]
        
        #Get the starting value from the tracked states
        initialVal = trackingData.getState(coordPath+'/value').to_numpy()[0]
        
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
    
########################## BOUNDING PELVIS_TX? ################################

# #Get coordinate absolute path
# coordPath = osimModel.updCoordinateSet().get('pelvis_tx').getAbsolutePathString()

# #Get how far pelvis travelled in previous speed
# initialVal = trackingData.getState(coordPath+'/value').to_numpy()[0]
# finalVal = trackingData.getState(coordPath+'/value').to_numpy()[-1]
# pelvisTranslation = finalVal - initialVal

# #Set previous speed and calculate factor
# previousSpeed = 10
# speedFactor = sprintSpeed / previousSpeed

# #Calculate desired new pelvis translation
# newPelvisTranslation = pelvisTranslation * speedFactor

# #Set bounds for pelvis translation
# problem.setStateInfo(coordPath+'/value',
#                      #Total bounds range
#                      [initialVal, initialVal + newPelvisTranslation],
#                      #Initial bounds
#                      initialVal,
#                      #Final bounds
#                      initialVal + newPelvisTranslation
#                      )

########################## BOUNDING PELVIS_TX? ################################

####################### BOUNDING PELVIS_TX_SPEED? #############################

#Get coordinate absolute path
coordPath = osimModel.updCoordinateSet().get('pelvis_tx').getAbsolutePathString()

#Set bounds for pelvis translation
problem.setStateInfo(coordPath+'/speed',
                     #Total bounds range
                     [],
                     #Initial bounds
                     sprintSpeed,
                     #Final bounds
                     sprintSpeed
                     )

####################### BOUNDING PELVIS_TX_SPEED? #############################

#Define the solver for the coarse first run through
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_optim_constraint_tolerance(0.01)
solver.set_optim_convergence_tolerance(0.01)
solver.set_multibody_dynamics_mode('explicit')
solver.set_num_mesh_intervals(25)
solver.resetProblem(problem)

#Set guess from marker tracking simulation
solver.setGuessFile('..\\flashTrackingSim\\flashTrackingSimSolution_coarse.sto')
solver.resetProblem(problem)

#Solve the tracking problem
solution = study.solve()
# study.visualize(solution)

#Save solution to file
solution.write(f'..\\flashSimSpeed{str(sprintSpeed)}\\flashSimSpeed{str(sprintSpeed)}Solution_tracking_timeTestVersionFiner.sto')

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
osim.STOFileAdapter().write(externalForcesTableFlat, f'..\\flashSimSpeed{str(sprintSpeed)}\\flashSimSpeed{str(sprintSpeed)}Solution_tracking_timeTestVersionFiner_grf.sto')

#Copy tracked states file over to main directory
shutil.move(f'flashSimSpeed{str(sprintSpeed)}_tracked_states.sto',
            f'..\\flashSimSpeed{str(sprintSpeed)}\\flashSimSpeed{str(sprintSpeed)}_tracked_states.sto')

# %% Run the Flash simulations - PREDICTIVE APPROACH...

"""

TODO:
    > Same as above - are activations being considered in cost function?

"""

#Set the baseline speed from the tracking simulations
trackingData = osim.MocoTrajectory('..\\flashTrackingSim\\flashTrackingSimSolution.sto')

#Get distance travelled based on pelvis translation alongside time taken
distTravelled = trackingData.getState('/jointset/ground_pelvis/pelvis_tx/value').to_numpy()[-1] - trackingData.getState('/jointset/ground_pelvis/pelvis_tx/value').to_numpy()[0]
timeTaken = trackingData.getTime().to_numpy()[-1] - trackingData.getTime().to_numpy()[0]

#Calculate average speed
baseSpeed = distTravelled / timeTaken

#Set a desired speed
sprintSpeed = 10

#Create folder to store data
try:
    os.mkdir(f'..\\flashSimSpeed{str(sprintSpeed)}')
except:
    print('Flash sim folder for current speed already detected...')
    
#Create the model to use in this simulation
osimModel = helper.createFlashModel(inputModelFile = '..\\..\\expSims\\data\\JA1_SCALED_Osim40_Muscles.osim',
                                    outputModelFile = f'..\\flashSimSpeed{str(sprintSpeed)}\\osimModel.osim',
                                    unilateralMuscles = False,
                                    jointsToWeld = ['subtalar_r', 'subtalar_l', 'radius_hand_r', 'radius_hand_l'],
                                    addMetabolicsModel = True)

#Define the predictive study problem
study = osim.MocoStudy()
study.setName(f'flashSimSpeed{str(sprintSpeed)}')

#Define the problem
problem = study.updProblem()

#Set model in problem
problem.setModelProcessor(osim.ModelProcessor(osimModel))

#Create a speed goal to match the desired speed
speedGoal = osim.MocoAverageSpeedGoal('speed')
speedGoal.setWeight(1)
speedGoal.set_desired_average_speed(sprintSpeed)

#Add to the problem
problem.addGoal(speedGoal)

#Create a peridocity goal similar to the tracking sim which will permit the simulation
#of a symmetrical half gait cycle

#Add periodicity constraint
periodicityGoal = osim.MocoPeriodicityGoal('symmetryGoal')

#Create symmetrical coordinate goals based on model coordinate set
for coordInd in range(osimModel.updCoordinateSet().getSize()):
    
    #Get coordinate name
    coordName = osimModel.updCoordinateSet().get(coordInd).getName()

    #Opposite coordinates are periodic except pelvis anterior-posterior translation    
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
            
#Create symmetrical muscle activation states based on state names
osimModel.initSystem()
modelStates = osimModel.getStateVariableNames()

# #Loop through states
# for stateInd in range(modelStates.getSize()):
    
#     #Get state name
#     stateName = modelStates.get(stateInd)
            
#     #Opposite muscle activations are periodic
    
#     #Right side muscles
#     if stateName.endswith('_r/activation'):
        
#         #Add state pair to goal for right side
#         periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair(stateName, stateName.replace('_r/','_l/')))
        
#     #Left side muscles
#     elif stateName.endswith('_l/activation'):
        
#         #Add state pair to goal for right side
#         periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair(stateName, stateName.replace('_l/','_r/')))
        
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
    # #Remaining controls
    # else:
    #     periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair(controlName))

#Add to problem
problem.addGoal(periodicityGoal)

#Create a goal to minimise muscle and actuator controls and set parameters
effortGoal = osim.MocoControlGoal('control_effort', 0.001) #set the weight when creating
effortGoal.setExponent(3)

#Update individual weights in control effort goal to be relative to
#actual muscle and reserve actuator names
#Set appropriate patterns in the weight set
#Muscles
effortGoal.setWeightForControlPattern('/forceset/.*/activation', 0.01)
#Idealised actuators
effortGoal.setWeightForControlPattern('/forceset/.*_actuator', 0.1)
#Reserves
effortGoal.setWeightForControlPattern('/forceset/.*_reserve', 10)
#Residuals
effortGoal.setWeightForControlPattern('/forceset/.*_residual', 10)

#Add muscle excitations to effort goal
for muscInd in range(osimModel.getMuscles().getSize()):
    #Get muscle name
    muscName = osimModel.getMuscles().get(muscInd).getName()
    #Set weight for control
    effortGoal.setWeightForControl(f'/forceset/{muscName}', 0.01)

#Add to problem
problem.addGoal(effortGoal)

#Add a low weighted states tracking for joint coordinates
#Exclude speeds and pelvis_tx as these will change with altered speed
kinematicGoal = osim.MocoStateTrackingGoal('kinematicTracking', 0.1)

#Set kinematic data to track
tableProcessor = osim.TableProcessor('..\\flashTrackingSim\\flashTrackingSimSolution.sto')
tableProcessor.append(osim.TabOpUseAbsoluteStateNames())
kinematicGoal.setReference(tableProcessor)
kinematicGoal.setAllowUnusedReferences(True)

#Create states weight set
stateWeights = osim.MocoWeightSet()

#Loop through model states and add to weight set
#Set a generic weight for kinematic states (except pelvis_tx)
#Set zero weight for everything else
for stateInd in range(modelStates.size()):
    
    #Get state name
    stateName = modelStates.get(stateInd)
    
    #Check for pelvis tx
    if 'pelvis_tx' not in stateName:
        
        #Check for joint coordinate value
        if stateName.endswith('/value'):
            
            #Append generic weight
            stateWeights.cloneAndAppend(osim.MocoWeight(stateName, 1))
            
        else:
            
            #Append zero weight
            stateWeights.cloneAndAppend(osim.MocoWeight(stateName, 0))
    else:
        
        #Append zero weight
        stateWeights.cloneAndAppend(osim.MocoWeight(stateName, 0))

#Set weight set in goal
kinematicGoal.setWeightSet(stateWeights)

#Add to problem
problem.addGoal(kinematicGoal)

###################### GRF TRACKING SHOULD BE SCALED! #######################

#Add contact tracking goal
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
trackLeftGRF = osim.MocoContactTrackingGoalGroup(forcesLeftFoot, 'LeftGRF1')
trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l')        

#Set right foot tracking parameters for no GRFs
forcesRightFoot = osim.StdVectorString()
for contactLabel in rightFootContacts:
    forcesRightFoot.append(contactLabel)
trackRightGRF = osim.MocoContactTrackingGoalGroup(forcesRightFoot, 'RightGRF1')
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r')

#Set list to create tracking goals in to avoid overwriting
contactTracking = []

#Create contact settings dictionary
contactTrackingSettings = {'vector': [osim.Vec3(1,0,0),
                                      osim.Vec3(0,1,0),
                                      osim.Vec3(0,0,1)],
                           'weight': [1e-3, 1e-3, 1e-3]
                           }

#Loop through contact tracking settings to create separate weighted goals
for ii in range(len(contactTrackingSettings['vector'])):
    
    #Create tracking goal
    contactTracking.append(osim.MocoContactTrackingGoal(f'contact{ii}',
                                                        contactTrackingSettings['weight'][ii]))
    
    #Set external loads
    contactTracking[ii].setExternalLoadsFile('..\\\..\\expSims\\data\\sprint_grf.xml')
    
    #Add the left and right tracking groups to the goal
    contactTracking[ii].addContactGroup(trackLeftGRF)
    contactTracking[ii].addContactGroup(trackRightGRF)
    
    #Set projection vector in problem
    contactTracking[ii].setProjection('vector')
    contactTracking[ii].setProjectionVector(contactTrackingSettings['vector'][ii])

    #Add to problem
    problem.addGoal(contactTracking[ii])
    
###################### GRF TRACKING SHOULD BE SCALED! #######################

# #Set the time bounds in the problem
# #The initial time is locked to the same as the tracking simulation
# #No final time bounds are set as this would limit the stride style (i.e. if running
# #a faster speed, locking the time does not permit shorter strides)
# problem.setTimeBounds(osim.MocoTrajectory('..\\flashTrackingSim\\flashTrackingSimSolution.sto').getInitialTime(),
#                       [])

#The alternative is to bound the final time to the relative increase in speed and
#the final time in the original solution (i.e. it shouldn't get any slower)

#Calculate percentage increase in speed
relSpeed = baseSpeed / sprintSpeed

#Get initial and final times to calculate duration
initialTime = osim.MocoTrajectory('..\\flashTrackingSim\\flashTrackingSimSolution.sto').getInitialTime()
finalTime = osim.MocoTrajectory('..\\flashTrackingSim\\flashTrackingSimSolution.sto').getFinalTime()
timeDur = finalTime - initialTime

#Get relative speed in the context of the final time
relTime = timeDur * relSpeed

#Set time bounds
problem.setTimeBounds(initialTime, [initialTime + relTime, finalTime])

#Set kinematic bounds in problem
#Here we allow a 25% window of the total range either side of the max and min value
#from the tracking simulation. We also set the initial bounds to be within 10% of
#the initial values from the tracking simulation. This needs to ignore pelis_tx,
#or at least provide wider bounds when speed is increased.
trackingData = osim.MocoTrajectory('..\\flashTrackingSim\\flashTrackingSimSolution.sto')

#Loop through coordinates
for coordInd in range(osimModel.updCoordinateSet().getSize()):
    
    #Get coordinate name and absolute path
    coordName = osimModel.updCoordinateSet().get(coordInd).getName()
    coordPath = osimModel.updCoordinateSet().get(coordName).getAbsolutePathString()
    
    #Check for pelvis_tx
    if coordName != 'pelvis_tx':
    
        #Get minimum and maximum value from the tracking states and calculate range
        minVal = trackingData.getState(coordPath+'/value').to_numpy().min()
        maxVal = trackingData.getState(coordPath+'/value').to_numpy().max()
        rangeVal = np.diff((minVal,maxVal))[0]
        
        #Get the starting value from the tracked states
        initialVal = trackingData.getState(coordPath+'/value').to_numpy()[0]
        
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

#Define the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_optim_constraint_tolerance(0.01)
solver.set_optim_convergence_tolerance(0.01)
solver.set_multibody_dynamics_mode('explicit')
solver.set_num_mesh_intervals(5)
solver.resetProblem(problem)

#Set guess from full tracking simulation
solver.setGuessFile('..\\flashTrackingSim\\flashTrackingSimSolution.sto')
solver.resetProblem(problem)

########################## MESSING WITH GUESS #################################

# # guess = solver.createGuess()
# guess = solver.getGuess()
# ### TODO: created guess time is all over the place for some reason?
# # guess.setTime(np.linspace(trackingData.getInitialTime(),trackingData.getFinalTime(),11))
# # solver.setGuess(guess)


# for stateName in guess.getStateNames():
#     if stateName.startswith('/forceset/'):
#         guess.getState(stateName).setTo(0)
        
# for controlName in guess.getControlNames():
#     if controlName.startswith('/forceset/'):
#         guess.getControl(controlName).setTo(0)

# guess.write('guess.sto')
# solver.setGuessFile('guess.sto')
# solver.resetProblem(problem)

########################## MESSING WITH GUESS #################################

#Solve the tracking problem
solution = study.solve()
# study.visualize(solution)

"""
For some reason a coarser mesh interval works better and doesn't blow out...?

How did the initial above start go?
    > Not great - it's possible that unnaturally strong muscles can't achieve low forces
        >> Perhaps strength stays the same, but activation max is increased...
    > Bounds causing problems?
        >> Removed all bounds    
        >> Joint coordinate speeds have specific bounds placed on them too...
    > Seems like the change in activation dynamics is what's causing the issues
        >> Can we achieve Flash speeds with this remaining "normal"?
        >> Turning off activation dynamics doesn't seem to have same problem as drastically changing values
    

Is using initial guess on first Flash sim appropriate? It's probably not that 
worthwhile from the muscle function perspective?
    > Coarse mesh start instead on first sim?
    > With the same max isometric force it may not be as problematic?

"""

#Save solution to file
solution.write(f'..\\flashSimSpeed{str(sprintSpeed)}\\flashSimSpeed{str(sprintSpeed)}Solution.sto')

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
osim.STOFileAdapter().write(externalForcesTableFlat, f'..\\flashSimSpeed{str(sprintSpeed)}\\flashSimSpeed{str(sprintSpeed)}Solution_grf.sto')


# %% ----- End of runSimulations.py ----- %% #