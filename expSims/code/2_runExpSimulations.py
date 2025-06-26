# -*- coding: utf-8 -*-
"""

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    Script that runs the simulations associated with the experimental Dorn et al
    dataset. These simulations will serve as a foundation for predictive simulations.
    
    NOTES:
        
        > Struggling to track the full gait cycle without bounding over 2nd strike
            >> Is it use of contact sphere parameters from Fallise et al.?
                >>> Testing this premise by using old tracking model...
                    >>>> Didn't solve issue...
            >> Is it the marker tracking sim data that's the problem?
                >>> Include contact tracking in this problem formulation?
                    >>>> Testing in experimental data processing script...
            >> This may have always been a problem given only half gait cycle was tracked...?
            >> Could be mesh interval related?
                >> Raising mesh interval started to help with torque driven marker-trcking simulation
                    >>> Testing even higher with mesh interval of 50...
                        >>>> Didn't help that much...
            >> Could also somewhat just be the elevated weighting on contact tracking
                >>> i.e. optimisation trying to achieve too much accuracy and not converging
                    >>>> Still got the same bounding result...
            >> States weight too high?


"""

# %% Import packages

import opensim as osim
import osimfunctions as helper
import numpy as np
import os
import shutil

# %% Modify any details in here

"""

This section of code can be modified to suit another computer set-up or to set
which simulations to run.

"""

#Set opensim install path directory (e.g. C:\\OpenSim 4.4)
opensimPath = os.path.join('C:', os.sep, 'OpenSim 4.4')

#Set simulations to run
runFullTrackingSim = True

# %% Set-up

#Add OpenSim geometry path (weird issues with this on new laptop)
osim.ModelVisualizer.addDirToGeometrySearchPaths(os.path.join(opensimPath, 'Geometry'))

#Set path to experimental data tracking sim results
expSimSolution = os.path.join('..', 'results', 'markerTracking', 'markerTracking_fullGaitCycle_solution.sto')

# %% Full states and GRF tracking simulation

"""

This section runs a simulation that generates a muscle-driven full gait cycle
by tracking the joint coordinate states estimated from the marker tracking simulation
while also tracking experimental GRFs from contact spheres attached to the feet.

"""

#Check whether to run simulation
if runFullTrackingSim:

    #Create folder to store data
    os.makedirs(os.path.join('..', 'results', 'fullTracking'), exist_ok = True)
    
    #Set-up logger for simulation
    osim.Logger.removeFileSink()
    osim.Logger.addFileSink(os.path.join('..', 'results', 'fullTracking', 'fullTrackingDataLog.log'))
    
    ##### SETTINGS #####
    
    #Set the states global tracking weight
    globalStatesTrackingWeight = 1
    
    #Set constant weight to scale tracking error speeds by
    speedsTrackingScale = 0.001
    
    #Create a dictionary that provides the kinematic task weights for function
    taskWeights = {'pelvis_tx': 5, 'pelvis_ty': 2.5, 'pelvis_tz': 5,
                   'pelvis_tilt': 10, 'pelvis_list': 10, 'pelvis_rotation': 10, 
                   'hip_flexion_r': 10, 'hip_adduction_r': 5, 'hip_rotation_r': 2.5, 
                   'knee_angle_r': 10, 'ankle_angle_r': 2.5,
                   'subtalar_angle_r': 1.0,
                   'mtp_angle_r': 0.5,
                   'hip_flexion_l': 10, 'hip_adduction_l': 5, 'hip_rotation_l': 2.5, 
                   'knee_angle_l': 10, 'ankle_angle_l': 2.5,
                   'subtalar_angle_l': 1.0,
                   'mtp_angle_l': 0.5,
                   'lumbar_extension': 15, 'lumbar_bending': 10, 'lumbar_rotation': 5,
                   'arm_flex_r': 2.5, 'arm_add_r': 1, 'arm_rot_r': 1,
                   'elbow_flex_r': 5, 'pro_sup_r': 1,
                   'arm_flex_l': 2.5, 'arm_add_l': 1, 'arm_rot_l': 1,
                   'elbow_flex_l': 5, 'pro_sup_l': 1}
    
    #Create contact settings dictionary
    contactTrackingSettings = {'vector': [osim.Vec3(1,0,0),
                                          osim.Vec3(0,1,0),
                                          osim.Vec3(0,0,1)],
                               'weight': [0.1, 0.1, 0.1]
                               }
    
    #Set mesh interval
    numMeshInterval = 25
    
    ##### ADJUST MODEL FOR PROBLEM #####
    
    #Set the optimal forces dictionary for appending torque actuators to model
    #This appends torque actuators to the non muscle driven side and reserves to the
    #muscle driven side
    optForces = {
        #upper body
        'lumbar_extension': [300, 'actuator'], 'lumbar_bending': [300, 'actuator'], 'lumbar_rotation': [300, 'actuator'],
        'arm_flex_r': [300, 'actuator'], 'arm_add_r': [300, 'actuator'], 'arm_rot_r': [300, 'actuator'],
        'elbow_flex_r': [100, 'actuator'], 'pro_sup_r': [100, 'actuator'],
        'arm_flex_l': [300, 'actuator'], 'arm_add_l': [300, 'actuator'], 'arm_rot_l': [300, 'actuator'],
        'elbow_flex_l': [100, 'actuator'], 'pro_sup_l': [100, 'actuator'],
        #left limb
        'hip_flexion_l': [300, 'actuator'], 'hip_adduction_l': [200, 'actuator'], 'hip_rotation_l': [100, 'actuator'],
        'knee_angle_l': [300, 'actuator'], 'ankle_angle_l': [300, 'actuator'],
        'subtalar_angle_l': [100, 'actuator'], 'mtp_angle_l': [50, 'actuator'],
        #right limb
        'hip_flexion_r': [2, 'reserve'], 'hip_adduction_r': [2, 'reserve'], 'hip_rotation_r': [2, 'reserve'],
        'knee_angle_r': [2, 'reserve'], 'ankle_angle_r': [2, 'reserve'],  
        'subtalar_angle_r': [2, 'reserve'], 'mtp_angle_r': [2, 'reserve'],
        #pelvis
        'pelvis_tx': [1, 'residual'], 'pelvis_ty': [1, 'residual'], 'pelvis_tz': [1, 'residual'],
        'pelvis_tilt': [1, 'residual'], 'pelvis_list': [1, 'residual'], 'pelvis_rotation': [1, 'residual']
        }
        
    #Create the model to use in this simulation
    #This model has muscles on the right side and torque actuators driving the left side
    #Contact spheres are added to simulate ground contact in lieu of external loads
    #A smooth Bhargava metabolics model is added to muscles
    osimModel = helper.createSimModel(inputModelFile = os.path.join('..', 'model', 'JA1_SCALED_Osim4_Muscles.osim'),
                                      outputModelFile = os.path.join('..', 'results', 'fullTracking', 'osimModel.osim'),
                                      optForces = optForces,
                                      unilateralMuscles = True,
                                      jointsToWeld = ['radius_hand_r', 'radius_hand_l'],
                                      addMetabolicsModel = True)
    
    ##### CREATE TRACKING PROBLEM #####
    
    #Create and name an instance of the MocoTrack tool.
    track = osim.MocoTrack()
    track.setName('fullTrackingSim')
    
    #Set model in tracking tool which requires some further edits
    
    #Edit muscle contraction velocity
    #Increase the maximum contraction velocity of muscles
    maxContractionVelocity = 30
    for muscleInd in range(osimModel.getMuscles().getSize()):
        osimModel.getMuscles().get(muscleInd).set_max_contraction_velocity(maxContractionVelocity)
    
    #Set model in tracking tool
    modelProcessor = osim.ModelProcessor(osimModel)
    
    #Append additional edits to the model using operators
    modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
    modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
    modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(2.0))
    modelProcessor.append(osim.ModOpScaleMaxIsometricForce(2))
    
    #Set model in tracking tool
    track.setModel(modelProcessor)
    
    #Set the coordinates reference in tracking tool
    
    #Convert the solution data from the tracking sim to coordinates
    trackingStates = osim.MocoTrajectory(expSimSolution).exportToStatesTable()
    osim.STOFileAdapter().write(trackingStates, os.path.join('..', 'results', 'fullTracking', 'coordinatesToTrack.sto'))
    
    #Set coordinates file as reference in tool
    track.setStatesReference(osim.TableProcessor(os.path.join('..', 'results', 'fullTracking', 'coordinatesToTrack.sto')))
    track.set_states_global_tracking_weight(globalStatesTrackingWeight)
     
    #Set allow unused references in case there are extra markers
    track.set_allow_unused_references(True)
    
    #Set speeds a derivatives of coordinate references
    track.set_track_reference_position_derivatives(True)
    
    #Set tracked states to guess
    track.set_apply_tracked_states_to_guess(True)
    
    #Create weight set for state tracking
    stateWeights = osim.MocoWeightSet()
    
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
                                                        speedsTrackingScale * taskWeights[coordName])) 
    
    #Add to tracking problem
    track.set_states_weight_set(stateWeights)
    
    #Set start and end times in problem based on experimental tracking data
    startTime = osim.MocoTrajectory(expSimSolution).getInitialTime()
    endTime = osim.MocoTrajectory(expSimSolution).getFinalTime()
    
    #Set in tracking problem
    track.set_initial_time(startTime)
    track.set_final_time(endTime)
    
    #Convert the tracking tool to a Moco study and problem
    study = track.initialize()
    problem = study.updProblem()
    
    ##### GOALS #####
            
    #Update the weight and exponent on the default control effort goal
    effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
    effort.setWeight(0.001)
    
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
        if 'contactForce_' in osimModel.updForceSet().get(forceInd).getAbsolutePathString():
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
    
    #Loop through contact tracking settings to create separate weighted goals
    for ii in range(len(contactTrackingSettings['vector'])):
        
        #Create tracking goal
        contactTracking.append(osim.MocoContactTrackingGoal(f'contact{ii}',
                                                            contactTrackingSettings['weight'][ii]))
        
        #Set external loads
        contactTracking[ii].setExternalLoadsFile(os.path.join('..', 'data', 'sprint_grf.xml'))
        
        #Add the left and right tracking groups to the goal
        contactTracking[ii].addContactGroup(trackLeftGRF)
        contactTracking[ii].addContactGroup(trackRightGRF)
        
        #Set projection vector in problem
        contactTracking[ii].setProjection('vector')
        contactTracking[ii].setProjectionVector(contactTrackingSettings['vector'][ii])
    
        #Add to problem
        problem.addGoal(contactTracking[ii])
        
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

    ##### BOUNDS #####    

    #Set kinematic bounds in problem
    #Bounds are set to model coordinate limits here
    #TODO: set initial bounds to help...?
    for coordInd in range(osimModel.updCoordinateSet().getSize()):

        #Set bounds to model max and min values
        problem.setStateInfo(f'{osimModel.updCoordinateSet().get(coordInd).getAbsolutePathString()}/value',
                             #Bounds set to specified ranges
                             [osimModel.updCoordinateSet().get(coordInd).getRangeMin(),
                              osimModel.updCoordinateSet().get(coordInd).getRangeMax()]
                             )
        
    ##### SOLVER #####

    #Define and reset the solver
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    solver.resetProblem(problem)
    
    #Set solver options
    solver.set_optim_max_iterations(1000)
    solver.set_num_mesh_intervals(10)
    # solver.set_minimize_implicit_multibody_accelerations(True)
    # solver.set_scale_variables_using_bounds(True) #https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=1815&t=17154&p=0&start=0&view=&sid=d92da2432f6fc11cc0ce2beb83ad7672
    solver.set_optim_constraint_tolerance(1e-4)
    solver.set_optim_convergence_tolerance(1e-2)   

    #Solve the tracking problem
    solution = study.solve()
    
    
    
    # study.visualize(solution)
    
    #Save solution to file
    solution.write('..\\fullTrackingSim\\fullTrackingSimSolution_coarse.sto')
    
    #Extract predicted GRFs
    externalForcesTableFlat = osim.createExternalLoadsTableForGait(trackModelProcessor.process(), solution,
                                                                   forcesRightFoot, forcesLeftFoot)
    
    #Write table to file
    osim.STOFileAdapter().write(externalForcesTableFlat, '..\\fullTrackingSim\\fullTrackingSimSolution_coarse_grf.sto')
    
    #Reset the solver to progress to a finer mesh
    
    #Reset solver intervals
    solver.set_num_mesh_intervals(25)
    
    #Set initial guess to coarse solution
    solver.setGuessFile('..\\fullTrackingSim\\fullTrackingSimSolution_coarse.sto')
    solver.resetProblem(problem)
    
    #Re-solve on finer grid
    solutionFiner = study.solve()
    # study.visualize(solutionFiner)
    
    #Copy tracked markers file over to main directory
    shutil.move('fullTrackingSim_tracked_states.sto',
                '..\\fullTrackingSim\\fullTrackingSim_tracked_states.sto')
    
    #Save solution to file
    solutionFiner.write('..\\fullTrackingSim\\fullTrackingSimSolution.sto')
    
    #Extract predicted GRFs
    externalForcesTableFlat = osim.createExternalLoadsTableForGait(trackModelProcessor.process(), solutionFiner,
                                                                   forcesRightFoot, forcesLeftFoot)
    
    #Write table to file
    osim.STOFileAdapter().write(externalForcesTableFlat, '..\\fullTrackingSim\\fullTrackingSimSolution_grf.sto')

# %% Predictive simulation

"""

This section runs a muscle-driven predictive simulation of the half gait cycle
that aims to replicate the tracking simulation results as a means to validate the
predictive approach used in subsequent simulations.

"""

#Create folder to store data
try:
    os.mkdir('..\\predictiveSim')
except:
    print('Predictive sim folder already detected...')
    
#Create the model to use in this simulation
osimModel = helper.createSimModel(inputModelFile = '..\\data\\JA1_SCALED_Osim40_Muscles.osim',
                                  outputModelFile = '..\\predictiveSim\\osimModel.osim',
                                  unilateralMuscles = False,
                                  jointsToWeld = ['subtalar_r', 'subtalar_l', 'radius_hand_r', 'radius_hand_l'])

#Define the predictive study problem
study = osim.MocoStudy()
study.setName('predictiveSim')

#Define the problem
problem = study.updProblem()

#Set model in problem
problem.setModelProcessor(osim.ModelProcessor(osimModel))

#Create a speed goal to match the experimental data sprint speed
#Set an average speed goal based on sprinting data
#Get the average sprint speed based on the pelvis translation from experimental data

#Get the tracking sim data
trackingData = osim.MocoTrajectory('..\\fullTrackingSim\\fullTrackingSimSolution.sto')

#Get distance travelled based on pelvis translation alongside time taken
distTravelled = trackingData.getState('/jointset/ground_pelvis/pelvis_tx/value').to_numpy()[-1] - trackingData.getState('/jointset/ground_pelvis/pelvis_tx/value').to_numpy()[0]
timeTaken = trackingData.getTime().to_numpy()[-1] - trackingData.getTime().to_numpy()[0]

#Calculate average speed
sprintSpeed = distTravelled / timeTaken

#Create the speed goal and set parameters
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

#Loop through states
for stateInd in range(modelStates.getSize()):
    
    #Get state name
    stateName = modelStates.get(stateInd)
            
    #Opposite muscle activations are periodic
    
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

#Create a goal to minimise muscle and actuator controls and set parameters
effortGoal = osim.MocoControlGoal('control_effort', 0.001) #set the weight when creating
effortGoal.setExponent(3)

#Update individual weights in control effort goal to be relative to
#actual muscle and reserve actuator names
#Set appropriate patterns in the weight set
#Muscles
effortGoal.setWeightForControlPattern('/forceset/.*/activation', 0.1)
#Idealised actuators
effortGoal.setWeightForControlPattern('/forceset/.*_actuator', 0.1)
#Reserves
effortGoal.setWeightForControlPattern('/forceset/.*_reserve', 10)
#Residuals
effortGoal.setWeightForControlPattern('/forceset/.*_residual', 10)

#Add to problem
problem.addGoal(effortGoal)

# #Add metabolics goal to guide muscle activation
# metGoal = osim.MocoOutputGoal('metCost', 0.1)
# metGoal.setOutputPath('/metabolicModel|total_metabolic_rate')
# metGoal.setDivideByDisplacement(True)
# metGoal.setDivideByMass(True)

# #Add to problem
# problem.addGoal(metGoal)

# #Add a low weighted states tracking for joint coordinates
# #Exclude speeds and pelvis_tx as these will change with altered speed
# kinematicGoal = osim.MocoStateTrackingGoal('kinematicTracking', 0.1)

# #Set kinematic data to track
# tableProcessor = osim.TableProcessor('..\\fullTrackingSim\\fullTrackingSimSolution.sto')
# tableProcessor.append(osim.TabOpUseAbsoluteStateNames())
# kinematicGoal.setReference(tableProcessor)
# kinematicGoal.setAllowUnusedReferences(True)

# #Create states weight set
# stateWeights = osim.MocoWeightSet()

# #Loop through model states and add to weight set
# #Set a generic weight for kinematic states (except pelvis_tx)
# #Set zero weight for everything else
# for stateInd in range(modelStates.size()):
    
#     #Get state name
#     stateName = modelStates.get(stateInd)
    
#     #Check for pelvis tx
#     if 'pelvis_tx' not in stateName:
        
#         #Check for joint coordinate value
#         if stateName.endswith('/value'):
            
#             #Append generic weight
#             stateWeights.cloneAndAppend(osim.MocoWeight(stateName, 1))
            
#         else:
            
#             #Append zero weight
#             stateWeights.cloneAndAppend(osim.MocoWeight(stateName, 0))
#     else:
        
#         #Append zero weight
#         stateWeights.cloneAndAppend(osim.MocoWeight(stateName, 0))

# #Set weight set in goal
# kinematicGoal.setWeightSet(stateWeights)

# #Add to problem
# problem.addGoal(kinematicGoal)

#Set the time bounds in the problem
#The initial time is locked to the same as the tracking simulation
#No final time bounds are set as this would limit the stride style (i.e. if running
#a faster speed, locking the time does not permit shorter strides)
problem.setTimeBounds(trackingData.getTime().to_numpy()[0], [])

#Set kinematic bounds in problem
#Here we allow a 25% window of the total range either side of the max and min value
#from the tracking simulation. We also set the initial bounds to be within 10% of
#the initial values from the tracking simulation. This needs to ignore pelis_tx,
#or at least provide wider bounds when speed is increased.

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
solver.set_num_mesh_intervals(25)
solver.resetProblem(problem)

#Set guess from full tracking simulation
solver.setGuessFile('..\\fullTrackingSim\\fullTrackingSimSolution.sto')
solver.resetProblem(problem)

#Solve the tracking problem
solution = study.solve()
# study.visualize(solution)

#Save solution to file
solution.write('..\\predictiveSim\\predictiveSimSolution.sto')

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
osim.STOFileAdapter().write(externalForcesTableFlat, '..\\predictiveSim\\predictiveSimSolution_grf.sto')

# %% ----- End of 2_runExpSimulations.py ----- %% #