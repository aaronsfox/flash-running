# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 15:23:54 2021

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This script goes through the process of running two-dimensional tracking sims
    of the experimental data. First, we run a torque-driven tracking simulation
    of experimental kinematics and GRF data. This produces kinematics that are
    consistent with the experimental data, but specifically with the GRFs being
    produced by the contact spheres on the model feet. Second, we run an identical
    muscle driven simulation - but with the kinematics from the torque driven tracking
    sim as the input. This produces a muscle-driven simulation of the running movement
    with GRFs tracked by the contact spheres. Third, we run a predictive simulation
    targeting the same running speed to examine the consistency of our predictive
    approach with the tracking sim.
    
"""

# %% Import packages

#import os
import opensim as osim
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# %% Define functions

#readSTO
def readSTO(fileName):
    
    #Read storage
    sto = osim.Storage(fileName)
        
    #Get column list
    labels = osimArrayToList(sto.getColumnLabels())

    #Get time
    time = osim.ArrayDouble()
    sto.getTimeColumn(time)  
    time = osimArrayToList(time)
    
    #Get data
    data = []
    for i in range(sto.getSize()):
        temp = osimArrayToList(sto.getStateVector(i).getData())
        temp.insert(0, time[i])
        data.append(temp)

    df = pd.DataFrame(data, columns=labels)
    df.index = df.time
    
    return df

#osimArrayToList
def osimArrayToList(array):
    temp = []
    for i in range(array.getSize()):
        temp.append(array.get(i))

    return temp

#addReserve
def addReserve(model, coord, optimal_force, max_control):
    actu = osim.ActivationCoordinateActuator()
    actu.set_coordinate(coord)
    if coord.startswith('pelvis'):
        prefix = 'residual_'
    elif (coord.startswith('lumbar') or
          coord.startswith('arm') or
          coord.startswith('elbow') or
          coord.startswith('pro_sup')):
        prefix = 'torque_'
    else:
        prefix = 'reserve_'
    actu.setName(prefix + coord)
    actu.setOptimalForce(optimal_force)
    actu.setMinControl(-max_control)
    actu.setMaxControl(max_control)
    model.addForce(actu)

# %% Set-up

#Some basic parameters to determine how the function runs.

#Specifically here we provide options to run or load the simulations
#Defaults here are False to not run as results can be imported
runTorqueTracking = False 
runMuscleTracking = False 
runMusclePrediction = False 

#Here we provide an option to compare and visualise the results of the function
#Defaults to true for loading results and comparing visually
compareResults = True

#Set matplotlib parameters
from matplotlib import rcParams
# rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Arial'
rcParams['font.weight'] = 'bold'
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 16
rcParams['axes.linewidth'] = 1.5
rcParams['axes.labelweight'] = 'bold'
rcParams['legend.fontsize'] = 10
rcParams['xtick.major.width'] = 1.5
rcParams['ytick.major.width'] = 1.5
rcParams['legend.framealpha'] = 0.0
rcParams['savefig.dpi'] = 300
rcParams['savefig.format'] = 'pdf'

# %% 2D basic tracking sim of experimental data (torque driven)

# This step is designed to generate relevant controls for the 2D model that
# appropriately track the sprint kinematics and GRF from the experimental data
# using torque actuators

#Check whether to run simulation
if runTorqueTracking:

    #Set simulation parameters
    
    #Timing and mesh data
    initialTime = osim.Storage('refQ_2D.sto').getFirstTime()
    finalTime = osim.Storage('refQ_2D.sto').getLastTime()
    duration = finalTime - initialTime
    meshNo = 50
    meshInterval = duration / meshNo
    
    #Weights
    trackingWeight = 10
    effortWeight = 0.1
    grfWeight = 1
    speedWeight = 1
    symmetryWeight = 1
    
    #Model parameters
    complexModel = False
    
    #General parameters
    visualiseSolution = False
    
    #Create the model processor for the tracking problem
    
    #Load the model and create processor
    if complexModel:
        osimModel = osim.Model('gaitModel2D_complex.osim')
    else:
        osimModel = osim.Model('gaitModel2D.osim')
    modelProcessor = osim.ModelProcessor(osimModel)
    
    #Add processing step to remove muscles
    modelProcessor.append(osim.ModOpRemoveMuscles())
    
    #Process model for additional parameters
    torqueModel = modelProcessor.process()
    torqueModel.initSystem()
    
    #Add torque actuators to support model
    #Set list to add reserves to
    if complexModel:
        raise ValueError('Complex model not supported here. Need to add model coordinates...')
    else:
        reserveList = ['pelvis_tilt',
                        'hip_flexion_l', 'hip_flexion_r',
                        'knee_angle_l','knee_angle_r',
                        'ankle_angle_l', 'ankle_angle_r']
        
    #Set optimal force and max torque
    optimalForce = 100
    maxTorque = np.inf
    #Add torque actuators
    for tt in range(len(reserveList)):
        addReserve(torqueModel, reserveList[tt], optimalForce, maxTorque)
    
    #Finalise connections
    torqueModel.finalizeConnections()
    
    #Get the track model as a processor object
    torqueModelProcessor = osim.ModelProcessor(torqueModel)
    
    #Set-up tracking problem
    
    #Process current model
    trackModel = torqueModelProcessor.process()
    trackModel.initSystem()
          
    #Construct the tracking object and set basic parameters
    track = osim.MocoTrack()
    track.setName('sprintTracking_torqueDriven')
    track.setModel(torqueModelProcessor)
    
    #Set kinematic data and parameters
    tableProcessor = osim.TableProcessor('refQ_2D.sto')
    # tableProcessor.append(osim.TabOpLowPassFilter(12)) ### data already filtered via RRA
    tableProcessor.append(osim.TabOpUseAbsoluteStateNames())
    track.setStatesReference(tableProcessor)
    track.set_states_global_tracking_weight(trackingWeight)
    #Set some coordinates not to be tracked to avoid poor motion
    #Similar process done in Dembia et al. Moco paper
    stateWeights = osim.MocoWeightSet()
    weightList = list()
    weightList.append(('/jointset/ground_pelvis/pelvis_ty', 0))
    weightList.append(('/jointset/ground_pelvis/pelvis_tilt', 0.1))
    for weight in weightList:
        stateWeights.cloneAndAppend(osim.MocoWeight(weight[0] + '/value', weight[1]))
        stateWeights.cloneAndAppend(osim.MocoWeight(weight[0] + '/speed', weight[1]))        
    track.set_states_weight_set(stateWeights)
    #Set tracked states to be used in guess
    track.set_apply_tracked_states_to_guess(True)
    #Set unused state references to be allowed in case of file errors
    track.set_allow_unused_references(True)
    #Set position derivatives to be used as speeds
    track.set_track_reference_position_derivatives(True)
    
    #Set control weights
    track.set_control_effort_weight(effortWeight)
    
    #Set timing parameters
    track.set_initial_time(initialTime)
    track.set_final_time(finalTime)
    track.set_mesh_interval(meshInterval)
    
    #Customise the base tracking problem with relevant goals
    study = track.initialize()
    problem = study.updProblem()
    
    #Adjust time bounds to allow for subtle fluctations in finish time
    problem.setTimeBounds(initialTime, [finalTime - 0.1, finalTime + 0.1])
    
    #Update the control effort goal to a cost of transport type cost
    effort = osim.MocoControlGoal().safeDownCast(problem.updGoal('control_effort'))
    effort.setDivideByDisplacement(True)
    
    #Set to not highly penalise the lumbar actuator
    if complexModel:
        effort.setWeightForControl('/forceset/tau_lumbar_ext', 0.001)
    else:
        effort.setWeightForControl('/forceset/lumbarAct', 0.001)
    
    #Set an average speed goal based on sprinting data
    df_kinematics = readSTO('refQ_2D.sto')
    startPos = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'].to_numpy()[0]
    endPos = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'].to_numpy()[-1]
    sprintSpeed = (endPos - startPos) / (finalTime - initialTime)
    speedGoal = osim.MocoAverageSpeedGoal('speed')
    speedGoal.set_desired_average_speed(sprintSpeed)
    speedGoal.setWeight(speedWeight)
    problem.addGoal(speedGoal)
    
    #Set symmetry constraint
    #Create symmetry goal
    symmetryGoal = osim.MocoPeriodicityGoal('symmetryGoal')
    #Set symmetric coordinate values (except for pelvis_tx) and speeds
    for cc in range(trackModel.getCoordinateSet().getSize()):
        #Get current coordinate name
        coordName = trackModel.getCoordinateSet().get(cc).getName()
        #Set symmetry parameters based on coordinate
        if coordName.endswith('_r'):
            #Joint angle value
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/value'))
            #Joint speed
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/speed'))
        elif coordName.endswith('_l'):
            #Joint angle value
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/value'))
            #Joint speed
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/speed'))
        elif coordName.endswith('_extension') or coordName.endswith('_tilt') or coordName.endswith('_ty'):
            #Joint angle value
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value'))
            #Joint speed
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed'))
    #Add a symmetry pair for pelvis_tx speed
    symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
        '/jointset/ground_pelvis/pelvis_tx/speed'))
    #Add symmetry goal
    symmetryGoal.setWeight(symmetryWeight)
    problem.addGoal(symmetryGoal)
    
    #Add contact tracking goal
    #Set force names
    if complexModel:
        forceNamesRightFoot = ['forceset/contactHeel_r',
                               'forceset/contactMH1_r',
                               'forceset/contactMH3_r',
                               'forceset/contactMH5_r',
                               'forceset/contactHallux_r',
                               'forceset/contactOtherToes_r']
        forceNamesLeftFoot = ['forceset/contactHeel_l',
                                'forceset/contactMH1_l',
                                'forceset/contactMH3_l',
                                'forceset/contactMH5_l',
                                'forceset/contactHallux_l',
                                'forceset/contactOtherToes_l']
    else:
        forceNamesRightFoot = ['forceset/contactHeel_r',
                               'forceset/contactMidfoot_r',
                               'forceset/contactToe_r']
        forceNamesLeftFoot = ['forceset/contactHeel_l',
                                'forceset/contactMidfoot_l',
                                'forceset/contactToe_l']
    
    #Create contact tracking goal
    contactGoal = osim.MocoContactTrackingGoal('contactGoal', grfWeight)
    #Set external loads
    contactGoal.setExternalLoadsFile('refGRF_2D.xml')
    #Set force name groups
    forceNames_r = osim.StdVectorString()
    forceNames_l = osim.StdVectorString()
    for ff in range(len(forceNamesRightFoot)):
        forceNames_r.append(forceNamesRightFoot[ff])
        forceNames_l.append(forceNamesLeftFoot[ff])
    #Create and add tracking groups
    trackRightGRF = osim.MocoContactTrackingGoalGroup(forceNames_r, 'RightGRF')
    trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r');
    contactGoal.addContactGroup(trackRightGRF)
    trackLeftGRF = osim.MocoContactTrackingGoalGroup(forceNames_l, 'LeftGRF')
    trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l')
    contactGoal.addContactGroup(trackLeftGRF)
    #Set parameters
    contactGoal.setProjection('plane')
    contactGoal.setProjectionVector(osim.Vec3(0, 0, 1))
    #Add contact tracking goal
    problem.addGoal(contactGoal)
    
    #Add state bounds
    #This also sets initial bounds on lumbar and pelvis coordinates that reflect
    #'normal' running motions - based off the experimental data
    problem.setStateInfo('/jointset/back/lumbar_extension/value', [np.deg2rad(-30),0],
                         [np.deg2rad(-10),0])
    problem.setStateInfo('/jointset/ground_pelvis/pelvis_tilt/value', [-10*np.pi/180, 0*np.pi/180],
                         [np.deg2rad(-5),0])
    problem.setStateInfo('/jointset/ground_pelvis/pelvis_tx/value', [0, 3])
    problem.setStateInfo('/jointset/ground_pelvis/pelvis_ty/value', [0.8, 1.25])
    problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [-25*np.pi/180, 75*np.pi/180])
    problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [-25*np.pi/180, 75*np.pi/180])
    if complexModel:
        problem.setStateInfo('/jointset/walker_knee_l/knee_angle_l/value', [0, 140*np.pi/180])
        problem.setStateInfo('/jointset/walker_knee_r/knee_angle_r/value', [0, 140*np.pi/180])
    else:
        problem.setStateInfo('/jointset/knee_l/knee_angle_l/value', [-140*np.pi/180, 0])
        problem.setStateInfo('/jointset/knee_r/knee_angle_r/value', [-140*np.pi/180, 0])
    problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [-20*np.pi/180, 30*np.pi/180])
    problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [-20*np.pi/180, 30*np.pi/180])
    
    #Configure the solver
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    solver.resetProblem(problem)
    solver.set_optim_constraint_tolerance(1e-2)
    solver.set_optim_convergence_tolerance(1e-2)
    
    #Solve
    solution = study.solve()
        
    #Option to visualise
    if visualiseSolution:
        study.visualize(solution)
    
    #Create a full gait cycle trajectory from the periodic solution.
    addPatterns = [".*pelvis_tx/value"]
    fullTraj = osim.createPeriodicTrajectory(solution, addPatterns)
    fullTraj.write('sprintTracking_torqueDriven_solution_fullTrajectory.sto')
    
    #Compute ground reaction forces generated by contact sphere from the 
    #full gait cycle trajectory
    externalLoads = osim.createExternalLoadsTableForGait(trackModel,
                                                         fullTraj,
                                                         forceNames_r,
                                                         forceNames_l)
    osim.STOFileAdapter.write(externalLoads,'trackedGRF_torqueDriven_2D.mot')

# %% 2D tracking sim of experimental data (muscle driven)

# This step is designed to generate relevant controls for the 2D model that
# appropriately track the sprint kinematics and GRF from the experimental data

#Check whether to run simulation
if runMuscleTracking:

    #Set simulation parameters
    
    #Timing and mesh data
    initialTime = osim.Storage('sprintTracking_torqueDriven_solution.sto').getFirstTime()
    finalTime = osim.Storage('sprintTracking_torqueDriven_solution.sto').getLastTime()
    duration = finalTime - initialTime
    meshNo = 100 #higher mesh here for muscle driven
    meshInterval = duration / meshNo
    
    #Model parameters
    complexModel = False
    
    #General parameters
    visualiseSolution = False
    
    #Muscle parameters
    passiveForces = False
    implicitTendonCompliance = True
    tendonDynamics = True
    contractVelScale = 2.5 #upscale contractile velocity for sprinting
    maxForceScale = 2 #upscale force for sprinting
    
    #Weights
    trackingWeight = 10
    effortWeight = 0.1
    grfWeight = 1
    speedWeight = 1
    symmetryWeight = 1
    
    #Create the model processor for the tracking problem
    
    #Load the model and create processor
    if complexModel:
        osimModel = osim.Model('gaitModel2D_complex.osim')
    else:
        osimModel = osim.Model('gaitModel2D.osim')  
    
    #Update contractile velocity maximum to assist with sprinting sim
    for mm in range(osimModel.getMuscles().getSize()):
        #Get the current muscle
        currMusc = osimModel.getMuscles().get(mm)
        #Get current max contraction velocity
        currVel = currMusc.getMaxContractionVelocity()
        #Set new max contraction velocity
        currMusc.set_max_contraction_velocity(currVel*contractVelScale)
        
    #Finalise connections
    osimModel.finalizeConnections()
    
    #Put the model in a processor
    modelProcessor = osim.ModelProcessor(osimModel)
    
    #Scale the maximum isometric force of muscles to assist with sprinting sim
    modelProcessor.append(osim.ModOpScaleMaxIsometricForce(maxForceScale))
    
    #Disable tendon compliance (will be re-enabled for certain muscles later)
    modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
    
    #Adjust fiber damping as per Dembia et al. Moco paper
    modelProcessor.append(osim.ModOpFiberDampingDGF(0.01))
    
    #Turn off passive forces if appropriate
    if not passiveForces:
        modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
        
    #Modify active force width
    #This appears most important for "de-stiffening" muscles to allow muscle
    #to function well over sprinting joint ranges
    modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
        
    #Process model for additional parameters
    trackModel = modelProcessor.process()
    trackModel.initSystem()
    
    #Enable tendon compliance in the gastrocnemius and soleus
    if tendonDynamics:
        muscles = trackModel.updMuscles()
        for mm in np.arange(muscles.getSize()):
            currMusc = osim.DeGrooteFregly2016Muscle.safeDownCast(muscles.get(int(mm)))
            muscName = currMusc.getName()
            #Enable tendon compliance dynamics in the plantarflexors
            if ('gastroc' in muscName) or ('soleus' in muscName):
                currMusc.set_ignore_tendon_compliance(False)
            
    #Get the track model as a processor object
    trackModelProcessor = osim.ModelProcessor(trackModel)
    
    #Set model to use implicit tendon compliance if appropriate
    if implicitTendonCompliance:
        trackModelProcessor.append(osim.ModOpUseImplicitTendonComplianceDynamicsDGF())
    
    #Set-up tracking problem
    
    #Process current model
    trackModel = trackModelProcessor.process()
    trackModel.initSystem()
    
    # #Count the number of Force objects in the model. This can be used to
    # #normalise control effort
    # numForces = 0
    # for ff in range(trackModel.getForceSet().getSize()):
    #     if (trackModel.getForceSet().get(ff).getConcreteClassName().endswith('Muscle') or
    #         trackModel.getForceSet().get(ff).getConcreteClassName().endswith('Actuator')):
    #         numForces += 1
            
    #Construct the tracking object and set basic parameters
    track = osim.MocoTrack()
    track.setName('sprintTracking_muscleDriven')
    track.setModel(trackModelProcessor)
    
    #Set kinematic data and parameters
    tableProcessor = osim.TableProcessor('sprintTracking_torqueDriven_solution.sto')
    tableProcessor.append(osim.TabOpUseAbsoluteStateNames())
    track.setStatesReference(tableProcessor)
    track.set_states_global_tracking_weight(trackingWeight)
    #Set tracked states to be used in guess
    track.set_apply_tracked_states_to_guess(True)
    #Set unused state references to be allowed in case of file errors
    track.set_allow_unused_references(True)
    #Set position derivatives to be used as speeds
    track.set_track_reference_position_derivatives(True)
    
    #Set control weights
    track.set_control_effort_weight(effortWeight)
    # track.set_control_effort_weight(effortWeight / numForces)
    
    #Set timing parameters
    track.set_initial_time(initialTime)
    track.set_final_time(finalTime)
    track.set_mesh_interval(meshInterval)
    
    #Customise the base tracking problem with relevant goals
    study = track.initialize()
    problem = study.updProblem()
    
    #Adjust time bounds to allow for subtle fluctations in finish time
    problem.setTimeBounds(initialTime, [finalTime - 0.1, finalTime + 0.1])
    
    #Update the control effort goal to a cost of transport type cost
    effort = osim.MocoControlGoal().safeDownCast(problem.updGoal('control_effort'))
    effort.setDivideByDisplacement(True)
    
    #Set to not highly penalise the lumbar actuator
    if complexModel:
        effort.setWeightForControl('/forceset/tau_lumbar_ext', 0.001)
    else:
        effort.setWeightForControl('/forceset/lumbarAct', 0.001)
    
    #Set an average speed goal based on sprinting data
    df_kinematics = readSTO('refQ_2D.sto')
    startPos = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'].to_numpy()[0]
    endPos = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'].to_numpy()[-1]
    sprintSpeed = (endPos - startPos) / (finalTime - initialTime)
    speedGoal = osim.MocoAverageSpeedGoal('speed')
    speedGoal.set_desired_average_speed(sprintSpeed)
    speedGoal.setWeight(speedWeight)
    problem.addGoal(speedGoal)
    
    #Set symmetry constraint
    #Create symmetry goal
    symmetryGoal = osim.MocoPeriodicityGoal('symmetryGoal')
    #Set symmetric coordinate values (except for pelvis_tx) and speeds
    for cc in range(trackModel.getCoordinateSet().getSize()):
        #Get current coordinate name
        coordName = trackModel.getCoordinateSet().get(cc).getName()
        #Set symmetry parameters based on coordinate
        if coordName.endswith('_r'):
            #Joint angle value
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/value'))
            #Joint speed
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/speed'))
        elif coordName.endswith('_l'):
            #Joint angle value
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/value'))
            #Joint speed
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/speed'))
        elif coordName.endswith('_extension') or coordName.endswith('_tilt') or coordName.endswith('_ty'):
            #Joint angle value
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value'))
            #Joint speed
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed'))
    #Add a symmetry pair for pelvis_tx speed
    symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
        '/jointset/ground_pelvis/pelvis_tx/speed'))
    #Set symmetric activations
    for ff in range(trackModel.getForceSet().getSize()):
        #Check for muscle or actuator
        if (trackModel.getForceSet().get(ff).getConcreteClassName().endswith('Muscle') or
            trackModel.getForceSet().get(ff).getConcreteClassName().endswith('Actuator')):
            #Get force name
            forceName = trackModel.getForceSet().get(ff).getName()
            #Set symmetry parameters based on coordinate
            if forceName.endswith('_r'):
                symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                    trackModel.getForceSet().get(ff).getAbsolutePathString()+'/activation',
                    trackModel.getForceSet().get(ff).getAbsolutePathString().replace('_r','_l')+'/activation'))
            elif forceName.endswith('_l'):
                symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                    trackModel.getForceSet().get(ff).getAbsolutePathString()+'/activation',
                    trackModel.getForceSet().get(ff).getAbsolutePathString().replace('_l','_r')+'/activation'))
            # elif forceName == 'lumbarAct':
            #     symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            #         trackModel.getForceSet().get(ff).getAbsolutePathString()))
    #Add symmetry goal
    symmetryGoal.setWeight(symmetryWeight)
    problem.addGoal(symmetryGoal)
    
    #Add contact tracking goal
    #Set force names
    if complexModel:
        forceNamesRightFoot = ['forceset/contactHeel_r',
                               'forceset/contactMH1_r',
                               'forceset/contactMH3_r',
                               'forceset/contactMH5_r',
                               'forceset/contactHallux_r',
                               'forceset/contactOtherToes_r']
        forceNamesLeftFoot = ['forceset/contactHeel_l',
                                'forceset/contactMH1_l',
                                'forceset/contactMH3_l',
                                'forceset/contactMH5_l',
                                'forceset/contactHallux_l',
                                'forceset/contactOtherToes_l']
    else:
        forceNamesRightFoot = ['forceset/contactHeel_r',
                               'forceset/contactMidfoot_r',
                               'forceset/contactToe_r']
        forceNamesLeftFoot = ['forceset/contactHeel_l',
                                'forceset/contactMidfoot_l',
                                'forceset/contactToe_l']
    
    #Create contact tracking goal
    contactGoal = osim.MocoContactTrackingGoal('contactGoal', grfWeight)
    #Set external loads
    contactGoal.setExternalLoadsFile('refGRF_2D.xml')
    #Set force name groups
    forceNames_r = osim.StdVectorString()
    forceNames_l = osim.StdVectorString()
    for ff in range(len(forceNamesRightFoot)):
        forceNames_r.append(forceNamesRightFoot[ff])
        forceNames_l.append(forceNamesLeftFoot[ff])
    #Create and add tracking groups
    trackRightGRF = osim.MocoContactTrackingGoalGroup(forceNames_r, 'RightGRF')
    trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r');
    contactGoal.addContactGroup(trackRightGRF)
    trackLeftGRF = osim.MocoContactTrackingGoalGroup(forceNames_l, 'LeftGRF')
    trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l')
    contactGoal.addContactGroup(trackLeftGRF)
    #Set parameters
    contactGoal.setProjection('plane')
    contactGoal.setProjectionVector(osim.Vec3(0, 0, 1))
    #Add contact tracking goal
    problem.addGoal(contactGoal)
    
    #Add state bounds
    #Set the initial values for the lumbar and pelvis coordinates that
    #reflect 'normal' running motions - bounds based off experimental data
    problem.setStateInfo('/jointset/back/lumbar_extension/value', [np.deg2rad(-30),0],
                         [np.deg2rad(-10),0])
    problem.setStateInfo('/jointset/ground_pelvis/pelvis_tilt/value', [-10*np.pi/180, 0*np.pi/180],
                         [np.deg2rad(-5),0])
    problem.setStateInfo('/jointset/ground_pelvis/pelvis_tx/value', [0, 3],
                         1.06261516975403) #tracking sim starting value
    problem.setStateInfo('/jointset/ground_pelvis/pelvis_ty/value', [0.8, 1.25],
                         [1.00, 1.02]) #likely reasonable starting value
    problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [-25*np.pi/180, 75*np.pi/180])
    problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [-25*np.pi/180, 75*np.pi/180])
    if complexModel:
        problem.setStateInfo('/jointset/walker_knee_l/knee_angle_l/value', [0, 140*np.pi/180])
        problem.setStateInfo('/jointset/walker_knee_r/knee_angle_r/value', [0, 140*np.pi/180])
    else:
        problem.setStateInfo('/jointset/knee_l/knee_angle_l/value', [-140*np.pi/180, 0])
        problem.setStateInfo('/jointset/knee_r/knee_angle_r/value', [-140*np.pi/180, 0])
    problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [-20*np.pi/180, 30*np.pi/180])
    problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [-20*np.pi/180, 30*np.pi/180])
    
    #Set reasonable bounds for tendon force
    problem.setStateInfoPattern('/forceset/.*/normalized_tendon_force', [0, 1.8], [], [])
    
    #Configure the solver
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    solver.set_optim_constraint_tolerance(1e-2)
    solver.set_optim_convergence_tolerance(1e-2)
    # solver.set_multibody_dynamics_mode('implicit')
    # solver.set_minimize_implicit_multibody_accelerations(True)
    # solver.set_implicit_multibody_accelerations_weight(1e-3 / trackModel.getNumCoordinates())
    # solver.set_implicit_multibody_acceleration_bounds(osim.MocoBounds(-500, 500))
    # solver.set_implicit_multibody_acceleration_bounds(osim.MocoBounds(-200, 200))
    # solver.set_minimize_implicit_auxiliary_derivatives(True)
    # solver.set_implicit_auxiliary_derivatives_weight(0.00001)
    solver.resetProblem(problem)

    #Set the normalized tendon forces if not loading initial guess from file
    if tendonDynamics:
        
        #Get the current default-ish guess
        guess = solver.getGuess()
        numRows = guess.getNumTimes()
        
        #Set the tendon force to reasonable values in initial guess to avoid a bad initial guess
        stateNames = trackModel.getStateVariableNames()
        for ii in range(trackModel.getNumStateVariables()):
            currState = stateNames.get(ii)
            if 'normalized_tendon_force' in currState:
                guess.setState(currState, np.linspace(0.2,0.2,numRows))
        
        #Set guess
        solver.setGuess(guess)
    
    #Solve
    solution = study.solve()
    
    #Option to visualise
    if visualiseSolution:
        study.visualize(solution)
    
    #Create a full gait cycle trajectory from the periodic solution.
    addPatterns = [".*pelvis_tx/value"]
    fullTraj = osim.createPeriodicTrajectory(solution, addPatterns)
    fullTraj.write('sprintTracking_muscleDriven_solution_fullTrajectory.sto')
    
    #Compute ground reaction forces generated by contact sphere from the 
    #full gait cycle trajectory
    externalLoads = osim.createExternalLoadsTableForGait(trackModel,
                                                         fullTraj,
                                                         forceNames_r,
                                                         forceNames_l)
    osim.STOFileAdapter.write(externalLoads,'trackedGRF_muscleDriven_2D.mot')

# %% 2D predictive at same speed (muscle driven)

# This section attempts to generate a predictive simulation of a half gait cycle
# at the same speed as tracking. Theoretically this should produce some similarities
# to the tracking sim as a pseudo-validation of the predictive approach

#Check whether to run simulation
if runMusclePrediction:

    #Set simulation parameters
    
    #Timing and mesh data
    initialTime = osim.Storage('sprintTracking_torqueDriven_solution.sto').getFirstTime()
    finalTime = osim.Storage('sprintTracking_torqueDriven_solution.sto').getLastTime()
    duration = finalTime - initialTime
    meshNo = 100 #higher mesh no for muscle driven
    meshInterval = duration / meshNo
    
    #Weights
    effortWeight = 0.1
    speedWeight = 1
    symmetryWeight = 1
    
    #Model parameters
    complexModel = False
    
    #General parameters
    visualiseSolution = False
    
    #Muscle parameters
    passiveForces = False
    implicitTendonCompliance = True
    tendonDynamics = True
    contractVelScale = 2.5 #upscale for sprinting
    maxForceScale = 2 #upscale for sprinting
    
    #Create the model processor for the tracking problem
    
    #Load the model and create processor
    if complexModel:
        osimModel = osim.Model('gaitModel2D_complex.osim')
    else:
        osimModel = osim.Model('gaitModel2D.osim')  
    
    #Update contractile velocity maximum to assist with sprinting sim
    for mm in range(osimModel.getMuscles().getSize()):
        #Get the current muscle
        currMusc = osimModel.getMuscles().get(mm)
        #Get current max contraction velocity
        currVel = currMusc.getMaxContractionVelocity()
        #Set new max contraction velocity
        currMusc.set_max_contraction_velocity(currVel*contractVelScale)
        
    #Finalise connections
    osimModel.finalizeConnections()
    
    #Put the model in a processor
    modelProcessor = osim.ModelProcessor(osimModel)
    
    #Scale the maximum isometric force of muscles to assist with sprinting sim
    modelProcessor.append(osim.ModOpScaleMaxIsometricForce(maxForceScale))
    
    #Disable tendon compliance (will be re-enabled for certain muscles later)
    modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
    
    #Adjust fiber damping as per Dembia et al. Moco paper
    modelProcessor.append(osim.ModOpFiberDampingDGF(0.01))
    
    #Turn off passive forces if appropriate
    if not passiveForces:
        modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
        
    #Modify active force width
    #This appears most important for "de-stiffening" muscles to allow appropriate movement
    modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
        
    #Process model for additional parameters
    studyModel = modelProcessor.process()
    studyModel.initSystem()
    
    #Enable tendon compliance in the gastrocnemius and soleus
    muscles = studyModel.updMuscles()
    for mm in np.arange(muscles.getSize()):
        currMusc = osim.DeGrooteFregly2016Muscle.safeDownCast(muscles.get(int(mm)))
        muscName = currMusc.getName()
        #Enable tendon compliance dynamics in the plantarflexors
        if ('gastroc' in muscName) or ('soleus' in muscName):
            currMusc.set_ignore_tendon_compliance(False)
    
    #Fianlise connections
    studyModel.finalizeConnections()
    
    #Get the model as a processor object
    studyModelProcessor = osim.ModelProcessor(studyModel)
    
    #Set model to use implicit tendon compliance if appropriate
    if implicitTendonCompliance:
        studyModelProcessor.append(osim.ModOpUseImplicitTendonComplianceDynamicsDGF())
     
    #Construct the study object and set basic parameters
    study = osim.MocoStudy()
    study.setName('sprintPrediction_matchedSpeed_muscleDriven')
    
    #Update problem
    problem = study.updProblem()
    
    #Set model
    studyModel = studyModelProcessor.process()
    studyModel.initSystem()
    problem.setModelCopy(studyModel)
    
    #Set the speed goal
    #Get the original speed
    df_kinematics = readSTO('refQ_2D.sto')
    startPos = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'].to_numpy()[0]
    endPos = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'].to_numpy()[-1]
    sprintSpeed = (endPos - startPos) / (finalTime - initialTime)
    #Add a speed goal to go at three times as fast
    problem.addGoal(osim.MocoAverageSpeedGoal('speed'))
    speedGoal = osim.MocoAverageSpeedGoal().safeDownCast(problem.updGoal('speed'))
    speedGoal.set_desired_average_speed(sprintSpeed)
    speedGoal.setWeight(speedWeight)
    
    #Set the time bounds
    problem.setTimeBounds(initialTime, [finalTime - 0.1, finalTime + 0.1])
    
    #Add and set the effort goal to the problem
    problem.addGoal(osim.MocoControlGoal('effort'))
    effort = osim.MocoControlGoal().safeDownCast(problem.updGoal('effort'))
    effort.setDivideByDisplacement(True)
    effort.setWeight(effortWeight)
    effort.setWeightForControl('/forceset/lumbarAct', 0.001)
    
    #Set the symmetry goal
    problem.addGoal(osim.MocoPeriodicityGoal('symmetry'))
    symmetryGoal = osim.MocoPeriodicityGoal().safeDownCast(problem.updGoal('symmetry'))
    #Set symmetric coordinate values (except for pelvis_tx) and speeds
    for cc in range(studyModel.getCoordinateSet().getSize()):
        #Get current coordinate name
        coordName = studyModel.getCoordinateSet().get(cc).getName()
        #Set symmetry parameters based on coordinate
        if coordName.endswith('_r'):
            #Joint angle value
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
                studyModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/value'))
            #Joint speed
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
                studyModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/speed'))
        elif coordName.endswith('_l'):
            #Joint angle value
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
                studyModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/value'))
            #Joint speed
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
                studyModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/speed'))
        elif coordName.endswith('_extension') or coordName.endswith('_tilt') or coordName.endswith('_ty'):
            #Joint angle value
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value'))
            #Joint speed
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed'))
    #Add a symmetry pair for pelvis_tx speed
    symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
        '/jointset/ground_pelvis/pelvis_tx/speed'))
    #Set symmetric activations
    for ff in range(studyModel.getForceSet().getSize()):
        #Check for muscle or actuator
        if (studyModel.getForceSet().get(ff).getConcreteClassName().endswith('Muscle') or
            studyModel.getForceSet().get(ff).getConcreteClassName().endswith('Actuator')):
            #Get force name
            forceName = studyModel.getForceSet().get(ff).getName()
            #Set symmetry parameters based on coordinate
            if forceName.endswith('_r'):
                symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                    studyModel.getForceSet().get(ff).getAbsolutePathString()+'/activation',
                    studyModel.getForceSet().get(ff).getAbsolutePathString().replace('_r','_l')+'/activation'))
            elif forceName.endswith('_l'):
                symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                    studyModel.getForceSet().get(ff).getAbsolutePathString()+'/activation',
                    studyModel.getForceSet().get(ff).getAbsolutePathString().replace('_l','_r')+'/activation'))
            # elif forceName == 'lumbarAct':
            #     symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            #         trackModel.getForceSet().get(ff).getAbsolutePathString()))
    #Add symmetry goal
    symmetryGoal.setWeight(symmetryWeight)
    
    #Add state bounds
    #Set the initial values for the lumbar and pelvis coordinates that
    #reflect 'normal' running motions - bounds based off experimental data
    problem.setStateInfo('/jointset/back/lumbar_extension/value', [np.deg2rad(-30),0],
                         [np.deg2rad(-10),0])
    problem.setStateInfo('/jointset/ground_pelvis/pelvis_tilt/value', [-10*np.pi/180, 0*np.pi/180],
                         [np.deg2rad(-5),0])
    problem.setStateInfo('/jointset/ground_pelvis/pelvis_tx/value', [0, 3],
                         1.06261516975403) #tracking sim starting value
    problem.setStateInfo('/jointset/ground_pelvis/pelvis_ty/value', [0.8, 1.25],
                         [1.00, 1.02]) #likely reasonable starting value
    problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [-25*np.pi/180, 75*np.pi/180])
    problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [-25*np.pi/180, 75*np.pi/180])
    if complexModel:
        problem.setStateInfo('/jointset/walker_knee_l/knee_angle_l/value', [0, 140*np.pi/180])
        problem.setStateInfo('/jointset/walker_knee_r/knee_angle_r/value', [0, 140*np.pi/180])
    else:
        problem.setStateInfo('/jointset/knee_l/knee_angle_l/value', [-140*np.pi/180, 0])
        problem.setStateInfo('/jointset/knee_r/knee_angle_r/value', [-140*np.pi/180, 0])
    problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [-20*np.pi/180, 30*np.pi/180])
    problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [-20*np.pi/180, 30*np.pi/180])
    
    #Set reasonable bounds for tendon force
    problem.setStateInfoPattern('/forceset/.*/normalized_tendon_force', [0, 1.8], [], [])
    
    #Configure the solver
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    solver.resetProblem(problem)
    solver.set_optim_constraint_tolerance(1e-2)
    solver.set_optim_convergence_tolerance(1e-2)
    solver.set_num_mesh_intervals(meshNo)
    
    #Set guess from tracking simulation
    solver.setGuessFile('sprintTracking_muscleDriven_solution.sto')
    
    #Solve
    solution = study.solve()
    
    #Option to visualise
    if visualiseSolution:
        study.visualize(solution)
    
    #Create a full gait cycle trajectory from the periodic solution.
    addPatterns = [".*pelvis_tx/value"]
    fullTraj = osim.createPeriodicTrajectory(solution, addPatterns)
    fullTraj.write('sprintPrediction_matchedSpeed_muscleDriven_fullTrajectory.sto')
    
    #Compute ground reaction forces generated by contact sphere from the 
    #full gait cycle trajectory
    #Set force names
    if complexModel:
        forceNamesRightFoot = ['forceset/contactHeel_r',
                               'forceset/contactMH1_r',
                               'forceset/contactMH3_r',
                               'forceset/contactMH5_r',
                               'forceset/contactHallux_r',
                               'forceset/contactOtherToes_r']
        forceNamesLeftFoot = ['forceset/contactHeel_l',
                                'forceset/contactMH1_l',
                                'forceset/contactMH3_l',
                                'forceset/contactMH5_l',
                                'forceset/contactHallux_l',
                                'forceset/contactOtherToes_l']
    else:
        forceNamesRightFoot = ['forceset/contactHeel_r',
                               'forceset/contactMidfoot_r',
                               'forceset/contactToe_r']
        forceNamesLeftFoot = ['forceset/contactHeel_l',
                                'forceset/contactMidfoot_l',
                                'forceset/contactToe_l']
    #Set force name groups
    forceNames_r = osim.StdVectorString()
    forceNames_l = osim.StdVectorString()
    for ff in range(len(forceNamesRightFoot)):
        forceNames_r.append(forceNamesRightFoot[ff])
        forceNames_l.append(forceNamesLeftFoot[ff])
    #Compute external loads
    externalLoads = osim.createExternalLoadsTableForGait(studyModel,
                                                         fullTraj,
                                                         forceNames_r,
                                                         forceNames_l)
    osim.STOFileAdapter.write(externalLoads,'predictedGRF_matchedSpeed_muscleDriven_2D.mot')

# %% Compare simulations

# This section loads in the simulation data and compares kinematics, muscle function
# (where appropriate), and the GRF data

#Check whether to compare
# if compareResults:
    
#Load in experimental data
df_expKinematics = readSTO('refQ_2D.sto')
df_expGRF = readSTO('refGRF_2D.mot')

#Load in torque driven tracking solution
df_torqueTracking = readSTO('sprintTracking_torqueDriven_solution.sto')
df_torqueTrackingGRF = readSTO('trackedGRF_torqueDriven_2D.mot')

#Load in muscle driven tracking solution
df_muscleTracking = readSTO('sprintTracking_muscleDriven_solution.sto')
df_muscleTrackingGRF = readSTO('trackedGRF_muscleDriven_2D.mot')

#Load in muscle driven predictive solution
df_musclePrediction = readSTO('sprintPrediction_matchedSpeed_muscleDriven_solution.sto')
df_musclePredictionGRF = readSTO('predictedGRF_matchedSpeed_muscleDriven_2D.mot')



#Visualisations...

### TODO: set a better colour palette...

#Set the dictionary colour palette
colourDict = {'expData': '#000000',
              'torqueTracking': '#4885ed',
              'muscleTracking': '#db3236',
              'musclePrediction': '#f4c20d'}

#Get the initial and final time from the experimental kinematics
expInitialTime = osim.Storage('refQ_2D.sto').getFirstTime()
expFinalTime = osim.Storage('refQ_2D.sto').getLastTime()

#GRF...

#Set figure and subplots
fig, ax = plt.subplots(figsize=(6, 3), nrows = 1, ncols = 2)

#Experimental data
#Extract the relevant section of data
vGRF = df_expGRF.loc[(df_expGRF['time'] >= expInitialTime) &
                     (df_expGRF['time'] <= expFinalTime),
                     ['ground_force_r_vy']].to_numpy().flatten()
apGRF = df_expGRF.loc[(df_expGRF['time'] >= expInitialTime) &
                      (df_expGRF['time'] <= expFinalTime),
                      ['ground_force_r_vx']].to_numpy().flatten()
timeVals = df_expGRF.loc[(df_expGRF['time'] >= expInitialTime) &
                         (df_expGRF['time'] <= expFinalTime),
                         ['time']].to_numpy().flatten()
#Normalise the data to 101-points
newTime = np.linspace(timeVals[0],timeVals[-1],101).flatten()
vGRF_norm = np.interp(newTime,timeVals,vGRF)
apGRF_norm = np.interp(newTime,timeVals,apGRF)
#Plot the data
ax[0].plot(np.linspace(0,100,101), vGRF_norm,
           linewidth = 2, color = colourDict['expData'])
ax[1].plot(np.linspace(0,100,101), apGRF_norm,
           linewidth = 2, color = colourDict['expData'])

#Torque tracking data (all data is relevant)
vGRF = df_torqueTrackingGRF['ground_force_r_vy'].to_numpy().flatten()
apGRF = df_torqueTrackingGRF['ground_force_r_vx'].to_numpy().flatten()
timeVals = df_torqueTrackingGRF['time'].to_numpy().flatten()
#Normalise the data to 101-points
newTime = np.linspace(timeVals[0],timeVals[-1],101).flatten()
vGRF_norm = np.interp(newTime,timeVals,vGRF)
apGRF_norm = np.interp(newTime,timeVals,apGRF)
#Plot the data
ax[0].plot(np.linspace(0,100,101), vGRF_norm,
           linewidth = 2, color = colourDict['torqueTracking'])
ax[1].plot(np.linspace(0,100,101), apGRF_norm,
           linewidth = 2, color = colourDict['torqueTracking'])

##### normalised data offset despite times matching???




###



