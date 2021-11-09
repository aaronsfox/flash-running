# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 15:23:54 2021

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This script goes through the process of running some tracking sims on the Dorn
    experimental data. Planned analysis is:
        - Torque driven tracking simulation of experimental kinematics and GRF data
        - Others...?
    
"""

# %% Import packages

#import os
import opensim as osim
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

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

#kinematicsToStates
def kinematicsToStates(kinematicsFileName = None, osimModelFileName = None,
                       outputFileName = 'coordinates.sto',
                       inDegrees = True, outDegrees = False):
    
    # Convenience function for converting IK results to a states storage.
    #
    # Input:    kinematicsFileName - file containing kinematic data. Header should only be coordinates name, rather than path to state
    #           osimModelFileName - opensim model filename that corresponds to kinematic data
    #           outputFileName - optional filename to output to (defaults to coordinates.sto)
    #           inDegrees - set to true if kinematics file is in degrees (defaults to True)
    #           outDegrees - set to true if desired output is in degrees (defaults to False)

    if kinematicsFileName is None:
        raise ValueError('Filename for kinematics is required')
    if osimModelFileName is None:
        raise ValueError('OpenSim model filename is required')

    #Load in the kinematic data
    kinematicsStorage = osim.Storage(kinematicsFileName)
    
    #Create a copy of the kinematics data to alter the column labels in
    statesStorage = osim.Storage(kinematicsFileName)
    
    #Resample the data points linearly to avoid any later issues with matching
    #time points. Use a time stamp for 250 Hz
    kinematicsStorage.resampleLinear(1/250)
    statesStorage.resampleLinear(1/250)
    
    #Get the column headers for the storage file
    angleNames = kinematicsStorage.getColumnLabels()
    
    #Get the corresponding full paths from the model to rename the
    #angles in the kinematics file
    kinematicModel = osim.Model(osimModelFileName)
    for ii in range(0,angleNames.getSize()):
        currAngle = angleNames.get(ii)
        if currAngle != 'time':
            #Get full path to coordinate
            fullPath = kinematicModel.updCoordinateSet().get(currAngle).getAbsolutePathString()+'/value'
            #Set angle name appropriately using full path
            angleNames.set(ii,fullPath)
    
    #Set the states storage object to have the updated column labels
    statesStorage.setColumnLabels(angleNames)
    
    #Appropriately set output in degrees or radians
    if inDegrees and not outDegrees:
        #Convert degrees values to radians for consistency with the current
        #file label (defaults back to inDegrees=no). Radians seem to work
        #better with the Moco process as well.
        kinematicModel.initSystem()
        kinematicModel.getSimbodyEngine().convertDegreesToRadians(statesStorage)
    elif inDegrees and outDegrees:
        #Change the storage label back to specifying indegrees=yes
        statesStorage.setInDegrees(True)
    elif not inDegrees and outDegrees:
        #Convert radians to degrees
        kinematicModel.initSystem()
        kinematicModel.getSimbodyEngine().convertRadiansToDegrees(statesStorage)
        #Reset labeling for degrees
        statesStorage.setInDegrees(True)
    
    #Write the states storage object to file
    statesStorage.printToXML(outputFileName)

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

# %% 3D basic tracking sim of experimental data (torque driven)

# This step is designed to generate relevant controls for the 2D model that
# appropriately track the sprint kinematics and GRF from the experimental data
# using torque actuators

#Check whether to run simulation
if runTorqueTracking:
    
    #Convert ik kinematics to states
    kinematicsToStates(kinematicsFileName = 'sprint_ik.mot',
                       osimModelFileName = 'JA1_SCALED_Osim40.osim',
                       outputFileName = 'sprint_ik_states.sto',
                       inDegrees = True, outDegrees = False)
    
    #Set simulation parameters
    
    #Timing and mesh data
    initialTime = osim.Storage('sprint_ik_states.sto').getFirstTime()
    finalTime = osim.Storage('sprint_ik_states.sto').getLastTime()
    duration = finalTime - initialTime
    meshNo = 50
    meshInterval = duration / meshNo
    
    #Weights
    trackingWeight = 10
    effortWeight = 0.1
    grfWeight = 1
    speedWeight = 1
    symmetryWeight = 1
    
    #General parameters
    visualiseSolution = False
    
    #Create the model processor for the tracking problem
    
    #Load the model and create processor
    osimModel = osim.Model('JA1_SCALED_Osim40.osim')
    #### NOTE: this model had some generic contact geometry added to it
    modelProcessor = osim.ModelProcessor(osimModel)
    
    ### Currently no muscles in this model
    # #Add processing step to remove muscles
    # modelProcessor.append(osim.ModOpRemoveMuscles())
    
    #Process model for additional parameters
    torqueModel = modelProcessor.process()
    torqueModel.initSystem()
    
    #Add torque actuators to support model
    #Set list to add reserves to
    #### Consider locking subtalar and mtp????
    reserveList = ['pelvis_tilt', 'pelvis_list', 'pelvis_rotation',
                   'hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r',
                   'knee_angle_r',
                   'ankle_angle_r', 'subtalar_angle_r', 'mtp_angle_r',
                   'hip_flexion_l', 'hip_adduction_l', 'hip_rotation_l',
                   'knee_angle_l',
                   'ankle_angle_l', 'subtalar_angle_l', 'mtp_angle_l',
                   'lumbar_extension', 'lumbar_bending', 'lumbar_rotation',
                   'arm_flex_r', 'arm_add_r', 'arm_rot_r',
                   'elbow_flex_r', 'pro_sup_r', 'wrist_flex_r', 'wrist_dev_r',
                   'arm_flex_l', 'arm_add_l', 'arm_rot_l',
                   'elbow_flex_l', 'pro_sup_l', 'wrist_flex_l', 'wrist_dev_l']
        
    #Set optimal force and max torque
    optimalForce = 100
    maxTorque = np.inf
    #Add torque actuators
    #### NOTE: this function needs updating with respect to labels for Dorn model
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
    tableProcessor = osim.TableProcessor('sprint_ik_states.sto')
    tableProcessor.append(osim.TabOpLowPassFilter(12))
    tableProcessor.append(osim.TabOpUseAbsoluteStateNames())
    track.setStatesReference(tableProcessor)
    track.set_states_global_tracking_weight(trackingWeight)
    #Set some coordinates not to be tracked to avoid poor motion
    #Similar process done in Dembia et al. Moco paper
    stateWeights = osim.MocoWeightSet()
    weightList = list()
    weightList.append(('/jointset/ground_pelvis/pelvis_ty', 0))
    weightList.append(('/jointset/ground_pelvis/pelvis_tilt', 0.1))
    # weightList.append(('/jointset/ground_pelvis/pelvis_list', 0.1))
    # weightList.append(('/jointset/ground_pelvis/pelvis_rotation', 0.1))
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
    
    #Set to not highly penalise the lumbar and arm actuators
    for forceNo in range(trackModel.updForceSet().getSize()):
        if 'lumbar' in trackModel.updForceSet().get(forceNo).getAbsolutePathString() or \
            'arm' in trackModel.updForceSet().get(forceNo).getAbsolutePathString() or \
                'wrist' in trackModel.updForceSet().get(forceNo).getAbsolutePathString() or \
                    'elbow' in trackModel.updForceSet().get(forceNo).getAbsolutePathString() or \
                        'pro_sup' in trackModel.updForceSet().get(forceNo).getAbsolutePathString():
                            effort.setWeightForControl(trackModel.updForceSet().get(forceNo).getAbsolutePathString(),
                                                       0.001)

    #Set an average speed goal based on sprinting data
    df_kinematics = readSTO('sprint_ik_states.sto')
    startPos = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'].to_numpy()[0]
    endPos = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'].to_numpy()[-1]
    sprintSpeed = (endPos - startPos) / (finalTime - initialTime)
    speedGoal = osim.MocoAverageSpeedGoal('speed')
    speedGoal.set_desired_average_speed(sprintSpeed)
    speedGoal.setWeight(speedWeight)
    problem.addGoal(speedGoal)
    
    #Set symmetry constraint
    #Create symmetry goal
    #### NOTE: this is a full gait cycle - so every coordinate should be matched with itself...
    symmetryGoal = osim.MocoPeriodicityGoal('symmetryGoal')
    #Set symmetric coordinate values (except for pelvis_tx) and speeds
    for cc in range(trackModel.getCoordinateSet().getSize()):
        #Get current coordinate name
        coordName = trackModel.getCoordinateSet().get(cc).getName()
        #Set symmetry for nearly all coordinates with self
        if '_tx' not in coordName:
            #Joint angle value
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value'))
            #Joint speed
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed'))        
        # #Set symmetry parameters based on coordinate
        # if coordName.endswith('_r'):
        #     #Joint angle value
        #     symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
        #         trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
        #         trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/value'))
        #     #Joint speed
        #     symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
        #         trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
        #         trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/speed')
        ### Doing this twice just doubles up right???
        # elif coordName.endswith('_l'):
        #     #Joint angle value
        #     symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
        #         trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
        #         trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/value'))
        #     #Joint speed
        #     symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
        #         trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
        #         trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/speed'))
        # elif coordName.endswith('_extension') or coordName.endswith('_tilt') or coordName.endswith('_ty'):
        #     #Joint angle value
        #     symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
        #         trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value'))
        #     #Joint speed
        #     symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
        #         trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed'))
    #Add a symmetry pair for pelvis_tx speed
    symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
        '/jointset/ground_pelvis/pelvis_tx/speed'))
    #Add symmetry goal
    symmetryGoal.setWeight(symmetryWeight)
    problem.addGoal(symmetryGoal)
    
    #Add contact tracking goal
    #Set force names
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
    
    #Create contact tracking goal
    contactGoal = osim.MocoContactTrackingGoal('contactGoal', grfWeight)
    #Set external loads
    contactGoal.setExternalLoadsFile('sprint_grf.xml')
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
    #Do this based off max and min coordinate values
    #Could be better (e.g. based off experimental data max/mins)
    for cc in range(trackModel.getCoordinateSet().getSize()):
        #Get max and min values
        mx = trackModel.getCoordinateSet().get(cc).getRangeMax()
        mn = trackModel.getCoordinateSet().get(cc).getRangeMin()
        #Set state info
        problem.setStateInfo(trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
                             [mn, mx])
    
    ###### Optimise contact parameters????? #####
    
    #Add parameter
    # heel_r_size = osim.MocoParameter()
    # heel_r_size.setName('heel_r_size')
    # heel_r_size.appendComponentPath('/contactgeometryset/heel_r')
    # heel_r_size.setPropertyName('radius')
    # heel_r_size.setBounds(osim.MocoBounds(0.025,0.045))
    # problem.addParameter(heel_r_size)
    
    
    ######
        
    # #Configure the solver
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    solver.resetProblem(problem)
    solver.set_optim_constraint_tolerance(1e-2) ### probably a bit high, but doing so for speedy solution
    solver.set_optim_convergence_tolerance(1e-2) ### probably a bit high, but doing so for speedy solution
    # solver.set_parameters_require_initsystem(False) #### for parameters --- don't think it's required
    
    #Solve
    solution = study.solve()
    
    ##### Stopped 3D solution after about 45 mins...kind of a disaster...
    
    #Print solution to file (doesn't seem to do this automatically anymore?)
    solution.write('sprintTracking_torqueDriven_solution.sto')
        
    #Option to visualise
    if visualiseSolution:
        study.visualize(solution)
    
    #####.....check trackingSims.py from Simulations for additional info

