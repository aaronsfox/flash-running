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
    
    TODO: Think there are some issues with using the same names for models across
    sections of the script --- doesn't update parameters?
    
"""

# %% Import packages

import os
import opensim as osim
import osimFunctions as helper
import numpy as np
import matplotlib.pyplot as plt

# %% Set-up

#Some basic parameters to determine how the function runs.

#Specifically here we provide options to run or load the simulations
#Defaults here are False to not run as results can be imported
runTorqueTracking = False 
runMuscleTracking = False 
runMusclePrediction = False 

#Specify whether to visualise solutions
#This is set to false by defulat as it can interupt analysis flow
visualiseSolution = False

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

#Add OpenSim geometry path (weird issues with this on new laptop)
osim.ModelVisualizer.addDirToGeometrySearchPaths('C:\\OpenSim 4.3\\Geometry')

#Set right and left contact sphere groups
rightFootContacts = ['/forceset/contactHeel_r', '/forceset/contactMidfoot_r', '/forceset/contactToe_r']
leftFootContacts = ['/forceset/contactHeel_l', '/forceset/contactMidfoot_l', '/forceset/contactToe_l']

# %% Run a 2D torque-driven tracking sim of experimental data

# This step is designed to generate a solution for the 2D model that appropriately
# track the sprint kinematics and GRF from the experimental data using torque actuators

"""

    > Symmetry goal seems appropriate but does change the state and contact tracking accuracy...
    > Messing with the weights too impacts convergence...
        > 0.1 for states and 1 for contact weights seems to get OK results
            > Although the motion looks a little wrong
            > Adding the symmetry goal to this converges but produces a bad motion

"""

#Check whether to run simulation
if runTorqueTracking:
    
    #Make directory for torque driven simulation results
    if not os.path.isdir('torqueDrivenTracking'):
        os.mkdir('torqueDrivenTracking')
    
    #Set-up logger for tracking simulation
    osim.Logger.removeFileSink()
    osim.Logger.addFileSink('torqueDrivenTracking\\torqueDrivenTrackingLog.log')
    
    #Load the model
    gaitModel2D = osim.Model('gaitModel2D.osim')
    
    #Add in the surface constraints from the dummy model
    #Load in the model
    dummyConstraintModel = osim.Model('..\\..\\ExpData\\Scaling\\gait9dof18musc_Ong_et_al_Moco_withSurfaceConstraints.osim')
    #Add the rolling on surface constraint
    #This can't be accessed through the API so needs to be worked around by
    #stealing it from a dummy model. The dummy model contains 4 rolling surface
    #constraints that we adopt across
    for conInd in range(dummyConstraintModel.updConstraintSet().getSize()):
        gaitModel2D.updConstraintSet().cloneAndAppend(dummyConstraintModel.updConstraintSet().get(conInd))
        
    #Finalise model connections
    gaitModel2D.finalizeConnections()
    
    #Create a model processor
    modelProcessor = osim.ModelProcessor(gaitModel2D)
    
    #Delete muscles
    modelProcessor.append(osim.ModOpRemoveMuscles())
    
    #Get model from the processor
    torqueModel = modelProcessor.process()
    
    #Create the dictionary of optimal forces to pass to function for torque actuators
    optForces = {'hip_flexion_r': 300, 'hip_flexion_l': 300,
                 'knee_angle_r': 300, 'knee_angle_l': 300,
                 'ankle_angle_r': 300, 'ankle_angle_l': 300,
                 'lumbar_extension': 300}
    
    #Add the actuators to the torque driven model
    torqueModel = helper.addTorqueActuators(osimModel = torqueModel,
                                            optForces = optForces,
                                            minControl = np.inf * -1,
                                            maxControl = np.inf)
    
    #Update friction parameters in contact model
    #This helps with smoother GRF predictions (I think...)
    #Settings for coefficients
    staticFriction = 0.5 #typical model = 0.8
    dynamicFriction = 0.5 #typical model = 0.8
    viscousFriction = 0.5 #typical model = 0.5
    #Update friction coefficients in forceset
    for forceInd in range(torqueModel.updForceSet().getSize()):
        #Check for contact geometry force
        if 'contact' in torqueModel.updForceSet().get(forceInd).getName():
            #Update the friction parameters
            osim.SmoothSphereHalfSpaceForce.safeDownCast(torqueModel.updForceSet().get(forceInd)).set_static_friction(staticFriction)
            osim.SmoothSphereHalfSpaceForce.safeDownCast(torqueModel.updForceSet().get(forceInd)).set_dynamic_friction(dynamicFriction)
            osim.SmoothSphereHalfSpaceForce.safeDownCast(torqueModel.updForceSet().get(forceInd)).set_viscous_friction(viscousFriction)
    
    #Finalise model connections
    torqueModel.finalizeConnections()
    
    #Save torque model to file for later use
    torqueModel.printToXML('torqueDrivenTracking\\torqueModel.osim')
    
    #Define the motion tracking problem
    track = osim.MocoTrack()
    track.setName('torqueDrivenTracking')
    
    #Set model
    track.setModel(osim.ModelProcessor(torqueModel))
    
    #Set kinematics
    tableProcessor = osim.TableProcessor('refQ_2D.sto')
    
    #Set states reference details
    track.setStatesReference(tableProcessor) # Apply the target data to the tracking problem
    track.set_states_global_tracking_weight(0.1) # Default tracking weight (is changed below) --- NOTE: have been using 1
    track.set_allow_unused_references(True) # Target data can include DoF not in this model
    track.set_track_reference_position_derivatives(True) # Track speed trajectories
    track.set_apply_tracked_states_to_guess(True) # Use target data in initial guess
    
    #Set times
    track.set_initial_time(osim.Storage('refQ_2D.sto').getFirstTime())
    track.set_final_time(osim.Storage('refQ_2D.sto').getLastTime())
    
    #Create weight set for state tracking
    stateWeights = osim.MocoWeightSet()
    
    #Create a dictionary that provides the kinematic task weights for function
    taskWeights = {'pelvis_tx': 500, 'pelvis_ty': 250,
                   'pelvis_tilt': 250, 
                   'hip_flexion_r': 100, 'knee_angle_r': 100,
                   'ankle_angle_r': 50,
                   'hip_flexion_l': 100, 'knee_angle_l': 100,
                   'ankle_angle_l': 50,
                   'lumbar_extension': 100}
    
    #Set constant weight to scale tracking error speeds by
    w = 0.001
    
    #Loop through coordinates to apply weights
    for coordInd in range(torqueModel.updCoordinateSet().getSize()):
        
        #Get name and absolute path to coordinate
        coordName = torqueModel.updCoordinateSet().get(coordInd).getName()
        coordPath = torqueModel.updCoordinateSet().get(coordInd).getAbsolutePathString()

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
    
    #Define the Moco study and problem
    study = track.initialize()
    problem = study.updProblem()
    
    ##### CREATE GOALS #####
    
    #Regularization term on MocoTrack problem (minimize squared muscle excitations)
    effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
    effort.setWeight(0.001)
    
    #Add contact tracking goal
    #Uses two separate vectors to apply different contact tracking weights    
    
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
                                          osim.Vec3(0,1,0)],
                               'weight': [2, 1]
                               }
    
    #Loop through contact tracking settings to create separate weighted goals
    for ii in range(len(contactTrackingSettings['vector'])):
        
        #Create tracking goal
        contactTracking.append(osim.MocoContactTrackingGoal(f'contact{ii}',
                                                            contactTrackingSettings['weight'][ii]))
        
        #Set external loads
        contactTracking[ii].setExternalLoadsFile('refGRF.xml')
        
        #Add the left and right tracking groups to the goal
        contactTracking[ii].addContactGroup(trackLeftGRF)
        contactTracking[ii].addContactGroup(trackRightGRF)
        
        #Set projection vector in problem
        contactTracking[ii].setProjection('vector')
        contactTracking[ii].setProjectionVector(contactTrackingSettings['vector'][ii])
    
        #Add to problem
        problem.addGoal(contactTracking[ii])
        
    #Create a symmetry goal to permit simulating one step
    symmetryGoal = osim.MocoPeriodicityGoal('symmetry')
    
    #Create a dictionary to prescribe state symmetry pairs or isolated states
    symmetryPairs = {
        #Joint coordinate
        'hip_angle': ['/jointset/hip_r/hip_flexion_r/value',
                      '/jointset/hip_l/hip_flexion_l/value'],
        'knee_angle': ['/jointset/knee_r/knee_angle_r/value',
                        '/jointset/knee_l/knee_angle_l/value'],
        'ankle_angle': ['/jointset/ankle_r/ankle_angle_r/value',
                        '/jointset/ankle_l/ankle_angle_l/value'],
        'lumbar_angle': ['/jointset/back/lumbar_extension/value'],
        'pelvis_angle': ['/jointset/ground_pelvis/pelvis_tilt/value'],
        'pelvis_ty': ['/jointset/ground_pelvis/pelvis_ty/value']
        }
    
    #Set symmetric pairings based on dictionary
    for symmetry in list(symmetryPairs.keys()):
        
        #Check for state pair or individual state
        if len(symmetryPairs[symmetry]) == 2:
            
            #Add state pair
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(symmetryPairs[symmetry][0],
                                                                    symmetryPairs[symmetry][1]))
            #Add reverse of state pair
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(symmetryPairs[symmetry][1],
                                                                    symmetryPairs[symmetry][0]))
            
        else:
            
            #Add individual state pairing
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(symmetryPairs[symmetry][0]))

    #Add to problem
    problem.addGoal(symmetryGoal)
    
    ##### SET BOUNDS #####
    
    #Create a dictionary to set kinematic bounds
    #Create this based on some generic values above/below the values in the
    #experimental data    
    
    #Load the experimental kinematics file as a table
    expKinematicsTable = osim.TimeSeriesTable('refQ_2D.sto')
    
    #Create the bounds dictionary
    kinematicBounds = {
        'pelvis_tilt': [expKinematicsTable.getDependentColumn('/jointset/ground_pelvis/pelvis_tilt/value').to_numpy().min() - np.deg2rad(10),
                        expKinematicsTable.getDependentColumn('/jointset/ground_pelvis/pelvis_tilt/value').to_numpy().max() + np.deg2rad(10)],
        'pelvis_tx': [expKinematicsTable.getDependentColumn('/jointset/ground_pelvis/pelvis_tx/value').to_numpy().min() - 0.2,
                      expKinematicsTable.getDependentColumn('/jointset/ground_pelvis/pelvis_tx/value').to_numpy().max() + 0.2],
        'pelvis_ty': [expKinematicsTable.getDependentColumn('/jointset/ground_pelvis/pelvis_ty/value').to_numpy().min() - 0.1,
                      expKinematicsTable.getDependentColumn('/jointset/ground_pelvis/pelvis_ty/value').to_numpy().max() + 0.1],
        'hip_flexion_l': [expKinematicsTable.getDependentColumn('/jointset/hip_l/hip_flexion_l/value').to_numpy().min() - np.deg2rad(15),
                          expKinematicsTable.getDependentColumn('/jointset/hip_l/hip_flexion_l/value').to_numpy().max() + np.deg2rad(15)],
        'knee_angle_l': [np.deg2rad(-140), np.deg2rad(0)],
        'ankle_angle_l': [expKinematicsTable.getDependentColumn('/jointset/ankle_l/ankle_angle_l/value').to_numpy().min() - np.deg2rad(15),
                          expKinematicsTable.getDependentColumn('/jointset/ankle_l/ankle_angle_l/value').to_numpy().max() + np.deg2rad(15)],
        'hip_flexion_r': [expKinematicsTable.getDependentColumn('/jointset/hip_r/hip_flexion_r/value').to_numpy().min() - np.deg2rad(15),
                          expKinematicsTable.getDependentColumn('/jointset/hip_r/hip_flexion_r/value').to_numpy().max() + np.deg2rad(15)],
        'knee_angle_r': [np.deg2rad(-140), np.deg2rad(0)],
        'ankle_angle_r': [expKinematicsTable.getDependentColumn('/jointset/ankle_r/ankle_angle_r/value').to_numpy().min() - np.deg2rad(15),
                          expKinematicsTable.getDependentColumn('/jointset/ankle_r/ankle_angle_r/value').to_numpy().max() + np.deg2rad(15)],
        'lumbar_extension': [expKinematicsTable.getDependentColumn('/jointset/back/lumbar_extension/value').to_numpy().min() - np.deg2rad(10),
                             5]
        }
    
    #Set the bounds using the dictionary values and experimental data
    for coordInd in range(torqueModel.updCoordinateSet().getSize()):
        
        #First check if coordinate is in kinematic bounds dictionary
        if torqueModel.updCoordinateSet().get(coordInd).getName() in list(kinematicBounds.keys()):
            
            #Get coordinate name for ease of use
            coordName = torqueModel.updCoordinateSet().get(coordInd).getName()

            #Set bounds in problem
            problem.setStateInfo(f'{torqueModel.updCoordinateSet().get(coordInd).getAbsolutePathString()}/value',
                                  #Bounds set to model ranges
                                  [kinematicBounds[coordName][0], kinematicBounds[coordName][1]]
                                  )
    
    ##### SOLVER SETTINGS #####
    
    #Define the solver
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    
    # #Get the initial guess to manipulate
    # initialGuess = solver.getGuess()
    
    #### TODO: is this actually needed for these data???
    #### Seems like it works just as well without and avoids a drop in the model
    #### to begin with...
    #### Might need with symmetry goal though...?
    #Push up pelvis_ty in starting guess to avoid contact spheres starting in ground
    #Set the pelvis_ty starting value to a small increment above where it is in coordinates data
    # initialGuess.setState('/jointset/ground_pelvis/pelvis_ty/value',
    #                       initialGuess.getState('/jointset/ground_pelvis/pelvis_ty/value').to_numpy() + 0.05)
    
    #Reset problem
    solver.resetProblem(problem)
    
    #Set updated guess in solver
    # solver.setGuess(initialGuess)
    
    #Set solver options
    #Run an initial pass with a coarse mesh to set as an initial guess in a
    #subsequent finer solution
    solver.set_optim_max_iterations(2000)
    # solver.set_num_mesh_intervals(15)
    solver.set_minimize_implicit_multibody_accelerations(True)
    solver.set_optim_constraint_tolerance(1e-2) 
    solver.set_optim_convergence_tolerance(1e-2)
    
    #Solve!
    #### TODO: add check if coarse solution exists? Or is this covered by the run boolean?
    # coarseSolution = study.solve()
    
    #Write coarse solution to file
    # coarseSolution.write('torqueDriven\\torqueDrivenSolution_coarse.sto')
    
    # #Set updated guess in solver for finer grid solution
    # solver.resetProblem(problem)
    # solver.setGuessFile('torqueDriven\\torqueDrivenSolution_coarse.sto')
    
    #Update solver to options to finer mesh
    solver.set_num_mesh_intervals(50)
        
    #Solve with finer mesh!
    solution = study.solve()

    #Write to file
    solution.write('torqueDrivenTracking\\torqueDrivenTrackingSolution.sto')

    #Visualise
    if visualiseSolution:
        study.visualize(solution)
    
    #Extract predicted GRFs
    
    #Add contact elements to extract from
    #Get right side
    contact_r = osim.StdVectorString()
    for contactLabel in rightFootContacts:
        contact_r.append(contactLabel)
    #Get left side
    contact_l = osim.StdVectorString()
    for contactLabel in leftFootContacts:
        contact_l.append(contactLabel)
            
    #Create forces table
    externalForcesTableFlat = osim.createExternalLoadsTableForGait(torqueModel, solution,
                                                                   contact_r, contact_l)
    
    #Write table to file
    osim.STOFileAdapter().write(externalForcesTableFlat,
                                'torqueDrivenTracking\\torqueDrivenTrackingSolution_grf.sto')
    
    #Remove tracked states file
    os.remove('torqueDrivenTracking_tracked_states.sto')
    
#Compare torque driven simulation to original experimental data
#### TODO: wrap both of these in a helper function that compares two solutions...

#Read in experimental coordinates data
expCoordinates = osim.TimeSeriesTable('refQ_2D.sto')

#Read in experimental GRF data
expGRF = osim.TimeSeriesTable('refGRF.mot')

#Read in torque diven solution data
torqueDrivenSolution = osim.TimeSeriesTable('torqueDrivenTracking\\torqueDrivenTrackingSolution.sto')

#Read in experimental GRF data
torqueDrivenSolutionGRF = osim.TimeSeriesTable('torqueDrivenTracking\\torqueDrivenTrackingSolution_grf.sto')

#Load torque model
torqueModel = osim.Model('torqueDrivenTracking\\torqueModel.osim')

#Set names for labelling data
expDataLabel = 'Experimental'
solDataLabel = 'Torque Driven Solution'

#Kinematics

#Set kinematic variables to plot
kinematicVars = ['lumbar_extension', 
                 'pelvis_tilt', 'pelvis_tx', 'pelvis_ty',
                 'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r',
                 'hip_flexion_l', 'knee_angle_l', 'ankle_angle_l']

#Create dictionary to match up kinematic variables to axes
kinematicAxes = {'lumbar_extension': [0,0],
                 'pelvis_tilt': [1,0], 'pelvis_tx': [1,1], 'pelvis_ty': [1,2],
                 'hip_flexion_r': [2,0], 'knee_angle_r': [2,1], 'ankle_angle_r': [2,2],
                 'hip_flexion_l': [3,0], 'knee_angle_l': [3,1], 'ankle_angle_l': [3,2]}

#Create subplot
fig, ax = plt.subplots(figsize = (9,15), nrows = 4, ncols = 3)

#Adjust subplots
plt.subplots_adjust(left = 0.075, right = 0.975, bottom = 0.05, top = 0.97,
                    hspace = 0.4, wspace = 0.4)

#Loop through kinematic variables
for var in kinematicVars:
    
    #Plot experimental data
    
    #Set current axes to plot on
    plotAx = ax[kinematicAxes[var][0], kinematicAxes[var][1]]
    
    #Create the full path string for current variable
    varPathStr = torqueModel.updCoordinateSet().get(var).getAbsolutePathString()+'/value'
    
    #Get the time data for experimental and solution
    expTime = np.array(expCoordinates.getIndependentColumn())
    solTime = np.array(torqueDrivenSolution.getIndependentColumn())
    
    #Get the joint coordinate data for experimental and solution
    expData = expCoordinates.getDependentColumn(varPathStr).to_numpy()
    solData = torqueDrivenSolution.getDependentColumn(varPathStr).to_numpy()
    
    #Check whether to convert to radians on the basis of whether the data is a
    #translation or not
    if var not in ['pelvis_tx', 'pelvis_ty']:
        #Convert data
        expData = np.rad2deg(expData)
        solData = np.rad2deg(solData)
        #Set y-label
        yLabel = 'Joint Angle (\u00b0)'
    else:
        #Just set y-label
        yLabel = 'Position (m)'
    
    #Plot the experimental data against time
    plotAx.plot(expTime, expData,
                color = '#000000', linewidth = 1.5,
                zorder = 2, label = expDataLabel)
    
    #Plot the solution data against time
    plotAx.plot(solTime, solData,
                color = '#db3236', linewidth = 1.5,
                zorder = 2, label = solDataLabel)
    
    #Set axes limits to time bounds
    plotAx.set_xlim([solTime[0], solTime[-1]])
    
    #Add labels
    #X-axis
    plotAx.set_xlabel('Time', fontsize = 10, fontweight = 'bold')
    
    #Y-axis
    plotAx.set_ylabel(yLabel, fontsize = 10, fontweight = 'bold')
    
    #Set title
    plotAx.set_title(var.replace('_',' ').title(),
                     pad = 3, fontsize = 12, fontweight = 'bold')
    
    #Add zero-dash line if necessary
    if plotAx.get_ylim()[0] < 0 < plotAx.get_ylim()[-1]:
        plotAx.axhline(y = 0, color = 'dimgrey',
                       linewidth = 0.5, ls = ':', 
                       zorder = 1)
        
    #Turn off top and right axes spines
    plotAx.spines['top'].set_visible(False)
    plotAx.spines['right'].set_visible(False)
    
    #Set ticks to inwards
    plotAx.tick_params(axis = 'both', direction = 'in', length = 3)
        
#Turn off two unused axes
ax.flatten()[1].axis('off')
ax.flatten()[2].axis('off')

#Save figure
plt.savefig('torqueDrivenTracking\\torqueDrivenTracking_experimentalData_comparisonKinematics.png',
            format = 'png', dpi = 300)

#Close figure
plt.close()

#Compare experimental and torque driven GRFs

#Set GRF variables to plot
grfVars = ['ground_force_r_vx', 'ground_force_r_vy']

#Set titles for variable
grfTitles = ['Ant. / Post. GRFs', 'Vertical GRFs']

#Create subplot
fig, ax = plt.subplots(figsize = (9,4), nrows = 1, ncols = 2)

#Adjust subplots
plt.subplots_adjust(left = 0.075, right = 0.975, bottom = 0.1, top = 0.9,
                    hspace = 0.4, wspace = 0.4)

#Loop through kinematic variables
for var in grfVars:
    
    #Get the time data for experimental and solution
    expTime = np.array(expGRF.getIndependentColumn())
    solTime = np.array(torqueDrivenSolutionGRF.getIndependentColumn())
    
    #Plot experimental data
    
    #Get all column labels that contain the force component
    #Set the force component based on the variable name
    forceComponent = var.split('_')[-1]
    #Get all values for current component from experimental data
    expForceList = [ii for ii in list(expGRF.getColumnLabels()) if forceComponent in ii]
    
    #Sum all of the experimental force components
    #Initialise array of zeros
    expData = np.zeros(len(expGRF.getIndependentColumn()))
    #Loop through and add forces
    for expForce in expForceList:
        expData += expGRF.getDependentColumn(expForce).to_numpy()
    
    #Get the GRF data from the solution
    solData = torqueDrivenSolutionGRF.getDependentColumn(var).to_numpy()
    
    #Plot the experimental data against time
    ax[grfVars.index(var)].plot(expTime, expData,
                                color = '#000000', linewidth = 1.5,
                                zorder = 2, label = expDataLabel)
    
    #Plot the solution data against time
    ax[grfVars.index(var)].plot(solTime, solData,
                                color = '#db3236', linewidth = 1.5,
                                zorder = 2, label = solDataLabel)
    
    #Set axes limits to time bounds
    ax[grfVars.index(var)].set_xlim([solTime[0], solTime[-1]])
    
    #Add labels
    #X-axis
    ax[grfVars.index(var)].set_xlabel('Time', fontsize = 10, fontweight = 'bold')
    
    #Y-axis
    ax[grfVars.index(var)].set_ylabel('Force (N)', fontsize = 10, fontweight = 'bold')
    
    #Set title
    ax[grfVars.index(var)].set_title(grfTitles[grfVars.index(var)],
                                     pad = 3, fontsize = 12, fontweight = 'bold')
    
    #Add zero-dash line if necessary
    if ax[grfVars.index(var)].get_ylim()[0] < 0 < ax[grfVars.index(var)].get_ylim()[-1]:
        ax[grfVars.index(var)].axhline(y = 0, color = 'dimgrey',
                                       linewidth = 0.5, ls = ':', 
                                       zorder = 1)
        
    #Turn off top and right axes spines
    ax[grfVars.index(var)].spines['top'].set_visible(False)
    ax[grfVars.index(var)].spines['right'].set_visible(False)
    
    #Set ticks to inwards
    ax[grfVars.index(var)].tick_params(axis = 'both', direction = 'in', length = 3)
        
#Save figure
plt.savefig('torqueDrivenTracking\\torqueDrivenTracking_experimentalData_comparisonGRFs.png',
            format = 'png', dpi = 300)

#Close figure
plt.close()

# %% TEST MUSCLE DRIVEN INVERSE -----TEST-----

#Make directory for torque driven simulation results
if not os.path.isdir('muscleDrivenInverse'):
    os.mkdir('muscleDrivenInverse')
    
#Set-up logger for tracking simulation
osim.Logger.removeFileSink()
osim.Logger.addFileSink('muscleDrivenInverse\\muscleDrivenInverseLog.log')

#Load the model
gaitModel2D = osim.Model('gaitModel2D.osim')

######## Replace wth larger muscle set???
    ##### Idea is to actuate right limb with muscles, left with torque actuators
        #### Can be done in experimental data function if this actually works....
    
#Remove original muscles from 2D model
#This needs to be done in 2 steps given the index of a force will reduce each time
#a muscle is taken out
removeInd = []
for forceInd in range(gaitModel2D.updForceSet().getSize()):
    if 'Muscle' in gaitModel2D.updForceSet().get(forceInd).getConcreteClassName():
        removeInd.append(forceInd)
#Remove muscles subtracting index value each time
for remove in removeInd:
    gaitModel2D.updForceSet().remove(remove-removeInd.index(remove))

#Replace right limb muscles with those from the 3D model
#Load the model
gaitModel3D = osim.Model('..\\..\\ExpData\\Scaling\\scaledModelAdjusted.osim')
#Get the muscles, only taking the right hand side
for forceInd in range(gaitModel3D.updForceSet().getSize()):
    if 'Muscle' in gaitModel3D.updForceSet().get(forceInd).getConcreteClassName() and \
        gaitModel3D.updForceSet().get(forceInd).getName().endswith('_r'):
            gaitModel2D.updForceSet().cloneAndAppend(gaitModel3D.updForceSet().get(forceInd))
            
#Update the pelvis and right hand side bodies from 3D model
#Once again this needs to be done in 2 steps given the change in index with removal
removeInd = []
removeName = []
for bodyInd in range(gaitModel2D.updBodySet().getSize()):
    if 'pelvis' in gaitModel2D.updBodySet().get(bodyInd).getName() or \
        gaitModel2D.updBodySet().get(bodyInd).getName().endswith('_r'):
            removeInd.append(bodyInd)
            removeName.append(gaitModel2D.updBodySet().get(bodyInd).getName())
#Remove bodies subtracting each time
for remove in removeInd:
    gaitModel2D.updBodySet().remove(remove-removeInd.index(remove))
#Add in the bodies that were removed from the 3D model
for addName in removeName:
    gaitModel2D.updBodySet().cloneAndAppend(gaitModel3D.updBodySet().get(addName))
            
#Clone and append the right patellofemoral body & joint
gaitModel2D.updBodySet().cloneAndAppend(gaitModel3D.updBodySet().get('patella_r'))
gaitModel2D.updJointSet().cloneAndAppend(gaitModel3D.updJointSet().get('patellofemoral_r'))

#### TODO: patella coordinate coupler constraints need to be added! Just right side!
    #### NOTE: negate the linear function coefficient to get knee angle beta correct!

#Remove lumbar actuator (will be replaced in next step)
removeActuator = []
for forceInd in range(gaitModel2D.updForceSet().getSize()):
    if gaitModel2D.updForceSet().get(forceInd).getName() == 'lumbarAct':
        removeActuator.append(forceInd)
gaitModel2D.updForceSet().remove(removeActuator[0])

#Add torque and reserve actuators to model
#Set actuators to add to model (reserves with low forces for muscle actuated side)
actuatorsToAdd = {'hip_flexion_r': 1, 'hip_flexion_l': 300,
                  'knee_angle_r': 1, 'knee_angle_l': 300,
                  'ankle_angle_r': 1, 'ankle_angle_l': 300,
                  'lumbar_extension': 300}

#Add torque actuators to model
gaitModel2D = helper.addTorqueActuators(osimModel = gaitModel2D,
                                        optForces = actuatorsToAdd,
                                        minControl = np.inf * -1,
                                        maxControl = np.inf)

#Finalise model connections
gaitModel2D.finalizeConnections()
gaitModel2D.finalizeFromProperties()
gaitModel2D.printToXML('gaitModel2D_extraMuscles.osim')

#Create a model processor
modelProcessor = osim.ModelProcessor('gaitModel2D_extraMuscles.osim')

#Add relevant model operators
modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
modelProcessor.append(osim.ModOpIgnoreTendonCompliance()) #will switch back on for individual muscles
modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(3.0)) #1.5 typically used in examples
modelProcessor.append(osim.ModOpScaleMaxIsometricForce(3.0))
modelProcessor.append(osim.ModOpTendonComplianceDynamicsModeDGF('implicit')) #Use muscle contractile dynamics

# #Add reserve actuators
# modelProcessor.append(osim.ModOpAddReserves(1))

#Process model to make further edits
muscleModel = modelProcessor.process()

#Update max contraction velocity of muscles
#Turn on tendon compliance for plantarflexors
for muscleInd in range(muscleModel.getMuscles().getSize()):
    #### TODO: how high does this need to be? 30 too low? 50?
    #Contraction velocity
    muscleModel.getMuscles().get(muscleInd).set_max_contraction_velocity(30)
    # #### TODO: check how changing tendon slack length by 10% effects convergence
    # #### Increasing by 10% didn't seem to help...
    # # #### Shortening TSL of all muscles causes big time issues...
    # # muscleModel.getMuscles().get(muscleInd).set_tendon_slack_length(0.95 * muscleModel.getMuscles().get(muscleInd).get_tendon_slack_length())
    #Tendon compliance
    # if 'gastroc' in muscleModel.getMuscles().get(muscleInd).getName() or 'soleus' in muscleModel.getMuscles().get(muscleInd).getName():
    if 'gas' in muscleModel.getMuscles().get(muscleInd).getName() or 'soleus' in muscleModel.getMuscles().get(muscleInd).getName():
        muscleModel.getMuscles().get(muscleInd).set_ignore_tendon_compliance(False)
        
#### TODO: any other muscle edits???
    #### I think adjusting tendon slack lengths to match Lai et al. works better
    #### It seems like there is buckling in the tendons (length < TSL) so it needs to be fixed
    #### More knee flexion needed in model too?
        
#Update friction parameters in contact model
#This helps with smoother GRF predictions (I think...)
#Settings for coefficients
staticFriction = 0.5 #typical model = 0.8
dynamicFriction = 0.5 #typical model = 0.8
viscousFriction = 0.5 #typical model = 0.5
#Update friction coefficients in forceset
for forceInd in range(muscleModel.updForceSet().getSize()):
    #Check for contact geometry force
    if 'contact' in muscleModel.updForceSet().get(forceInd).getName():
        #Update the friction parameters
        osim.SmoothSphereHalfSpaceForce.safeDownCast(muscleModel.updForceSet().get(forceInd)).set_static_friction(staticFriction)
        osim.SmoothSphereHalfSpaceForce.safeDownCast(muscleModel.updForceSet().get(forceInd)).set_dynamic_friction(dynamicFriction)
        osim.SmoothSphereHalfSpaceForce.safeDownCast(muscleModel.updForceSet().get(forceInd)).set_viscous_friction(viscousFriction)

#Finalise connections
muscleModel.finalizeConnections()

#Save model to file
muscleModel.printToXML('muscleDrivenInverse\\muscleModel.osim')

#Create Moco Inverse tool
inverse = osim.MocoInverse()

#Set model processor on tool
inverseModelProcessor = osim.ModelProcessor('muscleDrivenInverse\\muscleModel.osim')
inverse.setModel(inverseModelProcessor)

#Set kinematics in tool
tableProcessor = osim.TableProcessor(osim.MocoTrajectory('torqueDrivenTracking\\torqueDrivenTrackingSolution.sto').exportToStatesTable())
inverse.setKinematics(tableProcessor)
# inverse.setKinematics(osim.TableProcessor('refQ_2D.sto'))

#Set timings
inverse.set_initial_time(tableProcessor.process().getIndependentColumn()[0])
inverse.set_final_time(tableProcessor.process().getIndependentColumn()[-1])
inverse.set_mesh_interval(0.01)

#Set for extra columns
inverse.set_kinematics_allow_extra_columns(True)

#Convert to Moco study and problem to further manipulate aspects of problem
study = inverse.initialize()
problem = study.updProblem()

#Update the weight on the default control effort goal
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('excitation_effort'))
effort.setWeight(0.01)

#Update individual weights in control effort goal to be relative to actual muscle
#and reserve actuator names
#Set appropriate patterns in the weight set
#Muscles
effort.setWeightForControlPattern('/forceset/.*/activation', 0.1)
#Actuators
# effort.setWeightForControlPattern('/forceset/actuator_.*', 1)

#Set appropriate bounds on tendon force
problem.setStateInfoPattern('/forceset/.*/normalized_tendon_force', [0, 2.0], [], [])

#Define the solver and set its options
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_optim_constraint_tolerance(0.001)
solver.set_optim_convergence_tolerance(0.01)
solver.resetProblem(problem)

#Get the initial guess to manipulate
initialGuess = solver.createGuess()

#Set the normalised tendon force states to a reasonable value in the initial guess
muscleModel.initSystem()
for stateInd in range(muscleModel.getStateVariableNames().getSize()):    
    #Get current state name
    stateName = muscleModel.getStateVariableNames().get(stateInd)
    #Check for normalised tendon force
    if 'normalized_tendon_force' in stateName:
        #Set to a reasonable starting value in guess
        initialGuess.setState(stateName,
                              np.linspace(0.2, 0.2, initialGuess.getNumTimes()))
    
#Set updated guess in solver
solver.setGuess(initialGuess)

#Solve inverse problem
#### TODO: something wrong with this...using the updated model crashes MocoInverse...using original gaitModel2D.osim doesn't
    #### Looks like it's the rolling surface constraint that's the problem
        #### Which might not be needed given kinematics are held the same...
#### NOTE: solves quite quickly with original model though --- huge reserves though...
    #### Without these reserves the simulation massively struggles in restoration...it's infeasible...
        #### Worst offenders are the hip and knee flexion...
inverseSolution = study.solve()
# inverseSolution.write('test.sto')


# %% Run a 2D muscle-driven tracking sim of experimental data

# This step is designed to generate a solution for the 2D model that appropriately
# track the sprint kinematics and GRF using muscle actuators. The solution here
# is guided by data from the torque driven simulation (the kinematics in 
# particular)

#Check whether to run simulation
if runMuscleTracking:
    
    #Make directory for muscle driven simulation results
    if not os.path.isdir('muscleDrivenTracking'):
        os.mkdir('muscleDrivenTracking')
        
    #Set-up logger for tracking simulation
    osim.Logger.removeFileSink()
    osim.Logger.addFileSink('muscleDrivenTracking\\muscleDrivenTrackingLog.log')
    
    #Load the model
    osimModel = osim.Model('gaitModel2D.osim')
    # gaitModel2D = osim.Model('..\\..\\ExpData\\Scaling\\scaledModelAdjusted.osim')
    
    #Add in the surface constraints from the dummy model
    #Load in the model
    dummyConstraintModel = osim.Model('..\\..\\ExpData\\Scaling\\gait9dof18musc_Ong_et_al_Moco_withSurfaceConstraints.osim')
    #Add the rolling on surface constraint
    #This can't be accessed through the API so needs to be worked around by
    #stealing it from a dummy model. The dummy model contains 4 rolling surface
    #constraints that we adopt across
    for conInd in range(dummyConstraintModel.updConstraintSet().getSize()):
        osimModel.updConstraintSet().cloneAndAppend(dummyConstraintModel.updConstraintSet().get(conInd))

    #Add a metabolic cost model
    metabolics = osim.Bhargava2004SmoothedMuscleMetabolics()
    metabolics.setName('metabolicCost')
    metabolics.set_use_smoothing(True)
    
    #Add the specific muscles to the metabolics object
    for muscleInd in range(osimModel.getMuscles().getSize()):
        #Get muscle name    
        muscleName = osimModel.getMuscles().get(muscleInd).getName()
        #Add to metabolics object
        metabolics.addMuscle(muscleName, osimModel.getMuscles().get(muscleName))
        
    #Add metabolics to model
    osimModel.addComponent(metabolics)
    
    #Update lumbar actuator optimal force
    osim.CoordinateActuator.safeDownCast(osimModel.updForceSet().get('lumbarAct')).setOptimalForce(300)
    #### If adding more detailed muscles...
    # gaitModel2D = helper.addTorqueActuators(osimModel = osimModel,
    #                                         optForces = {'lumbar_extension': 300},
    #                                         minControl = np.inf * -1,
    #                                         maxControl = np.inf)
    
    #Finalise model connections
    osimModel.finalizeConnections()

    #Create a model processor
    modelProcessor = osim.ModelProcessor(osimModel)

    #Add relevant model operators
    modelProcessor.append(osim.ModOpIgnoreTendonCompliance()) #will switch back on for individual muscles
    modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
    modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(3.0)) #1.5 typically used in examples
    modelProcessor.append(osim.ModOpScaleMaxIsometricForce(2.0))
    # modelProcessor.append(osim.ModOpTendonComplianceDynamicsModeDGF('implicit')) #Use muscle contractile dynamics
    modelProcessor.append(osim.ModOpTendonComplianceDynamicsModeDGF('explicit')) #Seems to help with buckling issue?
    
    #Add reserve actuators
    # modelProcessor.append(osim.ModOpAddReserves(1))
    # modelProcessor.append(osim.ModOpAddReserves(100))
    
    #Process model to make further edits
    muscleModel = modelProcessor.process()

    #Update max contraction velocity of muscles
    #Turn on tendon compliance for plantarflexors
    for muscleInd in range(muscleModel.getMuscles().getSize()):
        #### TODO: how high does this need to be? 30 too low? 50?
        #Contraction velocity
        muscleModel.getMuscles().get(muscleInd).set_max_contraction_velocity(50)
        #### TODO: check how changing tendon slack length by 10% effects convergence
        #### Increasing by 10% didn't seem to help...
        # #### Shortening TSL of all muscles causes big time issues...
        # muscleModel.getMuscles().get(muscleInd).set_tendon_slack_length(0.95 * muscleModel.getMuscles().get(muscleInd).get_tendon_slack_length())
        #Tendon compliance
        if 'gastroc' in muscleModel.getMuscles().get(muscleInd).getName() or 'soleus' in muscleModel.getMuscles().get(muscleInd).getName():
        # if 'gas' in muscleModel.getMuscles().get(muscleInd).getName() or 'soleus' in muscleModel.getMuscles().get(muscleInd).getName():
            muscleModel.getMuscles().get(muscleInd).set_ignore_tendon_compliance(False)
            #### TODO: impact of stiffer tendon to avoid buckling???
            # osim.DeGrooteFregly2016Muscle.safeDownCast(muscleModel.getMuscles().get(muscleInd)).set_tendon_strain_at_one_norm_force(0.35)
            
    #### TODO: any other muscle edits???
        #### I think adjusting tendon slack lengths to match Lai et al. works better
        #### It seems like there is buckling in the tendons (length < TSL) so it needs to be fixed
        #### More knee flexion needed in model too?
            
    #Update friction parameters in contact model
    #This helps with smoother GRF predictions (I think...)
    #Settings for coefficients
    staticFriction = 0.5 #typical model = 0.8
    dynamicFriction = 0.5 #typical model = 0.8
    viscousFriction = 0.5 #typical model = 0.5
    #Update friction coefficients in forceset
    for forceInd in range(muscleModel.updForceSet().getSize()):
        #Check for contact geometry force
        if 'contact' in muscleModel.updForceSet().get(forceInd).getName():
            #Update the friction parameters
            osim.SmoothSphereHalfSpaceForce.safeDownCast(muscleModel.updForceSet().get(forceInd)).set_static_friction(staticFriction)
            osim.SmoothSphereHalfSpaceForce.safeDownCast(muscleModel.updForceSet().get(forceInd)).set_dynamic_friction(dynamicFriction)
            osim.SmoothSphereHalfSpaceForce.safeDownCast(muscleModel.updForceSet().get(forceInd)).set_viscous_friction(viscousFriction)

    #Finalise connections
    muscleModel.finalizeConnections()
    
    #Save model to file
    muscleModel.printToXML('muscleDrivenTracking\\muscleModel.osim')

    #Define the motion tracking problem
    track = osim.MocoTrack()
    track.setName('muscleDrivenTracking')
    
    #Set model
    track.setModel(osim.ModelProcessor(muscleModel))
    
    #Set kinematics
    #This uses the kinematics from the torque driven simulation
    tableProcessor = osim.TableProcessor('refQ_2D.sto')
    # tableProcessor = osim.TableProcessor(osim.MocoTrajectory('torqueDrivenTracking\\torqueDrivenTrackingSolution.sto').exportToStatesTable())

    #### TODO: do we need these variable weights now that we've got 'better' 
    #### kinematics from the torque driven solution? All edited to be balanced
    #### as 1 below

    #Set states reference details
    track.setStatesReference(tableProcessor) # Apply the target data to the tracking problem
    track.set_states_global_tracking_weight(0.1) # Default tracking weight (is changed below)
    track.set_allow_unused_references(True) # Target data can include DoF not in this model
    track.set_track_reference_position_derivatives(True) #Track speed trajectories
    track.set_apply_tracked_states_to_guess(True) # Use target data in initial guess
    
    #Set times
    track.set_initial_time(tableProcessor.process().getIndependentColumn()[0])
    track.set_final_time(tableProcessor.process().getIndependentColumn()[-1])
    
    #Create weight set for state tracking
    stateWeights = osim.MocoWeightSet()
    
    #Create a dictionary that provides the kinematic task weights for function
    taskWeights = {'pelvis_tx': 500, 'pelvis_ty': 250,
                    'pelvis_tilt': 250, 
                    'hip_flexion_r': 100, 'knee_angle_r': 100,
                    'ankle_angle_r': 50,
                    'hip_flexion_l': 100, 'knee_angle_l': 100,
                    'ankle_angle_l': 50,
                    'lumbar_extension': 100}
    # taskWeights = {'pelvis_tx': 1, 'pelvis_ty': 1,
    #                 'pelvis_tilt': 1, 
    #                 'hip_flexion_r': 1, 'knee_angle_r': 1,
    #                 'ankle_angle_r': 1,
    #                 'hip_flexion_l': 1, 'knee_angle_l': 1,
    #                 'ankle_angle_l': 1,
    #                 'lumbar_extension': 1}

    #Set constant weight to scale tracking error speeds by
    w = 0.001

    #Loop through coordinates to apply weights
    for coordInd in range(muscleModel.updCoordinateSet().getSize()):
        
        #Get name and absolute path to coordinate
        coordName = muscleModel.updCoordinateSet().get(coordInd).getName()
        coordPath = muscleModel.updCoordinateSet().get(coordInd).getAbsolutePathString()

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
    
    #Define the Moco study and problem
    study = track.initialize()
    problem = study.updProblem()
    
    ##### CREATE GOALS #####
    
    #Regularization term on MocoTrack problem (minimize squared muscle excitations)
    effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
    effort.setWeight(0.001)
    
    #Create the metabolic cost goal. The total metabolic rate includes activation
    #heat rate, shortening heat rate, mechanical work rate and basal metabolic
    #rate.
    metGoal = osim.MocoOutputGoal('metCost', 0.1)
    metGoal.setOutputPath('/metabolicCost|total_metabolic_rate')
    metGoal.setDivideByDisplacement(True)
    metGoal.setDivideByMass(True)
    
    #Add to problem
    problem.addGoal(metGoal)
    
    #Add contact tracking goal
    #Uses two separate vectors to apply different contact tracking weights    
    
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
                                          osim.Vec3(0,1,0)],
                               # 'weight': [2, 1],
                               'weight': [1, 0.5],
                               }
    
    #Loop through contact tracking settings to create separate weighted goals
    for ii in range(len(contactTrackingSettings['vector'])):
        
        #Create tracking goal
        contactTracking.append(osim.MocoContactTrackingGoal(f'contact{ii}',
                                                            contactTrackingSettings['weight'][ii]))
        
        #Set external loads
        contactTracking[ii].setExternalLoadsFile('refGRF.xml')
        
        #Add the left and right tracking groups to the goal
        contactTracking[ii].addContactGroup(trackLeftGRF)
        contactTracking[ii].addContactGroup(trackRightGRF)
        
        #Set projection vector in problem
        contactTracking[ii].setProjection('vector')
        contactTracking[ii].setProjectionVector(contactTrackingSettings['vector'][ii])
    
        #Add to problem
        problem.addGoal(contactTracking[ii])
        
    #Create a symmetry goal to permit simulating one step
    symmetryGoal = osim.MocoPeriodicityGoal('symmetry')
    
    #Create a dictionary to prescribe state symmetry pairs or isolated states
    symmetryPairs = {
        #Joint coordinate
        'hip_angle': ['/jointset/hip_r/hip_flexion_r/value',
                      '/jointset/hip_l/hip_flexion_l/value'],
        'knee_angle': ['/jointset/knee_r/knee_angle_r/value',
                        '/jointset/knee_l/knee_angle_l/value'],
        'ankle_angle': ['/jointset/ankle_r/ankle_angle_r/value',
                        '/jointset/ankle_l/ankle_angle_l/value'],
        'lumbar_angle': ['/jointset/back/lumbar_extension/value'],
        'pelvis_angle': ['/jointset/ground_pelvis/pelvis_tilt/value'],
        'pelvis_ty': ['/jointset/ground_pelvis/pelvis_ty/value']
        }
    
    #Set symmetric pairings based on dictionary
    for symmetry in list(symmetryPairs.keys()):
        
        #Check for state pair or individual state
        if len(symmetryPairs[symmetry]) == 2:
            
            #Add state pair
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(symmetryPairs[symmetry][0],
                                                                    symmetryPairs[symmetry][1]))
            #Add reverse of state pair
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(symmetryPairs[symmetry][1],
                                                                    symmetryPairs[symmetry][0]))
            
        else:
            
            #Add individual state pairing
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(symmetryPairs[symmetry][0]))

    #Add to problem
    problem.addGoal(symmetryGoal)
    
    ##### SET BOUNDS #####
    
    #Create a dictionary to set kinematic bounds
    #Create this based on some generic values above/below the values in the
    #experimental data    
    
    #Load the experimental kinematics file as a table
    expKinematicsTable = osim.TimeSeriesTable('refQ_2D.sto')
    # expKinematicsTable = osim.TimeSeriesTable('torqueDrivenTracking\\torqueDrivenTrackingSolution.sto')
    
    #Create the bounds dictionary
    kinematicBounds = {
        'pelvis_tilt': [expKinematicsTable.getDependentColumn('/jointset/ground_pelvis/pelvis_tilt/value').to_numpy().min() - np.deg2rad(10),
                        expKinematicsTable.getDependentColumn('/jointset/ground_pelvis/pelvis_tilt/value').to_numpy().max() + np.deg2rad(10)],
        'pelvis_tx': [expKinematicsTable.getDependentColumn('/jointset/ground_pelvis/pelvis_tx/value').to_numpy().min() - 0.2,
                      expKinematicsTable.getDependentColumn('/jointset/ground_pelvis/pelvis_tx/value').to_numpy().max() + 0.2],
        'pelvis_ty': [expKinematicsTable.getDependentColumn('/jointset/ground_pelvis/pelvis_ty/value').to_numpy().min() - 0.1,
                      expKinematicsTable.getDependentColumn('/jointset/ground_pelvis/pelvis_ty/value').to_numpy().max() + 0.1],
        'hip_flexion_l': [expKinematicsTable.getDependentColumn('/jointset/hip_l/hip_flexion_l/value').to_numpy().min() - np.deg2rad(15),
                          expKinematicsTable.getDependentColumn('/jointset/hip_l/hip_flexion_l/value').to_numpy().max() + np.deg2rad(15)],
        'knee_angle_l': [np.deg2rad(-140), np.deg2rad(0)],
        'ankle_angle_l': [expKinematicsTable.getDependentColumn('/jointset/ankle_l/ankle_angle_l/value').to_numpy().min() - np.deg2rad(15),
                          expKinematicsTable.getDependentColumn('/jointset/ankle_l/ankle_angle_l/value').to_numpy().max() + np.deg2rad(15)],
        'hip_flexion_r': [expKinematicsTable.getDependentColumn('/jointset/hip_r/hip_flexion_r/value').to_numpy().min() - np.deg2rad(15),
                          expKinematicsTable.getDependentColumn('/jointset/hip_r/hip_flexion_r/value').to_numpy().max() + np.deg2rad(15)],
        'knee_angle_r': [np.deg2rad(-140), np.deg2rad(0)],
        'ankle_angle_r': [expKinematicsTable.getDependentColumn('/jointset/ankle_r/ankle_angle_r/value').to_numpy().min() - np.deg2rad(15),
                          expKinematicsTable.getDependentColumn('/jointset/ankle_r/ankle_angle_r/value').to_numpy().max() + np.deg2rad(15)],
        'lumbar_extension': [expKinematicsTable.getDependentColumn('/jointset/back/lumbar_extension/value').to_numpy().min() - np.deg2rad(10),
                             5]
        }
    
    #Set the bounds using the dictionary values and experimental data
    for coordInd in range(muscleModel.updCoordinateSet().getSize()):
        
        #First check if coordinate is in kinematic bounds dictionary
        if muscleModel.updCoordinateSet().get(coordInd).getName() in list(kinematicBounds.keys()):
            
            #Get coordinate name for ease of use
            coordName = muscleModel.updCoordinateSet().get(coordInd).getName()

            #Set bounds in problem
            problem.setStateInfo(f'{muscleModel.updCoordinateSet().get(coordInd).getAbsolutePathString()}/value',
                                  #Bounds set to model ranges
                                  [kinematicBounds[coordName][0], kinematicBounds[coordName][1]]
                                  )
    
    #Set reasonable bounds on the tendon forces
    problem.setStateInfoPattern('/forceset/.*/normalized_tendon_force', [0, 2.0], [], [])

    #Define the solver
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())

    #Get the initial guess to manipulate
    initialGuess = solver.getGuess()
    
    #Set the normalised tendon force states to a reasonable value in the initial guess
    muscleModel.initSystem()
    for stateInd in range(muscleModel.getStateVariableNames().getSize()):    
        #Get current state name
        stateName = muscleModel.getStateVariableNames().get(stateInd)
        #Check for normalised tendon force
        if 'normalized_tendon_force' in stateName:
            #Set to a reasonable starting value in guess
            initialGuess.setState(stateName,
                                  np.linspace(0.2, 0.2, initialGuess.getNumTimes()))
            
    #Set updated guess in solver
    solver.resetProblem(problem)
    solver.setGuess(initialGuess)

    #Set solver options
    #Run an initial pass with a coarse mesh to set as an initial guess in a
    #subsequent finer solution
    solver.set_optim_max_iterations(2000) #higher for coarse solution...shouldn't be needed though...
    # solver.set_num_mesh_intervals(10)
    solver.set_minimize_implicit_multibody_accelerations(True)
    solver.set_optim_constraint_tolerance(1e-2) 
    solver.set_optim_convergence_tolerance(1e-2)
    solver.set_minimize_implicit_auxiliary_derivatives(True)
    solver.set_implicit_auxiliary_derivatives_weight(0.00001)

    #Solve!
    #### TODO: add check if coarse solution exists? Or is this covered by the run boolean?
    # coarseSolution = study.solve()
    
    #### Muscle driven solution having trouble with the inf_du parameter holding
    #### it back from converging, although it might be getting close. Need to determine
    #### what exactly is causing that --- tendon lengths? contraction velocities?
        #### Increasing TSL by 10% didn't seem to alleviate the problem, so it is
        #### potentially something else...
        
        #### [info] DeGrooteFregly2016Muscle 'gastroc_r' is buckling (length < tendon_slack_length) at time 0.342667 s.
        #### [info] DeGrooteFregly2016Muscle 'gastroc_r' is exceeding maximum contraction velocity at time 0.342667 s.
            #### this looks to be just before foot strike...
            #### if TSL > length, does TSL need to be reduced???
            
        #### What happens when you remove tendon compliance and set CV at 30?
            #### Still goes through a lot of iterations - seems to get over the
            #### issue of inf_du being inflated though, as it comes down more in
            #### this than other simulations...
                #### It might converge but probably still requires over 1000 iterations,
                #### which means it's a poorly designed problem
                    #### The muscle parameters are clearly the problem here...
            #### Eventually solves in ~1100 iterations
                #### Ankle kinematics are bad though...
                    #### Tendon likely too rigid to dorsiflex...
                    
        #### Implicit vs. explicit tendon dynamics?
            #### Does explicit tendon mode help?
                #### Most discussions on forum suggest implicit tendon with explicit
                #### multibody dynamics are best. Although some discussions verge
                #### on explicit tendon mode only struggling with more complex
                #### models...
                #### Explcit mode doesn't seem to help convergence, same sort
                #### of pattern across iterations appears to persist
                
        #### Passive forces?
            #### Nope...terrible problems here...
            
        #### Compliance for all tendons?
            #### Seems to slow down iterations...
            #### Seems to suffer from same issues with respect to inf_du, but may
            #### be slightly better???
                #### Problem very much wants to converge around ~1000 iterations
                #### This many iterations still therefore feels fairly ill conditioned
                    #### Wants to converge but refuses to...
                    #### It's obviously tendon buckling based on the info outputs 
                    #### when viewing these simulations (happens for most muscles)
                        #### [info] DeGrooteFregly2016Muscle 'hamstrings_r' is buckling (length < tendon_slack_length) at time 0.342667 s.
                        
        #### Stiff tendon (i.e. reduced strain at one norm force) help buckling?
            #### Solver succeeds in 1390 iterations
                #### Has some similar issues to no tendon compliance solution
                #### of ankle not dorsiflexing enough, likely due to stiff tendon
            
        #### I would say tendon slack length and/or optimal fibre length needs 
        #### to change...
        
        #### Muscles like gastroc and hamstrings are operating a long way short
        #### versus long, respectively, from their optimal length --- scale up the
        #### active force width to help this? Or change optimal fibre length?
            ##### Using a 3x scale of active force width seemed to promote convergence...
            ##### Still possibly some minor buckling/contraction velocity problems...
                ##### [info] DeGrooteFregly2016Muscle 'gastroc_r' is buckling (length < tendon_slack_length) at time 0.346 s.
                ##### [info] DeGrooteFregly2016Muscle 'gastroc_r' is exceeding maximum contraction velocity at time 0.346 s.
        
    # #Write coarse solution to file
    # coarseSolution.write('muscleDriven\\muscleDrivenSolution_coarse.sto')
    
    # #Set updated guess in solver for finer grid solution
    # solver.resetProblem(problem)
    # solver.setGuessFile('muscleDriven\\muscleDrivenSolution_coarse.sto')
    
    #Update solver to options to finer mesh
    solver.set_num_mesh_intervals(50) #NOTE: reduced mesh interval for finer muscle driven
        
    #Solve with finer mesh!
    solution = study.solve()

    #Write to file
    # solution.unseal()
    solution.write('muscleDrivenTracking\\Unsealed.sto')
    
    #Visualise
    if visualiseSolution:
        study.visualize(solution)
        
    #Extract predicted GRFs
    
    #Add contact elements to extract from
    #Get right side
    contact_r = osim.StdVectorString()
    for contactLabel in rightFootContacts:
        contact_r.append(contactLabel)
    #Get left side
    contact_l = osim.StdVectorString()
    for contactLabel in leftFootContacts:
        contact_l.append(contactLabel)
            
    #Create forces table
    externalForcesTableFlat = osim.createExternalLoadsTableForGait(osimModel, solution,
                                                                   contact_r, contact_l)
    
    #Write table to file
    osim.STOFileAdapter().write(externalForcesTableFlat,
                                'muscleDrivenTracking\\muscleDrivenTrackingSolution_grf.sto')
    
    #Remove tracked states file
    os.remove('muscleDrivenTracking_tracked_states.sto')
    
    #### Solution seems infeasible from a muscle force perspective whereby large
    #### reserve torques are required and also that the inf_du does not seem to
    #### be achievable either...
   


# %% 2D predictive sim at matched speed (muscle driven)

"""

NOTES:
    
    > Currently work in progress
    > Consider what goals are needed etc.
    > Consider what weights are needed for each goal...
    > Fiber damping on muscles?
    > Should there be a low GRF and kinematic weight?

"""


######## TEST THIS OUT...

#Make directory for torque driven simulation results
#### TODO: these need better labels for names
if not os.path.isdir('muscleDrivenPredictive'):
    os.mkdir('muscleDrivenPredictive')

#Set-up logger for tracking simulation
osim.Logger.removeFileSink()
osim.Logger.addFileSink('muscleDrivenPredictive\\muscleDrivenPredictiveLog.log')

#### TODO: create model in same way or grab from previous problem...
#### Just using existing for now

#Define the predictive study problem
study = osim.MocoStudy()
study.setName('muscleDrivenPredictive')

#Define the problem
problem = study.updProblem()

# #Load the model
# gaitModel2D = osim.Model('gaitModel2D.osim')

# #### TODO: do we need to update slack lengths etc.???

# #Create a model processor
# modelProcessor = osim.ModelProcessor(gaitModel2D)

# #Add relevant model operators
# modelProcessor.append(osim.ModOpIgnoreTendonCompliance()) #will switch back on for individual muscles
# modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF()) #TODO: ignore or keep?
# modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
# modelProcessor.append(osim.ModOpScaleMaxIsometricForce(2.0))
# modelProcessor.append(osim.ModOpTendonComplianceDynamicsModeDGF('implicit')); #Use muscle contractile dynamics
# modelProcessor.append(osim.ModOpFiberDampingDGF(0.01)) #Adjust fiber damping as per Dembia et al. Moco paper

# #Process model to make further edits
# muscleModel = modelProcessor.process()

# #Update max contraction velocity of muscles
# #Turn on tendon compliance for plantarflexors
# for muscleInd in range(muscleModel.getMuscles().getSize()):
#     #Contraction velocity
#     muscleModel.getMuscles().get(muscleInd).set_max_contraction_velocity(30)
#     #Tendon compliance
#     if 'gastroc' in muscleModel.getMuscles().get(muscleInd).getName() or 'soleus' in muscleModel.getMuscles().get(muscleInd).getName():
#         muscleModel.getMuscles().get(muscleInd).set_ignore_tendon_compliance(False)
    
# #### TODO: any other muscle edits???
#     #### I think adjusting tendon slack lengths to match Lai et al. works better?
#         #### Do this at scaling - alongside opening up more knee flexion to reflect Lai model

# #Update friction parameters in contact model
# #This helps with smoother GRF predictions (I think...)
# #Settings for coefficients
# staticFriction = 0.5 #typical model = 0.8
# dynamicFriction = 0.5 #typical model = 0.8
# viscousFriction = 0.5 #typical model = 0.5
# #Update friction coefficients in forceset
# for forceInd in range(muscleModel.updForceSet().getSize()):
#     #Check for contact geometry force
#     if 'contact' in muscleModel.updForceSet().get(forceInd).getName():
#         #Update the friction parameters
#         osim.SmoothSphereHalfSpaceForce.safeDownCast(muscleModel.updForceSet().get(forceInd)).set_static_friction(staticFriction)
#         osim.SmoothSphereHalfSpaceForce.safeDownCast(muscleModel.updForceSet().get(forceInd)).set_dynamic_friction(dynamicFriction)
#         osim.SmoothSphereHalfSpaceForce.safeDownCast(muscleModel.updForceSet().get(forceInd)).set_viscous_friction(viscousFriction)

# #Finalise connections
# muscleModel.finalizeConnections()

# #Save model to file
# muscleModel.printToXML('muscleDrivenPredictive\\muscleModel.osim')
# muscleModel.initSystem()

#Set model in problem
problem.setModelProcessor(osim.ModelProcessor(muscleModel))

#Create a speed goal to match the experimental data sprint speed

#Set an average speed goal based on sprinting data
#Get the average sprint speed based on the pelvis translation from experimental data
#Get distance travelled
distTravelled = osim.TimeSeriesTable('refQ_2D.sto').getDependentColumn('/jointset/ground_pelvis/pelvis_tx/value').to_numpy()[-1] - \
    osim.TimeSeriesTable('refQ_2D.sto').getDependentColumn('/jointset/ground_pelvis/pelvis_tx/value').to_numpy()[0]
#Get time taken
timeTaken = osim.TimeSeriesTable('refQ_2D.sto').getIndependentColumn()[-1] - osim.TimeSeriesTable('refQ_2D.sto').getIndependentColumn()[0]
#Calculate average speed
sprintSpeed = distTravelled / timeTaken
#Create the speed goal
speedGoal = osim.MocoAverageSpeedGoal('speed')
speedGoal.setWeight(1)
#Set the speed
speedGoal.set_desired_average_speed(sprintSpeed)
#Add to the problem
problem.addGoal(speedGoal)

#Create a symmetry goal to permit simulating one step
symmetryGoal = osim.MocoPeriodicityGoal('symmetry')
symmetryGoal.setWeight(1)
#Create a dictionary to prescribe state symmetry pairs or isolated states
symmetryPairs = {
    #Joint coordinate
    'hip_angle': ['/jointset/hip_r/hip_flexion_r/value',
                  '/jointset/hip_l/hip_flexion_l/value'],
    'knee_angle': ['/jointset/knee_r/knee_angle_r/value',
                    '/jointset/knee_l/knee_angle_l/value'],
    'ankle_angle': ['/jointset/ankle_r/ankle_angle_r/value',
                    '/jointset/ankle_l/ankle_angle_l/value'],
    'lumbar_angle': ['/jointset/back/lumbar_extension/value'],
    'pelvis_angle': ['/jointset/ground_pelvis/pelvis_tilt/value'],
    'pelvis_ty': ['/jointset/ground_pelvis/pelvis_ty/value'],
    #Muscle activations
    'hamstrings_act': ['/forceset/hamstrings_r/activation',
                        '/forceset/hamstrings_l/activation'],
    'bifemsh_act': ['/forceset/bifemsh_r/activation',
                    '/forceset/bifemsh_l/activation'],
    'glut_max_act': ['/forceset/glut_max_r/activation',
                      '/forceset/glut_max_l/activation'],
    'iliopsoas_act': ['/forceset/iliopsoas_r/activation',
                      '/forceset/iliopsoas_l/activation'],
    'rect_fem_act': ['/forceset/rect_fem_r/activation',
                      '/forceset/rect_fem_l/activation'],
    'vasti_act': ['/forceset/vasti_r/activation',
                  '/forceset/vasti_l/activation'],
    'gastroc_act': ['/forceset/gastroc_r/activation',
                    '/forceset/gastroc_l/activation'],
    'soleus_act': ['/forceset/soleus_r/activation',
                    '/forceset/soleus_l/activation'],
    'tib_ant_act': ['/forceset/tib_ant_r/activation',
                    '/forceset/tib_ant_l/activation'],
    #Tendon forces
    'gastroc_tendon': ['/forceset/gastroc_r/normalized_tendon_force',
                        '/forceset/gastroc_l/normalized_tendon_force'],
    'soleus_tendon': ['/forceset/soleus_r/normalized_tendon_force',
                      '/forceset/soleus_l/normalized_tendon_force'],
    }
#Set symmetric pairings based on dictionary
for symmetry in list(symmetryPairs.keys()):
    #Check for state pair or individual state
    if len(symmetryPairs[symmetry]) == 2:
        #Add state pair
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(symmetryPairs[symmetry][0],
                                                                symmetryPairs[symmetry][1]))
        #Add reverse of state pair
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(symmetryPairs[symmetry][1],
                                                                symmetryPairs[symmetry][0]))
        #Check whether joint coordinate and hence we add a joint speed pair too
        if '/value' in symmetryPairs[symmetry][0]:
            #Add the matched speed state pair
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(symmetryPairs[symmetry][0].replace('/value','/speed'),
                                                                   symmetryPairs[symmetry][1].replace('/value','/speed')))
            #Add the reverse of the matched state pair
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(symmetryPairs[symmetry][1].replace('/value','/speed'),
                                                                   symmetryPairs[symmetry][0].replace('/value','/speed')))
    else:
        #Add individual state
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(symmetryPairs[symmetry][0]))
        #Check whether joint coordinate and hence we add a joint speed state too
        if '/value' in symmetryPairs[symmetry][0]:
            #Add the matched speed state pair
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(symmetryPairs[symmetry][0].replace('/value','/speed')))
#Add the lumbar actuator control pair
symmetryGoal.addControlPair(osim.MocoPeriodicityGoalPair('/forceset/lumbarAct'))
#Add to problem
problem.addGoal(symmetryGoal)

#Add a control effort goal
### TODO: note that this is slightly different than earlier control goal...
effortGoal = osim.MocoControlGoal('effort', 0.001)
effortGoal.setExponent(3)
# effortGoal.setDivideByDisplacement(True)
effortGoal.setWeightForControl('/forceset/lumbarAct', 0.1) #Set to not highly penalise the lumbar actuator in control goal
problem.addGoal(effortGoal)

# #Set an upright torso goal
# #Grab the experimental data to fill with zeros that will represent the upright posture
# torsoRef = osim.MocoTrajectory('muscleDriven\\muscleDrivenSolution.sto').exportToStatesTable()
# #Loop through coordinates and set all to zero
# for coordinateInd in range(muscleModel.updCoordinateSet().getSize()):
#     #Set absolute path string to coordinate value
#     coordinateStr = muscleModel.updCoordinateSet().get(coordinateInd).getAbsolutePathString()+'/value'
#     #Set to zero in table (angles and speeds)
#     torsoRef.getDependentColumn(coordinateStr).setTo(0)
#     torsoRef.getDependentColumn(coordinateStr.replace('/value','/speed')).setTo(0)
# #Write torso reference coordinates to file
# osim.STOFileAdapter().write(torsoRef, 'muscleDrivenPredictive\\torsoRefQ.sto')
# #Create the upright torso goal
# torsoGoal = osim.MocoOrientationTrackingGoal('torso', 1)
# #Set the table processor reference in the goal
# torsoGoal.setStatesReference(osim.TableProcessor('muscleDrivenPredictive\\torsoRefQ.sto'))
# #Set the path to the torso frame
# bodyPaths = osim.StdVectorString()
# bodyPaths.append('/bodyset/torso')
# torsoGoal.setFramePaths(bodyPaths)
# #Add to problem
# problem.addGoal(torsoGoal)

#Add a low kinematic tracking weight to ensure running motion is replicated
coordinateGoal = osim.MocoStateTrackingGoal('kinematics', 0.1) #### TODO: consider weight
#Set kinematics reference
coordinateGoal.setReference(osim.TableProcessor(osim.MocoTrajectory('muscleDriven\\muscleDrivenSolution_2000iterUnsealed.sto').exportToStatesTable()))
#Set pattern to coordinate values
coordinateGoal.setPattern('/jointset/.*/value')
#Ignore unused references
coordinateGoal.setAllowUnusedReferences(True)
#Add to problem
problem.addGoal(coordinateGoal)

#Create the metabolic cost goal. The total metabolic rate includes activation
#heat rate, shortening heat rate, mechanical work rate and basal metabolic
#rate.
metGoal = osim.MocoOutputGoal('metCost', 1)
metGoal.setOutputPath('/metabolicCost|total_metabolic_rate')
metGoal.setDivideByDisplacement(True)
metGoal.setDivideByMass(True)

#Add to problem
problem.addGoal(metGoal)

#Add a low weight contact tracking goal
contactTracking = osim.MocoContactTrackingGoal('contact', 0.1)
#Set right foot tracking group for time period (i.e. relevant force name)
forcesRightFoot = osim.StdVectorString()
forcesRightFoot.append('/forceset/contactHeel_r')
forcesRightFoot.append('/forceset/contactMidfoot_r')
forcesRightFoot.append('/forceset/contactToe_r')
trackRightGRF = osim.MocoContactTrackingGoalGroup(forcesRightFoot, 'RightGRF1')
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r')
#Set left foot tracking group for time period (i.e. relevant force name)
forcesLeftFoot = osim.StdVectorString()
forcesLeftFoot.append('/forceset/contactHeel_l')
forcesLeftFoot.append('/forceset/contactMidfoot_l')
forcesLeftFoot.append('/forceset/contactToe_l')
trackLeftGRF = osim.MocoContactTrackingGoalGroup(forcesLeftFoot, 'LeftGRF1')
trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l')
#Set external loads
contactTracking.setExternalLoadsFile('refGRF.xml')
#Add the left and right tracking groups to the goal
contactTracking.addContactGroup(trackLeftGRF)
contactTracking.addContactGroup(trackRightGRF)
#Set projection plane to only track the XY forces (i.e. 2D simulation)
contactTracking.setProjection('plane')
contactTracking.setProjectionVector(osim.Vec3(0,0,1))

#Add to problem
problem.addGoal(contactTracking)

#Add joint coordinate bounds
#NOTE: currently reflective of allowable ranges in model --- could be more reasonable...
#NOTE: some random initial bounds provided
#Create dictionary with bounds on joint coordinate states
#Include an initial bounds state based off muscle driven solution too
#Get muscle driven solution as time series table
muscleDrivenSolutionTable = osim.TimeSeriesTable('muscleDriven\\muscleDrivenSolution_2000iterUnsealed.sto')
#Create a list of coordinate bounds to work through
coordinatesToBound = ['/jointset/ground_pelvis/pelvis_tilt/value',
                      '/jointset/ground_pelvis/pelvis_tx/value',
                      '/jointset/ground_pelvis/pelvis_ty/value',
                      '/jointset/hip_r/hip_flexion_r/value',
                      '/jointset/hip_l/hip_flexion_l/value',
                      '/jointset/knee_r/knee_angle_r/value',
                      '/jointset/knee_l/knee_angle_l/value',
                      '/jointset/ankle_r/ankle_angle_r/value',
                      '/jointset/ankle_l/ankle_angle_l/value',
                      '/jointset/back/lumbar_extension/value']
#Loop through coordinates and set them to be within reasonable max and min, and
#initial values on the basis of the solution
#Currently set to within +/- 10% of values from muscle driven solution
for coordinate in coordinatesToBound:
    problem.setStateInfo(coordinate,
                         [muscleDrivenSolutionTable.getDependentColumn(coordinate).to_numpy().min() - (np.abs(muscleDrivenSolutionTable.getDependentColumn(coordinate).to_numpy().min()) * 0.1),
                          muscleDrivenSolutionTable.getDependentColumn(coordinate).to_numpy().max() + (np.abs(muscleDrivenSolutionTable.getDependentColumn(coordinate).to_numpy().max()) * 0.1)],
                         [muscleDrivenSolutionTable.getDependentColumn(coordinate).to_numpy()[0] - (np.abs(muscleDrivenSolutionTable.getDependentColumn(coordinate).to_numpy()[0]) * 0.1),
                          muscleDrivenSolutionTable.getDependentColumn(coordinate).to_numpy()[0] + (np.abs(muscleDrivenSolutionTable.getDependentColumn(coordinate).to_numpy()[0]) * 0.1)]
                         )
    
#Set reasonable bounds on the tendon forces
problem.setStateInfoPattern('/forceset/.*/normalized_tendon_force', [0, 1.8], [], [])

#Set time bounds to allow for subtle fluctations in finish times with minor changes in speed
problem.setTimeBounds(osim.TimeSeriesTable('refQ_2D.sto').getIndependentColumn()[0],
                      [osim.TimeSeriesTable('refQ_2D.sto').getIndependentColumn()[-1] - 0.05,
                       osim.TimeSeriesTable('refQ_2D.sto').getIndependentColumn()[-1] + 0.05])

#Define the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())

#Set solver options
solver.set_optim_max_iterations(1000)
solver.set_num_mesh_intervals(10) ###TODO: set appropriately
solver.set_minimize_implicit_multibody_accelerations(True)
solver.set_optim_constraint_tolerance(1e-2) 
solver.set_optim_convergence_tolerance(1e-2)
solver.set_minimize_implicit_auxiliary_derivatives(True)
solver.set_implicit_auxiliary_derivatives_weight(0.00001)

#Set guess as tracking solution
solver.resetProblem(problem)
solver.setGuessFile('muscleDriven\\muscleDrivenSolution_2000iterUnsealed.sto')

#Solve!
solution = study.solve()

#### Solution converged quite quickly...
    #### But it's not correct...motion looks a little weird and trunk doesn't move...
    #### Need to figure out better goals for predictive sim...
#### Using updated code the solution seems to struggle with convergence
    #### Probably need to better guide the predictive problem?
    #### There's not much different with the problem exact for the periodicity
    #### constraints --- which could be causing the problems? Maybe these constraints
    #### need to come earlier with respect to the tracking simulations so that 
    #### they are embedded in earlier?

#Write to file
solution.write('muscleDrivenPredictive\\muscleDrivenPredictiveSolution.sto')

#Visualise
if visualiseSolution:
    study.visualize(solution)

#Extract predicted GRFs

#Add contact elements to extract from
contact_r = osim.StdVectorString()
contact_l = osim.StdVectorString()
contact_r.append('/forceset/contactHeel_r')
contact_r.append('/forceset/contactMidfoot_r')
contact_r.append('/forceset/contactToe_r')
contact_l.append('/forceset/contactHeel_l')
contact_l.append('/forceset/contactMidfoot_l')
contact_l.append('/forceset/contactToe_l')

#Create forces table
externalForcesTableFlat = osim.createExternalLoadsTableForGait(muscleModel, solution,
                                                               contact_r, contact_l)

#Write table to file
grfStoFile = osim.STOFileAdapter()
grfStoFile.write(externalForcesTableFlat, 'muscleDrivenPredictive\\muscleDrivenPredictiveSolution_grf.sto')


# %% NOTE: Some older stuff below here which might still be useful...

# %% 2D inverse sim of torque drive solution (muscle driven)


# This step is designed to generate relevant muscle controls for the 2D model that
# appropriately track the kinematics from the torque driven simulation. Effectively
# this will serve as a guess to the predictive muscle-driven simulation

#Check whether to run simulation
if runMuscleInverse:
    
    #Make directory for torque driven simulation results
    if not os.path.isdir('muscleInverse'):
        os.mkdir('muscleInverse')
    
    #Set-up logger for tracking simulation
    osim.Logger.removeFileSink()
    osim.Logger.addFileSink('muscleInverse\\muscleInverseLog.log')
    
    #Construct the inverse tool
    inverse = osim.MocoInverse()
    
    #Construct a ModelProcessor for use in the tool
    modelProcessor = osim.ModelProcessor('gaitModel2D.osim')
    
    #Add relevant model operators
    modelProcessor.append(osim.ModOpIgnoreTendonCompliance()) #will switch back on for individual muscles
    modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF()) #TODO: ignore or keep?
    modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
    modelProcessor.append(osim.ModOpScaleMaxIsometricForce(2.0)) ### TODO: fix back...
    modelProcessor.append(osim.ModOpAddReserves(1.0)) ### TODO: needed???
    modelProcessor.append(osim.ModOpAddExternalLoads('refGRF.xml')) ### TODO: this is 3D though...need estimates...spheres or exp forces?
    
    #Process model to make further edits
    inverseModel = modelProcessor.process()
    
    #Update max contraction velocity of muscles
    #Turn on tendon compliance for plantarflexors
    for muscleInd in range(inverseModel.getMuscles().getSize()):
        #Contraction velocity
        inverseModel.getMuscles().get(muscleInd).set_max_contraction_velocity(30)
        # #Tendon compliance
        # if 'gastroc' in inverseModel.getMuscles().get(muscleInd).getName() or 'soleus' in inverseModel.getMuscles().get(muscleInd).getName():
        #     inverseModel.getMuscles().get(muscleInd).set_ignore_tendon_compliance(False)
        
    #### TODO: any other muscle edits???
    
    #Update contact spheres in force set
    for forceInd in range(inverseModel.updForceSet().getSize()):
        #Check for contact geometry force
        if 'contact' in inverseModel.updForceSet().get(forceInd).getName():
            #Turn off force
            inverseModel.updForceSet().get(forceInd).set_appliesForce(False)

    #Finalise connections
    inverseModel.finalizeConnections()
    
    #Save model to file
    inverseModel.printToXML('muscleInverse\\inverseModel.osim')
    
    #Create a new processor with this updated model and set in tool
    inverseModelProcessor = osim.ModelProcessor(inverseModel)
    inverse.setModel(inverseModelProcessor)
    
    #Get the kinematic coordinates from the torque driven solution
    #Convert the solution to a states table storage
    torqueSolution = osim.MocoTrajectory('torqueDriven\\torqueDrivenSolution.sto')
    statesTable = torqueSolution.exportToStatesTable()
    osim.STOFileAdapter().write(statesTable, 'muscleInverse\\coordinates.sto')
    
    #Set the states to a table processor in the inverse problem
    inverse.setKinematics(osim.TableProcessor('muscleInverse\\coordinates.sto'))
    inverse.set_kinematics_allow_extra_columns(True)
    
    #Set the times and mesh intervals
    inverse.set_initial_time(osim.Storage('torqueDriven\\torqueDrivenSolution.sto').getFirstTime())
    inverse.set_final_time(osim.Storage('torqueDriven\\torqueDrivenSolution.sto').getLastTime())
    
    #Calculate mesh interval to match torque driven solution
    meshInterval = np.mean(np.diff(osim.TimeSeriesTable('torqueDriven\\torqueDrivenSolution.sto').getIndependentColumn()))
    
    #Set mesh interval
    inverse.set_mesh_interval(meshInterval) #### TODO: should this match torque solution?
    
    ####TODO: add tendon compliance in but set initial guess to have appropriate values
    ####TODO: consider tendon mode if including compliance
    ####TODO: consider that contact spheres are determining GRFs here...
    
    #Solve!
    solution = inverse.solve()
    
    #Write solution to file
    solution.getMocoSolution().write('muscleInverse\\muscleInverseSolution.sto')
    
    #### NOTES: Solution is shark-finning meaninng that mesh interval is probably
    #### too large. Some muscles saturating meaning that they need to be stronger
    #### or that passive forces/tendon compliance might be necessary...
    #### Some quite high reserves indicative of this too. Unsure what's happening
    #### with the contact spheres here too --- might be more appropriate to run
    #### this with external loads and those forces switched off???
        #### Still some similar problems with a solution with a smaller mesh interval
        #### and stronger muscles --- it's the plantarflexors that are saturating
        #### so tendon compliance or passive forces might help
            #### Solution struggles big time with passive forces
        #### Using external forces probably better --- need to create a 2D version
        #### from torque driven solution though
        
    #### Cranking muscle strength right up reduces saturation and minimises reserve
    #### use outside of residuals
        #### Obviously some issues with force production still at 2x strength
        #### Maybe due to muscle lengths? Could alter active force width more?
            #### Hamstrings still some high activation early could be indicative of this
        #### Tendon forces might also help here
    

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
    #### TODO: I think you can get around this by writing the solution yourself...
    mocoTraj = osim.MocoTrajectory('sprintTracking_torqueDriven_solution.sto')
    statesTab = mocoTraj.exportToStatesTable()
    tableProcessor = osim.TableProcessor(statesTab)    
    # tableProcessor = osim.TableProcessor('sprintTracking_torqueDriven_solution.sto')
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
    
    #Print solution to file (doesn't seem to do this automatically anymore?)
    solution.write('sprintTracking_muscleDriven_solution.sto')
    
    #Option to visualise
    if visualiseSolution:
        study.visualize(solution)
    
    #Create a full gait cycle trajectory from the periodic solution.
    addPatterns = [".*pelvis_tx/value"]
    fullTraj = osim.createPeriodicTrajectory(solution, addPatterns)
    fullTraj.write('sprintTracking_muscleDriven_solution_fullTrajectory.sto')
    
    #Compute ground reaction forces generated by contact sphere from the
    #original and full gait cycle trajectory
    #Half gait cycle
    externalLoads = osim.createExternalLoadsTableForGait(trackModel,
                                                         solution,
                                                         forceNames_r,
                                                         forceNames_l)
    osim.STOFileAdapter.write(externalLoads,'trackedGRF_muscleDriven_2D.mot')
    #Full gait cycle
    externalLoads = osim.createExternalLoadsTableForGait(trackModel,
                                                         fullTraj,
                                                         forceNames_r,
                                                         forceNames_l)
    osim.STOFileAdapter.write(externalLoads,'trackedGRF_fullTrajectory_muscleDriven_2D.mot')

# %% 2D predictive at same speed (muscle driven)

# This section attempts to generate a predictive simulation of a half gait cycle
# at the same speed as tracking. Theoretically this should produce some similarities
# to the tracking sim as a pseudo-validation of the predictive approach
#
# Note this has been adapted to a semi-predictive simulation as very low tracking
# weights for kinematics and GRF data

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
    #Additional low tracking weights for use
    motionTrackingWeight = 1e-3
    grfTrackingWeight = 1e-3
    
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
    
    #Finalise connections
    studyModel.finalizeConnections()
    
    #Get the model as a processor object
    studyModelProcessor = osim.ModelProcessor(studyModel)
    
    #Set model to use implicit tendon compliance if appropriate
    if implicitTendonCompliance:
        studyModelProcessor.append(osim.ModOpUseImplicitTendonComplianceDynamicsDGF())
     
    #Construct the study object and set basic parameters
    study = osim.MocoStudy()
    # study.setName('sprintPrediction_matchedSpeed_muscleDriven')
    study.setName('sprintPredictionLowTracking_matchedSpeed_muscleDriven')
    
    #Update problem
    problem = study.updProblem()
    
    #Set model
    studyModel = studyModelProcessor.process()
    studyModel.initSystem()
    problem.setModelAsCopy(studyModel)
    
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
    
    #Set low weight tracking goals
    
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
    contactGoal = osim.MocoContactTrackingGoal('contactGoal', grfTrackingWeight)
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
    
    #Create kinematic tracking goal
    kinematicGoal = osim.MocoStateTrackingGoal('kinematicTracking', motionTrackingWeight)
    #Set kinematic data and parameters
    tableProcessor = osim.TableProcessor('sprintTracking_muscleDriven_solution.sto')
    tableProcessor.append(osim.TabOpUseAbsoluteStateNames())
    kinematicGoal.setReference(tableProcessor)
    #Set unused state references to be allowed in case of file errors
    kinematicGoal.setAllowUnusedReferences(True)
    #Set pattern to track kinematics
    kinematicGoal.setPattern('/jointset/.*/value')
    #Add kineamtic tracking goal
    problem.addGoal(kinematicGoal)
        
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
    
    #Print solution to file (doesn't seem to do this automatically anymore?)
    solution.write('sprintPredictionLowTracking_matchedSpeed_muscleDriven_solution.sto')
    
    #Option to visualise
    if visualiseSolution:
        study.visualize(solution)
    
    #Create a full gait cycle trajectory from the periodic solution.
    addPatterns = [".*pelvis_tx/value"]
    fullTraj = osim.createPeriodicTrajectory(solution, addPatterns)
    # fullTraj.write('sprintPrediction_matchedSpeed_muscleDriven_fullTrajectory.sto')
    fullTraj.write('sprintPredictionLowTracking_matchedSpeed_muscleDriven_fullTrajectory.sto')
    
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
        
    #Compute ground reaction forces generated by contact sphere from the
    #original and full gait cycle trajectory
    #Half gait cycle
    externalLoads = osim.createExternalLoadsTableForGait(studyModel,
                                                         solution,
                                                         forceNames_r,
                                                         forceNames_l)
    osim.STOFileAdapter.write(externalLoads,'predictedLowTrackingGRF_matchedSpeed_muscleDriven_2D.mot')
    #Full gait cycle
    externalLoads = osim.createExternalLoadsTableForGait(studyModel,
                                                         fullTraj,
                                                         forceNames_r,
                                                         forceNames_l)
    osim.STOFileAdapter.write(externalLoads,'predictedLowTrackingGRF_matchedSpeed_fullTrajectory_muscleDriven_2D.mot')
    
# %% Compare simulations

# This section loads in the simulation data and compares kinematics, muscle function
# (where appropriate), and the GRF data

# %% Load data

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
df_musclePrediction = readSTO('sprintPredictionLowTracking_matchedSpeed_muscleDriven_solution.sto')
df_musclePredictionGRF = readSTO('predictedLowTrackingGRF_matchedSpeed_muscleDriven_2D.mot')



#Visualisations...

#Set the dictionary colour palette
colourDict = {'expData': '#000000',
              'torqueTracking': '#4885ed',
              'muscleTracking': '#db3236',
              'musclePrediction': '#f4c20d'}

#Get the initial and final time from the experimental kinematics
expInitialTime = osim.Storage('refQ_2D.sto').getFirstTime()
expFinalTime = osim.Storage('refQ_2D.sto').getLastTime()

# %% Compare GRFs

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
# ax[0].plot(np.linspace(0,100,101), vGRF_norm,
#            linewidth = 2, color = colourDict['expData'])
# ax[1].plot(np.linspace(0,100,101), apGRF_norm,
#            linewidth = 2, color = colourDict['expData'])
ax[0].plot(timeVals, vGRF, label = 'Exp. Data',
            linewidth = 2, color = colourDict['expData'])
ax[1].plot(timeVals, apGRF, label = 'Exp. Data',
            linewidth = 2, color = colourDict['expData'])

# #Identify the GRF section to plot based on exp. data
# #(i.e. contact phase)
# #Simple way to find the first zero in the exp. data
# startTime = timeVals[0]
# endTime = timeVals[np.where(vGRF == 0)[0][0]]

#Identify the GRF section to plot based on exp. data
#(i.e. half gait cycle)
#Simple way to find the first zero in the exp. data
startTime = timeVals[0]
endTime = timeVals[-1]

#Torque tracking data (all data is relevant)
vGRF = df_torqueTrackingGRF['ground_force_r_vy'].to_numpy().flatten()
apGRF = df_torqueTrackingGRF['ground_force_r_vx'].to_numpy().flatten()
timeVals = df_torqueTrackingGRF['time'].to_numpy().flatten()
#Normalise the data to 101-points
newTime = np.linspace(timeVals[0],timeVals[-1],101).flatten()
vGRF_norm = np.interp(newTime,timeVals,vGRF)
apGRF_norm = np.interp(newTime,timeVals,apGRF)
#Plot the data
# ax[0].plot(np.linspace(0,100,101), vGRF_norm,
#            linewidth = 2, color = colourDict['torqueTracking'])
# ax[1].plot(np.linspace(0,100,101), apGRF_norm,
#            linewidth = 2, color = colourDict['torqueTracking'])
ax[0].plot(timeVals, vGRF, label = 'Torque Driven Tracking',
           linewidth = 2, color = colourDict['torqueTracking'])
ax[1].plot(timeVals, apGRF, label = 'Torque Driven Tracking',
           linewidth = 2, color = colourDict['torqueTracking'])

#Muscle tracking data (all data is relevant)
vGRF = df_muscleTrackingGRF['ground_force_r_vy'].to_numpy().flatten()
apGRF = df_muscleTrackingGRF['ground_force_r_vx'].to_numpy().flatten()
timeVals = df_muscleTrackingGRF['time'].to_numpy().flatten()
#Normalise the data to 101-points
newTime = np.linspace(timeVals[0],timeVals[-1],101).flatten()
vGRF_norm = np.interp(newTime,timeVals,vGRF)
apGRF_norm = np.interp(newTime,timeVals,apGRF)
#Plot the data
# ax[0].plot(np.linspace(0,100,101), vGRF_norm,
#            linewidth = 2, color = colourDict['muscleTracking'])
# ax[1].plot(np.linspace(0,100,101), apGRF_norm,
#            linewidth = 2, color = colourDict['muscleTracking'])
ax[0].plot(timeVals, vGRF, label = 'Muscle Driven Tracking',
           linewidth = 2, color = colourDict['muscleTracking'])
ax[1].plot(timeVals, apGRF, label = 'Muscle Driven Tracking',
           linewidth = 2, color = colourDict['muscleTracking'])

#Muscle prediction data (all data is relevant)

#### TODO: Predictive data needs smoothing???
#### With the low tracking weight on predictive sim seems OK...

vGRF = df_musclePredictionGRF['ground_force_r_vy'].to_numpy().flatten()
apGRF = df_musclePredictionGRF['ground_force_r_vx'].to_numpy().flatten()
timeVals = df_muscleTrackingGRF['time'].to_numpy().flatten()

#Smooth data with 50Hz low-pass
#Define filter for forces data (50N Low-pass 4th Order Butterworth)
fs = 1 / np.mean(np.diff(timeVals))
nyq = 0.5 * fs
normCutoff = 50 / nyq
b, a = butter(4, normCutoff, btype = 'low', analog = False)
#Apply lowpass filter to force data (50N cut-off)
vGRF_f = lfilter(b,a,vGRF)
apGRF_f = lfilter(b,a,apGRF)

#Normalise the data to 101-points
newTime = np.linspace(timeVals[0],timeVals[-1],101).flatten()
vGRF_norm = np.interp(newTime,timeVals,vGRF_f)
apGRF_norm = np.interp(newTime,timeVals,apGRF_f)
#Plot the data
# ax[0].plot(np.linspace(0,100,101), vGRF_norm,
#            linewidth = 2, color = colourDict['musclePrediction'])
# ax[1].plot(np.linspace(0,100,101), apGRF_norm,
#            linewidth = 2, color = colourDict['musclePrediction'])
ax[0].plot(timeVals, vGRF, label = 'Muscle Driven Prediction',
           linewidth = 2, color = colourDict['musclePrediction'])
ax[1].plot(timeVals, apGRF, label = 'Muscle Driven Prediction',
           linewidth = 2, color = colourDict['musclePrediction'])

# #Add legend outside of second axes
# ax[1].legend(bbox_to_anchor = (1.85,0.4))

#Set each axes to time limits
ax[0].set_xlim([startTime,endTime])
ax[1].set_xlim([startTime,endTime])

#Tight layout
plt.tight_layout()

##### TODO: add labels etc.

# %% Compare kinematics

### TODO: set axes programatically for different parts rather than manually
### Use a dictionary to set the joint angle name to an axes subplot value

### TODO: currently just plotting one angle

### TODO: expand figure to include pelvis kinematics --- careful with rad2deg
### for pelvis translations

#Set kinematic variables to plot
kinVars = ['/jointset/ground_pelvis/pelvis_tilt/value',
           '/jointset/ground_pelvis/pelvis_tx/value',
           '/jointset/ground_pelvis/pelvis_ty/value',
           '/jointset/back/lumbar_extension/value',
           '/jointset/hip_r/hip_flexion_r/value',
           '/jointset/knee_r/knee_angle_r/value',
           '/jointset/ankle_r/ankle_angle_r/value',
           '/jointset/hip_l/hip_flexion_l/value',
           '/jointset/knee_l/knee_angle_l/value',
           '/jointset/ankle_l/ankle_angle_l/value']

#Set title labels
kinTitles = ['Pelvis Tilt',
             'Pelvis Tx',
             'Pelvis Ty',
             'Lumbar Ext.',
             'R. Hip Angle',
             'R. Knee Angle',
             'R. Ankle Angle',
             'L. Hip Angle',
             'L. Knee Angle',
             'L. Ankle Angle']

#Set the axes to plot the variables
whichAx = [[0,0], [0,1], [0,2],
           [1,1],
           [2,0],[2,1],[2,2],
           [3,0],[3,1],[3,2]]

#Set figure and subplots
fig, ax = plt.subplots(figsize=(9, 12), nrows = 4, ncols = 3)

#Loop through variables and plot different data
for kk in range(len(kinVars)):

    #Experimental data
    #Extract the relevant section of data
    angVals = df_expKinematics.loc[(df_expKinematics['time'] >= expInitialTime) &
                                   (df_expKinematics['time'] <= expFinalTime),
                                   [kinVars[kk]]].to_numpy().flatten()
    timeVals = df_expKinematics.loc[(df_expKinematics['time'] >= expInitialTime) &
                                    (df_expKinematics['time'] <= expFinalTime),
                                    ['time']].to_numpy().flatten()
    #Convert angular values if radians
    if 'pelvis_tx' not in kinVars[kk] and 'pelvis_ty' not in kinVars[kk]:
        angVals = np.rad2deg(angVals)
    #Plot data
    ax[whichAx[kk][0],whichAx[kk][1]].plot(timeVals, angVals, label = 'Exp. Data',
                                           linewidth = 2, color = colourDict['expData'])
    #Identify the experimental data section to cut axes down to
    #(i.e. half gait cycle)
    startTime = timeVals[0]
    endTime = timeVals[-1]
    
    #Torque tracking
    #Extract the relevant section of data
    angVals = df_torqueTracking[kinVars[kk]].to_numpy().flatten()
    timeVals = df_torqueTracking['time'].to_numpy().flatten()
    #Convert angular values if radians
    if 'pelvis_tx' not in kinVars[kk] and 'pelvis_ty' not in kinVars[kk]:
        angVals = np.rad2deg(angVals)
    #Plot data
    ax[whichAx[kk][0],whichAx[kk][1]].plot(timeVals, angVals, label = 'Torque Driven Tracking',
                                           linewidth = 2, color = colourDict['torqueTracking'])
    
    #Muscle tracking
    #Extract the relevant section of data
    angVals = df_muscleTracking[kinVars[kk]].to_numpy().flatten()
    timeVals = df_muscleTracking['time'].to_numpy().flatten()
    #Convert angular values if radians
    if 'pelvis_tx' not in kinVars[kk] and 'pelvis_ty' not in kinVars[kk]:
        angVals = np.rad2deg(angVals)
    #Plot data
    ax[whichAx[kk][0],whichAx[kk][1]].plot(timeVals, angVals, label = 'Muscle Driven Tracking',
                                           linewidth = 2, color = colourDict['muscleTracking'])
    
    #Muscle prediction
    #Extract the relevant section of data
    angVals = df_musclePrediction[kinVars[kk]].to_numpy().flatten()
    timeVals = df_musclePrediction['time'].to_numpy().flatten()
    #Convert angular values if radians
    if 'pelvis_tx' not in kinVars[kk] and 'pelvis_ty' not in kinVars[kk]:
        angVals = np.rad2deg(angVals)
    #Plot data
    ax[whichAx[kk][0],whichAx[kk][1]].plot(timeVals, angVals, label = 'Muscle Driven Prediction',
                                           linewidth = 2, color = colourDict['musclePrediction'])
    
    #Add title
    ax[whichAx[kk][0],whichAx[kk][1]].set_title(kinTitles[kk],
                                                fontsize = 12,
                                                fontweight = 'bold')
    
    #Set axes to time limits
    ax[whichAx[kk][0],whichAx[kk][1]].set_xlim([startTime,endTime])
    
#Tight layout
plt.tight_layout()

### TODO: add labels and legend

#Turn off axes not used
for rr in range(np.size(ax,0)):
    for cc in range(np.size(ax,1)):
        #Get current axes value
        currAx = [rr,cc]
        #Check if current axes is in plotting list
        if currAx not in whichAx:
            #Turn off axes
            ax[rr,cc].set_visible(False)

###



