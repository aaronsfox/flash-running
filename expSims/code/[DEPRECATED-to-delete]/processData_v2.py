# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 20:25:31 2020

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    [DEPRECATED - SUPERCEDED BY RUN SIMULATIONS I THINK...]
    
    This script processes the data provided by Dorn et al. (2012) using a tracking
    simulation approach to generate dynamically consistent sprinting data as a 
    feeder for further simulations.
    
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

# %% Experimental data marker tracking (WORKS)

"""

This section runs a partially (i.e. right limb) muscle-driven tracking simulation
of the experimental marker trajectories and the GRF data. No contact spheres or
rolling-on-surface constraints are considered, so the goal of this simulation is
to generate the model excitations that match the experimental data.

TODO:
    > Consider including tendon compliance for certain muscles?
        >> Nope turn it back off, not working well...
    > Consider removing floor body from model for visualisation purposes
    > Set guess from IK?

"""

#Create folder to store data
try:
    os.mkdir('..\\expTrackingSim')
except:
    print('Experimental tracking sim folder already detected...')
    
#Edit the model for use in the tracking sim tool
modelProcessor = osim.ModelProcessor('..\\data\\JA1_SCALED_Osim40_Muscles.osim')

#Append the necessary operators to the model
#Handles converting muscles to DeGrooteFregly and setting parameters
modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(2.0))
modelProcessor.append(osim.ModOpScaleMaxIsometricForce(2))
# modelProcessor.append(osim.ModOpTendonComplianceDynamicsModeDGF('implicit'))

#Process model for further edits
osimModel = modelProcessor.process()

#Remove left limb muscle and upper body forces to be replaced by torque actuators
#Set a list of forces to remove
removeForceInd = []
#Loop through forces and identify muscles to remove
for forceInd in range(osimModel.updForceSet().getSize()):
    #Check for muscle
    if osimModel.updForceSet().get(forceInd).getConcreteClassName().endswith('Muscle'):
        #Check for left hand side or upper body
        if osimModel.updForceSet().get(forceInd).getName().endswith('_l') or \
            osimModel.updForceSet().get(forceInd).getName().split('_')[0] in ['extobl', 'intobl', 'ercspn']:
            #Append index to list
            removeForceInd.append(forceInd)

#Remove the designated forces keeping in mind that the index reduces each time
#another force is removed
for removeInd in removeForceInd:
    osimModel.updForceSet().remove(removeInd - removeForceInd.index(removeInd))

#Set the coordinates that will need torque actuation
#Set an optimal force and label for each actuator
#Upper body and left side are idealised actuators
#Residuals and reserve actuators are lowly weighted for optimal force
optForces = {
    #upper body
    'lumbar_extension': [1000, 'actuator'], 'lumbar_bending': [1000, 'actuator'], 'lumbar_rotation': [1000, 'actuator'],
    'arm_flex_r': [300, 'actuator'], 'arm_add_r': [300, 'actuator'], 'arm_rot_r': [300, 'actuator'],
    'elbow_flex_r': [100, 'actuator'], 'pro_sup_r': [100, 'actuator'],
    'arm_flex_l': [300, 'actuator'], 'arm_add_l': [300, 'actuator'], 'arm_rot_l': [300, 'actuator'],
    'elbow_flex_l': [100, 'actuator'], 'pro_sup_l': [100, 'actuator'],
    #left limb
    'hip_flexion_l': [300, 'actuator'], 'hip_adduction_l': [300, 'actuator'], 'hip_rotation_l': [300, 'actuator'],
    'knee_angle_l': [300, 'actuator'], 'ankle_angle_l': [300, 'actuator'], 'mtp_angle_l': [300, 'actuator'],
    #right limb
    'hip_flexion_r': [2, 'reserve'], 'hip_adduction_r': [2, 'reserve'], 'hip_rotation_r': [2, 'reserve'],
    'knee_angle_r': [2, 'reserve'], 'ankle_angle_r': [2, 'reserve'],  'mtp_angle_r': [2, 'reserve'],
    #pelvis
    'pelvis_tx': [1, 'residual'], 'pelvis_ty': [1, 'residual'], 'pelvis_tz': [1, 'residual'],
    'pelvis_tilt': [1, 'residual'], 'pelvis_list': [1, 'residual'], 'pelvis_rotation': [1, 'residual']}

#Add torque actuators to model
for coordForce in optForces.keys():
    #Create the actuator
    actu = osim.CoordinateActuator()
    actu.setName(f'{coordForce}_{optForces[coordForce][1]}')
    actu.setCoordinate(osimModel.getCoordinateSet().get(coordForce))
    actu.setOptimalForce(optForces[coordForce][0])
    actu.setMinControl(np.inf*-1)
    actu.setMaxControl(np.inf)
    #Add to the models force set
    osimModel.updForceSet().cloneAndAppend(actu)

#Increase the maximum contraction velocity of muscles
maxContractionVelocity = 30
for muscleInd in range(osimModel.getMuscles().getSize()):
    #Contraction velocity
    osimModel.getMuscles().get(muscleInd).set_max_contraction_velocity(maxContractionVelocity)
    # #Tendon compliance
    # if '_gas_' in osimModel.getMuscles().get(muscleInd).getName() or 'soleus_' in osimModel.getMuscles().get(muscleInd).getName():
    #     osimModel.getMuscles().get(muscleInd).set_ignore_tendon_compliance(False)
        
#Finalise model connections
osimModel.finalizeConnections()

#Save model to file for later use if needed
osimModel.printToXML('..\\expTrackingSim\\expTrackingModel.osim')

#Create the tracking model processor for use in the tool
trackModelProcessor = osim.ModelProcessor(osimModel)

#Set joints to weld for simulation
jointsToWeld = ['subtalar_r', 'subtalar_l',
                # 'mtp_r', 'mtp_l',
                'radius_hand_r', 'radius_hand_l']
weldVectorStr = osim.StdVectorString()
[weldVectorStr.append(joint) for joint in jointsToWeld]

#Add remaining options to model processor
#Joints to weld and external loads
trackModelProcessor.append(osim.ModOpAddExternalLoads('..\\data\\sprint_grf.xml'))
trackModelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName('expTrackingSim')

#Set model in tracking tool
track.setModel(trackModelProcessor)

#Set the markers reference from TRC file
#Note that this needs to be done using a flattened table due to a bug
#See: https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=1815&t=14612&p=43550&start=0&view=
#### TODO: consider filtering?
markerData = osim.TimeSeriesTableVec3('..\\data\\sprint.trc')
markerTableProcessor = osim.TableProcessor(markerData.flatten())
track.setMarkersReference(markerTableProcessor)

#Set allow unused references in case there are extra markers
track.set_allow_unused_references(True)

#Set the global markers tracking weight
track.set_markers_global_tracking_weight(1)

#Create a dictionary to allocate marker weights
##### TODO: optimise this a little better
markerWeightDict = {'C7': 1, 'RSH': 1, 'LSH': 1, 'MAN': 1, 'T7': 1, 'LARM': 1,
                    'LELB': 1, 'LFOREARM': 1, 'LWR': 1, 'RARM': 1, 'RELB': 1, 
                    'RFOREARM': 1, 'RWR': 1, 'RASI': 1, 'LASI': 1, 'SACR': 1,
                    'LTHLP': 1, 'LTHLD': 1, 'LTHAP': 1, 'LTHAD': 1, 'LLEPI': 1,
                    'LPAT': 1, 'LTIAP': 1, 'LTIAD': 1, 'LTILAT': 1, 'LLMAL': 1,
                    'LHEEL': 1, 'LMFS': 1, 'LMFL': 1, 'LP1MT': 1, 'LTOE': 1,
                    'LP5MT': 1, 'RTHLP': 1, 'RTHLD': 1, 'RTHAP': 1, 'RTHAD': 1,
                    'RLEPI': 1, 'RPAT': 1, 'RTIAP': 1, 'RTIAD': 1, 'RTILAT': 1,
                    'RLMAL': 1, 'RHEEL': 1, 'RMFS': 1, 'RMFL': 1, 'RP1MT': 1,
                    'RTOE': 1, 'RP5MT': 1}

#Set a the tracking weights for markers in a weights object
markerWeights = osim.MocoWeightSet()
for marker in markerWeightDict.keys():
    markerWeights.cloneAndAppend(osim.MocoWeight(marker, markerWeightDict[marker]))
    
#Set marker weight set in problem
track.set_markers_weight_set(markerWeights)

#Set start and end times in problem
#Get the gait timings from helper function for a full gait cycle
startTime, endTime = helper.getGaitTimings(grfFile = '..\\data\\sprint_grf.mot',
                                           extLoads = '..\\data\\sprint_grf.xml',
                                           startForceName = 'RightGRF1',
                                           stopForceName = 'RightGRF2',
                                           forceThreshold = 20)

#Set in tracking problem
track.set_initial_time(startTime)
track.set_final_time(endTime)

#Set mesh interval
#### TODO: consider calculating this more objectively
#### Set as mesh intervals later
# track.set_mesh_interval(0.02)

#Initialise the Moco study and problem
study = track.initialize()
problem = study.updProblem()
        
#Update the weight & exponent on the default control effort goal
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
effort.setWeight(0.001)
effort.setExponent(3)

#Update individual weights in control effort goal to be relative to
#actual muscle and reserve actuator names
#Set appropriate patterns in the weight set
#Muscles
effort.setWeightForControlPattern('/forceset/.*/activation', 0.1)
#Idealised actuators
effort.setWeightForControlPattern('/forceset/.*_actuator', 0.01)
#Reserves
effort.setWeightForControlPattern('/forceset/.*_reserve', 10)
#Residuals
effort.setWeightForControlPattern('/forceset/.*_residual', 10)

#Set appropriate bounds on tendon force
problem.setStateInfoPattern('/forceset/.*/normalized_tendon_force', [0, 3.0], [], [])

#Define the solver and set its options
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_optim_constraint_tolerance(0.01)
solver.set_optim_convergence_tolerance(0.01)
solver.set_num_mesh_intervals(25) #### TODO: appropriate? too much? 50 seems to align with Falisse et al.?
solver.set_multibody_dynamics_mode('explicit')
solver.set_minimize_implicit_auxiliary_derivatives(True)
solver.set_implicit_auxiliary_derivatives_weight(0.001)
solver.resetProblem(problem)

#Set the normalized tendon forces if not loading initial guess from file
#Adapted from Ross Miller UMocoD code
guess = solver.getGuess()
numRows = guess.getNumTimes()
osimModel.initSystem()
stateNames = osimModel.getStateVariableNames()
for ii in range(osimModel.getNumStateVariables()):
    currentStateName = stateNames.getitem(ii)
    if 'normalized_tendon_force' in currentStateName:
        guess.setState(currentStateName, np.linspace(0.2,0.2,numRows))

#Set updated guess in solver
solver.setGuess(guess)

#Solve the tracking problem
solution = study.solve()
# study.visualize(solution)

#Copy tracked markers file over to main directory
shutil.move('expTrackingSim_tracked_markers.sto',
            '..\\expTrackingSim\\expTrackingSim_tracked_markers.sto')

#Save solution to file
solution.write('..\\expTrackingSim\\expTrackingSim_markerTrackingSolution.sto')

"""

Results from above simulation look pretty solid with respect to muscle activations
and minimal input from reserve/residual torques.

"""

# %% Kinematic & GRF tracking sim (SEEMS OK w/ LESS CONTACT SPHERES)

"""

This section runs a partially (i.e. right limb) muscle-driven tracking simulation
of the experimental kinematics and the GRF data. The difference in this 
simulation is that no external loads are considered and GRFs are tracked with
contact spheres attached to the feet. It tracks the coordinates from the previous
tracking sim generated from the marker data, alongside the GRFs.

TODO:
    > Consider including tendon compliance for certain muscles?
        >> Not really working well...
    > Consider a coarse solution as an initial guess for the more refined problem
      to help speed up convergence?
    > Consider periodicity goal?
        >> Phases to simulate the two strides with periodicity?
            >>> Or split into 2 for an initial sim?

"""

#Create folder to store data
try:
    os.mkdir('..\\fullTrackingSim')
except:
    print('Full tracking sim folder already detected...')
    
#Read in the previously created model for the tracking sim to edit
osimModel = osim.Model('..\\expTrackingSim\\expTrackingModel.osim')
    
#Add contact spheres at foot-ground contact model locations

#Get the markers that contain an fp reference for contact sphere locations
footGroundNames = []
for markerInd in range(osimModel.updMarkerSet().getSize()):
    #Check for fp indicator in marker
    if '_fp_' in osimModel.updMarkerSet().get(markerInd).getName():
        footGroundNames.append(osimModel.updMarkerSet().get(markerInd).getName())
        
#Get foot ground name locations in dictionary
footGroundLocs = {name: np.array((osimModel.updMarkerSet().get(name).get_location().get(0),
                                  osimModel.updMarkerSet().get(name).get_location().get(1),
                                  osimModel.updMarkerSet().get(name).get_location().get(2))) for name in footGroundNames}

#Convert to a singular heel and toe foot ground contact locations
#### NOTE: only heel and mtp in foot ground spheres currently...
footSphereLocs = {}
for name in footGroundLocs:
    if 'heel1' in name:
        #Set the other name of the heel location
        name2 = name.replace('1','2')
        #Get the midpoint of two locations
        mpLoc = np.array(((footGroundLocs[name][0]+footGroundLocs[name2][0])/2,
                          (footGroundLocs[name][1]+footGroundLocs[name2][1])/2,
                          (footGroundLocs[name][2]+footGroundLocs[name2][2])/2))
        #Set in dictionary
        footSphereLocs['heel_'+name[0]] = mpLoc
    # elif 'mt' in name or 'toe' in name:
    elif 'mt' in name:
        #Set the name and location as what is already done
        footSphereLocs[name.split('_')[-1]+'_'+name[0]] = footGroundLocs[name]
        
#Set foot ground contact sphere sizes
heelSphereRadius = 0.04
otherSphereRadius = 0.03

#Add sphere size into dictionary
footSphereSizes = {name: heelSphereRadius if 'heel' in name else otherSphereRadius for name in footSphereLocs.keys()}

#Create and connect the half space floor to ground
floorContact = osim.ContactHalfSpace() #create object
floorContact.setName('floor') #set name
floorContact.setOrientation(osim.Vec3(0, 0, -1.5707963267949001)) #set orientation
floorContact.connectSocket_frame(osim.PhysicalFrame.safeDownCast(osimModel.getGround())) #connect to ground
osimModel.addContactGeometry(floorContact) #add to model

#Iteratively create the contact geometry and spheres and attach to models force set
for contactPoint in footSphereLocs.keys():
    
    #Create the contact geometry and set parameters
    contactGeom = osim.ContactSphere() #create object
    contactGeom.setName(contactPoint) #set name
    
    #Conditional for heel or toe sphere
    if contactPoint.startswith('heel') and contactPoint.endswith('_r'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('calcn_r')) #connect to frame
    elif contactPoint.startswith('heel') and contactPoint.endswith('_l'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('calcn_l')) #connect to frame
    elif not contactPoint.startswith('heel') and contactPoint.endswith('_r'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('toes_r')) #connect to frame
    elif not contactPoint.startswith('heel') and contactPoint.endswith('_l'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('toes_l')) #connect to frame
        
    #Get sphere size to set and adjust location
    currSphereSize = footSphereSizes[contactPoint]    
        
    #Conditional whether to transform point location
    if 'heel' in contactPoint:
        contactGeom.setLocation(osim.Vec3(footSphereLocs[contactPoint][0],
                                          # footGroundLocs[contactPoint][1],
                                          footSphereLocs[contactPoint][1] + (currSphereSize/2), #adjust to bottom of foot height
                                          footSphereLocs[contactPoint][2])) #set location
    else:
        
        #Get translation from calcn to toe joint and subtract from location
        if contactPoint.endswith('r_'):
            jointTranslation = osimModel.updJointSet().get('mtp_r').get_frames(0).get_translation()
        else:
            jointTranslation = osimModel.updJointSet().get('mtp_l').get_frames(0).get_translation()
        contactGeom.setLocation(osim.Vec3(footSphereLocs[contactPoint][0] - jointTranslation.get(0),
                                          # footSphereLocs[contactPoint][1] - jointTranslation.get(1),
                                          footSphereLocs[contactPoint][1] - jointTranslation.get(1) + (currSphereSize/2), #adjust to bottom of foot height
                                          footSphereLocs[contactPoint][2] - jointTranslation.get(2))) #set location
    contactGeom.setRadius(currSphereSize) #set radius
    osimModel.addContactGeometry(contactGeom) #add to model
    
    #Create the sphere and set properties
    contactSphere = osim.SmoothSphereHalfSpaceForce() #create force
    contactSphere.setName('contact_'+contactPoint) #set name
    contactSphere.connectSocket_half_space(osimModel.updContactGeometrySet().get('floor')) #connect to floor
    contactSphere.connectSocket_sphere(osimModel.updContactGeometrySet().get(contactPoint)) #connect to sphere
    contactSphere.set_stiffness(3067776) #set stiffness
    contactSphere.set_dissipation(2) #set dissipation
    contactSphere.set_static_friction(0.8) #set static friction
    contactSphere.set_dynamic_friction(0.8) #set dynamic friction
    contactSphere.set_viscous_friction(0.5) #set viscous friction
    contactSphere.set_transition_velocity(0.2) #set transition velocity
    contactSphere.set_hertz_smoothing(300) #set hertz smoothing
    contactSphere.set_hunt_crossley_smoothing(50) #set hunt crossley smoothing
    osimModel.addForce(contactSphere) #add to model

#Finalise model connections
osimModel.finalizeConnections()

#Save model to file for later use if needed
osimModel.printToXML('..\\fullTrackingSim\\fullTrackingModel.osim')

#Create the tracking model processor for use in the tool
trackModelProcessor = osim.ModelProcessor(osimModel)

#Set joints to weld for simulation
jointsToWeld = ['subtalar_r', 'subtalar_l',
                # 'mtp_r', 'mtp_l',
                'radius_hand_r', 'radius_hand_l']
weldVectorStr = osim.StdVectorString()
[weldVectorStr.append(joint) for joint in jointsToWeld]

#Add remaining options to model processor
#Joints to weld
trackModelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName('fullTrackingSim')

#Set model in tracking tool
track.setModel(trackModelProcessor)

#Set the coordinates reference in tracking tool

#Convert the solution data from the tracking sim to coordinates
trackingTraj = osim.MocoTrajectory('..\\expTrackingSim\\expTrackingSim_markerTrackingSolution.sto')

#Export to states table and write to file
trackingStates = trackingTraj.exportToStatesTable()
osim.STOFileAdapter().write(trackingStates, '..\\fullTrackingSim\\coordinatesToTrack.sto')

#Set coordinates file as reference in tool
track.setStatesReference(osim.TableProcessor('..\\fullTrackingSim\\coordinatesToTrack.sto'))
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
# taskWeights = {'pelvis_tx': 1e-0, 'pelvis_ty': 1e-1, 'pelvis_tz': 1e-0,
#                'pelvis_tilt': 1e-0, 'pelvis_list': 1e-0, 'pelvis_rotation': 1e-0, 
#                'hip_flexion_r': 1e-0, 'hip_adduction_r': 1e-1, 'hip_rotation_r': 1e-2, 
#                'knee_angle_r': 1e-0, 'ankle_angle_r': 1e-2,
#                # 'subtalar_angle_r': 0,
#                'mtp_angle_r': 1e-2,
#                'hip_flexion_l': 1e-0, 'hip_adduction_l': 1e-1, 'hip_rotation_l': 1e-2, 
#                'knee_angle_l': 1e-0, 'ankle_angle_l': 1e-2,
#                # 'subtalar_angle_l': 0,
#                'mtp_angle_l': 1e-2,
#                'lumbar_extension': 1e-1, 'lumbar_bending': 1e-2, 'lumbar_rotation': 1e-2,
#                'arm_flex_r': 1e-1, 'arm_add_r': 1e-1, 'arm_rot_r': 1e-1,
#                'elbow_flex_r': 1e-1, 'pro_sup_r': 1e-1,
#                'arm_flex_l': 1e-1, 'arm_add_l': 1e-1, 'arm_rot_l': 1e-1,
#                'elbow_flex_l': 1e-1, 'pro_sup_l': 1e-1}

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
#Get the gait timings from helper function for a full gait cycle
startTime, endTime = helper.getGaitTimings(grfFile = '..\\data\\sprint_grf.mot',
                                            extLoads = '..\\data\\sprint_grf.xml',
                                            startForceName = 'RightGRF1',
                                            stopForceName = 'RightGRF2',
                                            forceThreshold = 20)
# startTime, endTime = helper.getGaitTimings(grfFile = '..\\data\\sprint_grf.mot',
#                                            extLoads = '..\\data\\sprint_grf.xml',
#                                            startForceName = 'RightGRF1',
#                                            stopForceName = 'LeftGRF1',
#                                            forceThreshold = 20)

#Set in tracking problem
track.set_initial_time(startTime)
track.set_final_time(endTime)

#Set mesh interval
#### TODO: consider calculating this more objectively
#### Set as mesh intervals later
# track.set_mesh_interval(0.02)

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
effort.setWeightForControlPattern('/forceset/.*_actuator', 0.01)
#Reserves
effort.setWeightForControlPattern('/forceset/.*_reserve', 10)
#Residuals
effort.setWeightForControlPattern('/forceset/.*_residual', 10) #lower than 100 previously...

#Add contact tracking goal
#Uses three separate vectors to apply different contact tracking weights

#Set right and left contact sphere groups
rightFootContacts = [f'/forceset/contact_{pointName}' for pointName in footSphereLocs.keys() if pointName.endswith('_r')]
leftFootContacts = [f'/forceset/contact_{pointName}' for pointName in footSphereLocs.keys() if pointName.endswith('_l')]

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
                            'weight': [1e-2, 1e-2, 1e-2],
                            # 'weight': [1, 0.5, 1],
                            # 'weight': [0.25, 0.1, 0.25],
                            }

#Loop through contact tracking settings to create separate weighted goals
for ii in range(len(contactTrackingSettings['vector'])):
    
    #Create tracking goal
    contactTracking.append(osim.MocoContactTrackingGoal(f'contact{ii}',
                                                        contactTrackingSettings['weight'][ii]))
    
    #Set external loads
    contactTracking[ii].setExternalLoadsFile('..\\data\\sprint_grf.xml')
    
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

#All states are periodic except pelvis anterior-posterior translation
for coordName in taskWeights.keys():
    if coordName != 'pelvis_tx':
        #Get full path to state
        stateName = osimModel.updCoordinateSet().get(coordName).getAbsolutePathString()+'/value'
        #Add state pair to goal
        periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair(stateName))

#Add to problem
problem.addGoal(periodicityGoal)
    
#Set kinematic bounds in problem

#### TODO: should do the same thing with final bounds...not if periodic constraint added

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

#Create a dictionary for kinematic bounds based on max and min values in experimental coordinates
#Here we allow a them to go 25% of the total range of motion below and above the min/max
kinematicBounds = {}

#Also create a dictionary which limits initial bounds to be relative to starting point
#Here we allow a 10% window either side of the starting value

#Loop through coordinates
for coordName in boundedCoords:
    
    #Get the full path to coordinate
    coordPath = osimModel.updCoordinateSet().get(coordName).getAbsolutePathString()
    
    #Get minimum and maximum value from the tracking states
    #Calculate range
    minVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy().min()
    maxVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy().max()
    rangeVal = np.diff((minVal,maxVal))[0]
    
    #Get the starting and final value from the tracked states
    initialVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy()[0]
    finalVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy()[-1]
    
    #Calculate total and initial bound ranges
    #Calculate and set the 20% range either side of the min and maximum
    #Set the initial value to be within 5% of the starting value
    totalBounds = [minVal - (rangeVal * 0.25), maxVal + (rangeVal * 0.25)]
    initialBounds = [initialVal - (np.abs(initialVal) * 0.1), initialVal + (np.abs(initialVal) * 0.1)]
    finalBounds = [finalVal - (np.abs(finalVal) * 0.1), finalVal + (np.abs(finalVal) * 0.1)]
    
    #Set bounds in problem
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
        
    #Set bounds in problem
    #Put check in place to correct initial bounds if outside of total ranges
    correctedFinalBounds = []
    #Lower bound
    if finalBounds[0] > totalBounds[0]:
        correctedFinalBounds.append(finalBounds[0])
    else:
        correctedFinalBounds.append(totalBounds[0])
    #Upper bound
    if finalBounds[1] < totalBounds[1]:
        correctedFinalBounds.append(finalBounds[1])
    else:
        correctedFinalBounds.append(totalBounds[1])
    
    #Set in problem
    problem.setStateInfo(coordPath+'/value',
                         #Total bounds range
                         totalBounds,
                         #Initial bounds range
                         correctedInitialBounds,
                         #Final bounds range
                         #correctedFinalBounds #### not if using a periodic constraint
                         )
    
#Define the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_optim_constraint_tolerance(0.01)
solver.set_optim_convergence_tolerance(0.01)
solver.set_multibody_dynamics_mode('explicit')
# solver.set_minimize_implicit_multibody_accelerations(True) #smoothness criterion?
# solver.set_implicit_multibody_accelerations_weight(0.0001)
solver.set_num_mesh_intervals(25) #### TODO: appropriate? too much? 50 seems to align with Falisse et al.? coarse to start with...?
solver.resetProblem(problem)

#Set guess from tracking simulation
# solver.setGuessFile('..\\expTrackingSim\\expTrackingSim_markerTrackingSolution.sto')
solver.setGuessFile('..\\fullTrackingSim\\fullTrackingSim_statesTrackingSolution_coarse.sto')
solver.resetProblem(problem)

# #Alter the initial guess to bump up pelvis_ty to ensure contact spheres don't
# #start by penetrating the ground
# initialGuess = solver.getGuess()
# initialGuess.setState('/jointset/ground_pelvis/pelvis_ty/value',
#                       initialGuess.getState('/jointset/ground_pelvis/pelvis_ty/value').to_numpy() + 0.1)
# solver.setGuess(initialGuess)

# #If using implicit mode then multibody accelerations need to be generated in guess
# initialGuess = osim.MocoTrajectory('..\\expTrackingSim\\expTrackingSim_markerTrackingSolution.sto')
# initialGuess.generateAccelerationsFromSpeeds()
# solver.setGuess(initialGuess)
# solver.resetProblem(problem)

#Solve the tracking problem
solution = study.solve()
# study.visualize(solution)

#Copy tracked markers file over to main directory
shutil.move('fullTrackingSim_tracked_states.sto',
            '..\\fullTrackingSim\\fullTrackingSim_tracked_states.sto')

#Save solution to file
solution.write('..\\fullTrackingSim\\fullTrackingSim_statesTrackingSolution.sto')

#Extract predicted GRFs

#Create forces table
externalForcesTableFlat = osim.createExternalLoadsTableForGait(trackModelProcessor.process(), solution,
                                                               forcesRightFoot, forcesLeftFoot)

#Write table to file
osim.STOFileAdapter().write(externalForcesTableFlat,
                            '..\\fullTrackingSim\\fullTrackingSim_statesTrackingSolution_grf.sto')

"""

NOTES: seems to work OK, but not perfect with the 3 spheres; torque tracking side doesn't
track as well - maybe because it's difficult to produce lower torques with the high
optimal forces?

Tested larger heel sphere with reduced state (0.1) and GRF (0.01) tracking weights, at a coarse mesh interval (10)...
    > Didn't seem to fix the issue
    
Maybe test out without contact traking to see how well state tracking can go...?
    > Can be achieved quite well
    > GRF on right foot is minimal but GRF on left foot is sort of OK, yet quite high...
    
Edited weights to scale back while increasing state tracking to 1
    > i.e. balanced specific coordinate weights but increased total focus relative to GRF and other goals
        >> States global weight of 1 with individual weights at 0-10 + GRF weights at 1e-2 seems pretty well scaled...
        >> Had to stop optimisation at function value of ~4.6 after 135 iters
            >> Visually was starting to look OK too...
    > Also widened the bounds (25% total and 10% initial) and increased mesh (50)
    > Seems to be a little glitchy though
        >> Minimise multibody accelerations?
        
Implicit dynamics and minimise accelerations
    > Minimising accelerations seems to make it smoother but not follow kinematics to achieve this
    > Having both of these makes iters increase and slower
    
Just using implicit mode?
    > Implicit not really helping

See: https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=1815&t=13756&p=0&start=0&view=&sid=0112f28ec76e50df8f39a72f6b4407e0
    > Suggest explicit multibody with implicit tendon dynamics while minimising derivatives
    > Also try with a coarse solution first to get a better initial guess...
    > Could also try scaling joint coordinate weights as per Miller based on values in data

"""

# %% Two-phase approach [part 1]

"""

This section runs a similar tracking simulation to above, but splits into two phases
to ensure symmetry and each step and hopefully improve the smoothness of the motion.

TODO:
    > ...

"""

#Create folder to store data
try:
    os.mkdir('..\\fullTrackingSim_rightLimb')
except:
    print('Full tracking sim (right limb) folder already detected...')
    
#Read in the previously created model for the tracking sim to edit
osimModel = osim.Model('..\\expTrackingSim\\expTrackingModel.osim')
    
#Add contact spheres at foot-ground contact model locations

#Get the markers that contain an fp reference for contact sphere locations
footGroundNames = []
for markerInd in range(osimModel.updMarkerSet().getSize()):
    #Check for fp indicator in marker
    if '_fp_' in osimModel.updMarkerSet().get(markerInd).getName():
        footGroundNames.append(osimModel.updMarkerSet().get(markerInd).getName())
        
#Get foot ground name locations in dictionary
footGroundLocs = {name: np.array((osimModel.updMarkerSet().get(name).get_location().get(0),
                                  osimModel.updMarkerSet().get(name).get_location().get(1),
                                  osimModel.updMarkerSet().get(name).get_location().get(2))) for name in footGroundNames}

#Convert to a singular heel and toe foot ground contact locations
#### NOTE: only heel and mtp in foot ground spheres currently...
footSphereLocs = {}
for name in footGroundLocs:
    if 'heel1' in name:
        #Set the other name of the heel location
        name2 = name.replace('1','2')
        #Get the midpoint of two locations
        mpLoc = np.array(((footGroundLocs[name][0]+footGroundLocs[name2][0])/2,
                          (footGroundLocs[name][1]+footGroundLocs[name2][1])/2,
                          (footGroundLocs[name][2]+footGroundLocs[name2][2])/2))
        #Set in dictionary
        footSphereLocs['heel_'+name[0]] = mpLoc
    # elif 'mt' in name or 'toe' in name:
    elif 'mt' in name:
        #Set the name and location as what is already done
        footSphereLocs[name.split('_')[-1]+'_'+name[0]] = footGroundLocs[name]
        
#Set foot ground contact sphere sizes
heelSphereRadius = 0.04
otherSphereRadius = 0.03

#Add sphere size into dictionary
footSphereSizes = {name: heelSphereRadius if 'heel' in name else otherSphereRadius for name in footSphereLocs.keys()}

#Create and connect the half space floor to ground
floorContact = osim.ContactHalfSpace() #create object
floorContact.setName('floor') #set name
floorContact.setOrientation(osim.Vec3(0, 0, -1.5707963267949001)) #set orientation
floorContact.connectSocket_frame(osim.PhysicalFrame.safeDownCast(osimModel.getGround())) #connect to ground
osimModel.addContactGeometry(floorContact) #add to model

#Iteratively create the contact geometry and spheres and attach to models force set
for contactPoint in footSphereLocs.keys():
    
    #Create the contact geometry and set parameters
    contactGeom = osim.ContactSphere() #create object
    contactGeom.setName(contactPoint) #set name
    
    #Conditional for heel or toe sphere
    if contactPoint.startswith('heel') and contactPoint.endswith('_r'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('calcn_r')) #connect to frame
    elif contactPoint.startswith('heel') and contactPoint.endswith('_l'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('calcn_l')) #connect to frame
    elif not contactPoint.startswith('heel') and contactPoint.endswith('_r'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('toes_r')) #connect to frame
    elif not contactPoint.startswith('heel') and contactPoint.endswith('_l'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('toes_l')) #connect to frame
        
    #Get sphere size to set and adjust location
    currSphereSize = footSphereSizes[contactPoint]    
        
    #Conditional whether to transform point location
    if 'heel' in contactPoint:
        contactGeom.setLocation(osim.Vec3(footSphereLocs[contactPoint][0],
                                          # footGroundLocs[contactPoint][1],
                                          footSphereLocs[contactPoint][1] + (currSphereSize/2), #adjust to bottom of foot height
                                          footSphereLocs[contactPoint][2])) #set location
    else:
        
        #Get translation from calcn to toe joint and subtract from location
        if contactPoint.endswith('r_'):
            jointTranslation = osimModel.updJointSet().get('mtp_r').get_frames(0).get_translation()
        else:
            jointTranslation = osimModel.updJointSet().get('mtp_l').get_frames(0).get_translation()
        contactGeom.setLocation(osim.Vec3(footSphereLocs[contactPoint][0] - jointTranslation.get(0),
                                          # footSphereLocs[contactPoint][1] - jointTranslation.get(1),
                                          footSphereLocs[contactPoint][1] - jointTranslation.get(1) + (currSphereSize/2), #adjust to bottom of foot height
                                          footSphereLocs[contactPoint][2] - jointTranslation.get(2))) #set location
    contactGeom.setRadius(currSphereSize) #set radius
    osimModel.addContactGeometry(contactGeom) #add to model
    
    #Create the sphere and set properties
    contactSphere = osim.SmoothSphereHalfSpaceForce() #create force
    contactSphere.setName('contact_'+contactPoint) #set name
    contactSphere.connectSocket_half_space(osimModel.updContactGeometrySet().get('floor')) #connect to floor
    contactSphere.connectSocket_sphere(osimModel.updContactGeometrySet().get(contactPoint)) #connect to sphere
    contactSphere.set_stiffness(3067776) #set stiffness
    contactSphere.set_dissipation(2) #set dissipation
    contactSphere.set_static_friction(0.8) #set static friction
    contactSphere.set_dynamic_friction(0.8) #set dynamic friction
    contactSphere.set_viscous_friction(0.5) #set viscous friction
    contactSphere.set_transition_velocity(0.2) #set transition velocity
    contactSphere.set_hertz_smoothing(300) #set hertz smoothing
    contactSphere.set_hunt_crossley_smoothing(50) #set hunt crossley smoothing
    osimModel.addForce(contactSphere) #add to model

#Finalise model connections
osimModel.finalizeConnections()

#Save model to file for later use if needed
osimModel.printToXML('..\\fullTrackingSim_rightLimb\\fullTrackingModel.osim')

#Create the tracking model processor for use in the tool
trackModelProcessor = osim.ModelProcessor(osimModel)

#Set joints to weld for simulation
jointsToWeld = ['subtalar_r', 'subtalar_l',
                # 'mtp_r', 'mtp_l',
                'radius_hand_r', 'radius_hand_l']
weldVectorStr = osim.StdVectorString()
[weldVectorStr.append(joint) for joint in jointsToWeld]

#Add remaining options to model processor
#Joints to weld
trackModelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName('fullTrackingSim_rightLimb')

#Set model in tracking tool
track.setModel(trackModelProcessor)

#Set the coordinates reference in tracking tool

#Convert the solution data from the tracking sim to coordinates
trackingTraj = osim.MocoTrajectory('..\\expTrackingSim\\expTrackingSim_markerTrackingSolution.sto')

#Export to states table and write to file
trackingStates = trackingTraj.exportToStatesTable()
osim.STOFileAdapter().write(trackingStates, '..\\fullTrackingSim_rightLimb\\coordinatesToTrack.sto')

#Set coordinates file as reference in tool
track.setStatesReference(osim.TableProcessor('..\\fullTrackingSim_rightLimb\\coordinatesToTrack.sto'))
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
# taskWeights = {'pelvis_tx': 1e-0, 'pelvis_ty': 1e-1, 'pelvis_tz': 1e-0,
#                'pelvis_tilt': 1e-0, 'pelvis_list': 1e-0, 'pelvis_rotation': 1e-0, 
#                'hip_flexion_r': 1e-0, 'hip_adduction_r': 1e-1, 'hip_rotation_r': 1e-2, 
#                'knee_angle_r': 1e-0, 'ankle_angle_r': 1e-2,
#                # 'subtalar_angle_r': 0,
#                'mtp_angle_r': 1e-2,
#                'hip_flexion_l': 1e-0, 'hip_adduction_l': 1e-1, 'hip_rotation_l': 1e-2, 
#                'knee_angle_l': 1e-0, 'ankle_angle_l': 1e-2,
#                # 'subtalar_angle_l': 0,
#                'mtp_angle_l': 1e-2,
#                'lumbar_extension': 1e-1, 'lumbar_bending': 1e-2, 'lumbar_rotation': 1e-2,
#                'arm_flex_r': 1e-1, 'arm_add_r': 1e-1, 'arm_rot_r': 1e-1,
#                'elbow_flex_r': 1e-1, 'pro_sup_r': 1e-1,
#                'arm_flex_l': 1e-1, 'arm_add_l': 1e-1, 'arm_rot_l': 1e-1,
#                'elbow_flex_l': 1e-1, 'pro_sup_l': 1e-1}

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
#Get the gait timings from helper function for the half gait cycle
startTime, endTime = helper.getGaitTimings(grfFile = '..\\data\\sprint_grf.mot',
                                            extLoads = '..\\data\\sprint_grf.xml',
                                            startForceName = 'RightGRF1',
                                            stopForceName = 'LeftGRF1',
                                            forceThreshold = 20)

#Set in tracking problem
track.set_initial_time(startTime)
track.set_final_time(endTime)

#Set mesh interval
#### TODO: consider calculating this more objectively
#### Set as mesh intervals later
# track.set_mesh_interval(0.02)

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
effort.setWeightForControlPattern('/forceset/.*_actuator', 0.01)
#Reserves
effort.setWeightForControlPattern('/forceset/.*_reserve', 10)
#Residuals
effort.setWeightForControlPattern('/forceset/.*_residual', 10) #lower than 100 previously...

#Add contact tracking goal
#Uses three separate vectors to apply different contact tracking weights

#Set right and left contact sphere groups
rightFootContacts = [f'/forceset/contact_{pointName}' for pointName in footSphereLocs.keys() if pointName.endswith('_r')]
leftFootContacts = [f'/forceset/contact_{pointName}' for pointName in footSphereLocs.keys() if pointName.endswith('_l')]

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
                            'weight': [1e-2, 1e-2, 1e-2],
                            # 'weight': [1, 0.5, 1],
                            # 'weight': [0.25, 0.1, 0.25],
                            }

#Loop through contact tracking settings to create separate weighted goals
for ii in range(len(contactTrackingSettings['vector'])):
    
    #Create tracking goal
    contactTracking.append(osim.MocoContactTrackingGoal(f'contact{ii}',
                                                        contactTrackingSettings['weight'][ii]))
    
    #Set external loads
    contactTracking[ii].setExternalLoadsFile('..\\data\\sprint_grf.xml')
    
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

#Opposite states are periodic except pelvis anterior-posterior translation
for coordName in taskWeights.keys():
    
    if coordName != 'pelvis_tx':
        
        #Get full path to state
        stateName = osimModel.updCoordinateSet().get(coordName).getAbsolutePathString()+'/value'
        
        #Check for singular periodicity conditions
        if coordName in ['pelvis_ty', 'pelvis_tz', 'lumbar_extension']:
            
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

#Add to problem
problem.addGoal(periodicityGoal)
    
#Set kinematic bounds in problem

#### TODO: should do the same thing with final bounds...not if periodic constraint added

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
#Here we allow a 10% window either side of the starting value

#Loop through coordinates
for coordName in boundedCoords:
    
    #Get the full path to coordinate
    coordPath = osimModel.updCoordinateSet().get(coordName).getAbsolutePathString()
    
    #Get minimum and maximum value from the tracking states
    #Calculate range
    minVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy().min()
    maxVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy().max()
    rangeVal = np.diff((minVal,maxVal))[0]
    
    #Get the starting and final value from the tracked states
    initialVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy()[0]
    finalVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy()[-1]
    
    #Calculate total and initial bound ranges
    #Calculate and set the 20% range either side of the min and maximum
    #Set the initial value to be within 5% of the starting value
    totalBounds = [minVal - (rangeVal * 0.25), maxVal + (rangeVal * 0.25)]
    initialBounds = [initialVal - (np.abs(initialVal) * 0.1), initialVal + (np.abs(initialVal) * 0.1)]
    finalBounds = [finalVal - (np.abs(finalVal) * 0.1), finalVal + (np.abs(finalVal) * 0.1)]
    
    #Set bounds in problem
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
        
    #Set bounds in problem
    #Put check in place to correct initial bounds if outside of total ranges
    correctedFinalBounds = []
    #Lower bound
    if finalBounds[0] > totalBounds[0]:
        correctedFinalBounds.append(finalBounds[0])
    else:
        correctedFinalBounds.append(totalBounds[0])
    #Upper bound
    if finalBounds[1] < totalBounds[1]:
        correctedFinalBounds.append(finalBounds[1])
    else:
        correctedFinalBounds.append(totalBounds[1])
    
    #Set in problem
    problem.setStateInfo(coordPath+'/value',
                         #Total bounds range
                         totalBounds,
                         #Initial bounds range
                         #correctedInitialBounds, ### not if using periodic constraint?
                         #Final bounds range
                         #correctedFinalBounds #### not if using a periodic constraint
                         )
    
#Define the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_optim_constraint_tolerance(0.01)
solver.set_optim_convergence_tolerance(0.01)
solver.set_multibody_dynamics_mode('explicit')
# solver.set_minimize_implicit_multibody_accelerations(True) #smoothness criterion?
# solver.set_implicit_multibody_accelerations_weight(0.0001)
solver.set_num_mesh_intervals(10) #### TODO: coarse to start with?
solver.resetProblem(problem)

# #Set guess from tracking simulation
# # solver.setGuessFile('..\\expTrackingSim\\expTrackingSim_markerTrackingSolution.sto')
# solver.setGuessFile('..\\fullTrackingSim\\fullTrackingSim_statesTrackingSolution_coarse.sto')
# solver.resetProblem(problem)

#Solve the tracking problem
solution = study.solve()
# study.visualize(solution)

#Copy tracked states file over to main directory
shutil.move('fullTrackingSim_rightLimb_tracked_states.sto',
            '..\\fullTrackingSim_rightLimb\\fullTrackingSim_rightLimb_tracked_states.sto')

#Save solution to file
solution.write('..\\fullTrackingSim_rightLimb\\fullTrackingSim_rightLimb_statesTrackingSolution.sto')

#Extract predicted GRFs

#Create forces table
externalForcesTableFlat = osim.createExternalLoadsTableForGait(trackModelProcessor.process(), solution,
                                                                forcesRightFoot, forcesLeftFoot)

#Write table to file
osim.STOFileAdapter().write(externalForcesTableFlat,
                            '..\\fullTrackingSim_rightLimb\\fullTrackingSim_rightLimb_statesTrackingSolution_grf.sto')

"""

NOTES: 
    
    Seems to be working well...

"""

# %% Two-phase approach [part 1 - with all muscles]

"""

This section runs a similar tracking simulation to above, but splits into two phases
to ensure symmetry and each step and hopefully improve the smoothness of the motion.

TODO:
    > Could add periodicity to muscle activations so that the activations are symmetrical
      and hence useable in the initial guess
    > Link up initial experimental marker tracking sim to this approach with all
      muscles and stance phase so that it is consistent?

"""

#Create folder to store data
try:
    os.mkdir('..\\fullTrackingSim_rightLimb_allMuscles')
except:
    print('Full tracking sim (right limb all muscles) folder already detected...')
    
#Edit the model for use in the tracking sim tool
modelProcessor = osim.ModelProcessor('..\\data\\JA1_SCALED_Osim40_Muscles.osim')

#Append the necessary operators to the model
#Handles converting muscles to DeGrooteFregly and setting parameters
modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(2.0))
modelProcessor.append(osim.ModOpScaleMaxIsometricForce(2))
# modelProcessor.append(osim.ModOpTendonComplianceDynamicsModeDGF('implicit'))

#Process model for further edits
osimModel = modelProcessor.process()

#Remove left limb muscle and upper body forces to be replaced by torque actuators
#Set a list of forces to remove
removeForceInd = []
#Loop through forces and identify muscles to remove
for forceInd in range(osimModel.updForceSet().getSize()):
    #Check for muscle
    if osimModel.updForceSet().get(forceInd).getConcreteClassName().endswith('Muscle'):
        #Check for upper body
        if osimModel.updForceSet().get(forceInd).getName().split('_')[0] in ['extobl', 'intobl', 'ercspn']:
            #Append index to list
            removeForceInd.append(forceInd)

#Remove the designated forces keeping in mind that the index reduces each time
#another force is removed
for removeInd in removeForceInd:
    osimModel.updForceSet().remove(removeInd - removeForceInd.index(removeInd))

#Set the coordinates that will need torque actuation
#Set an optimal force and label for each actuator
#Upper body and left side are idealised actuators
#Residuals and reserve actuators are lowly weighted for optimal force
optForces = {
    #upper body
    'lumbar_extension': [1000, 'actuator'], 'lumbar_bending': [1000, 'actuator'], 'lumbar_rotation': [1000, 'actuator'],
    'arm_flex_r': [300, 'actuator'], 'arm_add_r': [300, 'actuator'], 'arm_rot_r': [300, 'actuator'],
    'elbow_flex_r': [100, 'actuator'], 'pro_sup_r': [100, 'actuator'],
    'arm_flex_l': [300, 'actuator'], 'arm_add_l': [300, 'actuator'], 'arm_rot_l': [300, 'actuator'],
    'elbow_flex_l': [100, 'actuator'], 'pro_sup_l': [100, 'actuator'],
    #left limb
    'hip_flexion_l': [2, 'reserve'], 'hip_adduction_l': [2, 'reserve'], 'hip_rotation_l': [2, 'reserve'],
    'knee_angle_l': [2, 'reserve'], 'ankle_angle_l': [2, 'reserve'],  'mtp_angle_l': [2, 'reserve'],
    #right limb
    'hip_flexion_r': [2, 'reserve'], 'hip_adduction_r': [2, 'reserve'], 'hip_rotation_r': [2, 'reserve'],
    'knee_angle_r': [2, 'reserve'], 'ankle_angle_r': [2, 'reserve'],  'mtp_angle_r': [2, 'reserve'],
    #pelvis
    'pelvis_tx': [1, 'residual'], 'pelvis_ty': [1, 'residual'], 'pelvis_tz': [1, 'residual'],
    'pelvis_tilt': [1, 'residual'], 'pelvis_list': [1, 'residual'], 'pelvis_rotation': [1, 'residual']}

#Add torque actuators to model
for coordForce in optForces.keys():
    #Create the actuator
    actu = osim.CoordinateActuator()
    actu.setName(f'{coordForce}_{optForces[coordForce][1]}')
    actu.setCoordinate(osimModel.getCoordinateSet().get(coordForce))
    actu.setOptimalForce(optForces[coordForce][0])
    actu.setMinControl(np.inf*-1)
    actu.setMaxControl(np.inf)
    #Add to the models force set
    osimModel.updForceSet().cloneAndAppend(actu)

#Increase the maximum contraction velocity of muscles
maxContractionVelocity = 30
for muscleInd in range(osimModel.getMuscles().getSize()):
    #Contraction velocity
    osimModel.getMuscles().get(muscleInd).set_max_contraction_velocity(maxContractionVelocity)
    # #Tendon compliance
    # if '_gas_' in osimModel.getMuscles().get(muscleInd).getName() or 'soleus_' in osimModel.getMuscles().get(muscleInd).getName():
    #     osimModel.getMuscles().get(muscleInd).set_ignore_tendon_compliance(False)
        
#Add contact spheres at foot-ground contact model locations

#Get the markers that contain an fp reference for contact sphere locations
footGroundNames = []
for markerInd in range(osimModel.updMarkerSet().getSize()):
    #Check for fp indicator in marker
    if '_fp_' in osimModel.updMarkerSet().get(markerInd).getName():
        footGroundNames.append(osimModel.updMarkerSet().get(markerInd).getName())
        
#Get foot ground name locations in dictionary
footGroundLocs = {name: np.array((osimModel.updMarkerSet().get(name).get_location().get(0),
                                  osimModel.updMarkerSet().get(name).get_location().get(1),
                                  osimModel.updMarkerSet().get(name).get_location().get(2))) for name in footGroundNames}

#Convert to a singular heel and toe foot ground contact locations
#### NOTE: only heel and mtp in foot ground spheres currently...
footSphereLocs = {}
for name in footGroundLocs:
    if 'heel1' in name:
        #Set the other name of the heel location
        name2 = name.replace('1','2')
        #Get the midpoint of two locations
        mpLoc = np.array(((footGroundLocs[name][0]+footGroundLocs[name2][0])/2,
                          (footGroundLocs[name][1]+footGroundLocs[name2][1])/2,
                          (footGroundLocs[name][2]+footGroundLocs[name2][2])/2))
        #Set in dictionary
        footSphereLocs['heel_'+name[0]] = mpLoc
    # elif 'mt' in name or 'toe' in name:
    elif 'mt' in name:
        #Set the name and location as what is already done
        footSphereLocs[name.split('_')[-1]+'_'+name[0]] = footGroundLocs[name]
        
#Set foot ground contact sphere sizes
heelSphereRadius = 0.04
otherSphereRadius = 0.03

#Add sphere size into dictionary
footSphereSizes = {name: heelSphereRadius if 'heel' in name else otherSphereRadius for name in footSphereLocs.keys()}

#Create and connect the half space floor to ground
floorContact = osim.ContactHalfSpace() #create object
floorContact.setName('floor') #set name
floorContact.setOrientation(osim.Vec3(0, 0, -1.5707963267949001)) #set orientation
floorContact.connectSocket_frame(osim.PhysicalFrame.safeDownCast(osimModel.getGround())) #connect to ground
osimModel.addContactGeometry(floorContact) #add to model

#Iteratively create the contact geometry and spheres and attach to models force set
for contactPoint in footSphereLocs.keys():
    
    #Create the contact geometry and set parameters
    contactGeom = osim.ContactSphere() #create object
    contactGeom.setName(contactPoint) #set name
    
    #Conditional for heel or toe sphere
    if contactPoint.startswith('heel') and contactPoint.endswith('_r'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('calcn_r')) #connect to frame
    elif contactPoint.startswith('heel') and contactPoint.endswith('_l'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('calcn_l')) #connect to frame
    elif not contactPoint.startswith('heel') and contactPoint.endswith('_r'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('toes_r')) #connect to frame
    elif not contactPoint.startswith('heel') and contactPoint.endswith('_l'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('toes_l')) #connect to frame
        
    #Get sphere size to set and adjust location
    currSphereSize = footSphereSizes[contactPoint]    
        
    #Conditional whether to transform point location
    if 'heel' in contactPoint:
        contactGeom.setLocation(osim.Vec3(footSphereLocs[contactPoint][0],
                                          # footGroundLocs[contactPoint][1],
                                          footSphereLocs[contactPoint][1] + (currSphereSize/2), #adjust to bottom of foot height
                                          footSphereLocs[contactPoint][2])) #set location
    else:
        
        #Get translation from calcn to toe joint and subtract from location
        if contactPoint.endswith('r_'):
            jointTranslation = osimModel.updJointSet().get('mtp_r').get_frames(0).get_translation()
        else:
            jointTranslation = osimModel.updJointSet().get('mtp_l').get_frames(0).get_translation()
        contactGeom.setLocation(osim.Vec3(footSphereLocs[contactPoint][0] - jointTranslation.get(0),
                                          # footSphereLocs[contactPoint][1] - jointTranslation.get(1),
                                          footSphereLocs[contactPoint][1] - jointTranslation.get(1) + (currSphereSize/2), #adjust to bottom of foot height
                                          footSphereLocs[contactPoint][2] - jointTranslation.get(2))) #set location
    contactGeom.setRadius(currSphereSize) #set radius
    osimModel.addContactGeometry(contactGeom) #add to model
    
    #Create the sphere and set properties
    contactSphere = osim.SmoothSphereHalfSpaceForce() #create force
    contactSphere.setName('contact_'+contactPoint) #set name
    contactSphere.connectSocket_half_space(osimModel.updContactGeometrySet().get('floor')) #connect to floor
    contactSphere.connectSocket_sphere(osimModel.updContactGeometrySet().get(contactPoint)) #connect to sphere
    contactSphere.set_stiffness(3067776) #set stiffness
    contactSphere.set_dissipation(2) #set dissipation
    contactSphere.set_static_friction(0.8) #set static friction
    contactSphere.set_dynamic_friction(0.8) #set dynamic friction
    contactSphere.set_viscous_friction(0.5) #set viscous friction
    contactSphere.set_transition_velocity(0.2) #set transition velocity
    contactSphere.set_hertz_smoothing(300) #set hertz smoothing
    contactSphere.set_hunt_crossley_smoothing(50) #set hunt crossley smoothing
    osimModel.addForce(contactSphere) #add to model    

#Finalise model connections
osimModel.finalizeConnections()

#Save model to file for later use if needed
osimModel.printToXML('..\\fullTrackingSim_rightLimb_allMuscles\\fullTrackingSim_rightLimb_allMuscles.osim')

#Create the tracking model processor for use in the tool
trackModelProcessor = osim.ModelProcessor(osimModel)

#Set joints to weld for simulation
jointsToWeld = ['subtalar_r', 'subtalar_l',
                # 'mtp_r', 'mtp_l',
                'radius_hand_r', 'radius_hand_l']
weldVectorStr = osim.StdVectorString()
[weldVectorStr.append(joint) for joint in jointsToWeld]

#Add remaining options to model processor
#Joints to weld
trackModelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName('fullTrackingSim_rightLimb_allMuscles')

#Set model in tracking tool
track.setModel(trackModelProcessor)

#Set the coordinates reference in tracking tool

#Convert the solution data from the tracking sim to coordinates
trackingTraj = osim.MocoTrajectory('..\\expTrackingSim\\expTrackingSim_markerTrackingSolution.sto')

#Export to states table and write to file
trackingStates = trackingTraj.exportToStatesTable()
osim.STOFileAdapter().write(trackingStates, '..\\fullTrackingSim_rightLimb\\coordinatesToTrack.sto')

#Set coordinates file as reference in tool
track.setStatesReference(osim.TableProcessor('..\\fullTrackingSim_rightLimb\\coordinatesToTrack.sto'))
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
# taskWeights = {'pelvis_tx': 1e-0, 'pelvis_ty': 1e-1, 'pelvis_tz': 1e-0,
#                'pelvis_tilt': 1e-0, 'pelvis_list': 1e-0, 'pelvis_rotation': 1e-0, 
#                'hip_flexion_r': 1e-0, 'hip_adduction_r': 1e-1, 'hip_rotation_r': 1e-2, 
#                'knee_angle_r': 1e-0, 'ankle_angle_r': 1e-2,
#                # 'subtalar_angle_r': 0,
#                'mtp_angle_r': 1e-2,
#                'hip_flexion_l': 1e-0, 'hip_adduction_l': 1e-1, 'hip_rotation_l': 1e-2, 
#                'knee_angle_l': 1e-0, 'ankle_angle_l': 1e-2,
#                # 'subtalar_angle_l': 0,
#                'mtp_angle_l': 1e-2,
#                'lumbar_extension': 1e-1, 'lumbar_bending': 1e-2, 'lumbar_rotation': 1e-2,
#                'arm_flex_r': 1e-1, 'arm_add_r': 1e-1, 'arm_rot_r': 1e-1,
#                'elbow_flex_r': 1e-1, 'pro_sup_r': 1e-1,
#                'arm_flex_l': 1e-1, 'arm_add_l': 1e-1, 'arm_rot_l': 1e-1,
#                'elbow_flex_l': 1e-1, 'pro_sup_l': 1e-1}

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
#Get the gait timings from helper function for the half gait cycle
startTime, endTime = helper.getGaitTimings(grfFile = '..\\data\\sprint_grf.mot',
                                            extLoads = '..\\data\\sprint_grf.xml',
                                            startForceName = 'RightGRF1',
                                            stopForceName = 'LeftGRF1',
                                            forceThreshold = 20)

#Set in tracking problem
track.set_initial_time(startTime)
track.set_final_time(endTime)

#Set mesh interval
#### TODO: consider calculating this more objectively
#### Set as mesh intervals later
# track.set_mesh_interval(0.02)

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
effort.setWeightForControlPattern('/forceset/.*_actuator', 0.01)
#Reserves
effort.setWeightForControlPattern('/forceset/.*_reserve', 10)
#Residuals
effort.setWeightForControlPattern('/forceset/.*_residual', 10) #lower than 100 previously...

#Add contact tracking goal
#Uses three separate vectors to apply different contact tracking weights

#Set right and left contact sphere groups
rightFootContacts = [f'/forceset/contact_{pointName}' for pointName in footSphereLocs.keys() if pointName.endswith('_r')]
leftFootContacts = [f'/forceset/contact_{pointName}' for pointName in footSphereLocs.keys() if pointName.endswith('_l')]

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
                            'weight': [1e-2, 1e-2, 1e-2],
                            # 'weight': [1, 0.5, 1],
                            # 'weight': [0.25, 0.1, 0.25],
                            }

#Loop through contact tracking settings to create separate weighted goals
for ii in range(len(contactTrackingSettings['vector'])):
    
    #Create tracking goal
    contactTracking.append(osim.MocoContactTrackingGoal(f'contact{ii}',
                                                        contactTrackingSettings['weight'][ii]))
    
    #Set external loads
    contactTracking[ii].setExternalLoadsFile('..\\data\\sprint_grf.xml')
    
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

#Opposite states are periodic except pelvis anterior-posterior translation
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

#Add to problem
problem.addGoal(periodicityGoal)
    
#Set kinematic bounds in problem

#### TODO: should do the same thing with final bounds...not if periodic constraint added

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
#Here we allow a 10% window either side of the starting value

#Loop through coordinates
for coordName in boundedCoords:
    
    #Get the full path to coordinate
    coordPath = osimModel.updCoordinateSet().get(coordName).getAbsolutePathString()
    
    #Get minimum and maximum value from the tracking states
    #Calculate range
    minVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy().min()
    maxVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy().max()
    rangeVal = np.diff((minVal,maxVal))[0]
    
    #Get the starting and final value from the tracked states
    initialVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy()[0]
    finalVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy()[-1]
    
    #Calculate total and initial bound ranges
    #Calculate and set the 20% range either side of the min and maximum
    #Set the initial value to be within 5% of the starting value
    totalBounds = [minVal - (rangeVal * 0.25), maxVal + (rangeVal * 0.25)]
    initialBounds = [initialVal - (np.abs(initialVal) * 0.1), initialVal + (np.abs(initialVal) * 0.1)]
    finalBounds = [finalVal - (np.abs(finalVal) * 0.1), finalVal + (np.abs(finalVal) * 0.1)]
    
    #Set bounds in problem
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
        
    #Set bounds in problem
    #Put check in place to correct initial bounds if outside of total ranges
    correctedFinalBounds = []
    #Lower bound
    if finalBounds[0] > totalBounds[0]:
        correctedFinalBounds.append(finalBounds[0])
    else:
        correctedFinalBounds.append(totalBounds[0])
    #Upper bound
    if finalBounds[1] < totalBounds[1]:
        correctedFinalBounds.append(finalBounds[1])
    else:
        correctedFinalBounds.append(totalBounds[1])
    
    #Set in problem
    problem.setStateInfo(coordPath+'/value',
                         #Total bounds range
                         totalBounds,
                         #Initial bounds range
                         #correctedInitialBounds, ### not if using periodic constraint?
                         #Final bounds range
                         #correctedFinalBounds #### not if using a periodic constraint
                         )
    
#Define the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_optim_constraint_tolerance(0.01)
solver.set_optim_convergence_tolerance(0.01)
solver.set_multibody_dynamics_mode('explicit')
# solver.set_minimize_implicit_multibody_accelerations(True) #smoothness criterion?
# solver.set_implicit_multibody_accelerations_weight(0.0001)
solver.set_num_mesh_intervals(10) #### TODO: coarse to start with?
solver.resetProblem(problem)

# #Set guess from tracking simulation
# # solver.setGuessFile('..\\expTrackingSim\\expTrackingSim_markerTrackingSolution.sto')
# solver.setGuessFile('..\\fullTrackingSim\\fullTrackingSim_statesTrackingSolution_coarse.sto')
# solver.resetProblem(problem)

#Solve the tracking problem
solution = study.solve()
# study.visualize(solution)

#Copy tracked states file over to main directory
shutil.move('fullTrackingSim_rightLimb_allMuscles_tracked_states.sto',
            '..\\fullTrackingSim_rightLimb_allMuscles\\fullTrackingSim_rightLimb_allMuscles_tracked_states.sto')

#Save solution to file
solution.write('..\\fullTrackingSim_rightLimb_allMuscles\\fullTrackingSim_rightLimb_allMuscles_statesTrackingSolution.sto')

#Extract predicted GRFs

#Create forces table
externalForcesTableFlat = osim.createExternalLoadsTableForGait(trackModelProcessor.process(), solution,
                                                               forcesRightFoot, forcesLeftFoot)

#Write table to file
osim.STOFileAdapter().write(externalForcesTableFlat,
                            '..\\fullTrackingSim_rightLimb_allMuscles\\fullTrackingSim_rightLimb_allMuscles_statesTrackingSolution_grf.sto')

"""

NOTES: 
    
    A solution with all muscles can be obtained
    Follow-up step is to develop a mirrored full cycle replicant
    
    If we want to set muscle states then would need periodicity on activations?

"""

#Create folder to store data
try:
    os.mkdir('..\\fullTrackingSim_symmetrical')
except:
    print('Full tracking sim (symmetrical) folder already detected...')
    
#Convert solution into a full gait cycle that mirrors the first stride
rightLimbSolution = osim.MocoTrajectory('..\\fullTrackingSim_rightLimb_allMuscles\\fullTrackingSim_rightLimb_allMuscles_statesTrackingSolution.sto')

#Get the time stamps from the right limb solution
rightLimbTime = rightLimbSolution.getTime().to_numpy()

#Create a new time stamp that doubles the length and time
fullTime = np.linspace(rightLimbTime[0], rightLimbTime[-1] + (rightLimbTime[-1] - rightLimbTime[0]), len(rightLimbTime)*2 - 1)

#Double number of time stamps to fill inverted states
rightLimbSolution.resampleWithNumTimes(len(fullTime))

#Set new time in trajectory
rightLimbSolution.setTime(fullTime)

#Create dictionary to store adjusted data in so that it isn't edited progressively
newCoordinateVals = {}

#Adjust coordinate values in solution
#This includes adding the opposite right/left side together, adding symmetrical
#coordinates together, and inverting negated coordinates
stateNames = rightLimbSolution.exportToStatesTable().getColumnLabels()
for currState in stateNames:
    
    #Right hand side
    if currState.endswith('_r/value'): ### or currState.endswith('_r/activation'):
        
        #Get the state values
        stateVals = rightLimbSolution.getState(currState).to_numpy()
        
        #Get the alternate sided state values
        stateVals_opp = rightLimbSolution.getState(currState.replace('_r/', '_l/')).to_numpy()
        
        #Join together while ignoring the first/end shared value
        joinedStateVals = np.concatenate((stateVals, stateVals_opp[1:]))
        
        #Resample to number of times in solution
        resampledStateVals = np.interp(fullTime, 
                                       np.linspace(rightLimbTime[0], rightLimbTime[-1] + (rightLimbTime[-1] - rightLimbTime[0]), len(joinedStateVals)),
                                       joinedStateVals)
        
        #Set new values in dictionary
        newCoordinateVals[currState] = resampledStateVals
        
    #Left hand side
    elif currState.endswith('_l/value'): ### or currState.endswith('_r/activation'):
        
        #Get the state values
        stateVals = rightLimbSolution.getState(currState).to_numpy()
        
        #Get the alternate sided state values
        stateVals_opp = rightLimbSolution.getState(currState.replace('_l/', '_r/')).to_numpy()
        
        #Join together while ignoring the first/end shared value
        joinedStateVals = np.concatenate((stateVals, stateVals_opp[1:]))
        
        #Resample to number of times in solution
        resampledStateVals = np.interp(fullTime, 
                                       np.linspace(rightLimbTime[0], rightLimbTime[-1] + (rightLimbTime[-1] - rightLimbTime[0]), len(joinedStateVals)),
                                       joinedStateVals)
        
        #Set new values in dictionary
        newCoordinateVals[currState] = resampledStateVals
        
    #Symmetrical coordinates
    elif currState.endswith('/value') and currState.split('/')[3] in ['pelvis_ty', 'pelvis_tz', 'pelvis_tilt', 'lumbar_extension']:
        
        #Get the state values
        stateVals = rightLimbSolution.getState(currState).to_numpy()
        
        #Join together while ignoring the first/end shared value
        joinedStateVals = np.concatenate((stateVals, stateVals[1:]))
        
        #Resample to number of times in solution
        resampledStateVals = np.interp(fullTime, 
                                       np.linspace(rightLimbTime[0], rightLimbTime[-1] + (rightLimbTime[-1] - rightLimbTime[0]), len(joinedStateVals)),
                                       joinedStateVals)
        
        #Set new values in dictionary
        newCoordinateVals[currState] = resampledStateVals
        
    #Negated coordinates
    elif currState.endswith('/value') and currState.split('/')[3] in ['lumbar_bending', 'lumbar_rotation', 'pelvis_list', 'pelvis_rotation']:
        
        #Get the state values
        stateVals = rightLimbSolution.getState(currState).to_numpy()
        
        #Get the flipped version
        stateVals_opp = stateVals[::-1]
        
        #Join together while ignoring the first/end shared value
        joinedStateVals = np.concatenate((stateVals, stateVals_opp[1:]))
        
        #Resample to number of times in solution
        resampledStateVals = np.interp(fullTime, 
                                       np.linspace(rightLimbTime[0], rightLimbTime[-1] + (rightLimbTime[-1] - rightLimbTime[0]), len(joinedStateVals)),
                                       joinedStateVals)
        
        #Set new values in dictionary
        newCoordinateVals[currState] = resampledStateVals
        
    #Special case for pelvis translation
    elif currState.endswith('_tx/value'):        
        
        #Get state values
        stateVals = rightLimbSolution.getState(currState).to_numpy()
        
        #Get additional values
        stateVals_t = np.zeros(len(stateVals))
        for ii in range(len(stateVals)):
            stateVals_t[ii] = stateVals[-1] + (stateVals[ii] - stateVals[0])
        
        #Join together while ignoring the first/end shared value
        joinedStateVals = np.concatenate((stateVals, stateVals_t[1:]))
        
        #Resample to number of times in solution
        resampledStateVals = np.interp(fullTime, 
                                       np.linspace(rightLimbTime[0], rightLimbTime[-1] + (rightLimbTime[-1] - rightLimbTime[0]), len(joinedStateVals)),
                                       joinedStateVals)
        
        #Set new values in dictionary
        newCoordinateVals[currState] = resampledStateVals

#Set edited state values in solution
for currState in list(newCoordinateVals.keys()):
    rightLimbSolution.setState(currState, newCoordinateVals[currState])

#Export to states table for selecting columns
statesTable = rightLimbSolution.exportToStatesTable()

#Remove any unwanted columns (forces, speeds)
for colLabel in statesTable.getColumnLabels():
    #Check for whether to remove column
    if colLabel.startswith('/forceset') or colLabel.endswith('/speed'):
        statesTable.removeColumn(colLabel)
        
#Write statescoordinates to file
osim.STOFileAdapter().write(statesTable, '..\\fullTrackingSim_symmetrical\\coordinatesToTrack.sto')

#Edit GRF file to create a symmetrical right-left force profile
grfTable = osim.TimeSeriesTable('..\\fullTrackingSim_rightLimb_allMuscles\\fullTrackingSim_rightLimb_allMuscles_statesTrackingSolution_grf.sto')

#Create new table to fill
newGrfTable = osim.TimeSeriesTable()

#Set the column labels
#Create a simplified set being just the force values
tableLabels = osim.StdVectorString()
for currLabel in grfTable.getColumnLabels():
    if '_force_' in currLabel and '_v' in currLabel:
        tableLabels.append(currLabel)
newGrfTable.setColumnLabels(tableLabels)

#Create an array that stores the column data
newForceVals = np.zeros((len(fullTime),len(newGrfTable.getColumnLabels())))
for currForce in newGrfTable.getColumnLabels():
    
    #Check for force that requires replication
    if '_force_' in currForce and currForce.split('_')[-1] in ['vx', 'vy']:
        
        #Get the force values
        forceVals = grfTable.getDependentColumn(currForce).to_numpy()
        
        #Get the alternate force values
        if '_force_r' in currForce:
            forceVals_opp = grfTable.getDependentColumn(currForce.replace('_r_', '_l_')).to_numpy()
        else:
            forceVals_opp = grfTable.getDependentColumn(currForce.replace('_l_', '_r_')).to_numpy()
        
        #Join together while ignoring the first/end shared value
        joinedForceVals = np.concatenate((forceVals, forceVals_opp[1:]))
        
        #Set new values in array
        newForceVals[:,list(newGrfTable.getColumnLabels()).index(currForce)] = joinedForceVals
        
    #Check for force that requires inverting
    if '_force_' in currForce and currForce.split('_')[-1] in ['vz']:
        
        #Get the force values
        forceVals = grfTable.getDependentColumn(currForce).to_numpy()
        
        #Get the alternate force values
        if '_force_r' in currForce:
            forceVals_opp = grfTable.getDependentColumn(currForce.replace('_r_', '_l_')).to_numpy()
        else:
            forceVals_opp = grfTable.getDependentColumn(currForce.replace('_l_', '_r_')).to_numpy()
        
        #Join together while ignoring the first/end shared value
        joinedForceVals = np.concatenate((forceVals, forceVals_opp[1:]*-1))
        
        #Set new values in array
        newForceVals[:,list(newGrfTable.getColumnLabels()).index(currForce)] = joinedForceVals

#Fill the table with rows by looping through the rows and columns
for iRow in range(len(fullTime)):
    #Create the current row
    osimRow = osim.RowVector.createFromMat(newForceVals[iRow,:])
    #Append to table
    newGrfTable.appendRow(iRow, osimRow)
    #Set time for current row
    newGrfTable.setIndependentValueAtIndex(iRow, fullTime[iRow])
    
#Write to mot file format
osim.STOFileAdapter().write(newGrfTable, '..\\fullTrackingSim_symmetrical\\grfToTrack.mot')

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
forceXML.printToXML('..\\fullTrackingSim_symmetrical\\grfToTrack.xml')

#Generate tracking simulation for symmetrical data

#Read in the previously created model for the tracking sim to edit
osimModel = osim.Model('..\\expTrackingSim\\expTrackingModel.osim')
    
#Add contact spheres at foot-ground contact model locations

#Get the markers that contain an fp reference for contact sphere locations
footGroundNames = []
for markerInd in range(osimModel.updMarkerSet().getSize()):
    #Check for fp indicator in marker
    if '_fp_' in osimModel.updMarkerSet().get(markerInd).getName():
        footGroundNames.append(osimModel.updMarkerSet().get(markerInd).getName())
        
#Get foot ground name locations in dictionary
footGroundLocs = {name: np.array((osimModel.updMarkerSet().get(name).get_location().get(0),
                                  osimModel.updMarkerSet().get(name).get_location().get(1),
                                  osimModel.updMarkerSet().get(name).get_location().get(2))) for name in footGroundNames}

#Convert to a singular heel and toe foot ground contact locations
#### NOTE: only heel and mtp in foot ground spheres currently...
footSphereLocs = {}
for name in footGroundLocs:
    if 'heel1' in name:
        #Set the other name of the heel location
        name2 = name.replace('1','2')
        #Get the midpoint of two locations
        mpLoc = np.array(((footGroundLocs[name][0]+footGroundLocs[name2][0])/2,
                          (footGroundLocs[name][1]+footGroundLocs[name2][1])/2,
                          (footGroundLocs[name][2]+footGroundLocs[name2][2])/2))
        #Set in dictionary
        footSphereLocs['heel_'+name[0]] = mpLoc
    # elif 'mt' in name or 'toe' in name:
    elif 'mt' in name:
        #Set the name and location as what is already done
        footSphereLocs[name.split('_')[-1]+'_'+name[0]] = footGroundLocs[name]
        
#Set foot ground contact sphere sizes
heelSphereRadius = 0.04
otherSphereRadius = 0.03

#Add sphere size into dictionary
footSphereSizes = {name: heelSphereRadius if 'heel' in name else otherSphereRadius for name in footSphereLocs.keys()}

#Create and connect the half space floor to ground
floorContact = osim.ContactHalfSpace() #create object
floorContact.setName('floor') #set name
floorContact.setOrientation(osim.Vec3(0, 0, -1.5707963267949001)) #set orientation
floorContact.connectSocket_frame(osim.PhysicalFrame.safeDownCast(osimModel.getGround())) #connect to ground
osimModel.addContactGeometry(floorContact) #add to model

#Iteratively create the contact geometry and spheres and attach to models force set
for contactPoint in footSphereLocs.keys():
    
    #Create the contact geometry and set parameters
    contactGeom = osim.ContactSphere() #create object
    contactGeom.setName(contactPoint) #set name
    
    #Conditional for heel or toe sphere
    if contactPoint.startswith('heel') and contactPoint.endswith('_r'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('calcn_r')) #connect to frame
    elif contactPoint.startswith('heel') and contactPoint.endswith('_l'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('calcn_l')) #connect to frame
    elif not contactPoint.startswith('heel') and contactPoint.endswith('_r'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('toes_r')) #connect to frame
    elif not contactPoint.startswith('heel') and contactPoint.endswith('_l'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('toes_l')) #connect to frame
        
    #Get sphere size to set and adjust location
    currSphereSize = footSphereSizes[contactPoint]    
        
    #Conditional whether to transform point location
    if 'heel' in contactPoint:
        contactGeom.setLocation(osim.Vec3(footSphereLocs[contactPoint][0],
                                          # footGroundLocs[contactPoint][1],
                                          footSphereLocs[contactPoint][1] + (currSphereSize/2), #adjust to bottom of foot height
                                          footSphereLocs[contactPoint][2])) #set location
    else:
        
        #Get translation from calcn to toe joint and subtract from location
        if contactPoint.endswith('r_'):
            jointTranslation = osimModel.updJointSet().get('mtp_r').get_frames(0).get_translation()
        else:
            jointTranslation = osimModel.updJointSet().get('mtp_l').get_frames(0).get_translation()
        contactGeom.setLocation(osim.Vec3(footSphereLocs[contactPoint][0] - jointTranslation.get(0),
                                          # footSphereLocs[contactPoint][1] - jointTranslation.get(1),
                                          footSphereLocs[contactPoint][1] - jointTranslation.get(1) + (currSphereSize/2), #adjust to bottom of foot height
                                          footSphereLocs[contactPoint][2] - jointTranslation.get(2))) #set location
    contactGeom.setRadius(currSphereSize) #set radius
    osimModel.addContactGeometry(contactGeom) #add to model
    
    #Create the sphere and set properties
    contactSphere = osim.SmoothSphereHalfSpaceForce() #create force
    contactSphere.setName('contact_'+contactPoint) #set name
    contactSphere.connectSocket_half_space(osimModel.updContactGeometrySet().get('floor')) #connect to floor
    contactSphere.connectSocket_sphere(osimModel.updContactGeometrySet().get(contactPoint)) #connect to sphere
    contactSphere.set_stiffness(3067776) #set stiffness
    contactSphere.set_dissipation(2) #set dissipation
    contactSphere.set_static_friction(0.8) #set static friction
    contactSphere.set_dynamic_friction(0.8) #set dynamic friction
    contactSphere.set_viscous_friction(0.5) #set viscous friction
    contactSphere.set_transition_velocity(0.2) #set transition velocity
    contactSphere.set_hertz_smoothing(300) #set hertz smoothing
    contactSphere.set_hunt_crossley_smoothing(50) #set hunt crossley smoothing
    osimModel.addForce(contactSphere) #add to model

#Finalise model connections
osimModel.finalizeConnections()

#Save model to file for later use if needed
osimModel.printToXML('..\\fullTrackingSim_symmetrical\\fullTrackingSim_symmetrical.osim')

#Create the tracking model processor for use in the tool
trackModelProcessor = osim.ModelProcessor(osimModel)

#Set joints to weld for simulation
jointsToWeld = ['subtalar_r', 'subtalar_l',
                # 'mtp_r', 'mtp_l',
                'radius_hand_r', 'radius_hand_l']
weldVectorStr = osim.StdVectorString()
[weldVectorStr.append(joint) for joint in jointsToWeld]

#Add remaining options to model processor
#Joints to weld
trackModelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName('fullTrackingSim_symmetrical')

#Set model in tracking tool
track.setModel(trackModelProcessor)

#Set the coordinates reference in tracking tool
track.setStatesReference(osim.TableProcessor('..\\fullTrackingSim_symmetrical\\coordinatesToTrack.sto'))
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
# taskWeights = {'pelvis_tx': 1e-0, 'pelvis_ty': 1e-1, 'pelvis_tz': 1e-0,
#                'pelvis_tilt': 1e-0, 'pelvis_list': 1e-0, 'pelvis_rotation': 1e-0, 
#                'hip_flexion_r': 1e-0, 'hip_adduction_r': 1e-1, 'hip_rotation_r': 1e-2, 
#                'knee_angle_r': 1e-0, 'ankle_angle_r': 1e-2,
#                # 'subtalar_angle_r': 0,
#                'mtp_angle_r': 1e-2,
#                'hip_flexion_l': 1e-0, 'hip_adduction_l': 1e-1, 'hip_rotation_l': 1e-2, 
#                'knee_angle_l': 1e-0, 'ankle_angle_l': 1e-2,
#                # 'subtalar_angle_l': 0,
#                'mtp_angle_l': 1e-2,
#                'lumbar_extension': 1e-1, 'lumbar_bending': 1e-2, 'lumbar_rotation': 1e-2,
#                'arm_flex_r': 1e-1, 'arm_add_r': 1e-1, 'arm_rot_r': 1e-1,
#                'elbow_flex_r': 1e-1, 'pro_sup_r': 1e-1,
#                'arm_flex_l': 1e-1, 'arm_add_l': 1e-1, 'arm_rot_l': 1e-1,
#                'elbow_flex_l': 1e-1, 'pro_sup_l': 1e-1}

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

#Set in tracking problem
track.set_initial_time(fullTime[0])
track.set_final_time(fullTime[-1])

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
effort.setWeightForControlPattern('/forceset/.*_actuator', 0.01)
#Reserves
effort.setWeightForControlPattern('/forceset/.*_reserve', 10)
#Residuals
effort.setWeightForControlPattern('/forceset/.*_residual', 10) #lower than 100 previously...

#Add contact tracking goal
#Uses three separate vectors to apply different contact tracking weights

#Set right and left contact sphere groups
rightFootContacts = [f'/forceset/contact_{pointName}' for pointName in footSphereLocs.keys() if pointName.endswith('_r')]
leftFootContacts = [f'/forceset/contact_{pointName}' for pointName in footSphereLocs.keys() if pointName.endswith('_l')]

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
                            'weight': [1e-2, 1e-2, 1e-2],
                            # 'weight': [1, 0.5, 1],
                            # 'weight': [0.25, 0.1, 0.25],
                            }

#Loop through contact tracking settings to create separate weighted goals
for ii in range(len(contactTrackingSettings['vector'])):
    
    #Create tracking goal
    contactTracking.append(osim.MocoContactTrackingGoal(f'contact{ii}',
                                                        contactTrackingSettings['weight'][ii]))
    
    #Set external loads
    contactTracking[ii].setExternalLoadsFile('..\\fullTrackingSim_symmetrical\\grfToTrack.xml')
    
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

#Opposite states are periodic except pelvis anterior-posterior translation
for coordName in taskWeights.keys():
    
    if coordName != 'pelvis_tx':
        
        #Get full path to state
        stateName = osimModel.updCoordinateSet().get(coordName).getAbsolutePathString()+'/value'
        
        #Set periodicity from start to finish for state
        periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair(stateName))

#Add to problem
problem.addGoal(periodicityGoal)
    
#Set kinematic bounds in problem

#### TODO: should do the same thing with final bounds...not if periodic constraint added

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
#Here we allow a 10% window either side of the starting value

#Loop through coordinates
for coordName in boundedCoords:
    
    #Get the full path to coordinate
    coordPath = osimModel.updCoordinateSet().get(coordName).getAbsolutePathString()
    
    #Get minimum and maximum value from the tracking states
    #Calculate range
    minVal = statesTable.getDependentColumn(coordPath+'/value').to_numpy().min()
    maxVal = statesTable.getDependentColumn(coordPath+'/value').to_numpy().max()
    rangeVal = np.diff((minVal,maxVal))[0]
    
    #Get the starting and final value from the tracked states
    initialVal = statesTable.getDependentColumn(coordPath+'/value').to_numpy()[0]
    
    #Calculate total and initial bound ranges
    #Calculate and set the 25% range either side of the min and maximum
    totalBounds = [minVal - (rangeVal * 0.25), maxVal + (rangeVal * 0.25)]
        
    #Set in problem
    problem.setStateInfo(coordPath+'/value',
                         #Total bounds range
                         totalBounds,
                         #Initial bounds range
                         initialVal
                         )
    
#Define the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_optim_constraint_tolerance(0.01)
solver.set_optim_convergence_tolerance(0.01)
solver.set_multibody_dynamics_mode('explicit')
# solver.set_minimize_implicit_multibody_accelerations(True) #smoothness criterion?
# solver.set_implicit_multibody_accelerations_weight(0.0001)
solver.set_num_mesh_intervals(20) #### TODO: coarse to start with?
solver.resetProblem(problem)

#Solve the tracking problem
solution = study.solve()
# study.visualize(solution)

"""

NOTES: need to finish up saving this solution and clean-up files...

Still think optimal force needs to be reduced on the torque actuators to allow
a smoother solution value...? Or an even lower value on it's weighting in function...

"""


# %% Two-phase approach [part 2]

"""

This section runs a similar tracking simulation to above, but splits into two phases
to ensure symmetry and each step and hopefully improve the smoothness of the motion.

TODO:
    > Should it be tracking the original trajectory or an inverted right limb solution?
        >> Attempts in here to do this as well
        >> If it's doing this does it make sense to track an inverted GRF too?
    > Or should it be just like the other solution and then stitch together to find a holistic solution?
        >> Currently testing this approach

"""

#Create folder to store data
try:
    os.mkdir('..\\fullTrackingSim_leftLimb')
except:
    print('Full tracking sim (left limb) folder already detected...')
    
#Read in the previously created model for the tracking sim to edit
osimModel = osim.Model('..\\expTrackingSim\\expTrackingModel.osim')
    
#Add contact spheres at foot-ground contact model locations

#Get the markers that contain an fp reference for contact sphere locations
footGroundNames = []
for markerInd in range(osimModel.updMarkerSet().getSize()):
    #Check for fp indicator in marker
    if '_fp_' in osimModel.updMarkerSet().get(markerInd).getName():
        footGroundNames.append(osimModel.updMarkerSet().get(markerInd).getName())
        
#Get foot ground name locations in dictionary
footGroundLocs = {name: np.array((osimModel.updMarkerSet().get(name).get_location().get(0),
                                  osimModel.updMarkerSet().get(name).get_location().get(1),
                                  osimModel.updMarkerSet().get(name).get_location().get(2))) for name in footGroundNames}

#Convert to a singular heel and toe foot ground contact locations
#### NOTE: only heel and mtp in foot ground spheres currently...
footSphereLocs = {}
for name in footGroundLocs:
    if 'heel1' in name:
        #Set the other name of the heel location
        name2 = name.replace('1','2')
        #Get the midpoint of two locations
        mpLoc = np.array(((footGroundLocs[name][0]+footGroundLocs[name2][0])/2,
                          (footGroundLocs[name][1]+footGroundLocs[name2][1])/2,
                          (footGroundLocs[name][2]+footGroundLocs[name2][2])/2))
        #Set in dictionary
        footSphereLocs['heel_'+name[0]] = mpLoc
    # elif 'mt' in name or 'toe' in name:
    elif 'mt' in name:
        #Set the name and location as what is already done
        footSphereLocs[name.split('_')[-1]+'_'+name[0]] = footGroundLocs[name]
        
#Set foot ground contact sphere sizes
heelSphereRadius = 0.04
otherSphereRadius = 0.03

#Add sphere size into dictionary
footSphereSizes = {name: heelSphereRadius if 'heel' in name else otherSphereRadius for name in footSphereLocs.keys()}

#Create and connect the half space floor to ground
floorContact = osim.ContactHalfSpace() #create object
floorContact.setName('floor') #set name
floorContact.setOrientation(osim.Vec3(0, 0, -1.5707963267949001)) #set orientation
floorContact.connectSocket_frame(osim.PhysicalFrame.safeDownCast(osimModel.getGround())) #connect to ground
osimModel.addContactGeometry(floorContact) #add to model

#Iteratively create the contact geometry and spheres and attach to models force set
for contactPoint in footSphereLocs.keys():
    
    #Create the contact geometry and set parameters
    contactGeom = osim.ContactSphere() #create object
    contactGeom.setName(contactPoint) #set name
    
    #Conditional for heel or toe sphere
    if contactPoint.startswith('heel') and contactPoint.endswith('_r'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('calcn_r')) #connect to frame
    elif contactPoint.startswith('heel') and contactPoint.endswith('_l'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('calcn_l')) #connect to frame
    elif not contactPoint.startswith('heel') and contactPoint.endswith('_r'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('toes_r')) #connect to frame
    elif not contactPoint.startswith('heel') and contactPoint.endswith('_l'):
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('toes_l')) #connect to frame
        
    #Get sphere size to set and adjust location
    currSphereSize = footSphereSizes[contactPoint]    
        
    #Conditional whether to transform point location
    if 'heel' in contactPoint:
        contactGeom.setLocation(osim.Vec3(footSphereLocs[contactPoint][0],
                                          # footGroundLocs[contactPoint][1],
                                          footSphereLocs[contactPoint][1] + (currSphereSize/2), #adjust to bottom of foot height
                                          footSphereLocs[contactPoint][2])) #set location
    else:
        
        #Get translation from calcn to toe joint and subtract from location
        if contactPoint.endswith('r_'):
            jointTranslation = osimModel.updJointSet().get('mtp_r').get_frames(0).get_translation()
        else:
            jointTranslation = osimModel.updJointSet().get('mtp_l').get_frames(0).get_translation()
        contactGeom.setLocation(osim.Vec3(footSphereLocs[contactPoint][0] - jointTranslation.get(0),
                                          # footSphereLocs[contactPoint][1] - jointTranslation.get(1),
                                          footSphereLocs[contactPoint][1] - jointTranslation.get(1) + (currSphereSize/2), #adjust to bottom of foot height
                                          footSphereLocs[contactPoint][2] - jointTranslation.get(2))) #set location
    contactGeom.setRadius(currSphereSize) #set radius
    osimModel.addContactGeometry(contactGeom) #add to model
    
    #Create the sphere and set properties
    contactSphere = osim.SmoothSphereHalfSpaceForce() #create force
    contactSphere.setName('contact_'+contactPoint) #set name
    contactSphere.connectSocket_half_space(osimModel.updContactGeometrySet().get('floor')) #connect to floor
    contactSphere.connectSocket_sphere(osimModel.updContactGeometrySet().get(contactPoint)) #connect to sphere
    contactSphere.set_stiffness(3067776) #set stiffness
    contactSphere.set_dissipation(2) #set dissipation
    contactSphere.set_static_friction(0.8) #set static friction
    contactSphere.set_dynamic_friction(0.8) #set dynamic friction
    contactSphere.set_viscous_friction(0.5) #set viscous friction
    contactSphere.set_transition_velocity(0.2) #set transition velocity
    contactSphere.set_hertz_smoothing(300) #set hertz smoothing
    contactSphere.set_hunt_crossley_smoothing(50) #set hunt crossley smoothing
    osimModel.addForce(contactSphere) #add to model

#Finalise model connections
osimModel.finalizeConnections()

#Save model to file for later use if needed
osimModel.printToXML('..\\fullTrackingSim_leftLimb\\fullTrackingModel.osim')

#Create the tracking model processor for use in the tool
trackModelProcessor = osim.ModelProcessor(osimModel)

#Set joints to weld for simulation
jointsToWeld = ['subtalar_r', 'subtalar_l',
                # 'mtp_r', 'mtp_l',
                'radius_hand_r', 'radius_hand_l']
weldVectorStr = osim.StdVectorString()
[weldVectorStr.append(joint) for joint in jointsToWeld]

#Add remaining options to model processor
#Joints to weld
trackModelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName('fullTrackingSim_leftLimb')

#Set model in tracking tool
track.setModel(trackModelProcessor)

#Set the coordinates reference in tracking tool

#### OPTION 1: ORIGINAL EXPERIMENTAL DATA TRACKING

# #Convert the solution data from the tracking sim to coordinates
# trackingTraj = osim.MocoTrajectory('..\\expTrackingSim\\expTrackingSim_markerTrackingSolution.sto')

# #Export to states table and write to file
# trackingStates = trackingTraj.exportToStatesTable()
# osim.STOFileAdapter().write(trackingStates, '..\\fullTrackingSim_leftLimb\\coordinatesToTrack.sto')

#### OPTION 2: INVERTED RIGHT LIMB SOLUTION

#Convert the solution data from the tracking sim to coordinates
rightLimbSolution = osim.MocoTrajectory('..\\fullTrackingSim_rightLimb\\fullTrackingSim_rightLimb_statesTrackingSolution.sto')

#Adjust pelvis_tx in trajectory for new starting point and consistent translation
pelvisTx = rightLimbSolution.getState('/jointset/ground_pelvis/pelvis_tx/value').to_numpy()
pelvisTx_new = np.zeros(len(pelvisTx))
for ii in range(len(pelvisTx)):
    if ii == 0:
        pelvisTx_new[ii] = pelvisTx[-1]
    else:
        pelvisTx_new[ii] = pelvisTx[-1] + (pelvisTx[ii] - pelvisTx[0])

#Set in trajectory
rightLimbSolution.setState('/jointset/ground_pelvis/pelvis_tx/value', pelvisTx_new)

#Export to states table for further editing
trackingStates = rightLimbSolution.exportToStatesTable()

#Remove any unwanted columns (forces, speeds)
for colLabel in trackingStates.getColumnLabels():
    #Check for whether to remove column
    if colLabel.startswith('/forceset') or colLabel.endswith('/speed'):
        trackingStates.removeColumn(colLabel)
        
#Get column labels to flip left to right sides
trackingColumnLabels = trackingStates.getColumnLabels()

#Create inverted column labels
invertedColumnLabels = osim.StdVectorString()
for colLabel in trackingColumnLabels:
    if colLabel.startswith('/jointset') and '_r/' in colLabel:
        #Invert the left to right and append
        invertedColumnLabels.append(colLabel.replace('_r/', '_l/'))
    elif colLabel.startswith('/jointset') and '_l/' in colLabel:
        #Invert the left to right and append
        invertedColumnLabels.append(colLabel.replace('_l/', '_r/'))
    else:
        #Append the original column label
        invertedColumnLabels.append(colLabel)
    
#Set the new column labels
trackingStates.setColumnLabels(invertedColumnLabels)

#Adjusted negated periodic joint coordinates
for colLabel in trackingColumnLabels:
    if colLabel.split('/')[3] in ['lumbar_bending', 'lumbar_rotation', 'pelvis_list', 'pelvis_rotation']:
        #Negate the values
        trackingStates.updDependentColumn(colLabel).negateInPlace()
        
#Update the time values to reflect the left limb start and end times

#Set start and end times in problem
#Get the gait timings from helper function for the half gait cycle
startTime, endTime = helper.getGaitTimings(grfFile = '..\\data\\sprint_grf.mot',
                                            extLoads = '..\\data\\sprint_grf.xml',
                                            startForceName = 'LeftGRF1',
                                            stopForceName = 'RightGRF2',
                                            forceThreshold = 20)

#Set appropriately spaced time values in updated states
newTime = np.linspace(startTime, endTime, len(trackingStates.getIndependentColumn()))
for timeStamp in newTime[::-1]:
    trackingStates.setIndependentValueAtIndex(list(newTime).index(timeStamp), timeStamp)

#Write updated coordinates to file
osim.STOFileAdapter().write(trackingStates, '..\\fullTrackingSim_leftLimb\\coordinatesToTrack.sto')

#Set coordinates file as reference in tool
track.setStatesReference(osim.TableProcessor('..\\fullTrackingSim_leftLimb\\coordinatesToTrack.sto'))
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
# taskWeights = {'pelvis_tx': 1e-0, 'pelvis_ty': 1e-1, 'pelvis_tz': 1e-0,
#                'pelvis_tilt': 1e-0, 'pelvis_list': 1e-0, 'pelvis_rotation': 1e-0, 
#                'hip_flexion_r': 1e-0, 'hip_adduction_r': 1e-1, 'hip_rotation_r': 1e-2, 
#                'knee_angle_r': 1e-0, 'ankle_angle_r': 1e-2,
#                # 'subtalar_angle_r': 0,
#                'mtp_angle_r': 1e-2,
#                'hip_flexion_l': 1e-0, 'hip_adduction_l': 1e-1, 'hip_rotation_l': 1e-2, 
#                'knee_angle_l': 1e-0, 'ankle_angle_l': 1e-2,
#                # 'subtalar_angle_l': 0,
#                'mtp_angle_l': 1e-2,
#                'lumbar_extension': 1e-1, 'lumbar_bending': 1e-2, 'lumbar_rotation': 1e-2,
#                'arm_flex_r': 1e-1, 'arm_add_r': 1e-1, 'arm_rot_r': 1e-1,
#                'elbow_flex_r': 1e-1, 'pro_sup_r': 1e-1,
#                'arm_flex_l': 1e-1, 'arm_add_l': 1e-1, 'arm_rot_l': 1e-1,
#                'elbow_flex_l': 1e-1, 'pro_sup_l': 1e-1}

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

#Set in tracking problem
track.set_initial_time(startTime)
track.set_final_time(endTime)

#Set mesh interval
#### TODO: consider calculating this more objectively
#### Set as mesh intervals later
# track.set_mesh_interval(0.02)

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
effort.setWeightForControlPattern('/forceset/.*_actuator', 0.01)
#Reserves
effort.setWeightForControlPattern('/forceset/.*_reserve', 10)
#Residuals
effort.setWeightForControlPattern('/forceset/.*_residual', 10) #lower than 100 previously...

#Add contact tracking goal
#Uses three separate vectors to apply different contact tracking weights

#Set right and left contact sphere groups
rightFootContacts = [f'/forceset/contact_{pointName}' for pointName in footSphereLocs.keys() if pointName.endswith('_r')]
leftFootContacts = [f'/forceset/contact_{pointName}' for pointName in footSphereLocs.keys() if pointName.endswith('_l')]

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
                            'weight': [1e-2, 1e-2, 1e-2],
                            # 'weight': [1, 0.5, 1],
                            # 'weight': [0.25, 0.1, 0.25],
                            }

#Loop through contact tracking settings to create separate weighted goals
for ii in range(len(contactTrackingSettings['vector'])):
    
    #Create tracking goal
    contactTracking.append(osim.MocoContactTrackingGoal(f'contact{ii}',
                                                        contactTrackingSettings['weight'][ii]))
    
    #Set external loads
    contactTracking[ii].setExternalLoadsFile('..\\data\\sprint_grf.xml')
    
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

#Opposite states are periodic except pelvis anterior-posterior translation
for coordName in taskWeights.keys():
    
    if coordName != 'pelvis_tx':
        
        #Get full path to state
        stateName = osimModel.updCoordinateSet().get(coordName).getAbsolutePathString()+'/value'
        
        #Check for singular periodicity conditions
        if coordName in ['pelvis_ty', 'pelvis_tz', 'lumbar_extension']:
            
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

#Add to problem
problem.addGoal(periodicityGoal)
    
#Set kinematic bounds in problem

#### OPTION 1: STANDARD LIKE RIGHT LIMB TRACKING ORIGINAL SOLUTION

#### TODO: should do the same thing with final bounds...not if periodic constraint added

# #Set the joint coordinates to place bounds on
# boundedCoords = ['pelvis_tx', 'pelvis_ty', 'pelvis_tz', 'pelvis_tilt', 'pelvis_list', 'pelvis_rotation', 
#                  'hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r', 
#                  'knee_angle_r', 'ankle_angle_r', 
#                  'mtp_angle_r', # 'subtalar_angle_r',
#                  'hip_flexion_l', 'hip_adduction_l', 'hip_rotation_l', 
#                  'knee_angle_l', 'ankle_angle_l',
#                  'mtp_angle_l', # 'subtalar_angle_l',
#                  'lumbar_extension', 'lumbar_bending', 'lumbar_rotation',
#                  'arm_flex_r', 'arm_add_r', 'arm_rot_r',
#                  'elbow_flex_r', 'pro_sup_r',
#                  'arm_flex_l', 'arm_add_l', 'arm_rot_l',
#                  'elbow_flex_l', 'pro_sup_l']

# #Also create a dictionary which limits initial bounds to be relative to starting point
# #Here we allow a 10% window either side of the starting value

# #Loop through coordinates
# for coordName in boundedCoords:
    
#     #Get the full path to coordinate
#     coordPath = osimModel.updCoordinateSet().get(coordName).getAbsolutePathString()
    
#     #Get minimum and maximum value from the tracking states
#     #Calculate range
#     minVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy().min()
#     maxVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy().max()
#     rangeVal = np.diff((minVal,maxVal))[0]
    
#     #Get the starting and final value from the tracked states
#     initialVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy()[0]
#     finalVal = trackingStates.getDependentColumn(coordPath+'/value').to_numpy()[-1]
    
#     #Calculate total and initial bound ranges
#     #Calculate and set the 20% range either side of the min and maximum
#     #Set the initial value to be within 5% of the starting value
#     totalBounds = [minVal - (rangeVal * 0.25), maxVal + (rangeVal * 0.25)]
#     initialBounds = [initialVal - (np.abs(initialVal) * 0.1), initialVal + (np.abs(initialVal) * 0.1)]
#     finalBounds = [finalVal - (np.abs(finalVal) * 0.1), finalVal + (np.abs(finalVal) * 0.1)]
    
#     #Set bounds in problem
#     #Put check in place to correct initial bounds if outside of total ranges
#     correctedInitialBounds = []
#     #Lower bound
#     if initialBounds[0] > totalBounds[0]:
#         correctedInitialBounds.append(initialBounds[0])
#     else:
#         correctedInitialBounds.append(totalBounds[0])
#     #Upper bound
#     if initialBounds[1] < totalBounds[1]:
#         correctedInitialBounds.append(initialBounds[1])
#     else:
#         correctedInitialBounds.append(totalBounds[1])
        
#     #Set bounds in problem
#     #Put check in place to correct initial bounds if outside of total ranges
#     correctedFinalBounds = []
#     #Lower bound
#     if finalBounds[0] > totalBounds[0]:
#         correctedFinalBounds.append(finalBounds[0])
#     else:
#         correctedFinalBounds.append(totalBounds[0])
#     #Upper bound
#     if finalBounds[1] < totalBounds[1]:
#         correctedFinalBounds.append(finalBounds[1])
#     else:
#         correctedFinalBounds.append(totalBounds[1])
    
#     #Set in problem
#     problem.setStateInfo(coordPath+'/value',
#                          #Total bounds range
#                          totalBounds,
#                          #Initial bounds range
#                          #correctedInitialBounds, ### not if using periodic constraint?
#                          #Final bounds range
#                          #correctedFinalBounds #### not if using a periodic constraint
#                          )
    
#### OPTION 2: BASED ON END VALUES OF RIGHT LIMB TRACKING ORIGINAL SOLUTION

# #Get the right limb solution
# rightLimbSolution = osim.MocoTrajectory('..\\fullTrackingSim_rightLimb\\fullTrackingSim_rightLimb_statesTrackingSolution.sto')

# #Get state names
# stateNames = rightLimbSolution.getStateNames()

# #Loop through states and set initial bounds
# for currState in stateNames:
    
#     #Check if a kinematic state that we want to set bounds on
#     if currState.endswith('/value'):
        
#         #Get minimum and maximum value from the tracking states
#         #Calculate range
#         minVal = trackingStates.getDependentColumn(currState).to_numpy().min()
#         maxVal = trackingStates.getDependentColumn(currState).to_numpy().max()
#         rangeVal = np.diff((minVal,maxVal))[0]
        
#         #Calculate total bound ranges
#         totalBounds = [minVal - (rangeVal * 0.25), maxVal + (rangeVal * 0.25)]
        
#         #Get the starting value as the end value from the previous solution
#         #Set 25% bounds around this
#         initialVal = rightLimbSolution.getState(currState).to_numpy()[-1]
#         initialBounds = [initialVal - (np.abs(initialVal) * 0.25), initialVal + (np.abs(initialVal) * 0.25)]
        
#         #Put check in place to correct initial bounds if outside of total ranges
#         correctedInitialBounds = []
#         #Lower bound
#         if initialBounds[0] > totalBounds[0]:
#             correctedInitialBounds.append(initialBounds[0])
#         else:
#             correctedInitialBounds.append(totalBounds[0])
#         #Upper bound
#         if initialBounds[1] < totalBounds[1]:
#             correctedInitialBounds.append(initialBounds[1])
#         else:
#             correctedInitialBounds.append(totalBounds[1])
        
#         #Set in problem
#         problem.setStateInfo(currState,
#                              #Total bounds range
#                              totalBounds,
#                              #Initial bound value
#                              correctedInitialBounds)
    
#     else:
        
#         #Just set initial bound value
#         #Set 25% bounds around this
#         initialVal = rightLimbSolution.getState(currState).to_numpy()[-1]
#         initialBounds = [initialVal - (np.abs(initialVal) * 0.25), initialVal + (np.abs(initialVal) * 0.25)]
        
#         #Set in problem
#         problem.setStateInfo(currState,
#                              #Total bounds range
#                              [],
#                              #Initial bound value
#                              initialBounds)
        
# #Get control names        
# controlNames = rightLimbSolution.getControlNames()

# #Loop through controls and set bounds
# for currControl in controlNames:
    
#     #Get initial value as final value from previous solution
#     #Set 25% bounds around this
#     initialVal = rightLimbSolution.getControl(currControl).to_numpy()[-1]
#     initialBounds = [initialVal - (np.abs(initialVal) * 0.25), initialVal + (np.abs(initialVal) * 0.25)]
    
#     #Set in problem
#     problem.setControlInfo(currControl,
#                            #Total bounds range
#                            [],
#                            #Initial bound value
#                            initialBounds)
    
#### OPTION 3: TRACKING INVERTED RIGHT LIMB SOLUTION

#Get state names
stateNames = rightLimbSolution.getStateNames()

#Loop through states and set initial bounds
for currState in stateNames:
    
    #Check if a kinematic state that we want to set bounds on
    if currState.endswith('/value'):
        
        #Get minimum and maximum value from the tracking states
        #Calculate range
        minVal = trackingStates.getDependentColumn(currState).to_numpy().min()
        maxVal = trackingStates.getDependentColumn(currState).to_numpy().max()
        rangeVal = np.diff((minVal,maxVal))[0]
        
        #Calculate total bound ranges
        totalBounds = [minVal - (rangeVal * 0.25), maxVal + (rangeVal * 0.25)]
        
        #Get the starting value as the end value from the previous solution
        initialVal = rightLimbSolution.getState(currState).to_numpy()[-1]
        
        #Set in problem
        problem.setStateInfo(currState,
                              #Total bounds range
                              totalBounds,
                              #Initial bound value
                              initialVal)
    
    else:
        
        #Just set initial bound value
        initialVal = rightLimbSolution.getState(currState).to_numpy()[-1]
        
        #Set in problem
        problem.setStateInfo(currState,
                              #Total bounds range
                              [],
                              #Initial bound value
                              initialVal)
        
# #Get control names        
# controlNames = rightLimbSolution.getControlNames()

# #Loop through controls and set bounds
# for currControl in controlNames:
    
#     #Get initial value as final value from previous solution
#     initialVal = rightLimbSolution.getControl(currControl).to_numpy()[-1]
    
#     #Set in problem
#     problem.setControlInfo(currControl,
#                             #Total bounds range
#                             [],
#                             #Initial bound value
#                             initialVal)

#Define the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_optim_constraint_tolerance(0.01)
solver.set_optim_convergence_tolerance(0.01)
solver.set_multibody_dynamics_mode('explicit')
# solver.set_minimize_implicit_multibody_accelerations(True) #smoothness criterion?
# solver.set_implicit_multibody_accelerations_weight(0.0001)
solver.set_num_mesh_intervals(10) #### TODO: coarse to start with?
solver.resetProblem(problem)

# #Set guess from tracking simulation
# # solver.setGuessFile('..\\expTrackingSim\\expTrackingSim_markerTrackingSolution.sto')
# solver.setGuessFile('..\\fullTrackingSim\\fullTrackingSim_statesTrackingSolution_coarse.sto')
# solver.resetProblem(problem)

#Solve the tracking problem
solution = study.solve()
# study.visualize(solution)

#Copy tracked states file over to main directory
shutil.move('fullTrackingSim_leftLimb_tracked_states.sto',
            '..\\fullTrackingSim_leftLimb\\fullTrackingSim_leftLimb_tracked_states.sto')

#Save solution to file
solution.write('..\\fullTrackingSim_leftLimb\\fullTrackingSim_leftLimb_statesTrackingSolution222.sto')

#Extract predicted GRFs

#Create forces table
externalForcesTableFlat = osim.createExternalLoadsTableForGait(trackModelProcessor.process(), solution,
                                                                forcesRightFoot, forcesLeftFoot)

#Write table to file
osim.STOFileAdapter().write(externalForcesTableFlat,
                            '..\\fullTrackingSim_leftLimb\\fullTrackingSim_leftLimb_statesTrackingSolution_grf.sto')

"""

NOTES: 
    
    Inverted simulation not going well for some reason?
    
    Do we need to actually track this, rather than just moving to predictive that
    combines the inverted solution to the original as an initial guess?
        > This won't have a swing phase though, unless we keep the muscles on the
          right hand side for this...?

"""

# %% ----- End of processData.py ----- %% #