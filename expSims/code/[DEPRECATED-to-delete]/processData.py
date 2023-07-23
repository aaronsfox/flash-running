# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 20:25:31 2020

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This script processes the data provided by Dorn et al. (2012) using a tracking
    simulation approach to generate dynamically consistent sprinting data as a 
    feeder for further simulations.
    
"""

# %% Import packages

import opensim as osim
import osimHelper as helper
import numpy as np
import os

# %% Set-up

#Add OpenSim geometry path (weird issues with this on new laptop)
osim.ModelVisualizer.addDirToGeometrySearchPaths('C:\\OpenSim 4.3\\Geometry')

#Set-up initial logger for processing experimental data
osim.Logger.removeFileSink()
osim.Logger.addFileSink('processDataLog.log')

# %% Run a muscle-driven marker tracking problem

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName('muscleDrivenMarkerTracking')

#Construct a model processor for the tool
modelProcessor = osim.ModelProcessor('..\\data\\JA1_SCALED_Osim40_Muscles.osim')

#Append the necessary operators to the model
#Handles converting muscles to DeGrooteFregly and setting parameters
modelProcessor.append(osim.ModOpAddExternalLoads('..\\data\\sprint_grf.xml'))
modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
modelProcessor.append(osim.ModOpScaleMaxIsometricForce(2))

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
#Upper body and left side
torqueDrivenCoords = ['lumbar_extension', 'lumbar_bending', 'lumbar_rotation',
                      'arm_flex_r', 'arm_add_r', 'arm_rot_r',
                      'elbow_flex_r', 'pro_sup_r',
                      'arm_flex_l', 'arm_add_l', 'arm_rot_l',
                      'elbow_flex_l', 'pro_sup_l',
                      'hip_flexion_l', 'hip_adduction_l', 'hip_rotation_l',
                      'knee_angle_l', 'ankle_angle_l', 'mtp_angle_l']

#Create the optimal forces dictionary with a set torque value
torqueVal = 1000
optForces = {coord: torqueVal for coord in torqueDrivenCoords}

#Add torque actuators to model
osimModel = helper.addTorqueActuators(osimModel = osimModel, optForces = optForces,
                                      minControl = np.inf*-1, maxControl = np.inf)

#Increase the maximum contraction velocity of muscles
maxContractionVelocity = 30
for muscleInd in range(osimModel.getMuscles().getSize()):
    osimModel.getMuscles().get(muscleInd).set_max_contraction_velocity(maxContractionVelocity)
    
#Finalise model connections
osimModel.finalizeConnections()

#Set the list of joints to be welded and create array string
jointsToWeld = ['subtalar_r', 'subtalar_l',
                #'mtp_r', 'mtp_l',
                'radius_hand_r', 'radius_hand_l']
weldVectorStr = osim.StdVectorString()
[weldVectorStr.append(joint) for joint in jointsToWeld]

#### TODO: other individual changes to model?
    ##### > tendon compliance in plantarflexors?
    ##### > rolling surface constraint to eliminate ground projection?

#Create the tracking model processor
#Add the weld operator
trackModelProcessor = osim.ModelProcessor(osimModel)
trackModelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))

#Set model in tracking tool
track.setModel(modelProcessor)

#Set the markers reference from TRC file
#Note that this needs to be done using a flattened table due to a bug
#See: https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=1815&t=14612&p=43550&start=0&view=
#### TODO: consider filtering?
markerData = osim.TimeSeriesTableVec3('..\\data\\sprint.trc')
tableProcessor = osim.TableProcessor(markerData.flatten())
track.setMarkersReference(tableProcessor)

#Set allow unused references in case there are extra markers
track.set_allow_unused_references(True)

#Set the global markers tracking weight
track.set_markers_global_tracking_weight(10)

#Set a default tracking weight for now on all markers
##### TODO: optimise this a little better
markerWeights = osim.MocoWeightSet()
for marker in osim.TimeSeriesTableVec3('..\\data\\sprint.trc').getColumnLabels():
    markerWeights.cloneAndAppend(osim.MocoWeight(marker, 1))
    
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
track.set_mesh_interval(0.02)

#Initialise the Moco study and problem
study = track.initialize()
problem = study.updProblem()

#### TODO: any updates on goals in problem?
    ### > metabolics?

#Update the weight on the default control effort goal
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
effort.setWeight(0.001)

#Update individual weights in control effort goal to be relative to
#actual muscle and reserve actuator names
#Set appropriate patterns in the weight set
#Muscles
effort.setWeightForControlPattern('/forceset/.*/activation', 0.1)
#Actuators
effort.setWeightForControlPattern('/forceset/.*_actuator', 0.01)

#Solve the tracking problem
solution = study.solve()
study.visualize(solution)

#Remove tracked markers file
os.remove('muscleDrivenMarkerTracking_tracked_markers.sto')

#Save solution to file
solution.write('..\\data\\sprint_markerTrackingSolution_muscleDriven.sto')


# %% Run a torque-driven marker tracking problem

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName('torqueDrivenMarkerTracking')

#Read in model
osimModel = osim.Model('..\\data\\JA1_SCALED_Osim40.osim')

#Create a dictionary for adding coordinate actuators to the model
optForces = {}
#Skip over residual pelvis coordinates
for coordInd in range(osimModel.updCoordinateSet().getSize()):
    if not osimModel.updCoordinateSet().get(coordInd).getName().startswith('pelvis'):
        optForces[osimModel.updCoordinateSet().get(coordInd).getName()] = 300

#Add torque actuators to model
osimModel = helper.addTorqueActuators(osimModel = osimModel, optForces = optForces,
                                      minControl = np.inf*-1, maxControl = np.inf)

#Construct a model processor for the tool
modelProcessor = osim.ModelProcessor(osimModel)

#Add the external loads to the model processor
modelProcessor.append(osim.ModOpAddExternalLoads('..\\data\\sprint_grf.xml'))

#Set model in tracking tool
track.setModel(modelProcessor)

#Set the markers reference from TRC file
#Note that this needs to be done using a flattened table due to an error
#See: https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=1815&t=14612&p=43550&start=0&view=
#### TODO: consider filtering?
markerData = osim.TimeSeriesTableVec3('..\\data\\sprint.trc')
tableProcessor = osim.TableProcessor(markerData.flatten())
track.setMarkersReference(tableProcessor)

#Set allow unused references in case there are extra markers
track.set_allow_unused_references(True)

#Set the global markers tracking weight
track.set_markers_global_tracking_weight(10)

#Set a default tracking weight for now on all markers
##### TODO: optimise this a little better
markerWeights = osim.MocoWeightSet()
for marker in osim.TimeSeriesTableVec3('..\\data\\sprint.trc').getColumnLabels():
    markerWeights.cloneAndAppend(osim.MocoWeight(marker, 1))
    
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
track.set_mesh_interval(0.01)

#### TODO: or consider shifting to a MocoStudy to edit the solver with a bit more finer detail?

#Solve the tracking problem
solution = track.solve()

#Remove tracked markers file
os.remove('torqueDrivenMarkerTracking_tracked_markers.sto')

#Save solution to file
solution.write('..\\data\\sprint_markerTrackingSolution.sto')

# %% Run a torque-driven tracking simulation of markers & GRF data

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName('torqueDrivenMarkerGrfTracking')

#Read in model
osimModel = osim.Model('..\\data\\JA1_SCALED_Osim40.osim')

#Create a dictionary for adding coordinate actuators to the model
optForces = {}
#Skip over residual pelvis coordinates and welded wrist coordinates
for coordInd in range(osimModel.updCoordinateSet().getSize()):
    if osimModel.updCoordinateSet().get(coordInd).getName().startswith('pelvis') or osimModel.updCoordinateSet().get(coordInd).getName().startswith('wrist'):
        print(f'Skipping coordinate {osimModel.updCoordinateSet().get(coordInd).getName()}')
    else:
        optForces[osimModel.updCoordinateSet().get(coordInd).getName()] = 300

#Add torque actuators to model
osimModel = helper.addTorqueActuators(osimModel = osimModel, optForces = optForces,
                                      minControl = np.inf*-1, maxControl = np.inf)

#Append rolling surface constraints to model
#Load the dummy model with surface constraints (not accessible through Python API)
dummyModel = osim.Model('..\\data\\dummyModel_withSurfaceConstraints.osim')
#Clone and append surface constraints over to model
for constraintInd in range(dummyModel.updConstraintSet().getSize()):
    osimModel.updConstraintSet().cloneAndAppend(dummyModel.updConstraintSet().get(constraintInd))

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

#Create and connect the half space floor to ground
floorContact = osim.ContactHalfSpace() #create object
floorContact.setName('floor') #set name
floorContact.setOrientation(osim.Vec3(0, 0, -1.5707963267949001)) #set orientation
floorContact.connectSocket_frame(osim.PhysicalFrame.safeDownCast(osimModel.getGround())) #connect to ground
osimModel.addContactGeometry(floorContact) #add to model

#Iteratively create the contact geometry and spheres and attach to models force set
sphereRadius = 0.025
for contactPoint in footGroundNames:
    
    #Create the contact geometry and set parameters
    contactGeom = osim.ContactSphere() #create object
    contactGeom.setName(contactPoint) #set name
    #Conditional for heel or toe sphere
    if contactPoint.startswith('r_') and 'heel' in contactPoint:
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('calcn_r')) #connect to frame
    elif contactPoint.startswith('l_') and 'heel' in contactPoint:
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('calcn_l')) #connect to frame
    elif contactPoint.startswith('r_') and 'heel' not in contactPoint:
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('toes_r')) #connect to frame
    elif contactPoint.startswith('l_') and 'heel' not in contactPoint:
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('toes_l')) #connect to frame
    #Conditional whether to transform point location
    if 'heel' in contactPoint:
        contactGeom.setLocation(osim.Vec3(footGroundLocs[contactPoint][0],
                                          footGroundLocs[contactPoint][1],
                                          footGroundLocs[contactPoint][2])) #set location
    else:
        #Get translation from calcn to toe joint and subtract from location
        if contactPoint.startswith('r_'):
            jointTranslation = osimModel.updJointSet().get('mtp_r').get_frames(0).get_translation()
        else:
            jointTranslation = osimModel.updJointSet().get('mtp_l').get_frames(0).get_translation()
        contactGeom.setLocation(osim.Vec3(footGroundLocs[contactPoint][0] - jointTranslation.get(0),
                                          footGroundLocs[contactPoint][1] - jointTranslation.get(1),
                                          footGroundLocs[contactPoint][2] - jointTranslation.get(2))) #set location
    contactGeom.setRadius(sphereRadius) #set radius
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

#Set the list of joints to be welded and create array string
jointsToWeld = ['radius_hand_r', 'radius_hand_l']
weldVectorStr = osim.StdVectorString()
[weldVectorStr.append(joint) for joint in jointsToWeld]

#Construct a model processor for the tool
trackModelProcessor = osim.ModelProcessor(osimModel)

#Add the external loads to the model processor
trackModelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))

#Set model in tracking tool
track.setModel(trackModelProcessor)

#Set the markers reference from TRC file
#Note that this needs to be done using a flattened table due to an error
#See: https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=1815&t=14612&p=43550&start=0&view=
#### TODO: consider filtering?
markerData = osim.TimeSeriesTableVec3('..\\data\\sprint.trc')
tableProcessor = osim.TableProcessor(markerData.flatten())
track.setMarkersReference(tableProcessor)

#Set allow unused references in case there are extra markers
track.set_allow_unused_references(True)

#Set the global markers tracking weight
track.set_markers_global_tracking_weight(10)

#Set a default tracking weight for now on all markers
##### TODO: optimise this a little better
markerWeights = osim.MocoWeightSet()
for marker in osim.TimeSeriesTableVec3('..\\data\\sprint.trc').getColumnLabels():
    markerWeights.cloneAndAppend(osim.MocoWeight(marker, 1))
    
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
#### TODO: consider calculating this more objectively --- actually done with solver later...
# track.set_mesh_interval(0.01)

#Define the Moco study and problem
study = track.initialize()
problem = study.updProblem()

#Set the effort regularisation goal
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
effort.setWeight(0.001)

#Add contact tracking goal
#Uses three separate vectors to apply different contact tracking weights

#Set right and left contact sphere groups
rightFootContacts = [f'/forceset/contact_{pointName}' for pointName in footGroundNames if pointName.startswith('r_')]
leftFootContacts = [f'/forceset/contact_{pointName}' for pointName in footGroundNames if pointName.startswith('l_')]

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
                           'weight': [2, 1, 2],
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

#Define the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())

# #Set guess from earlier simulation
# solver.setGuessFile('..\\data\\sprint_markerTrackingSolution.sto')
# solver.resetProblem(problem)

#Alter the initial guess to bump up pelvis_ty to ensure contact spheres don't
#start by penetrating the ground
initialGuess = solver.getGuess()
initialGuess.setState('/jointset/ground_pelvis/pelvis_ty/value',
                      initialGuess.getState('/jointset/ground_pelvis/pelvis_ty/value').to_numpy() + 0.05)
solver.resetProblem(problem)
solver.setGuess(initialGuess)

#Set solver options
solver.set_optim_max_iterations(1000)
solver.set_num_mesh_intervals(50)
solver.set_minimize_implicit_multibody_accelerations(True)
solver.set_optim_constraint_tolerance(1e-2) #default
solver.set_optim_convergence_tolerance(1e-2) #default
solver.set_minimize_implicit_auxiliary_derivatives(True)
solver.set_implicit_auxiliary_derivatives_weight(0.00001)

#Solve!
solution = study.solve()
# study.visualize(solution)

#Remove tracked markers file
os.remove('torqueDrivenMarkerGrfTracking_tracked_markers.sto')

#Save solution to file (TODO: GRF data too...)
# solution.write('..\\data\\sprint_markerTrackingSolution.sto')

##### Above original plan not really working for some reason...
    ##### Remove contact goal to check if this is the problem...
        #### Didn't really help --- could be constraints?

# %% Run a torque-driven tracking simulation of IK & GRF data

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName('torqueDrivenIkGrfTracking')

#Read in model
osimModel = osim.Model('..\\data\\JA1_SCALED_Osim40.osim')

#Create a dictionary for adding coordinate actuators to the model
optForces = {}
#Skip over residual pelvis coordinates and welded wrist coordinates
for coordInd in range(osimModel.updCoordinateSet().getSize()):
    # if osimModel.updCoordinateSet().get(coordInd).getName().startswith('pelvis') or osimModel.updCoordinateSet().get(coordInd).getName().startswith('wrist'):
    if osimModel.updCoordinateSet().get(coordInd).getName().startswith('wrist'):
        print(f'Skipping coordinate {osimModel.updCoordinateSet().get(coordInd).getName()}')
    else:
        optForces[osimModel.updCoordinateSet().get(coordInd).getName()] = 300

#Add torque actuators to model
osimModel = helper.addTorqueActuators(osimModel = osimModel, optForces = optForces,
                                      minControl = np.inf*-1, maxControl = np.inf)

#Append rolling surface constraints to model
#Load the dummy model with surface constraints (not accessible through Python API)
dummyModel = osim.Model('..\\data\\dummyModel_withSurfaceConstraints.osim')
#Clone and append surface constraints over to model
for constraintInd in range(dummyModel.updConstraintSet().getSize()):
    osimModel.updConstraintSet().cloneAndAppend(dummyModel.updConstraintSet().get(constraintInd))

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

#Create and connect the half space floor to ground
floorContact = osim.ContactHalfSpace() #create object
floorContact.setName('floor') #set name
floorContact.setOrientation(osim.Vec3(0, 0, -1.5707963267949001)) #set orientation
floorContact.connectSocket_frame(osim.PhysicalFrame.safeDownCast(osimModel.getGround())) #connect to ground
osimModel.addContactGeometry(floorContact) #add to model

#Iteratively create the contact geometry and spheres and attach to models force set
sphereRadius = 0.025
for contactPoint in footGroundNames:
    
    #Create the contact geometry and set parameters
    contactGeom = osim.ContactSphere() #create object
    contactGeom.setName(contactPoint) #set name
    #Conditional for heel or toe sphere
    if contactPoint.startswith('r_') and 'heel' in contactPoint:
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('calcn_r')) #connect to frame
    elif contactPoint.startswith('l_') and 'heel' in contactPoint:
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('calcn_l')) #connect to frame
    elif contactPoint.startswith('r_') and 'heel' not in contactPoint:
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('toes_r')) #connect to frame
    elif contactPoint.startswith('l_') and 'heel' not in contactPoint:
        contactGeom.connectSocket_frame(osimModel.updBodySet().get('toes_l')) #connect to frame
    #Conditional whether to transform point location
    if 'heel' in contactPoint:
        contactGeom.setLocation(osim.Vec3(footGroundLocs[contactPoint][0],
                                          footGroundLocs[contactPoint][1] + sphereRadius,
                                          footGroundLocs[contactPoint][2])) #set location
    else:
        #Get translation from calcn to toe joint and subtract from location
        if contactPoint.startswith('r_'):
            jointTranslation = osimModel.updJointSet().get('mtp_r').get_frames(0).get_translation()
        else:
            jointTranslation = osimModel.updJointSet().get('mtp_l').get_frames(0).get_translation()
        contactGeom.setLocation(osim.Vec3(footGroundLocs[contactPoint][0] - jointTranslation.get(0),
                                          (footGroundLocs[contactPoint][1] - jointTranslation.get(1)) + sphereRadius,
                                          footGroundLocs[contactPoint][2] - jointTranslation.get(2))) #set location
    contactGeom.setRadius(sphereRadius) #set radius
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

#Set the list of joints to be welded and create array string
jointsToWeld = ['radius_hand_r', 'radius_hand_l']
weldVectorStr = osim.StdVectorString()
[weldVectorStr.append(joint) for joint in jointsToWeld]

#Construct a model processor for the tool
trackModelProcessor = osim.ModelProcessor(osimModel)

#Add the external loads to the model processor
trackModelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))

#Set model in tracking tool
track.setModel(trackModelProcessor)
torqueModel = trackModelProcessor.process()

#Set the coordinates reference in tracking tool

#Convert the ik data to a coordinates file
helper.kinematicsToStates(kinematicsFileName = '..\\data\\ik.mot',
                          osimModelFileName = '..\\data\\JA1_SCALED_Osim40.osim',
                          outputFileName = '..\\data\\coordinates.sto',
                          inDegrees = True, outDegrees = False)

#Set coordinates file as reference in tool
track.setStatesReference(osim.TableProcessor('..\\data\\coordinates.sto'))
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
taskWeights = {'pelvis_tx': 500, 'pelvis_ty': 1, 'pelvis_tz': 250,
               'pelvis_tilt': 250, 'pelvis_list': 250, 'pelvis_rotation': 250, 
               'hip_flexion_r': 100, 'hip_adduction_r': 50, 'hip_rotation_r': 25, 
               'knee_angle_r': 100, 'ankle_angle_r': 50,
               'subtalar_angle_r': 10, 'mtp_angle_r': 10,
               'hip_flexion_l': 100, 'hip_adduction_l': 50, 'hip_rotation_l': 25, 
               'knee_angle_l': 100, 'ankle_angle_l': 50,
               'subtalar_angle_l': 10, 'mtp_angle_l': 10,
               'lumbar_extension': 100, 'lumbar_bending': 50, 'lumbar_rotation': 50,
               'arm_flex_r': 25, 'arm_add_r': 10, 'arm_rot_r': 10,
               'elbow_flex_r': 50, 'pro_sup_r': 10,
               'arm_flex_l': 25, 'arm_add_l': 10, 'arm_rot_l': 10,
               'elbow_flex_l': 50, 'pro_sup_l': 10}

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
#### TODO: consider calculating this more objectively --- actually done with solver later...
# track.set_mesh_interval(0.01)

#Define the Moco study and problem
study = track.initialize()
problem = study.updProblem()

#Set the effort regularisation goal
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
effort.setWeight(0.001)

#Add contact tracking goal
#Uses three separate vectors to apply different contact tracking weights

#Set right and left contact sphere groups
rightFootContacts = [f'/forceset/contact_{pointName}' for pointName in footGroundNames if pointName.startswith('r_')]
leftFootContacts = [f'/forceset/contact_{pointName}' for pointName in footGroundNames if pointName.startswith('l_')]

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
                           'weight': [2, 1, 2],
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

#Define the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())

# #Set guess from earlier simulation
# solver.setGuessFile('..\\data\\sprint_markerTrackingSolution.sto')
# solver.resetProblem(problem)

#Alter the initial guess to bump up pelvis_ty to ensure contact spheres don't
#start by penetrating the ground
initialGuess = solver.getGuess()
initialGuess.setState('/jointset/ground_pelvis/pelvis_ty/value',
                      initialGuess.getState('/jointset/ground_pelvis/pelvis_ty/value').to_numpy() + 0.05)
solver.resetProblem(problem)
solver.setGuess(initialGuess)

#Set solver options
solver.set_optim_max_iterations(1000)
solver.set_num_mesh_intervals(50)
solver.set_minimize_implicit_multibody_accelerations(True)
solver.set_optim_constraint_tolerance(1e-2) #default
solver.set_optim_convergence_tolerance(1e-2) #default
solver.set_minimize_implicit_auxiliary_derivatives(True)
solver.set_implicit_auxiliary_derivatives_weight(0.00001)

#Solve!
solution = study.solve()
# study.visualize(solution)

#Remove tracked markers file
os.remove('torqueDrivenIkGrfTracking_tracked_states.sto')

#Save solution to file (TODO: GRF data too...)
solution.write('..\\data\\sprint_markerTrackingSolution.sto')

##### Seem to be ground contact issues
    ##### perhaps bump sphere's up so that bottom is in line with floor markers?
        #### Or maybe half radius --- entire radius seems to high?
        #### Did it solve the problem though?
            #### Might not have been this and might just have been impatient with solver?
    ##### Removing constraints didn't help
    ##### Residual actuators? Or is there too many contact spheres?

# %% ----- End of processData.py ----- %% #