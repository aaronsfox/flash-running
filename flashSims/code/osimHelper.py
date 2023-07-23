# Opensim helper functions for processing data

# %% Import packages

import opensim as osim
import numpy as np

# %% Set-up

#Add OpenSim geometry path (weird issues with this on new laptop)
osim.ModelVisualizer.addDirToGeometrySearchPaths('C:\\OpenSim 4.3\\Geometry')

# %% Function for creating models for use in various simulations

def createFlashModel(inputModelFile = None, outputModelFile = None,
                     unilateralMuscles = False, jointsToWeld = [],
                     addMetabolicsModel = False):

    """
    
    Convenience function for converting IK results to a states storage.
    
    Input:    inputModelFile - file name of model to edit
              outputModelFile - file name for saving model
              unilateralMuscles - whether to remove left side muscles (default = False)
              jointsToWeld - list of joints to weld in model
              addMetabolicsModel - flag whether to add Bhargava metabolics to model (default = False)
    
    """
    
    #Check inputs
    if inputModelFile is None or outputModelFile is None:
        raise ValueError('Both an input and output model file is required.')
    
    #Edit the model for use in the tracking sim tool
    modelProcessor = osim.ModelProcessor(inputModelFile)

    #Append the necessary operators to the model
    #Handles converting muscles to DeGrooteFregly and setting parameters
    modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
    modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
    modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
    modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(5.0))
    modelProcessor.append(osim.ModOpScaleMaxIsometricForce(2))
    # modelProcessor.append(osim.ModOpIgnoreActivationDynamics()) ### added for Flash model

    #Process model for further edits
    osimModel = modelProcessor.process()

    #Remove left limb muscles (if desired) and upper body forces
    #Set a list of forces to remove
    removeForceInd = []
    #Loop through forces and identify muscles to remove
    for forceInd in range(osimModel.updForceSet().getSize()):
        #Check for muscle
        if osimModel.updForceSet().get(forceInd).getConcreteClassName().endswith('Muscle'):
            #Check for left hand side or upper body
            if unilateralMuscles and osimModel.updForceSet().get(forceInd).getName().endswith('_l'):
                #Append index to list
                removeForceInd.append(forceInd)
            elif osimModel.updForceSet().get(forceInd).getName().split('_')[0] in ['extobl', 'intobl', 'ercspn']:
                #Append index to list
                removeForceInd.append(forceInd)

    #Remove the designated forces keeping in mind that the index reduces each time
    #another force is removed
    for removeInd in removeForceInd:
        osimModel.updForceSet().remove(removeInd - removeForceInd.index(removeInd))

    #Set the coordinates that will need torque actuation
    #Set an optimal force and label for each actuator
    #Upper body and left side are idealised actuators (when muscles are removed)
    #Residuals and reserve actuators are lowly weighted for optimal force
    if unilateralMuscles:
        optForces = {
            #upper body
            'lumbar_extension': [1000, 'actuator'], 'lumbar_bending': [1000, 'actuator'], 'lumbar_rotation': [1000, 'actuator'],
            'arm_flex_r': [300, 'actuator'], 'arm_add_r': [300, 'actuator'], 'arm_rot_r': [300, 'actuator'],
            'elbow_flex_r': [100, 'actuator'], 'pro_sup_r': [100, 'actuator'],
            'arm_flex_l': [300, 'actuator'], 'arm_add_l': [300, 'actuator'], 'arm_rot_l': [300, 'actuator'],
            'elbow_flex_l': [100, 'actuator'], 'pro_sup_l': [100, 'actuator'],
            #left limb
            'hip_flexion_l': [100, 'actuator'], 'hip_adduction_l': [50, 'actuator'], 'hip_rotation_l': [50, 'actuator'],
            'knee_angle_l': [100, 'actuator'], 'ankle_angle_l': [100, 'actuator'],
            # 'mtp_angle_l': [50, 'actuator'],
            #right limb
            'hip_flexion_r': [2, 'reserve'], 'hip_adduction_r': [2, 'reserve'], 'hip_rotation_r': [2, 'reserve'],
            'knee_angle_r': [2, 'reserve'], 'ankle_angle_r': [2, 'reserve'], 
            # 'mtp_angle_r': [2, 'reserve'],
            # #pelvis
            # 'pelvis_tx': [1, 'residual'], 'pelvis_ty': [1, 'residual'], 'pelvis_tz': [1, 'residual'],
            # 'pelvis_tilt': [1, 'residual'], 'pelvis_list': [1, 'residual'], 'pelvis_rotation': [1, 'residual']
            }
    else:
        optForces = {
            #upper body
            'lumbar_extension': [300, 'actuator'], 'lumbar_bending': [300, 'actuator'], 'lumbar_rotation': [300, 'actuator'],
            'arm_flex_r': [300, 'actuator'], 'arm_add_r': [300, 'actuator'], 'arm_rot_r': [300, 'actuator'],
            'elbow_flex_r': [100, 'actuator'], 'pro_sup_r': [100, 'actuator'],
            'arm_flex_l': [300, 'actuator'], 'arm_add_l': [300, 'actuator'], 'arm_rot_l': [300, 'actuator'],
            'elbow_flex_l': [100, 'actuator'], 'pro_sup_l': [100, 'actuator'],
            #left limb
            'hip_flexion_l': [2, 'reserve'], 'hip_adduction_l': [2, 'reserve'], 'hip_rotation_l': [2, 'reserve'],
            'knee_angle_l': [2, 'reserve'], 'ankle_angle_l': [2, 'reserve'],
            # 'mtp_angle_l': [2, 'reserve'],
            #right limb
            'hip_flexion_r': [2, 'reserve'], 'hip_adduction_r': [2, 'reserve'], 'hip_rotation_r': [2, 'reserve'],
            'knee_angle_r': [2, 'reserve'], 'ankle_angle_r': [2, 'reserve'], 
            # 'mtp_angle_r': [2, 'reserve'],
            # #pelvis
            # 'pelvis_tx': [1, 'residual'], 'pelvis_ty': [1, 'residual'], 'pelvis_tz': [1, 'residual'],
            # 'pelvis_tilt': [1, 'residual'], 'pelvis_list': [1, 'residual'], 'pelvis_rotation': [1, 'residual']
            }

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

    #Increase the maximum contraction velocity of muscles to extreme
    #Reduce activation and deactivation time constants
    maxContractionVelocity = 1e5
    # activationConstant = 1e-4
    # deactivationConstant = 4e-4
    
    #Loop through muscles    
    for muscleInd in range(osimModel.getMuscles().getSize()):
        #Contraction velocity
        osimModel.getMuscles().get(muscleInd).set_max_contraction_velocity(maxContractionVelocity)
        # #Activation dynamics
        # osim.DeGrooteFregly2016Muscle().safeDownCast(osimModel.getMuscles().get(muscleInd)).set_activation_time_constant(activationConstant)
        # osim.DeGrooteFregly2016Muscle().safeDownCast(osimModel.getMuscles().get(muscleInd)).set_deactivation_time_constant(deactivationConstant)
        #Increase max control of muscle to infitinity
        osimModel.getMuscles().get(muscleInd).set_max_control(np.inf)
            
    #Reduce the moment of inertia across bodies
    for bodyInd in range(osimModel.updBodySet().getSize()):
        #Get current bodies inertia
        currInertia = osimModel.updBodySet().get(bodyInd).get_inertia()
        #Reduce moments of inertia by a factor of 100 and set in model
        osimModel.updBodySet().get(bodyInd).set_inertia(osim.Vec6(currInertia.get(0) / 1000,
                                                                  currInertia.get(1) / 1000,
                                                                  currInertia.get(2) / 1000,
                                                                  currInertia.get(3), currInertia.get(4), currInertia.get(5)))

    #Finalise model connections
    osimModel.finalizeConnections()
    
    #Create the tracking model processor for use in the tool
    modelProcessor = osim.ModelProcessor(osimModel)

    #Set joints to weld for simulation
    if len(jointsToWeld) > 0:
        #Create vector string object
        weldVectorStr = osim.StdVectorString()
        [weldVectorStr.append(joint) for joint in jointsToWeld]
        #Append model processor
        modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))

    #Process to model form before adding contact spheres
    finalModel = modelProcessor.process()
    
    #Add contact spheres at foot-ground contact model locations

    #Get the markers that contain an fp reference for contact sphere locations
    footGroundNames = []
    for markerInd in range(finalModel.updMarkerSet().getSize()):
        #Check for fp indicator in marker
        if '_fp_' in finalModel.updMarkerSet().get(markerInd).getName():
            footGroundNames.append(finalModel.updMarkerSet().get(markerInd).getName())
            
    #Get foot ground name locations in dictionary
    footGroundLocs = {name: np.array((finalModel.updMarkerSet().get(name).get_location().get(0),
                                      finalModel.updMarkerSet().get(name).get_location().get(1),
                                      finalModel.updMarkerSet().get(name).get_location().get(2))) for name in footGroundNames}

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
        elif 'mt' in name: ### or 'toe' in name:
        # elif 'mt' in name:
            #Set the name and location as what is already done
            footSphereLocs[name.split('_')[-1]+'_'+name[0]] = footGroundLocs[name]
    
    #Add the additional sphere at the mid metatarsal location
    footSphereLocs['mtMid_r'] = np.array(((footSphereLocs['mt1_r'][0]+footSphereLocs['mt2_r'][0])/2,
                                          (footSphereLocs['mt1_r'][1]+footSphereLocs['mt2_r'][1])/2,
                                          (footSphereLocs['mt1_r'][2]+footSphereLocs['mt2_r'][2])/2))
    footSphereLocs['mtMid_l'] = np.array(((footSphereLocs['mt1_l'][0]+footSphereLocs['mt2_l'][0])/2,
                                          (footSphereLocs['mt1_l'][1]+footSphereLocs['mt2_l'][1])/2,
                                          (footSphereLocs['mt1_l'][2]+footSphereLocs['mt2_l'][2])/2))
            
    #Set foot ground contact sphere sizes
    heelSphereRadius = 0.04
    otherSphereRadius = 0.03

    #Add sphere size into dictionary
    footSphereSizes = {name: heelSphereRadius if 'heel' in name else otherSphereRadius for name in footSphereLocs.keys()}

    #Create and connect the half space floor to ground
    floorContact = osim.ContactHalfSpace() #create object
    floorContact.setName('floor') #set name
    floorContact.setOrientation(osim.Vec3(0, 0, -1.5707963267949001)) #set orientation
    floorContact.connectSocket_frame(osim.PhysicalFrame.safeDownCast(finalModel.getGround())) #connect to ground
    finalModel.addContactGeometry(floorContact) #add to model

    #Iteratively create the contact geometry and spheres and attach to models force set
    for contactPoint in footSphereLocs.keys():
        
        #Create the contact geometry and set parameters
        contactGeom = osim.ContactSphere() #create object
        contactGeom.setName(contactPoint) #set name
        
        #Conditional for heel or toe sphere
        if contactPoint.startswith('heel') and contactPoint.endswith('_r'):
            contactGeom.connectSocket_frame(finalModel.updBodySet().get('calcn_r')) #connect to frame
        elif contactPoint.startswith('heel') and contactPoint.endswith('_l'):
            contactGeom.connectSocket_frame(finalModel.updBodySet().get('calcn_l')) #connect to frame
        elif not contactPoint.startswith('heel') and contactPoint.endswith('_r'):
            contactGeom.connectSocket_frame(finalModel.updBodySet().get('toes_r')) #connect to frame
        elif not contactPoint.startswith('heel') and contactPoint.endswith('_l'):
            contactGeom.connectSocket_frame(finalModel.updBodySet().get('toes_l')) #connect to frame
            
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
            if contactPoint.endswith('_r'):
                jointTranslation = finalModel.updJointSet().get('mtp_r').get_frames(0).get_translation()
            else:
                jointTranslation = finalModel.updJointSet().get('mtp_l').get_frames(0).get_translation()
            contactGeom.setLocation(osim.Vec3(footSphereLocs[contactPoint][0] - jointTranslation.get(0),
                                              # footSphereLocs[contactPoint][1] - jointTranslation.get(1),
                                              footSphereLocs[contactPoint][1] - jointTranslation.get(1) + (currSphereSize/2), #adjust to bottom of foot height
                                              footSphereLocs[contactPoint][2] - jointTranslation.get(2))) #set location
        contactGeom.setRadius(currSphereSize) #set radius
        finalModel.addContactGeometry(contactGeom) #add to model
        
        #Create the sphere and set properties
        contactSphere = osim.SmoothSphereHalfSpaceForce() #create force
        contactSphere.setName('contact_'+contactPoint) #set name
        contactSphere.connectSocket_half_space(finalModel.updContactGeometrySet().get('floor')) #connect to floor
        contactSphere.connectSocket_sphere(finalModel.updContactGeometrySet().get(contactPoint)) #connect to sphere
        contactSphere.set_stiffness(3067776) #set stiffness
        contactSphere.set_dissipation(2) #set dissipation
        contactSphere.set_static_friction(0.8) #set static friction
        contactSphere.set_dynamic_friction(0.8) #set dynamic friction
        contactSphere.set_viscous_friction(0.5) #set viscous friction
        contactSphere.set_transition_velocity(0.2) #set transition velocity
        contactSphere.set_hertz_smoothing(300) #set hertz smoothing
        contactSphere.set_hunt_crossley_smoothing(50) #set hunt crossley smoothing
        finalModel.addForce(contactSphere) #add to model
            
    #Check whether to add metabolics model
    if addMetabolicsModel:
        
        #Create the metabolics model
        metabolics = osim.Bhargava2004SmoothedMuscleMetabolics()
        metabolics.setName('metabolicModel')
        metabolics.set_use_smoothing(True)
        
        #Add the specific muscles to the metabolics object
        for muscleInd in range(finalModel.getMuscles().getSize()):
            #Get muscle name    
            muscleName = finalModel.getMuscles().get(muscleInd).getName()
            #Add to metabolics object
            metabolics.addMuscle(muscleName, finalModel.getMuscles().get(muscleName))
        
        #Add metabolics to model
        finalModel.addComponent(metabolics)

    #Finalise model connections
    finalModel.finalizeConnections()

    #Save model to file for later use if needed
    finalModel.printToXML(outputModelFile)
    
    return finalModel

# %% getGaitTimings

def getGaitTimings(grfFile = None, extLoads = None,
                   startForceName = 'RightGRF1',
                   stopForceName = 'LeftGRF1',
                   forceThreshold = 20):
    
    # Convenience function for getting gait cycle timings based on force data.
    # Note that this function is only applicable to external loads where foot strikes
    # are labelled as individual forces within an XML file (i.e. individual strikes
    # on force plates during overground running is the scenario this is used in).
    # Timings for a half gait cycle can be achieved by using a right then left, 
    # or left then right force as the start and stop. Timings for a full gait cycle
    #
    # Input:    grfFile - .mot file containing GRF time history
    #           extLoads - .xml file for external loads linked to GRF data
    #           startForceName - the force name and hence force identifier to start the timings with
    #           stopForceName - the force name and hence force identifier to start the timings with
    #           forceThreshold - force in Newtons for detecting contact
    
    #Check for input
    if grfFile is None:
        raise ValueError('A GRF file in .mot format is required')
    if extLoads is None:
        raise ValueError('An external loads file in .xml format is required')
        
    #Load the GRF data
    grfTable = osim.TimeSeriesTable(grfFile)
    
    #Read in the external loads file
    externalLoads = osim.ExternalLoads(extLoads, True)
    
    #Get the force identifiers for the desired force names
    startForceIdentifier = externalLoads.get(startForceName).getForceIdentifier()
    stopForceIdentifier = externalLoads.get(stopForceName).getForceIdentifier()
    
    #Get the vertical force from the GRF data as a numpy array
    startVGRF = grfTable.getDependentColumn(f'{startForceIdentifier}y').to_numpy()
    stopVGRF = grfTable.getDependentColumn(f'{stopForceIdentifier}y').to_numpy()
    
    #Find the first instance where vGRF > 50N in forces
    startOnInd = np.argmax(startVGRF > forceThreshold)
    stopOnInd = np.argmax(stopVGRF > forceThreshold)
    
    #Check whether stop is after start and raise an error if not
    if stopOnInd < startOnInd:
        raise ValueError('End time before start time. Check names and order of forces used.')
    
    #Get the times for the gait timings
    startTime = grfTable.getIndependentColumn()[startOnInd]
    endTime = grfTable.getIndependentColumn()[stopOnInd]
    
    #Print outputs
    print(f'Start Time: {startTime}')
    print(f'End Time: {endTime}')
    
    return startTime,endTime

# %%
    
