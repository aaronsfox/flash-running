# Opensim helper functions for processing data

# %% Import packages

import opensim as osim
import numpy as np

# %% Set-up

#Add OpenSim geometry path (weird issues with this on new laptop)
osim.ModelVisualizer.addDirToGeometrySearchPaths('C:\\OpenSim 4.3\\Geometry')

# %% Function for creating models for use in various simulations

def createSimModel(inputModelFile = None, outputModelFile = None,
                   unilateralMuscles = False, jointsToWeld = [],
                   externalLoadsFile = None,
                   addMetabolicsModel = False):

    """
    
    Convenience function for converting IK results to a states storage.
    
    Input:    inputModelFile - file name of model to edit
              outputModelFile - file name for saving model
              unilateralMuscles - whether to remove left side muscles (default = False)
              jointsToWeld - list of joints to weld in model
              externalLoadsFile - provide external loads to append or else add contact spheres
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
    modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(2.0))
    modelProcessor.append(osim.ModOpScaleMaxIsometricForce(2))
    # modelProcessor.append(osim.ModOpTendonComplianceDynamicsModeDGF('implicit'))

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
            'knee_angle_l': [100, 'actuator'], 'ankle_angle_l': [100, 'actuator'], 'mtp_angle_l': [50, 'actuator'],
            #right limb
            'hip_flexion_r': [2, 'reserve'], 'hip_adduction_r': [2, 'reserve'], 'hip_rotation_r': [2, 'reserve'],
            'knee_angle_r': [2, 'reserve'], 'ankle_angle_r': [2, 'reserve'],  'mtp_angle_r': [2, 'reserve'],
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
            'knee_angle_l': [2, 'reserve'], 'ankle_angle_l': [2, 'reserve'], 'mtp_angle_l': [2, 'reserve'],
            #right limb
            'hip_flexion_r': [2, 'reserve'], 'hip_adduction_r': [2, 'reserve'], 'hip_rotation_r': [2, 'reserve'],
            'knee_angle_r': [2, 'reserve'], 'ankle_angle_r': [2, 'reserve'],  'mtp_angle_r': [2, 'reserve'],
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
    
    #Create the tracking model processor for use in the tool
    modelProcessor = osim.ModelProcessor(osimModel)

    #Set joints to weld for simulation
    if len(jointsToWeld) > 0:
        #Create vector string object
        weldVectorStr = osim.StdVectorString()
        [weldVectorStr.append(joint) for joint in jointsToWeld]
        #Append model processor
        modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))
        
    #Append external loadds if provided
    #Otherwise add contact spheres
    if externalLoadsFile is not None:
        
        #Append external loads to processor
        modelProcessor.append(osim.ModOpAddExternalLoads(externalLoadsFile))
        
        #Process the model as no more operators are left
        finalModel = modelProcessor.process()
        
    else:
        
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
                if contactPoint.endswith('r_'):
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

# %% Function for converting IK results to states storage file

def kinematicsToStates(kinematicsFileName = None, osimModelFileName = None,
                       outputFileName = 'coordinates.sto',
                       inDegrees = True, outDegrees = False):
    
    """
    
    Convenience function for converting IK results to a states storage.
    
    Input:    kinematicsFileName - file containing kinematic data. Header should only be coordinates name, rather than path to state
              osimModelFileName - opensim model filename that corresponds to kinematic data
              outputFileName - optional filename to output to (defaults to coordinates.sto)
              inDegrees - set to true if kinematics file is in degrees (defaults to True)
              outDegrees - set to true if desired output is in degrees (defaults to False)
    
    """
    
    #Check inputs
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
    
# %% getMassOfModel

def getMassOfModel(osimModelFileName = None):
    
    # Convenience function for getting total mass of model.
    #
    # Input:    osimModelFileName - opensim model filename that corresponds to kinematic data
    
    if osimModelFileName is None:
        raise ValueError('OpenSim model filename is required')
        
    #Set starting mass
    totalMass = 0
    
    #Load in model
    osimModel = osim.Model(osimModelFileName)
    
    #Get bodies
    allBodies = osimModel.getBodySet()
    
    #Loop through bodies and get mass
    for ii in range(allBodies.getSize()):
        totalMass = totalMass + allBodies.get(ii).getMass()
        
    #Return total mass
    return totalMass

# %% Function to add set of torque actuators to model

def addTorqueActuators(osimModel = None,
                       optForces = None,
					   minControl = -1,
					   maxControl = +1):
    
    """
    
    Convenience function for adding series of torque actuators to model
    
    Input:    osimModel - OpenSim model object for use
              optForces - dict of coordinates and their associated optimal forces to add
			  minControl - minimum control signal value for actuators (default = -1)
			  maxControl - maximum control signal value for actuators (default = +1)
              
    Output:   osimModel - updated torque driven model
                  
    """
    
    #Check inputs
    if osimModel is None or optForces is None:
            raise ValueError('All inputs for this function are required!')
    
    #Remove the original lumbar actuators to not apply force
    #Get the force indices to remove from the model
    forceRemove = []
    for forceInd in range(osimModel.updForceSet().getSize()):
        if 'lumbar' in osimModel.updForceSet().get(forceInd).getName():
            forceRemove.append(forceInd)
    #Each time the force is removed the indices reduce, so need to account for this
    for ind in range(len(forceRemove)):
        osimModel.updForceSet().remove(forceRemove[ind]-ind)
    
    #Intialise model system
    osimModel.initSystem()
    
    #Get coordinate list
    coordinatesList = list(optForces.keys())
    
    #Get coordinate set
    coordSet = osimModel.getCoordinateSet()
    
    #Loop through coordinates and add actuators
    for coordinate in coordinatesList:
        #Create actuator
        actu = osim.CoordinateActuator()
        #Set name
        actu.setName(f'{coordinate}_actuator')
        #Set coordinate
        actu.setCoordinate(coordSet.get(coordinate))
        #Set optimal force
        actu.setOptimalForce(optForces[coordinate])
        #Set min and max control
        actu.setMinControl(minControl)
        actu.setMaxControl(maxControl)
        #Append to model force set
        osimModel.updForceSet().cloneAndAppend(actu)
        # #Append to model components
        # osimModel.addComponent(actu)
    
    #Finalise model connections
    osimModel.finalizeConnections()
    
    #Return model
    return osimModel

# %% addCoordinateActuator

def addCoordinateActuator(osimModel = None,
                          coordName = None,
                          optForce = 1000,
                          controlVals = None):
    
    # Convenience function for adding torque coordinate actuators to model
    #
    # Input:    osimModel - model object to add actuator to
    #           coordName - name of coordinate to actuate
    #           optForce - optimal force value for coordinate actuator
    #           controlVals - min and maximum control values for actuator (if None Inf is used)
    
    #Input checks
    if osimModel is None:
        raise ValueError('A model object is required.')
        
    if coordName is None:
        raise ValueError('A coordinate name is required.')
    
    #Get coordinate set from model
    coordSet = osimModel.updCoordinateSet()
    
    #Create actuator
    actu = osim.CoordinateActuator()
    
    #Set actuator name
    actu.setName('tau_'+coordName)
    
    #Set coordinate for actuator
    actu.setCoordinate(coordSet.get(coordName))
    
    #Set optimal force for actuator
    actu.setOptimalForce(optForce)
    
    #Set max and min controls
    if controlVals is None:
        actu.setMinControl(np.inf*-1)
        actu.setMaxControl(np.inf)
    else:
        actu.setMinControl(controlVals[0])
        actu.setMaxControl(controlVals[1])
        
    #Add actuator to model
    osimModel.updForceSet().cloneAndAppend(actu)
    
    
