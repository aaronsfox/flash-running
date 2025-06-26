# Opensim helper functions for processing data

# %% Import packages

import opensim as osim
import numpy as np

# %% Function for creating models for use in various simulations

def createSimModel(inputModelFile = None,
                   outputModelFile = None,
                   optForces = None,
                   unilateralMuscles = False,
                   jointsToWeld = [],
                   externalLoadsFile = None,
                   addMetabolicsModel = False):

    """
    
    Convenience function for converting IK results to a states storage.
    
    Input:    inputModelFile - file name of model to edit
              outputModelFile - file name for saving model
              optForces - a dictionary that specifies the torque actuators to append to the model with labels
                          each key in the dictionary is a coordinate with it's value a list of the optimal force and label
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

    #Convert muscles to DeGrooteFregly type
    modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
    # modelProcessor.append(osim.ModOpTendonComplianceDynamicsModeDGF('implicit'))
    
    #Set joints to weld for simulation
    if len(jointsToWeld) > 0:
        
        #Create vector string object
        weldVectorStr = osim.StdVectorString()
        [weldVectorStr.append(joint) for joint in jointsToWeld]
        
        #Append model processor
        modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))
        
    #Append external loads if provided
    #Otherwise add contact spheres
    if externalLoadsFile is not None:
        
        #Append external loads to processor
        modelProcessor.append(osim.ModOpAddExternalLoads(externalLoadsFile))

    #Process model for further edits
    osimModel = modelProcessor.process()

    #Remove left limb muscles (if desired) and upper body forces (torque actuated)
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

    #Add torque actuators to model based on optimal forces dictionary
    if optForces is not None:
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
        
    #Add contact spheres if no external loads are provided
    if externalLoadsFile is None:
                
        #Set the reference contact sphere locations based on Faliise et al. OpenSim AD process
        
        #### TODO: do we need more participant-specific locations?
        
        #Set the contact half space details
        referenceHalfSpace = {'name': 'floor', 'location': np.array([0, 0, 0]),
                              'orientation': np.array([0, 0, -np.pi/2]),
                              'frame': 'ground'}
        
        #Set contact spheres dictionary
        refContactSpheres = {
            's1_r': {'radius': 0.032, 'location': np.array([0.0019011578840796601,   -0.01,  -0.00382630379623308]), 'orientation': np.array([0, 0, 0]), 'socket_frame': 'calcn_r'},
            's2_r': {'radius': 0.032, 'location': np.array([0.14838639994206301,     -0.01,  -0.028713422052654002]), 'orientation': np.array([0, 0, 0]), 'socket_frame': 'calcn_r'},
            's3_r': {'radius': 0.032, 'location': np.array([0.13300117060705099,     -0.01,  0.051636247344956601]), 'orientation': np.array([0, 0, 0]), 'socket_frame': 'calcn_r'},
            # 's4_r': {'radius': 0.032, 'location': np.array([0.066234666199163503,    -0.01,  0.026364160674169801]), 'orientation': np.array([0, 0, 0]), 'socket_frame': 'calcn_r'},
            's5_r': {'radius': 0.032, 'location': np.array([0.059999999999999998,    -0.01,  -0.018760308461917698]), 'orientation': np.array([0, 0, 0]), 'socket_frame': 'toes_r' },
            # 's6_r': {'radius': 0.032, 'location': np.array([0.044999999999999998,    -0.01,  0.061856956754965199]), 'orientation': np.array([0, 0, 0]), 'socket_frame': 'toes_r' },
            's1_l': {'radius': 0.032, 'location': np.array([0.0019011578840796601,   -0.01,  0.00382630379623308]), 'orientation': np.array([0, 0, 0]), 'socket_frame': 'calcn_l'},
            's2_l': {'radius': 0.032, 'location': np.array([0.14838639994206301,     -0.01,  0.028713422052654002]), 'orientation': np.array([0, 0, 0]), 'socket_frame': 'calcn_l'},
            's3_l': {'radius': 0.032, 'location': np.array([0.13300117060705099,     -0.01,  -0.051636247344956601]), 'orientation': np.array([0, 0, 0]), 'socket_frame': 'calcn_l'},
            # 's4_l': {'radius': 0.032, 'location': np.array([0.066234666199163503,    -0.01,  -0.026364160674169801]), 'orientation': np.array([0, 0, 0]), 'socket_frame': 'calcn_l'},
            's5_l': {'radius': 0.032, 'location': np.array([0.059999999999999998,    -0.01,  0.018760308461917698]), 'orientation': np.array([0, 0, 0]), 'socket_frame': 'toes_l' },
            # 's6_l': {'radius': 0.032, 'location': np.array([0.044999999999999998,    -0.01,  -0.061856956754965199]), 'orientation': np.array([0, 0, 0]), 'socket_frame': 'toes_l' }
            }  
        
        #Set the reference scale factors
        refScaleFactors = {'calcn_r': np.array([0.91392399999999996, 0.91392399999999996, 0.91392399999999996]),
                           'toes_r':  np.array([0.91392399999999996, 0.91392399999999996, 0.91392399999999996]),
                           'calcn_l': np.array([0.91392399999999996, 0.91392399999999996, 0.91392399999999996]),
                           'toes_l':  np.array([0.91392399999999996, 0.91392399999999996, 0.91392399999999996])}
        
        #Set contact sphere parameters
        stiffness = 3067776
        dissipation = 2.0
        staticFriction = 0.8
        dynamicFriction = 0.8
        viscousFriction = 0.5
        transitionVelocity = 0.2
        
        #Add the contact half space
        
        #Create the half space
        contactHalfSpace = osim.ContactHalfSpace(
            osim.Vec3(referenceHalfSpace['location']),
            osim.Vec3(referenceHalfSpace['orientation']),
            osimModel.get_ground(), referenceHalfSpace['name'])
        
        #Add the half space to model geometry
        osimModel.addContactGeometry(contactHalfSpace)
        
        #Add the contact spheres and forces
        
        #Loop through contact spheres
        for refSphere in refContactSpheres:
            
            #Create the contact spheres
            
            #Get the body to attach to
            body = osimModel.updBodySet().get(refContactSpheres[refSphere]['socket_frame'])
            
            #Calculate the scale factors
            bodyScaleFactors = body.get_attached_geometry(0).get_scale_factors().to_numpy() 
            setScaleFactors = refScaleFactors[refContactSpheres[refSphere]['socket_frame']]
            finalScaleFactors = setScaleFactors / bodyScaleFactors
            finalScaleFactors[1] = 1
            
            #Set the scaled location
            scaledLocation = refContactSpheres[refSphere]['location'] / finalScaleFactors
            
            #Create the contact sphere object
            contactSphere = osim.ContactSphere(
                refContactSpheres[refSphere]['radius'],
                osim.Vec3(scaledLocation), body, refSphere)
            
            #Connect the sphere to the body
            contactSphere.connectSocket_frame(body)
            
            #Add to the models contact geometry
            osimModel.addContactGeometry(contactSphere)
            
            #Create the smooth sphere half space force
            
            #Create the force object
            sphereForce = osim.SmoothSphereHalfSpaceForce(
                'contactForce_' + refSphere, 
                contactSphere, contactHalfSpace)
            
            #Set the parameters in the sphere force
            sphereForce.set_stiffness(stiffness)
            sphereForce.set_dissipation(dissipation)
            sphereForce.set_static_friction(staticFriction)
            sphereForce.set_dynamic_friction(dynamicFriction)
            sphereForce.set_viscous_friction(viscousFriction)
            sphereForce.set_transition_velocity(transitionVelocity)        
            sphereForce.connectSocket_half_space(contactHalfSpace)
            sphereForce.connectSocket_sphere(contactSphere)
            
            #Add the force to the model
            osimModel.addForce(sphereForce)

    #Check whether to add metabolics model
    if addMetabolicsModel:
        
        #Create the metabolics model
        metabolics = osim.Bhargava2004SmoothedMuscleMetabolics()
        metabolics.setName('metabolicModel')
        metabolics.set_use_smoothing(True)
        
        #Add the specific muscles to the metabolics object
        for muscleInd in range(osimModel.getMuscles().getSize()):
            #Get muscle name    
            muscleName = osimModel.getMuscles().get(muscleInd).getName()
            #Add to metabolics object
            metabolics.addMuscle(muscleName, osimModel.getMuscles().get(muscleName))
        
        #Add metabolics to model
        osimModel.addComponent(metabolics)

    #Finalise model connections
    osimModel.finalizeConnections()

    #Save model to file for later use if needed
    osimModel.printToXML(outputModelFile)
    
    return osimModel

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
    
    
