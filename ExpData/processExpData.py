# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 11:31:10 2020

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This code processes the experimental max sprinting trial contained in this
    folder for use as an initial guess in future predictive simulations.
    
"""

# %% Import packages

import os
import btk
import math
import opensim as osim
import osimHelper
import numpy as np
import pandas as pd
from scipy.signal import butter, lfilter
import matplotlib.pyplot as plt
import xml.etree.ElementTree as et

# %% Participant details

#Get participant mass from static trial

#Load in the c3d data
staticC3D = btk.btkAcquisitionFileReader()
staticC3D.SetFilename('static.c3d')
staticC3D.Update()
staticAcq = staticC3D.GetOutput()

#Get force data
#Participant force data is on analog channel 4, we can get the vertical force
staticForce = staticAcq.GetAnalog('Force.Fz4').GetValues()

#Calculate mass
massKg = (np.mean(staticForce) * -1) / 9.80665

# %% Convert c3d data to OpenSim formats

#Convert static c3d to .trc

#Set up data adapters
staticC3DFile = osim.C3DFileAdapter()
staticTRCFile = osim.TRCFileAdapter()

#Get markers
staticData = staticC3DFile.read('static.c3d')
staticMarkers = staticC3DFile.getMarkersTable(staticData)

#Rotate marker data
#First rotation
rotMat1 = osim.Rotation(np.deg2rad(-90),
                        osim.Vec3(0, 0, 1))
for ii in range(staticMarkers.getNumRows()):
    vec = staticMarkers.getRowAtIndex(ii)
    vec_rotated = rotMat1.multiply(vec)
    staticMarkers.setRowAtIndex(ii, vec_rotated)
#Second rotation
rotMat2 = osim.Rotation(np.deg2rad(-90),
                        osim.Vec3(1, 0, 0))
for ii in range(staticMarkers.getNumRows()):
    vec = staticMarkers.getRowAtIndex(ii)
    vec_rotated = rotMat2.multiply(vec)
    staticMarkers.setRowAtIndex(ii, vec_rotated)
    
#Write static data to file
staticTRCFile.write(staticMarkers, 'static.trc')

#Convert dynamic c3d to .trc

#Set up data adapters
dynamicC3DFile = osim.C3DFileAdapter()
dynamicTRCFile = osim.TRCFileAdapter()

#Get markers
dynamicData = dynamicC3DFile.read('maxSprint.c3d')
dynamicMarkers = dynamicC3DFile.getMarkersTable(dynamicData)

#Rotate marker data
#First rotation
for ii in range(dynamicMarkers.getNumRows()):
    vec = dynamicMarkers.getRowAtIndex(ii)
    vec_rotated = rotMat1.multiply(vec)
    dynamicMarkers.setRowAtIndex(ii, vec_rotated)
#Second rotation
for ii in range(dynamicMarkers.getNumRows()):
    vec = dynamicMarkers.getRowAtIndex(ii)
    vec_rotated = rotMat2.multiply(vec)
    dynamicMarkers.setRowAtIndex(ii, vec_rotated)
    
#Write static data to file
dynamicTRCFile.write(dynamicMarkers, 'maxSprint.trc')

#Convert c3d force data to .mot and .xml files

#Get forces data
dynamicForces = dynamicC3DFile.getForcesTable(dynamicData)

#Rotate forces
#First rotation
rotMat1 = osim.Rotation(np.deg2rad(-90),
                        osim.Vec3(0, 0, 1))
for ii in range(dynamicForces.getNumRows()):
    vec = dynamicForces.getRowAtIndex(ii)
    vec_rotated = rotMat1.multiply(vec)
    dynamicForces.setRowAtIndex(ii, vec_rotated)
#Second rotation
rotMat2 = osim.Rotation(np.deg2rad(-90),
                        osim.Vec3(1, 0, 0))
for ii in range(dynamicForces.getNumRows()):
    vec = dynamicForces.getRowAtIndex(ii)
    vec_rotated = rotMat2.multiply(vec)
    dynamicForces.setRowAtIndex(ii, vec_rotated)
    
#Convert mm units to m
pList = ['p1','p2','p3','p4','p5','p6','p7','p8']
mList = ['m1','m2','m3','m4','m5','m6','m7','m8']
for ff in range(0,len(pList)):
    osimHelper.mm_to_m(dynamicForces, pList[ff])
    osimHelper.mm_to_m(dynamicForces, mList[ff])

#Define filter for forces data (50N Low-pass 4th Order Butterworth)
fs = 2000
nyq = 0.5 * fs
normCutoff = 50 / nyq
b, a = butter(4, normCutoff, btype = 'low', analog = False)

#Extract the forces data
forcesTime = dynamicForces.getIndependentColumn()
forcesLabels = dynamicForces.getColumnLabels()
forcesMat = np.zeros((len(forcesTime), 3 * len(forcesLabels)))
for ll in range(0,len(forcesLabels)):
    
    #Set current label
    label = forcesLabels[ll]
    
    #Get Vec3 data
    f = dynamicForces.getDependentColumn(label)
    
    #Convert Vec3 into numpy array
    v_np = np.zeros((len(forcesTime), 3))
    for ii in range(0,len(forcesTime)):
        v_np[ii, 0] = f[ii][0]
        v_np[ii, 1] = f[ii][1]
        v_np[ii, 2] = f[ii][2]
        
    #Interpolate to replace NaNs
    v_int = np.zeros((len(forcesTime), 3))
    for jj in range(3):
        v_pd = pd.Series(v_np[:, jj])
        v_pd = v_pd.interpolate(limit_direction = "both", kind = "cubic")
        v_int[:, jj] = v_pd.to_numpy()
        
    #Apply lowpass filter to force data (50N cut-off)
    v_filt = np.zeros((len(forcesTime), 3))
    for jj in range(3):
        v_filt[:,jj] = lfilter(b, a, v_int[:,jj])
        
    #Append to array
    forcesMat[:,ll*3:ll*3+3] = v_filt
        
#Convert array to dataframe
allForceLabels = list()
for ff in range(0,len(forcesLabels)):
    allForceLabels.append(forcesLabels[ff]+'_x')
    allForceLabels.append(forcesLabels[ff]+'_y')
    allForceLabels.append(forcesLabels[ff]+'_z')
df_forces = pd.DataFrame(forcesMat, columns = allForceLabels)

#Extract the forces relating to each foot
#This is a manual process given the singular trial and knowledge of which plates
#have which contacts
#Right foot = fp7,fp3 
#Left foot = fp5,fp1
#Right forces
ground_force_r_vx = df_forces['f7_x'].values + df_forces['f3_x'].values
ground_force_r_vy = df_forces['f7_y'].values + df_forces['f3_y'].values
ground_force_r_vz = df_forces['f7_z'].values + df_forces['f3_z'].values
#Left forces
ground_force_l_vx = df_forces['f5_x'].values + df_forces['f1_x'].values
ground_force_l_vy = df_forces['f5_y'].values + df_forces['f1_y'].values
ground_force_l_vz = df_forces['f5_z'].values + df_forces['f1_z'].values
#Right torques
ground_torque_r_x = df_forces['m7_x'].values + df_forces['m3_x'].values
ground_torque_r_y = df_forces['m7_y'].values + df_forces['m3_y'].values
ground_torque_r_z = df_forces['m7_z'].values + df_forces['m3_z'].values
#Left torques
ground_torque_l_x = df_forces['m5_x'].values + df_forces['m1_x'].values
ground_torque_l_y = df_forces['m5_y'].values + df_forces['m1_y'].values
ground_torque_l_z = df_forces['m5_z'].values + df_forces['m1_z'].values
#Right position
ground_force_r_px = df_forces['p7_x'].values + df_forces['p3_x'].values
ground_force_r_py = df_forces['p7_y'].values + df_forces['p3_y'].values
ground_force_r_pz = df_forces['p7_z'].values + df_forces['p3_z'].values
#Left position
ground_force_l_px = df_forces['p5_x'].values + df_forces['p1_x'].values
ground_force_l_py = df_forces['p5_y'].values + df_forces['p1_y'].values
ground_force_l_pz = df_forces['p5_z'].values + df_forces['p1_z'].values

#Export forces to .mot file
#Set labels for file
motLabels = ['time',
             'ground_force_r_vx', 'ground_force_r_vy', 'ground_force_r_vz',
             'ground_force_l_vx', 'ground_force_l_vy', 'ground_force_l_vz',
             'ground_force_r_px', 'ground_force_r_py', 'ground_force_r_pz',
             'ground_force_l_px', 'ground_force_l_py', 'ground_force_l_pz',
             'ground_torque_r_x', 'ground_torque_r_y', 'ground_torque_r_z',
             'ground_torque_l_x', 'ground_torque_l_y', 'ground_torque_l_z']
#Create storage object
forcesStorage = osim.Storage()
#Set storage file labels
colLabels = osim.ArrayStr()
for cc in range(0,len(motLabels)):
    colLabels.append(motLabels[cc])
forcesStorage.setColumnLabels(colLabels)
#Create data array
forceData = np.transpose(np.array([ground_force_r_vx, ground_force_r_vy, ground_force_r_vz,
                                   ground_force_l_vx, ground_force_l_vy, ground_force_l_vz,
                                   ground_force_r_px, ground_force_r_py, ground_force_r_pz,
                                   ground_force_l_px, ground_force_l_py, ground_force_l_pz,
                                   ground_torque_r_x, ground_torque_r_y, ground_torque_r_z,
                                   ground_torque_l_x, ground_torque_l_y, ground_torque_l_z]))
#Append data to storage object
nrow, ncol = forceData.shape
for ii in range(nrow):
    row = osim.ArrayDouble()
    for jj in range(ncol):
        row.append(forceData[ii,jj])
    #Add data to storage
    forcesStorage.append(forcesTime[ii], row)
#Set name for storage object
forcesStorage.setName('refGRF')
#Print to file
forcesStorage.printResult(forcesStorage, 'refGRF', os.getcwd(), 0.001, '.mot')

#Create a 2D version of the forces with no medial/lateral force and point/torque data
#Create storage object
forcesStorage2D = osim.Storage()
#Set storage file labels
forcesStorage2D.setColumnLabels(colLabels)
#Create dataset
forceData2D = np.transpose(np.array([ground_force_r_vx, ground_force_r_vy, np.zeros([nrow]),
                                     ground_force_l_vx, ground_force_l_vy, np.zeros([nrow]),
                                     np.zeros([nrow]), np.zeros([nrow]), np.zeros([nrow]),
                                     np.zeros([nrow]), np.zeros([nrow]), np.zeros([nrow]),
                                     np.zeros([nrow]), np.zeros([nrow]), np.zeros([nrow]),
                                     np.zeros([nrow]), np.zeros([nrow]), np.zeros([nrow])]))
#Append data to storage object
for ii in range(nrow):
    row = osim.ArrayDouble()
    for jj in range(ncol):
        row.append(forceData2D[ii,jj])
    #Add data to storage
    forcesStorage2D.append(forcesTime[ii], row)
#Set name for storage object
forcesStorage2D.setName('refGRF_2D')
#Print to file
forcesStorage2D.printResult(forcesStorage2D, 'refGRF_2D', os.getcwd(), 0.001, '.mot')

#Create external forces .xml file
forceXML = osim.ExternalLoads()
#Create and append the right GRF external force
rightGRF = osim.ExternalForce()
rightGRF.setName('RightGRF')
rightGRF.setAppliedToBodyName('calcn_r')
rightGRF.setForceExpressedInBodyName('ground')
rightGRF.setPointExpressedInBodyName('ground')
rightGRF.setForceIdentifier('ground_force_r_v')
rightGRF.setPointIdentifier('ground_force_r_p')
rightGRF.setTorqueIdentifier('ground_torque_r_')
forceXML.cloneAndAppend(rightGRF)
#Create and append the left GRF external force
leftGRF = osim.ExternalForce()
leftGRF.setName('LeftGRF')
leftGRF.setAppliedToBodyName('calcn_l')
leftGRF.setForceExpressedInBodyName('ground')
leftGRF.setPointExpressedInBodyName('ground')
leftGRF.setForceIdentifier('ground_force_l_v')
leftGRF.setPointIdentifier('ground_force_l_p')
leftGRF.setTorqueIdentifier('ground_torque_l_')
forceXML.cloneAndAppend(leftGRF)
#Set GRF datafile
forceXML.setDataFileName('refGRF.mot')
#Set filtering for kinematics
forceXML.setLowpassCutoffFrequencyForLoadKinematics(12)
#Write to file
forceXML.printToXML('refGRF.xml')

#Adapt for the 2D version 
forceXML.setDataFileName('refGRF_2D.mot')
forceXML.printToXML('refGRF_2D.xml')

# %% Model scaling

#Set up the scale tool
scaleTool = osim.ScaleTool()

#Set participant mass
scaleTool.setSubjectMass(massKg)

#Set generic model file
genModelFileName = 'genericModel_LaiArnold.osim'
scaleTool.getGenericModelMaker().setModelFileName(genModelFileName)

#Set the measurement set
measurementSetObject = osim.OpenSimObject.makeObjectFromFile('scaleMeasurementSet.xml')
measurementSet = osim.MeasurementSet.safeDownCast(measurementSetObject)
scaleTool.getModelScaler().setMeasurementSet(measurementSet)

#Set scale tasks
taskSet = osim.IKTaskSet('scaleTasks.xml')
for k in range(0,taskSet.getSize()-1):
    scaleTool.getMarkerPlacer().getIKTaskSet().adoptAndAppend(taskSet.get(k))

#Set marker file

#Add the virtual markers and create a new .trc file to use in scaling
osimHelper.addVirtualMarkersStatic('static.trc','staticVirtualMarkers.trc')

#Place in scale tool
scaleTool.getMarkerPlacer().setMarkerFileName('staticVirtualMarkers.trc')
scaleTool.getModelScaler().setMarkerFileName('staticVirtualMarkers.trc')

#Set options
scaleTool.getModelScaler().setPreserveMassDist(True)
scaleOrder = osim.ArrayStr(); scaleOrder.set(0,'measurements')
scaleTool.getModelScaler().setScalingOrder(scaleOrder)

#Set time ranges
timeRange = osim.ArrayDouble()
timeRange.set(0,0.5); timeRange.set(1,1.5)
scaleTool.getMarkerPlacer().setTimeRange(timeRange)
scaleTool.getModelScaler().setTimeRange(timeRange)

#Set output files
scaleTool.getModelScaler().setOutputModelFileName('scaledModel.osim')
scaleTool.getModelScaler().setOutputScaleFileName('scaleSet.xml')

#Set marker adjuster parameters
scaleTool.getMarkerPlacer().setOutputMotionFileName('static_motion.mot')
scaleTool.getMarkerPlacer().setOutputModelFileName('scaledModelAdjusted.osim')

#Save and run scale tool
scaleTool.printToXML('scaleSetup.xml')
scaleTool.run()

#Load in scaled model
scaledModel = osim.Model('scaledModelAdjusted.osim')

#The locations and size of the contact sphere parameters go unchanged with
#standard model scaling, so these need to be edited to ensure they are in
#an appropriate place. This can be done based on the scale set parameters
#for the representative bodies.

#Load in the scale set, parsed from the XML tree
xmlTree = et.parse('scaleSet.xml')
xmlRoot = xmlTree.getroot()

#Create a dictionary with segments and scale factors to access for calculations
scaleFactors = {}

#Loop through segments and place  in dictionary
for segment in range(0,len(xmlRoot.findall('./ScaleSet/objects/Scale/segment'))):
    #Get current segment name
    currSegment = xmlRoot.findall('./ScaleSet/objects/Scale/segment')[segment].text
    #Get current scale name and parse to 0 x 3 array
    currScale = xmlRoot.findall('./ScaleSet/objects/Scale/scales')[segment].text
    currScale = str.split(currScale)
    scaleFactors[currSegment] = [float(currScale[0]),float(currScale[1]),float(currScale[2])]
    
#Get the 3D scale factors for the relevant bodies and average to scale the
#sphere radii. Note that these scale factors for each foot will be the same
#for heel and toes as it looks like same scale factors are applied.
heelSphereScale_r = sum(scaleFactors['calcn_r']) / 3
heelSphereScale_l = sum(scaleFactors['calcn_l']) / 3
toesSphereScale_r = sum(scaleFactors['toes_r']) / 3
toesSphereScale_l = sum(scaleFactors['toes_l']) / 3

#Scale the radii to each of their respective factors
#While accessing the spheres, also adjust their position based on the scale
#factor for the respective axes
#Create a list of the sheres to loop through and edit
sphereList = ['heel_r','mh1_r','mh3_r','mh5_r','hallux_r','othertoes_r',
              'heel_l','mh1_l','mh3_l','mh5_l','hallux_l','othertoes_l']
#Loop through sphere list
for ss in range(0,len(sphereList)):
    #Get the current sphere
    currSphere = osim.ContactSphere.safeDownCast(scaledModel.getContactGeometrySet().get(sphereList[ss]))
    #Set the current scaling factor based on sphere name
    if '_r' in sphereList[ss]:
        if 'hallux' in sphereList[ss] or 'othertoes' in sphereList[ss]:
            scaleFac = toesSphereScale_r
            scaleBod = 'toes_r'
        else:
            scaleFac = heelSphereScale_r
            scaleBod = 'calcn_r'
    elif '_l' in sphereList[ss]:
        if 'hallux' in sphereList[ss] or 'othertoes' in sphereList[ss]:
            scaleFac = toesSphereScale_l
            scaleBod = 'toes_l'
        else:
            scaleFac = heelSphereScale_l
            scaleBod = 'calcn_l'
    #Rescale the radius
    currSphere.setRadius(currSphere.getRadius() * scaleFac)
    #Set new location
    newLoc = []
    newLoc.append(currSphere.getLocation().get(0) * scaleFactors[scaleBod][0])
    newLoc.append(currSphere.getLocation().get(1) * scaleFactors[scaleBod][1])
    newLoc.append(currSphere.getLocation().get(2) * scaleFactors[scaleBod][2])
    currSphere.setLocation(osim.Vec3(newLoc[0],newLoc[1],newLoc[2]))    
    
#Reset model name
scaledModel.setName('scaledModel')

#Finalise connections and update scaled model
scaledModel.finalizeConnections()
scaledModel.printToXML('scaledModelAdjusted.osim')

#Scale muscle strength based on linear function presented in Handsfield
#et al. (2014). This uses some convenience functions that are packaged
#with the Rajagopal et al. (2016) gait model. Note that the height of
#the generic model is 1.700; the current participant height is 1.79
osimHelper.scaleOptimalForceSubjectSpecific(genericModelFileName = 'genericModel_LaiArnold.osim',
                                            scaledModelFileName = 'scaledModelAdjusted.osim',
                                            genericHeight = 1.700, scaledHeight = 1.79,
                                            outputModelFileName = 'scaledModelMuscle.osim')
#Load in new scaled model
scaledModelMuscle = osim.Model('scaledModelMuscle.osim')

# %% Set simulation parameters

#Identify time parameters for simulation
#This will be based off a half gait cycle of the right limb, that being from
#right foot heel strike to left foot heel strike
[startTime,endTime] = osimHelper.getHalfGaitCycle('refGRF.mot')

#Add the virtual torso, pelvis and hip joint markers to the .trc file
osimHelper.addVirtualMarkersDynamic(staticTRC = 'staticVirtualMarkers.trc',
                                    dynamicTRC = 'maxSprint.trc',
                                    outputTRC = 'maxSprintVirtualMarkers.trc')

# %% Inverse kinematics

#Initialise IK tool
ikTool = osim.InverseKinematicsTool()

#Set model
ikTool.setModel(scaledModelMuscle)

#Set task set
ikTool.set_IKTaskSet(osim.IKTaskSet('ikTasks.xml'))

#Set marker file
ikTool.set_marker_file('maxSprintVirtualMarkers.trc')

#Set times
ikTool.setStartTime(startTime)
ikTool.setEndTime(endTime)

#Set output filename
ikTool.set_output_motion_file('ikResults.mot')

#Run IK
ikTool.run()

#Convert IK results to a states file
osimHelper.kinematicsToStates(kinematicsFileName = 'ikResults.mot',
                              osimModelFileName = 'scaledModelMuscle.osim',
                              outputFileName = 'refQ.sto',
                              inDegrees = True, outDegrees = False)

###### TODO: check IK errors - some look high???

# %% Run a tracking simulation on experimental data

# This section contains adapted code from Ross Miller's UMocoD walking project

# At this point we convert the model and data to a 2D problem, so there are elements
# here that are necessary to do this

#Lock the frontal and transverse coordinates of the model
lockList = ['pelvis_list', 'pelvis_rotation', 'pelvis_tz',
            'hip_adduction_r', 'hip_rotation_r', 'hip_adduction_l', 'hip_rotation_l',
            'lumbar_bending', 'lumbar_rotation']
for kk in range(0,len(lockList)):
    scaledModelMuscle.getCoordinateSet().get(lockList[kk]).set_locked(True)
    
#Set lumbar bending and rotational torques to non-existent values
#Removing or disabling them generates errors from problem to solver
osim.CoordinateActuator.safeDownCast(scaledModelMuscle.getForceSet().get('tau_lumbar_bend')).setMinControl(-1e-10)
osim.CoordinateActuator.safeDownCast(scaledModelMuscle.getForceSet().get('tau_lumbar_bend')).setMaxControl(1e-10)
osim.CoordinateActuator.safeDownCast(scaledModelMuscle.getForceSet().get('tau_lumbar_bend')).setOptimalForce(1)
osim.CoordinateActuator.safeDownCast(scaledModelMuscle.getForceSet().get('tau_lumbar_rot')).setMinControl(-1e-10)
osim.CoordinateActuator.safeDownCast(scaledModelMuscle.getForceSet().get('tau_lumbar_rot')).setMaxControl(1e-10)
osim.CoordinateActuator.safeDownCast(scaledModelMuscle.getForceSet().get('tau_lumbar_rot')).setOptimalForce(1)
      
#Export model for processing
scaledModelMuscle.finalizeConnections()
scaledModelMuscle.printToXML('model2D.osim')
    
#Define the motion tracking problem
track = osim.MocoTrack()
track.setName('sprintTracking')

#Set kinematics
tableProcessor = osim.TableProcessor('refQ.sto')

#Set model and parameters
modelProcessor = osim.ModelProcessor('model2D.osim')
modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
modelProcessor.append(osim.ModOpTendonComplianceDynamicsModeDGF('implicit')) # Use muscle contractile dynamics
#modelProcessor.append(ModOpIgnorePassiveFiberForcesDGF()); # Set passive muscle fiber forces to zero
track.setModel(modelProcessor)

#Set states reference details
track.setStatesReference(tableProcessor) # Apply the target data to the tracking problem
track.set_states_global_tracking_weight(1.0) # Default tracking weight (is changed below)
track.set_allow_unused_references(True) # Target data can include DoF not in this model
track.set_track_reference_position_derivatives(True) # Track speed trajectories
track.set_apply_tracked_states_to_guess(True) # Use target data in initial guess

#Set times
track.set_initial_time(startTime)
track.set_final_time(endTime)

#Specify tracking weights as mean standard deviations from Miller et al. (2014)
#Pelvis and lumbar targets have arbitrarily large weights
stateWeights = osim.MocoWeightSet()
#Joint values
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tx/value',(1/(1*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_ty/value',(1/(2*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tilt/value',(1/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/value',(1/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/value',(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/walker_knee_r/knee_angle_r/value',(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/value',(1/(1*0.0574))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/value',(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/walker_knee_l/knee_angle_l/value',(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/value',(1/(1*0.0574))**2))
#Joint speeds
w = 0.001 #Scale the generalized speed tracking errors by this constant
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tx/speed',w*(1/(1*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_ty/speed',w*(1/(2*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tilt/speed',w*(1/(1*0.0585))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/speed',w*(1/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/speed',w*(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/walker_knee_r/knee_angle_r/speed',w*(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/speed',w*(1/(1*0.0574))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/speed',w*(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/walker_knee_l/knee_angle_l/speed',w*(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/speed',w*(1/(1*0.0574))**2))
#Add to tracking problem
track.set_states_weight_set(stateWeights)

#Define the Moco study and problem
study = track.initialize()
problem = study.updProblem()

#Define the periodicity goal
#Keeping in mind that the simulation is of a half gait cycle, similar to that
#proposed in Umberger's Moco examples
periodicityGoal = osim.MocoPeriodicityGoal('symmetryGoal')
problem.addGoal(periodicityGoal)
model = modelProcessor.process()
model.initSystem()

#Set symmetric coordinate values
#Set symmetry pairs
symPairs = [['hip_flexion_r','hip_flexion_l'],
            ['knee_angle_r','knee_angle_l'],
            ['ankle_angle_r','ankle_angle_l']]
for ii in range(0,len(symPairs)):
    #Set the jointsets depending on current pair
    if 'hip' in symPairs[ii][0]:
        jointName = ['/jointset/hip_r/','/jointset/hip_l/']
    elif 'knee' in symPairs[ii][0]:
        jointName = ['/jointset/walker_knee_r/','/jointset/walker_knee_l/']
    elif 'ankle' in symPairs[ii][0]:
        jointName = ['/jointset/ankle_r/','/jointset/ankle_l/']
    
    #Set the pair for coordinate value
    pair = osim.MocoPeriodicityGoalPair(jointName[0]+symPairs[ii][0]+'/value',
                                        jointName[1]+symPairs[ii][1]+'/value')
    #Add to the goal
    periodicityGoal.addStatePair(pair)
    
    #Set the pair for coordinate speed
    pair = osim.MocoPeriodicityGoalPair(jointName[0]+symPairs[ii][0]+'/speed',
                                        jointName[1]+symPairs[ii][1]+'/speed')
    #Add to the goal
    periodicityGoal.addStatePair(pair)
    
#Symmetric muscle activations
#Get muscle names
muscNames = list()
for mm in range(0,model.getMuscles().getSize()):
    muscNames.append(model.getMuscles().get(mm).getName()[0:-2])
#Get only unique names (i.e. remove double ups for both sides)
muscNames = list(np.unique(np.array(muscNames)))
#Loop through muscles and set periodicity pairs
for mm in range(0,len(muscNames)):
    #Set the pair for activations
    pair = osim.MocoPeriodicityGoalPair('/forceset/'+muscNames[mm]+'_r/activation',
                                        '/forceset/'+muscNames[mm]+'_l/activation')
    #Add to the goal
    periodicityGoal.addStatePair(pair)
    
    
    #Only set tendon force pair if tendon compliance is considered
    if model.getMuscles().get(muscNames[mm]+'_r').get_ignore_tendon_compliance() is False:
        #Set the pair for normalised tendon force
        pair = osim.MocoPeriodicityGoalPair('/forceset/'+muscNames[mm]+'_r/normalized_tendon_force',
                                            '/forceset/'+muscNames[mm]+'_l/normalized_tendon_force')
        #Add to the goal
        periodicityGoal.addStatePair(pair)
    
#Symmetric coordinate actuator controls
periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair('/forceset/tau_lumbar_ext'));

#Prescribed average sprint speed
#Calculate speed from kinematic data
df_kinematics = osimHelper.readSTO('refQ.sto')
startX = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'][df_kinematics.index[0]]
endX = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'][df_kinematics.index[-1]]
sprintSpeed = (endX - startX) / (df_kinematics.index[-1] - df_kinematics.index[0])

#Add the speed goal
speedGoal = osim.MocoAverageSpeedGoal('speed')
problem.addGoal(speedGoal)
speedGoal.set_desired_average_speed(sprintSpeed)

#Regularization term on MocoTrack problem (minimize squared muscle excitations)
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
effort.setWeight(0.1)

#Create GRF tracking goal
contactTracking = osim.MocoContactTrackingGoal('contact',1.0)
#Set external loads
contactTracking.setExternalLoadsFile('refGRF_2D.xml')
#Set right foot tracking
forcesRightFoot = osim.StdVectorString()
forcesRightFoot.append('/forceset/contactHeel_r')
forcesRightFoot.append('/forceset/contactMH1_r')
forcesRightFoot.append('/forceset/contactMH3_r')
forcesRightFoot.append('/forceset/contactMH5_r')
forcesRightFoot.append('/forceset/contactHallux_r')
forcesRightFoot.append('/forceset/contactOtherToes_r')
trackRightGRF = osim.MocoContactTrackingGoalGroup(forcesRightFoot,'RightGRF')
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r')
#Set left foot tracking
forcesLeftFoot = osim.StdVectorString()
forcesLeftFoot.append('/forceset/contactHeel_l')
forcesLeftFoot.append('/forceset/contactMH1_l')
forcesLeftFoot.append('/forceset/contactMH3_l')
forcesLeftFoot.append('/forceset/contactMH5_l')
forcesLeftFoot.append('/forceset/contactHallux_l')
forcesLeftFoot.append('/forceset/contactOtherToes_l')
trackLeftGRF = osim.MocoContactTrackingGoalGroup(forcesLeftFoot,'LeftGRF')
trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l')
#ADd to problem
contactTracking.addContactGroup(trackRightGRF)
contactTracking.addContactGroup(trackLeftGRF)
#Set parameters in problem
contactTracking.setProjection('plane')
contactTracking.setProjectionVector(osim.Vec3(0, 0, 1))
#Add to problem
problem.addGoal(contactTracking)

#Set joint coordinate bounds
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tx/value', [0, 5])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_ty/value', [0.75, 1.25])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tilt/value', [math.radians(-20), math.radians(10)])
problem.setStateInfo('/jointset/back/lumbar_extension/value', [math.radians(-30), math.radians(5)])
problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [math.radians(-30), math.radians(90)])
problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [math.radians(-30), math.radians(90)])
problem.setStateInfo('/jointset/walker_knee_r/knee_angle_r/value', [math.radians(0), math.radians(140)])
problem.setStateInfo('/jointset/walker_knee_l/knee_angle_l/value', [math.radians(0), math.radians(140)])
problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [math.radians(-40), math.radians(30)])
problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [math.radians(-40), math.radians(30)])

#Set muscle limits in problem
#### TODO: these may need to be adjusted for sprinting (from walking data)
problem.setStateInfoPattern('/forceset/.*/normalized_tendon_force', [0, 1.8], [], [])
problem.setStateInfoPattern('/forceset/.*/activation',   [0.001, 1.0], [], [])

#Define the solver and set its options
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
#solver.set_multibody_dynamics_mode('implicit')
#solver.set_minimize_implicit_multibody_accelerations(true)
#solver.set_implicit_multibody_accelerations_weight(0.00001)
solver.set_optim_max_iterations(1000)
solver.set_num_mesh_intervals(25) ####TODO: might need to be upped, low for speed right now
solver.set_optim_constraint_tolerance(1e-2) ####TODO: might need to be lower
solver.set_optim_convergence_tolerance(1e-2) ####TODO: might need to be lower
solver.set_minimize_implicit_auxiliary_derivatives(True)
solver.set_implicit_auxiliary_derivatives_weight(0.001)
solver.resetProblem(problem)

#Set the normalized tendon forces if not loading initial guess from file
guess = solver.getGuess()
numRows = guess.getNumTimes()
stateNames = model.getStateVariableNames()
for ii in range(0,model.getNumStateVariables()):
    #Get current state name
    currentStateName = stateNames.get(ii)
    if 'normalized_tendon_force' in currentStateName:
        guess.setState(currentStateName, np.linspace(0.2,0.2,numRows))

# %% Solve the tracking problem

#Solve!
sprintTrackingSolution = study.solve()

study.printToXML('test.omoco')

# % Write the solution to a file
# sprintTrackingSolution.write('sprintTracking_solution_halfStride.sto')

# % Visualize the solution
# study.visualize(gaitTrackingSolution);


#### TODO: must unlock model coordinates to avoid NaN's in force table...

forcesLeftFoot.append('/forceset/contactHeel_l')
forcesLeftFoot.append('/forceset/contactMH1_l')
forcesLeftFoot.append('/forceset/contactMH3_l')
forcesLeftFoot.append('/forceset/contactMH5_l')
forcesLeftFoot.append('/forceset/contactHallux_l')
forcesLeftFoot.append('/forceset/contactOtherToes_l')

#Write solution's GRF to a file
contact_r = osim.StdVectorString()
contact_l = osim.StdVectorString()
contact_r.append('/forceset/contactHeel_r')
contact_r.append('/forceset/contactMH1_r')
contact_r.append('/forceset/contactMH3_r')
contact_r.append('/forceset/contactMH5_r')
contact_r.append('/forceset/contactHallux_r')
contact_r.append('/forceset/contactOtherToes_r')
contact_l.append('/forceset/contactHeel_l')
contact_l.append('/forceset/contactMH1_l')
contact_l.append('/forceset/contactMH3_l')
contact_l.append('/forceset/contactMH5_l')
contact_l.append('/forceset/contactHallux_l')
contact_l.append('/forceset/contactOtherToes_l')
externalForcesTableFlat = osim.createExternalLoadsTableForGait(model,sprintTrackingSolution,contact_r,contact_l)
osim.writeTableToFile(externalForcesTableFlat,'sprintTracking_solution_halfStride_GRF.sto')


# %% DON'T THINK RRA IS NECESSARY GIVEN CHANGES TO KINEMATICS WILL OCCUR DURING TRACKING SIM...?

# # %% Residual reduction algorithm

# #Initialise RRA tool
# rraTool = osim.RRATool()

# #Set general settings outside of loop

# #Set tool to replace model force set
# rraTool.setReplaceForceSet(True)

# #Set actuators file
# forceSetFiles = osim.ArrayStr()
# forceSetFiles.set(0,'..\\rraActuators.xml')
# rraTool.setForceSetFiles(forceSetFiles)

# #Set tracking tasks file
# rraTool.setTaskSetFileName('..\\rraTasks.xml')

# #Set a low pass filter frequency on the kinematics data
# rraTool.setLowpassCutoffFrequency(12)

# #Set to adjust the COM to reduce residuals
# rraTool.setAdjustCOMToReduceResiduals(True)

# #Set the torso body COM to be the one that gets adjusted
# rraTool.setAdjustedCOMBody('torso')

# #Set external loads file
# rraTool.setExternalLoadsFileName('..\\refGRF.xml')

# #Loop through and perform rra three times
# #Each time we'll adjust the mass specified by the RRA tool, and then re-use
# #this adjusted model in the next iteration. On the first iteration, we'll
# #use the IK motion data - but on subsequent iterations we'll use the
# #adjusted rra kinematics. Similarly, on the first iteration we'll set the
# #model to our scaled model - but then use the adjusted version on
# #subsequent iterations.
# for rr in range(3):
    
#     #Make directory for current iteration
#     os.mkdir('rra'+str(rr+1))
#     os.chdir('rra'+str(rr+1))

#     #Set tool name based on iteration
#     rraTool.setName('rra'+str(rr+1))
    
#     #Set desired kinematics file
#     if rr == 0:

#         #Use IK data
#         rraTool.setDesiredKinematicsFileName('..\\ikResults.mot')

#         #Set initial and final time
#         #You need to use the IK first and last times here as I don't think the tool
#         #likes if the IK data doesn't have the times you put in
#         rraTool.setStartTime(osim.Storage('..\\ikResults.mot').getFirstTime())
#         rraTool.setFinalTime(osim.Storage('..\\ikResults.mot').getLastTime())

#         #Set to original scaled model
#         rraTool.setModelFilename('..\\scaledModelMuscle.osim')
        
#     #Set output model file
#     rraTool.setOutputModelFileName('rraAdjustedModel_'+str(rr+1)+'.osim')

#     #Print out the tool
#     rraTool.printToXML('setupRRA'+str(rr+1)+'.xml')

#     #Run RRA
#     #Reloading the setup tool seems to help with running the tool smoothly
#     runRRA = osim.RRATool('setupRRA'+str(rr+1)+'.xml')
#     runRRA.run()
    








# %%
