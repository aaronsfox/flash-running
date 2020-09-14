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
from scipy.spatial.transform import Rotation as R
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

#Load in the c3d data via btk
dynamicC3D = btk.btkAcquisitionFileReader()
dynamicC3D.SetFilename('maxSprint.c3d')
dynamicC3D.Update()
dynamicAcq = dynamicC3D.GetOutput()

#Extract the force platforms data
forcePlatforms = btk.btkForcePlatformsExtractor()
forcePlatforms.SetInput(dynamicAcq)
forcePlatforms.Update()

#Get the wrenchs for position and force data
groundReactionWrenches = btk.btkGroundReactionWrenchFilter()
groundReactionWrenches.SetInput(forcePlatforms.GetOutput())
groundReactionWrenches.Update()

#Loop through the 8 plates and get their info
grf = list()
grm = list()
cop = list()
for fp in range(0,8):
    #Get the different sets of values and append to lists
    grf.append(groundReactionWrenches.GetOutput().GetItem(fp).GetForce().GetValues())
    grm.append(groundReactionWrenches.GetOutput().GetItem(fp).GetMoment().GetValues())
    cop.append(groundReactionWrenches.GetOutput().GetItem(fp).GetPosition().GetValues())

#Rotate force plate data
#First rotation
rotMag = np.radians(-90)
rotAx = np.array([0,0,1])
rotVec = rotMag * rotAx
rotMat1 = R.from_rotvec(rotVec)
for fp in range(0,8):
    grf[fp] = rotMat1.apply(grf[fp])
    grm[fp] = rotMat1.apply(grm[fp])
    cop[fp] = rotMat1.apply(cop[fp])
#Second rotation
rotMag = np.radians(-90)
rotAx = np.array([1,0,0])
rotVec = rotMag * rotAx
rotMat2 = R.from_rotvec(rotVec)
for fp in range(0,8):
    grf[fp] = rotMat2.apply(grf[fp])
    grm[fp] = rotMat2.apply(grm[fp])
    cop[fp] = rotMat2.apply(cop[fp])

#Convert mm to m units
for fp in range(0,8):
    grm[fp] = grm[fp] / 1000
    cop[fp] = cop[fp] / 1000

#Extract data into matrix
forcesMat = np.zeros((len(grf[0]), 3 * 3 * 8))
colNo = 0
for fp in range(0,8):
    forcesMat[:,colNo:colNo+3] = grf[fp]
    forcesMat[:,colNo+3:colNo+3+3] = cop[fp]
    forcesMat[:,colNo+6:colNo+6+3] = grm[fp]
    colNo = colNo+9

#Define filter for forces data (50N Low-pass 4th Order Butterworth)
fs = 2000
nyq = 0.5 * fs
normCutoff = 50 / nyq
b, a = butter(4, normCutoff, btype = 'low', analog = False)

#Apply lowpass filter to force data (50N cut-off)
forcesFilt = np.zeros((np.size(forcesMat,0), np.size(forcesMat,1)))
for jj in range(0,np.size(forcesMat,1)):
    forcesFilt[:,jj] = lfilter(b, a, forcesMat[:,jj])
    
#Cast array to dataframe
#Set force labels
forceLabels = ['Fx','Fy','Fz','Px','Py','Pz','Mx','My','Mz']
forceNames = list()
for ff in range(0,8):
    for pp in range(0,len(forceLabels)):
        forceNames.append(forceLabels[pp]+str(ff+1))
#Place in dataframe
df_forces = pd.DataFrame(forcesFilt, columns = forceNames)

#Convert all force data below a 20N threshold to zeros

#Set vertical ground reaction force names
vgrfNames = list()
for ff in range(0,8):
    vgrfNames.append('Fy'+str(ff+1))
    
#Loop through and correct data
for vv in range(0,len(vgrfNames)):
    
    #Extract force data
    currVGRF = df_forces[vgrfNames[vv]].values
    
    #Get boolean of force threshold
    offForce = currVGRF < 20
    
    #Set any relevant column data to zero
    #Forces
    df_forces.loc[offForce, 'Fx'+str(vv+1)] = 0
    df_forces.loc[offForce, 'Fy'+str(vv+1)] = 0
    df_forces.loc[offForce, 'Fz'+str(vv+1)] = 0
    #COP
    df_forces.loc[offForce, 'Px'+str(vv+1)] = 0
    df_forces.loc[offForce, 'Py'+str(vv+1)] = 0
    df_forces.loc[offForce, 'Pz'+str(vv+1)] = 0
    #Moments
    df_forces.loc[offForce, 'Mx'+str(vv+1)] = 0
    df_forces.loc[offForce, 'My'+str(vv+1)] = 0
    df_forces.loc[offForce, 'Mz'+str(vv+1)] = 0
    
#Extract the forces relating to each foot
#This is a manual process given the singular trial and knowledge of which plates
#have which contacts
#Right foot = fp7,fp3 
#Left foot = fp5,fp1
#Right forces
ground_force_r_vx = df_forces['Fx7'].values + df_forces['Fx3'].values
ground_force_r_vy = df_forces['Fy7'].values + df_forces['Fy3'].values
ground_force_r_vz = df_forces['Fz7'].values + df_forces['Fz3'].values
#Left forces
ground_force_l_vx = df_forces['Fx5'].values + df_forces['Fx1'].values
ground_force_l_vy = df_forces['Fy5'].values + df_forces['Fy1'].values
ground_force_l_vz = df_forces['Fz5'].values + df_forces['Fz1'].values
#Right torques
ground_torque_r_x = df_forces['Mx7'].values + df_forces['Mx3'].values
ground_torque_r_y = df_forces['My7'].values + df_forces['My3'].values
ground_torque_r_z = df_forces['Mz7'].values + df_forces['Mz3'].values
#Left torques
ground_torque_l_x = df_forces['Mx5'].values + df_forces['Mx1'].values
ground_torque_l_y = df_forces['My5'].values + df_forces['My1'].values
ground_torque_l_z = df_forces['Mz5'].values + df_forces['Mz1'].values
#Right position
ground_force_r_px = df_forces['Px7'].values + df_forces['Px3'].values
ground_force_r_py = df_forces['Py7'].values + df_forces['Py3'].values
ground_force_r_pz = df_forces['Pz7'].values + df_forces['Pz3'].values
#Left position
ground_force_l_px = df_forces['Px5'].values + df_forces['Px1'].values
ground_force_l_py = df_forces['Py5'].values + df_forces['Py1'].values
ground_force_l_pz = df_forces['Pz5'].values + df_forces['Pz1'].values

#Export forces to .mot file
#Set labels for file
motLabels = ['time',
             'ground_force_r_vx', 'ground_force_r_vy', 'ground_force_r_vz',
             'ground_force_r_px', 'ground_force_r_py', 'ground_force_r_pz',
             'ground_torque_r_x', 'ground_torque_r_y', 'ground_torque_r_z',
             'ground_force_l_vx', 'ground_force_l_vy', 'ground_force_l_vz',
             'ground_force_l_px', 'ground_force_l_py', 'ground_force_l_pz',
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
                                   ground_force_r_px, ground_force_r_py, ground_force_r_pz,
                                   ground_torque_r_x, ground_torque_r_y, ground_torque_r_z,
                                   ground_force_l_vx, ground_force_l_vy, ground_force_l_vz,
                                   ground_force_l_px, ground_force_l_py, ground_force_l_pz,
                                   ground_torque_l_x, ground_torque_l_y, ground_torque_l_z]))

#Append data to storage object
nrow, ncol = forceData.shape
#Create timestamp based on sampling rate and length
forcesTime = np.linspace(0,nrow/fs,num = nrow)
#Add data
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
                                     ground_force_r_px, ground_force_r_py, np.zeros([nrow]),
                                     ground_torque_r_x, ground_torque_r_y, np.zeros([nrow]),
                                     ground_force_l_vx, ground_force_l_vy, np.zeros([nrow]),
                                     ground_force_l_px, ground_force_l_py, np.zeros([nrow]),
                                     ground_torque_l_x, ground_torque_l_y, np.zeros([nrow])]))
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

# %% Run an inverse simulation to on experimental data

# # At this point we convert the model and data to a 2D problem, so there are elements
# # here that are necessary to do this

# #Lock the frontal and transverse coordinates of the model
# lockList = ['pelvis_list', 'pelvis_rotation', 'pelvis_tz',
#             'hip_adduction_r', 'hip_rotation_r', 'hip_adduction_l', 'hip_rotation_l',
#             'lumbar_bending', 'lumbar_rotation']
# for kk in range(0,len(lockList)):
#     scaledModelMuscle.getCoordinateSet().get(lockList[kk]).set_locked(True)
    
# #Set lumbar bending and rotational torques to non-existent values
# #Removing or disabling them generates errors from problem to solver
# osim.CoordinateActuator.safeDownCast(scaledModelMuscle.getForceSet().get('tau_lumbar_bend')).setMinControl(-1e-10)
# osim.CoordinateActuator.safeDownCast(scaledModelMuscle.getForceSet().get('tau_lumbar_bend')).setMaxControl(1e-10)
# osim.CoordinateActuator.safeDownCast(scaledModelMuscle.getForceSet().get('tau_lumbar_bend')).setOptimalForce(1)
# osim.CoordinateActuator.safeDownCast(scaledModelMuscle.getForceSet().get('tau_lumbar_rot')).setMinControl(-1e-10)
# osim.CoordinateActuator.safeDownCast(scaledModelMuscle.getForceSet().get('tau_lumbar_rot')).setMaxControl(1e-10)
# osim.CoordinateActuator.safeDownCast(scaledModelMuscle.getForceSet().get('tau_lumbar_rot')).setOptimalForce(1)
      
# #Export model for processing
# scaledModelMuscle.finalizeConnections()
# scaledModelMuscle.printToXML('model2D.osim')

# #Construct the MocoInverse tool.
# inverse = osim.MocoInverse()

# #Set model and parameters
# modelProcessor = osim.ModelProcessor('model2D.osim')
# modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
# modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
# modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
# # modelProcessor.append(osim.ModOpTendonComplianceDynamicsModeDGF('implicit')) # Use muscle contractile dynamics
# #modelProcessor.append(ModOpIgnorePassiveFiberForcesDGF()); # Set passive muscle fiber forces to zero
# modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
# modelProcessor.append(osim.ModOpAddReserves(1.0))
# inverse.setModel(modelProcessor)

# #Construct a TableProcessor of the coordinate data and pass it to the
# #inverse tool. TableProcessors can be used in the same way as
# #ModelProcessors by appending TableOperators to modify the base table.
# #A TableProcessor with no operators, as we have here, simply returns the
# #base table.
# inverse.setKinematics(osim.TableProcessor('refQ.sto'))

# # Initial time, final time, and mesh interval.
# inverse.set_initial_time(startTime)
# inverse.set_final_time(endTime)
# inverse.set_mesh_interval(0.02) #### currently just adopted from example

# #By default, Moco gives an error if the kinematics contains extra columns.
# #Here, we tell Moco to allow (and ignore) those extra columns.
# inverse.set_kinematics_allow_extra_columns(True)

# #Solve the problem and write the solution to a Storage file.
# inverseSolution = inverse.solve()
# # inverseSolution.getMocoSolution().unseal()
# inverseSolution.getMocoSolution().write('stopped_MocoInverse_solution.sto')

# %% Run a 3D tracking simulation without GRF tracking

#Turn off the contact force elements
osim.SmoothSphereHalfSpaceForce.safeDownCast(scaledModelMuscle.getForceSet().get('contactHeel_r')).set_appliesForce(False)
osim.SmoothSphereHalfSpaceForce.safeDownCast(scaledModelMuscle.getForceSet().get('contactMH1_r')).set_appliesForce(False)
osim.SmoothSphereHalfSpaceForce.safeDownCast(scaledModelMuscle.getForceSet().get('contactMH3_r')).set_appliesForce(False)
osim.SmoothSphereHalfSpaceForce.safeDownCast(scaledModelMuscle.getForceSet().get('contactMH5_r')).set_appliesForce(False)
osim.SmoothSphereHalfSpaceForce.safeDownCast(scaledModelMuscle.getForceSet().get('contactHallux_r')).set_appliesForce(False)
osim.SmoothSphereHalfSpaceForce.safeDownCast(scaledModelMuscle.getForceSet().get('contactOtherToes_r')).set_appliesForce(False)
osim.SmoothSphereHalfSpaceForce.safeDownCast(scaledModelMuscle.getForceSet().get('contactHeel_l')).set_appliesForce(False)
osim.SmoothSphereHalfSpaceForce.safeDownCast(scaledModelMuscle.getForceSet().get('contactMH1_l')).set_appliesForce(False)
osim.SmoothSphereHalfSpaceForce.safeDownCast(scaledModelMuscle.getForceSet().get('contactMH3_l')).set_appliesForce(False)
osim.SmoothSphereHalfSpaceForce.safeDownCast(scaledModelMuscle.getForceSet().get('contactMH5_l')).set_appliesForce(False)
osim.SmoothSphereHalfSpaceForce.safeDownCast(scaledModelMuscle.getForceSet().get('contactHallux_l')).set_appliesForce(False)
osim.SmoothSphereHalfSpaceForce.safeDownCast(scaledModelMuscle.getForceSet().get('contactOtherToes_l')).set_appliesForce(False)

#Export model for processing
scaledModelMuscle.finalizeConnections()
scaledModelMuscle.printToXML('model3D_noContact.osim')

#Define the motion tracking problem
track = osim.MocoTrack()
track.setName('sprintTracking3D_noContact')

#Set kinematics
tableProcessor = osim.TableProcessor('refQ.sto')

#Set model and parameters
modelProcessor = osim.ModelProcessor('model3D_noContact.osim')
modelProcessor.append(osim.ModOpAddExternalLoads('refGRF.xml'))
modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
modelProcessor.append(osim.ModOpAddReserves(2))
# modelProcessor.append(osim.ModOpTendonComplianceDynamicsModeDGF('implicit')) # Use muscle contractile dynamics
#modelProcessor.append(ModOpIgnorePassiveFiberForcesDGF()); # Set passive muscle fiber forces to zero
track.setModel(modelProcessor)

#Set states reference details
track.setStatesReference(tableProcessor) # Apply the target data to the tracking problem
track.set_states_global_tracking_weight(5) # Default tracking weight (is changed below)
track.set_allow_unused_references(True) # Target data can include DoF not in this model
track.set_track_reference_position_derivatives(True) # Track speed trajectories
track.set_apply_tracked_states_to_guess(True) # Use target data in initial guess

#Set times
track.set_initial_time(osim.Storage('refQ.sto').getFirstTime())
track.set_final_time(osim.Storage('refQ.sto').getLastTime())

#Specify tracking weights as mean standard deviations from Miller et al. (2014)
#Pelvis and lumbar targets have arbitrarily large weights
stateWeights = osim.MocoWeightSet()
#Joint values
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tx/value',10))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_ty/value',10))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tz/value',10))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tilt/value',50))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_list/value',50))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_rotation/value',50))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/value',5))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_bending/value',5))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_rotation/value',5))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/value',25))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_adduction_r/value',10))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_rotation_r/value',5))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/walker_knee_r/knee_angle_r/value',50))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/value',25))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/value',25))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_adduction_l/value',10))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_rotation_l/value',5))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/walker_knee_l/knee_angle_l/value',50))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/value',25))
#Add to tracking problem
track.set_states_weight_set(stateWeights)

#Define the Moco study and problem
study = track.initialize()
problem = study.updProblem()

#Regularization term on MocoTrack problem (minimize squared muscle excitations)
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
effort.setWeight(1)

#Set joint coordinate bounds
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tx/value', [0, 5])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_ty/value', [0.75, 1.25])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tz/value', [0, 0.3])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tilt/value', [math.radians(-20), math.radians(10)])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_list/value', [math.radians(-20), math.radians(20)])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_rotation/value', [math.radians(-20), math.radians(20)])
problem.setStateInfo('/jointset/back/lumbar_extension/value', [math.radians(-30), math.radians(5)])
problem.setStateInfo('/jointset/back/lumbar_bending/value', [math.radians(-20), math.radians(20)])
problem.setStateInfo('/jointset/back/lumbar_rotation/value', [math.radians(-10), math.radians(25)])
problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [math.radians(-30), math.radians(90)])
problem.setStateInfo('/jointset/hip_r/hip_adduction_r/value', [math.radians(-20), math.radians(15)])
problem.setStateInfo('/jointset/hip_r/hip_rotation_r/value', [math.radians(-30), math.radians(15)])
problem.setStateInfo('/jointset/hip_l/hip_adduction_l/value', [math.radians(-20), math.radians(15)])
problem.setStateInfo('/jointset/hip_l/hip_rotation_l/value', [math.radians(-30), math.radians(15)])
problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [math.radians(-30), math.radians(90)])
problem.setStateInfo('/jointset/walker_knee_r/knee_angle_r/value', [math.radians(0), math.radians(140)])
problem.setStateInfo('/jointset/walker_knee_l/knee_angle_l/value', [math.radians(0), math.radians(140)])
problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [math.radians(-40), math.radians(30)])
problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [math.radians(-40), math.radians(30)])

#Set muscle limits in problem
#### TODO: these may need to be adjusted for sprinting (from walking data)
# problem.setStateInfoPattern('/forceset/.*/normalized_tendon_force', [0, 1.8], [], [])
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
# solver.set_minimize_implicit_auxiliary_derivatives(True)
# solver.set_implicit_auxiliary_derivatives_weight(0.001)
solver.resetProblem(problem)

# %% Solve the tracking problem without contact

#Solve!
sprintTrackingSolution3D_noContact = study.solve()

# %% Run a tracking simulation on experimental data

# This section contains adapted code from Ross Miller's UMocoD walking project

# At this point we convert the model and data to a 2D problem, so there are elements
# here that are necessary to do this

#Re-import model
model2D = osim.Model('scaledModelMuscle.osim')

#Lock the frontal and transverse coordinates of the model
lockList = ['pelvis_list', 'pelvis_rotation', 'pelvis_tz',
            'hip_adduction_r', 'hip_rotation_r', 'hip_adduction_l', 'hip_rotation_l',
            'lumbar_bending', 'lumbar_rotation']
# for kk in range(0,len(lockList)):
#     model2D.getCoordinateSet().get(lockList[kk]).set_locked(True)

#Unlock everything
for kk in range(0,model2D.getCoordinateSet().getSize()):
    model2D.getCoordinateSet().get(kk).set_locked(False)    
    
# #Set lumbar bending and rotational torques to non-existent values
# #Removing or disabling them generates errors from problem to solver
# osim.CoordinateActuator.safeDownCast(model2D.getForceSet().get('tau_lumbar_bend')).setMinControl(-1e-10)
# osim.CoordinateActuator.safeDownCast(model2D.getForceSet().get('tau_lumbar_bend')).setMaxControl(1e-10)
# osim.CoordinateActuator.safeDownCast(model2D.getForceSet().get('tau_lumbar_bend')).setOptimalForce(1)
# osim.CoordinateActuator.safeDownCast(model2D.getForceSet().get('tau_lumbar_rot')).setMinControl(-1e-10)
# osim.CoordinateActuator.safeDownCast(model2D.getForceSet().get('tau_lumbar_rot')).setMaxControl(1e-10)
# osim.CoordinateActuator.safeDownCast(model2D.getForceSet().get('tau_lumbar_rot')).setOptimalForce(1)
      
#Export model for processing
model2D.finalizeConnections()
model2D.printToXML('model2D.osim')
    
#Define the motion tracking problem
track = osim.MocoTrack()
track.setName('sprintTracking_2D_withContact')

#Set kinematics
tableProcessor = osim.TableProcessor('refQ.sto')

#Set model and parameters
modelProcessor = osim.ModelProcessor('model2D.osim')
modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
# modelProcessor.append(osim.ModOpTendonComplianceDynamicsModeDGF('implicit')) # Use muscle contractile dynamics
#modelProcessor.append(ModOpIgnorePassiveFiberForcesDGF()); # Set passive muscle fiber forces to zero
track.setModel(modelProcessor)

#Set states reference details
track.setStatesReference(tableProcessor) # Apply the target data to the tracking problem
track.set_states_global_tracking_weight(1.0) # Default tracking weight (is changed below)
track.set_allow_unused_references(True) # Target data can include DoF not in this model
track.set_track_reference_position_derivatives(True) # Track speed trajectories
track.set_apply_tracked_states_to_guess(True) # Use target data in initial guess

#Set times
track.set_initial_time(osim.Storage('refQ.sto').getFirstTime())
track.set_final_time(osim.Storage('refQ.sto').getLastTime())

#Specify tracking weights as mean standard deviations from Miller et al. (2014)
#Pelvis and lumbar targets have arbitrarily large weights
#Set 'locked' coordinates to zero weight
stateWeights = osim.MocoWeightSet()
#Joint values
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tx/value',(1/(1*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_ty/value',(1/(2*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tz/value',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tilt/value',(1/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_list/value',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_rotation/value',0))                    
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/value',(1/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_bending/value',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_rotation/value',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/value',(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_adduction_r/value',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_rotation_r/value',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/walker_knee_r/knee_angle_r/value',(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/value',(1/(1*0.0574))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/value',(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_adduction_l/value',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_rotation_l/value',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/walker_knee_l/knee_angle_l/value',(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/value',(1/(1*0.0574))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/subtalar_r/subtalar_angle_r/value',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/subtalar_l/subtalar_angle_l/value',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/mtp_r/mtp_angle_r/value',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/mtp_l/mtp_angle_l/value',0))
#Joint speeds
w = 0.001 #Scale the generalized speed tracking errors by this constant
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tx/speed',w*(1/(1*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_ty/speed',w*(1/(2*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tz/speed',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tilt/speed',w*(1/(1*0.0585))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_list/speed',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_rotation/speed',0))    
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/speed',w*(1/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_bending/speed',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_rotation/speed',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/speed',w*(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_adduction_r/speed',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_rotation_r/speed',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/walker_knee_r/knee_angle_r/speed',w*(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/speed',w*(1/(1*0.0574))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/speed',w*(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_adduction_l/speed',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_rotation_l/speed',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/walker_knee_l/knee_angle_l/speed',w*(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/speed',w*(1/(1*0.0574))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/subtalar_r/subtalar_angle_r/speed',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/subtalar_l/subtalar_angle_l/speed',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/mtp_r/mtp_angle_r/speed',0))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/mtp_l/mtp_angle_l/speed',0))
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
    
    
    # #Only set tendon force pair if tendon compliance is considered
    # if model.getMuscles().get(muscNames[mm]+'_r').get_ignore_tendon_compliance() is False:
    #     #Set the pair for normalised tendon force
    #     pair = osim.MocoPeriodicityGoalPair('/forceset/'+muscNames[mm]+'_r/normalized_tendon_force',
    #                                         '/forceset/'+muscNames[mm]+'_l/normalized_tendon_force')
    #     #Add to the goal
    #     periodicityGoal.addStatePair(pair)
    
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
#Add to problem
contactTracking.addContactGroup(trackRightGRF)
contactTracking.addContactGroup(trackLeftGRF)
#Set parameters in problem
contactTracking.setProjection('plane')
contactTracking.setProjectionVector(osim.Vec3(0, 0, 1))
#Add to problem
problem.addGoal(contactTracking)

#Set joint coordinate bounds
#This includes artificially 'locking' certain joints
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tx/value', [0, 5])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_ty/value', [0.75, 1.25])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tz/value', [0])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tilt/value', [math.radians(-15), math.radians(5)])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_list/value', [0])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_rotation/value', [0])
problem.setStateInfo('/jointset/back/lumbar_extension/value', [math.radians(-25), math.radians(5)])
problem.setStateInfo('/jointset/back/lumbar_bending/value', [0])
problem.setStateInfo('/jointset/back/lumbar_rotation/value', [0])
problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [math.radians(-30), math.radians(90)])
problem.setStateInfo('/jointset/hip_r/hip_adduction_r/value', [0])
problem.setStateInfo('/jointset/hip_r/hip_rotation_r/value', [0])
problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [math.radians(-30), math.radians(90)])
problem.setStateInfo('/jointset/hip_l/hip_adduction_l/value', [0])
problem.setStateInfo('/jointset/hip_l/hip_rotation_l/value', [0])
problem.setStateInfo('/jointset/walker_knee_r/knee_angle_r/value', [math.radians(0), math.radians(140)])
problem.setStateInfo('/jointset/walker_knee_l/knee_angle_l/value', [math.radians(0), math.radians(140)])
problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [math.radians(-40), math.radians(30)])
problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [math.radians(-40), math.radians(30)])
problem.setStateInfo('/jointset/subtalar_r/subtalar_angle_r/value', [0])
problem.setStateInfo('/jointset/mtp_r/mtp_angle_r/value', [0])
problem.setStateInfo('/jointset/subtalar_l/subtalar_angle_l/value', [0])
problem.setStateInfo('/jointset/mtp_l/mtp_angle_l/value', [0])

#Set muscle limits in problem
#### TODO: these may need to be adjusted for sprinting (from walking data)
# problem.setStateInfoPattern('/forceset/.*/normalized_tendon_force', [0, 1.8], [], [])
problem.setStateInfoPattern('/forceset/.*/activation',   [0.001, 1.0], [], [])


######

# #Create rotation time series table
# rotTable = osim.TimeSeriesTableRotation()

# #Set column labels for pelvis bodyset
# rotLabels = osim.StdVectorString()
# rotLabels.append('/bodyset/pelvis')
# rotTable.setColumnLabels(rotLabels)

# #Create a zero rotation value
# zeroRot = osim.Rotation()
# zeroRot.setToZero()

# #Set the data to zeros
# #Use length of reference kinematic data as a reference
# nrows = len(osim.TimeSeriesTable('refQ.sto').getIndependentColumn())
# for ii in range(nrows):
#     #Get an opensim row vector rotation with a zero value
#     row = osim.RowVectorRotation(1,zeroRot)
#     #Append to the table
#     rotTable.appendRow(ii,row)
    
# #Set the time data
# for ii in range(nrows):
#     rotTable.setIndependentValueAtIndex(ii,osim.TimeSeriesTable('refQ.sto').getIndependentColumn()[ii])

# # Add rotation table of zeros to goal
# pelvisGoal.setRotationReference(rotTable)   ###TABLE WON'T WORK!!!! )


######


#Add low weight goal to keep pelvis orientation vertical
#Create goal
pelvisGoal = osim.MocoOrientationTrackingGoal('pelvisOrientation',0.1)
#Set path to pelvis frame
framePaths = osim.StdVectorString()
framePaths.append('/bodyset/pelvis')
pelvisGoal.setFramePaths(framePaths)

#Create a states rotation reference table of default values for neutral pelvis orientation
#Start with the previous solution base
pelvisO = osim.MocoTrajectory('sprintTracking3D_noContact_solution.sto')
#Set all states data to zero
numRows = pelvisO.getNumTimes()
stateNames = model.getStateVariableNames()
for ii in range(0,model.getNumStateVariables()):
    currState = stateNames.get(ii)
    if 'jointset' in currState:
        pelvisO.setState(currState, np.linspace(0.0,0.0,numRows))
    elif 'forceset' in currState:
        pelvisO.setState(currState, np.linspace(0.02,0.02,numRows))
#Export to states file and set as reference file
pelvisSTO = osim.STOFileAdapter()
pelvisSTO.write(pelvisO.exportToStatesTable(),'neutralOrientationStates.sto')
pelvisGoal.setStatesReference(osim.TableProcessor('neutralOrientationStates.sto'))
#Add to problem
problem.addGoal(pelvisGoal)

#Define the solver and set its options
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
#solver.set_multibody_dynamics_mode('implicit')
#solver.set_minimize_implicit_multibody_accelerations(true)
#solver.set_implicit_multibody_accelerations_weight(0.00001)
solver.set_optim_max_iterations(1000)
solver.set_num_mesh_intervals(50) ####TODO: upped from 25
solver.set_optim_constraint_tolerance(1e-2) ####TODO: might need to be lower
solver.set_optim_convergence_tolerance(1e-2) ####TODO: might need to be lower
solver.set_minimize_implicit_auxiliary_derivatives(True)
solver.set_implicit_auxiliary_derivatives_weight(0.001)
solver.resetProblem(problem)

# ##### Use existing try guess for now...
# solver.setGuessFile('sprintTracking_2D_withContact_solution_try1.sto')

##### guess won't work with 'unlocked' approach...

#Create a blank guess to fill with relevant 3D tracking info
guess = solver.getGuess()

#Grab the previous solution as a fresh MocoTrajectory
##### A bit wonky --- but should do OK
traj = osim.MocoTrajectory('sprintTracking3D_noContact_solution.sto')

#Convert the blank guess to the same number of nodes as the solution trajectory
guess.resampleWithNumTimes(traj.getNumTimes())

#Get the state names of the guess to work through
guessStates = guess.getStateNames()

#Loop through and fill states from solution
#Ensure that altered/locked coordinates remain as zero
for ss in range(0,len(guessStates)):
    #Check for newly locked coordinates
    if any(checkList in guessStates[ss] for checkList in lockList):
        #Set as zeros
        guess.setState(guessStates[ss], np.linspace(0.0,0.0,traj.getNumTimes()))
    else:
        #Copy the state from the trajectory
        guess.setState(guessStates[ss], traj.getStateMat(guessStates[ss]))

#Get the control names of the guess to work through
guessControls = guess.getControlNames()

#Loop through and fill controls from solution
#Ensure that locked lumbar controls drop to zero
for cc in range(0,len(guessControls)):
    #Copy the control from the trajectory
    guess.setControl(guessControls[cc], traj.getControlMat(guessControls[cc]))
    # #Check for lumbar controls
    # if 'lumbar_bend' in guessControls[cc] or 'lumbar_rot' in guessControls[cc]:
    #     #Set as zeros
    #     guess.setControl(guessControls[cc], np.linspace(0.0,0.0,traj.getNumTimes()))
    # else:
    #     #Copy the control from the trajectory
    #     guess.setControl(guessControls[cc], traj.getControlMat(guessControls[cc]))

#Save the edited guess to file
guess.write('initialGuess_2D_tracking.sto')

#Set the guess file in the solver
solver.setGuessFile('initialGuess_2D_tracking.sto')

# #Set the normalized tendon forces if not loading initial guess from file
# guess = solver.getGuess()
# numRows = guess.getNumTimes()
# stateNames = model.getStateVariableNames()
# for ii in range(0,model.getNumStateVariables()):
#     #Get current state name
#     currentStateName = stateNames.get(ii)
#     if 'normalized_tendon_force' in currentStateName:
#         guess.setState(currentStateName, np.linspace(0.2,0.2,numRows))

# %% Solve the tracking problem

#Solve!
sprintTrackingSolution = study.solve()

###### Iterations exceeded...
sprintTrackingSolution.unseal()

# study.printToXML('test.omoco')

# % Write the solution to a file
# sprintTrackingSolution.write('sprintTracking_solution_halfStride.sto')

# % Visualize the solution
# study.visualize(gaitTrackingSolution);


#### TODO: must unlock model coordinates to avoid NaN's in force table...

#Unlock a model to get GRFs

# # grfModel = osim.Model('model2D.osim')
# for kk in range(0,model.getCoordinateSet().getSize()):
#     model.getCoordinateSet().get(kk).set_locked(False)
# model.finalizeConnections()

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
osim.writeTableToFile(externalForcesTableFlat,'sprintTracking_2D_withContact_solution_GRF.sto')

##### still get some z axis forces (check UMocoD paramteters)

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
