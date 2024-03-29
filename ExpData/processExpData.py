# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 11:31:10 2020

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This code processes the experimental max sprinting trial contained in this
    folder for use in future optimal control simulations.
    
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
# import matplotlib.pyplot as plt
import xml.etree.ElementTree as et
import shutil
from distutils.dir_util import copy_tree

# %% Participant details

#Get participant mass from static trial

#Navigate to data directory
os.chdir('Data')

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
# #Set filtering for kinematics
# forceXML.setLowpassCutoffFrequencyForLoadKinematics(12)
#Write to file
forceXML.printToXML('refGRF.xml')

#Adapt for the 2D version 
forceXML.setDataFileName('refGRF_2D.mot')
forceXML.printToXML('refGRF_2D.xml')

# %% Model scaling (3D Model)

#Navigate to scaling directory
os.chdir('..\\Scaling')

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
osimHelper.addVirtualMarkersStatic('..\\Data\\static.trc',
                                   '..\\Data\\staticVirtualMarkers.trc')

#Place in scale tool
scaleTool.getMarkerPlacer().setMarkerFileName('..\\Data\\staticVirtualMarkers.trc')
scaleTool.getModelScaler().setMarkerFileName('..\\Data\\staticVirtualMarkers.trc')

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

#Loop through segments and place in dictionary
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
    #Option to scale size here...
    # currSphere.setRadius(currSphere.getRadius() * scaleFac)
    #Set new location
    newLoc = []
    newLoc.append(currSphere.getLocation().get(0) * scaleFactors[scaleBod][0])
    newLoc.append(currSphere.getLocation().get(1) * scaleFactors[scaleBod][1])
    newLoc.append(currSphere.getLocation().get(2) * scaleFactors[scaleBod][2])
    currSphere.setLocation(osim.Vec3(newLoc[0],newLoc[1],newLoc[2]))    

#Ensure all spheres lie on a flat plane with the foot in a neutral position (i.e. default pose)
#This will done by finding minimum base of sphere and ensuring each aligns with that
#on the Y-axis
#Loop through sphere list to identify minimum Y-value
minY = 999 #default unnecessarily high starting value
for ss in range(0,len(sphereList)):
    #Get the current sphere
    currSphere = osim.ContactSphere.safeDownCast(scaledModel.getContactGeometrySet().get(sphereList[ss]))
    #Check the y-value minimum location taking centre position and radius
    currY = currSphere.getLocation().get(1) - currSphere.getRadius()
    #Replace min-Y if appropriate
    if currY < minY:
        minY = currSphere.getLocation().get(1) - currSphere.getRadius()
#Re-loop through spheres and adjust Y-location if it is higher than the min
for ss in range(0,len(sphereList)): 
    #Get the current sphere
    currSphere = osim.ContactSphere.safeDownCast(scaledModel.getContactGeometrySet().get(sphereList[ss]))
    #Check min location and adjust if necessary
    #Need to do absolutes as Y-value is negative
    if abs((currSphere.getLocation().get(1) - currSphere.getRadius())) < abs(minY):
        #Drop the y-value on the sphere
        #Calculate new y-loc
        diffY = abs(minY) - abs((currSphere.getLocation().get(1) - currSphere.getRadius()))
        newY = currSphere.getLocation().get(1) - diffY
        currSphere.setLocation(osim.Vec3(
            currSphere.getLocation().get(0),
            newY,
            currSphere.getLocation().get(2)))
    
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

# %% Model scaling (2D Model)

#Set up the scale tool
scaleTool2D = osim.ScaleTool()

#Set participant mass
scaleTool2D.setSubjectMass(massKg)

#Set generic model file
genModelFileName2D = 'gait9dof18musc_Ong_et_al_Moco.osim'
scaleTool2D.getGenericModelMaker().setModelFileName(genModelFileName2D)

#Set the measurement set
measurementSetObject = osim.OpenSimObject.makeObjectFromFile('scaleMeasurementSet.xml')
measurementSet = osim.MeasurementSet.safeDownCast(measurementSetObject)
scaleTool2D.getModelScaler().setMeasurementSet(measurementSet)

#Set scale tasks
taskSet = osim.IKTaskSet('scaleTasks.xml')
for k in range(0,taskSet.getSize()-1):
    scaleTool2D.getMarkerPlacer().getIKTaskSet().adoptAndAppend(taskSet.get(k))

#Set marker file
scaleTool2D.getMarkerPlacer().setMarkerFileName('..\\Data\\staticVirtualMarkers.trc')
scaleTool2D.getModelScaler().setMarkerFileName('..\\Data\\staticVirtualMarkers.trc')

#Set options
scaleTool2D.getModelScaler().setPreserveMassDist(True)
scaleOrder = osim.ArrayStr(); scaleOrder.set(0,'measurements')
scaleTool2D.getModelScaler().setScalingOrder(scaleOrder)

#Set time ranges
timeRange = osim.ArrayDouble()
timeRange.set(0,0.5); timeRange.set(1,1.5)
scaleTool2D.getMarkerPlacer().setTimeRange(timeRange)
scaleTool2D.getModelScaler().setTimeRange(timeRange)

#Set output files
scaleTool2D.getModelScaler().setOutputModelFileName('scaledModel2D.osim')
scaleTool2D.getModelScaler().setOutputScaleFileName('scaleSet2D.xml')

#Set marker adjuster parameters
scaleTool2D.getMarkerPlacer().setOutputMotionFileName('static_motion_2D.mot')
scaleTool2D.getMarkerPlacer().setOutputModelFileName('scaledModelAdjusted2D.osim')

#Save and run scale tool
scaleTool2D.printToXML('scaleSetup2D.xml')
scaleTool2D.run()

#Load in scaled model
scaledModel2D = osim.Model('scaledModelAdjusted2D.osim')

#The locations and size of the contact sphere parameters go unchanged with
#standard model scaling, so these need to be edited to ensure they are in
#an appropriate place. This can be done based on the scale set parameters
#for the representative bodies.

#Load in the scale set, parsed from the XML tree
xmlTree2D = et.parse('scaleSet2D.xml')
xmlRoot2D = xmlTree2D.getroot()

#Create a dictionary with segments and scale factors to access for calculations
scaleFactors2D = {}

#Loop through segments and place in dictionary
for segment in range(0,len(xmlRoot2D.findall('./ScaleSet/objects/Scale/segment'))):
    #Get current segment name
    currSegment = xmlRoot2D.findall('./ScaleSet/objects/Scale/segment')[segment].text
    #Get current scale name and parse to 0 x 3 array
    currScale = xmlRoot2D.findall('./ScaleSet/objects/Scale/scales')[segment].text
    currScale = str.split(currScale)
    scaleFactors2D[currSegment] = [float(currScale[0]),float(currScale[1]),float(currScale[2])]
    
#Get the 3D scale factors for the relevant bodies and average to scale the
#sphere radii. Note that these scale factors for each foot will be the same
#for heel and toes as it looks like same scale factors are applied.
heelSphereScale_r_2D = sum(scaleFactors2D['calcn_r']) / 3
heelSphereScale_l_2D = sum(scaleFactors2D['calcn_l']) / 3
toesSphereScale_r_2D = sum(scaleFactors2D['toes_r']) / 3
toesSphereScale_l_2D = sum(scaleFactors2D['toes_l']) / 3

#Scale the radii to each of their respective factors
#While accessing the spheres, also adjust their position based on the scale
#factor for the respective axes
#Create a list of the sheres to loop through and edit
sphereList2D = ['heel_r','midfoot_r', 'toe_r',
                'heel_r','midfoot_l', 'toe_l']
#Loop through sphere list
for ss in range(len(sphereList2D)):
    #Get the current sphere
    currSphere = osim.ContactSphere.safeDownCast(scaledModel2D.getContactGeometrySet().get(sphereList2D[ss]))
    #Set the current scaling factor based on sphere name
    if '_r' in sphereList2D[ss]:
        if 'toe' in sphereList2D[ss] or 'midfoot' in sphereList2D[ss]:
            scaleFac = toesSphereScale_r
            scaleBod = 'toes_r'
        else:
            scaleFac = heelSphereScale_r
            scaleBod = 'calcn_r'
    elif '_l' in sphereList2D[ss]:
        if 'toe' in sphereList2D[ss] or 'midfoot' in sphereList2D[ss]:
            scaleFac = toesSphereScale_l
            scaleBod = 'toes_l'
        else:
            scaleFac = heelSphereScale_l
            scaleBod = 'calcn_l'
    #Rescale the radius
    #Option to scale size here
    # currSphere.setRadius(currSphere.getRadius() * scaleFac)
    #Set new location
    newLoc = []
    newLoc.append(currSphere.getLocation().get(0) * scaleFactors[scaleBod][0])
    newLoc.append(currSphere.getLocation().get(1) * scaleFactors[scaleBod][1])
    newLoc.append(currSphere.getLocation().get(2) * scaleFactors[scaleBod][2])
    currSphere.setLocation(osim.Vec3(newLoc[0],newLoc[1],newLoc[2]))    

#Ensure all spheres lie on a flat plane with the foot in a neutral position (i.e. default pose)
#This will done by finding minimum base of sphere and ensuring each aligns with that
#on the Y-axis
#Loop through sphere list to identify minimum Y-value
minY = 999 #default unnecessarily high starting value
for ss in range(0,len(sphereList2D)):
    #Get the current sphere
    currSphere = osim.ContactSphere.safeDownCast(scaledModel2D.getContactGeometrySet().get(sphereList2D[ss]))
    #Check the y-value minimum location taking centre position and radius
    currY = currSphere.getLocation().get(1) - currSphere.getRadius()
    #Replace min-Y if appropriate
    if currY < minY:
        minY = currSphere.getLocation().get(1) - currSphere.getRadius()
#Re-loop through spheres and adjust Y-location if it is higher than the min
for ss in range(0,len(sphereList2D)): 
    #Get the current sphere
    currSphere = osim.ContactSphere.safeDownCast(scaledModel2D.getContactGeometrySet().get(sphereList2D[ss]))
    #Check min location and adjust if necessary
    #Need to do absolutes as Y-value is negative
    if abs((currSphere.getLocation().get(1) - currSphere.getRadius())) < abs(minY):
        #Drop the y-value on the sphere
        #Calculate new y-loc
        diffY = abs(minY) - abs((currSphere.getLocation().get(1) - currSphere.getRadius()))
        newY = currSphere.getLocation().get(1) - diffY
        currSphere.setLocation(osim.Vec3(
            currSphere.getLocation().get(0),
            newY,
            currSphere.getLocation().get(2)))

#Reset model name
scaledModel2D.setName('scaledModel2D')

#Remove the marker set as it's no use in this 2D context
scaledModel2D.updMarkerSet().clearAndDestroy()

#Update scaled model
scaledModel2D.printToXML('scaledModelAdjusted2D.osim')

#Scale muscle strength based on linear function presented in Handsfield
#et al. (2014). This uses some convenience functions that are packaged
#with the Rajagopal et al. (2016) gait model. Note that the height of
#the generic model is 1.700; the current participant height is 1.79
osimHelper.scaleOptimalForceSubjectSpecific(genericModelFileName = 'gait9dof18musc_Ong_et_al_Moco.osim',
                                            scaledModelFileName = 'scaledModelAdjusted2D.osim',
                                            genericHeight = 1.700, scaledHeight = 1.79,
                                            outputModelFileName = 'scaledModelMuscle2D.osim')
#Load in new scaled model
scaledModelMuscle2D = osim.Model('scaledModelMuscle2D.osim')

# %% Adapt 3D model to create a complex 2D version
    
### TODO: does this get used anymore???

#Set a new copy of the model
model2D = osim.Model('scaledModelMuscle.osim')

#Get hip joints from scaled model to work from
origHip_r = scaledModelMuscle.getJointSet().get('hip_r')
origHip_l = scaledModelMuscle.getJointSet().get('hip_l')
    
#Create pin joint for right hip
#Set parent and child locations
locationInParent = origHip_r.get_frames(0).get_translation()
orientationInParent = origHip_r.get_frames(0).get_orientation()
locationInChild = origHip_r.get_frames(1).get_translation()
orientationInChild = origHip_r.get_frames(1).get_orientation()
#Get bodies for joint
pelvisBody = model2D.getBodySet().get('pelvis')
femurBody = model2D.getBodySet().get('femur_r')
#Create joint
hipJoint_r = osim.PinJoint('hip_r', pelvisBody, locationInParent,
                            orientationInParent, femurBody, locationInChild,
                            orientationInChild)
#Update additional joint parameters
hipJoint_r.getCoordinate().setName('hip_flexion_r') #set coordinate name
hipJoint_r.getCoordinate().setRangeMin(origHip_r.get_coordinates(0).getRangeMin()) #get min from old joint
hipJoint_r.getCoordinate().setRangeMax(origHip_r.get_coordinates(0).getRangeMax()) #get max from old joint
hipJoint_r.getCoordinate().setDefaultValue(0) #set default value to zero
hipJoint_r.getCoordinate().setDefaultSpeedValue(0) #set default speed to zero
hipJoint_r.getCoordinate().set_clamped(True) #set joint to clamped
hipJoint_r.getCoordinate().set_locked(False) #set locked to false
hipJoint_r.getCoordinate().set_prescribed(False) #set prescribed to false

#Create pin joint for left hip
#Set parent and child locations
locationInParent = origHip_l.get_frames(0).get_translation()
orientationInParent = origHip_l.get_frames(0).get_orientation()
locationInChild = origHip_l.get_frames(1).get_translation()
orientationInChild = origHip_l.get_frames(1).get_orientation()
#Get bodies for joint
pelvisBody = model2D.getBodySet().get('pelvis')
femurBody = model2D.getBodySet().get('femur_l')
#Create joint
hipJoint_l = osim.PinJoint('hip_l', pelvisBody, locationInParent,
                            orientationInParent, femurBody, locationInChild,
                            orientationInChild)
#Update additional joint parameters
hipJoint_l.getCoordinate().setName('hip_flexion_l') #set coordinate name
hipJoint_l.getCoordinate().setRangeMin(origHip_l.get_coordinates(0).getRangeMin()) #get min from old joint
hipJoint_l.getCoordinate().setRangeMax(origHip_l.get_coordinates(0).getRangeMax()) #get max from old joint
hipJoint_l.getCoordinate().setDefaultValue(0) #set default value to zero
hipJoint_l.getCoordinate().setDefaultSpeedValue(0) #set default speed to zero
hipJoint_l.getCoordinate().set_clamped(True) #set joint to clamped
hipJoint_l.getCoordinate().set_locked(False) #set locked to false
hipJoint_l.getCoordinate().set_prescribed(False) #set prescribed to false

#Remove existing hip joints from model
model2D.getJointSet().remove(model2D.getJointSet().get('hip_r'))
model2D.getJointSet().remove(model2D.getJointSet().get('hip_l'))

#Add new hip joints
model2D.addJoint(hipJoint_r)
model2D.addJoint(hipJoint_l)

#Get knee joints from scaled model to work from
origKnee_r = scaledModelMuscle.getJointSet().get('walker_knee_r')
origKnee_l = scaledModelMuscle.getJointSet().get('walker_knee_l')
    
#Create pin joint for right knee
#Set parent and child locations
locationInParent = origKnee_r.get_frames(0).get_translation()
orientationInParent = osim.Vec3(0,math.radians(180),0)
locationInChild = origKnee_r.get_frames(1).get_translation()
orientationInChild = osim.Vec3(0,math.radians(180),0)
#Get bodies for joint
femurBody = model2D.getBodySet().get('femur_r')
tibiaBody = model2D.getBodySet().get('tibia_r')
#Create joint
kneeJoint_r = osim.PinJoint('walker_knee_r', femurBody, locationInParent,
                            orientationInParent, tibiaBody, locationInChild,
                            orientationInChild)
#Update additional joint parameters
kneeJoint_r.getCoordinate().setName('knee_angle_r') #set coordinate name
kneeJoint_r.getCoordinate().setRangeMin(origKnee_r.get_coordinates(0).getRangeMin()) #get min from old joint
kneeJoint_r.getCoordinate().setRangeMax(origKnee_r.get_coordinates(0).getRangeMax()) #get max from old joint
kneeJoint_r.getCoordinate().setDefaultValue(0) #set default value to zero
kneeJoint_r.getCoordinate().setDefaultSpeedValue(0) #set default speed to zero
kneeJoint_r.getCoordinate().set_clamped(True) #set joint to clamped
kneeJoint_r.getCoordinate().set_locked(False) #set locked to false
kneeJoint_r.getCoordinate().set_prescribed(False) #set prescribed to false

#Create pin joint for left knee
#Set parent and child locations
locationInParent = origKnee_l.get_frames(0).get_translation()
orientationInParent = osim.Vec3(0,math.radians(180),0)
locationInChild = origKnee_l.get_frames(1).get_translation()
orientationInChild = osim.Vec3(0,math.radians(180),0)
#Get bodies for joint
femurBody = model2D.getBodySet().get('femur_l')
tibiaBody = model2D.getBodySet().get('tibia_l')
#Create joint
kneeJoint_l = osim.PinJoint('walker_knee_l', femurBody, locationInParent,
                            orientationInParent, tibiaBody, locationInChild,
                            orientationInChild)
#Update additional joint parameters
kneeJoint_l.getCoordinate().setName('knee_angle_l') #set coordinate name
kneeJoint_l.getCoordinate().setRangeMin(origKnee_l.get_coordinates(0).getRangeMin()) #get min from old joint
kneeJoint_l.getCoordinate().setRangeMax(origKnee_l.get_coordinates(0).getRangeMax()) #get max from old joint
kneeJoint_l.getCoordinate().setDefaultValue(0) #set default value to zero
kneeJoint_l.getCoordinate().setDefaultSpeedValue(0) #set default speed to zero
kneeJoint_l.getCoordinate().set_clamped(True) #set joint to clamped
kneeJoint_l.getCoordinate().set_locked(False) #set locked to false
kneeJoint_l.getCoordinate().set_prescribed(False) #set prescribed to false

#Remove existing hip joints from model
model2D.getJointSet().remove(model2D.getJointSet().get('walker_knee_r'))
model2D.getJointSet().remove(model2D.getJointSet().get('walker_knee_l'))

#Add new hip joints
model2D.addJoint(kneeJoint_r)
model2D.addJoint(kneeJoint_l)

#Get back joint from scaled model to work from
origLumbar = scaledModelMuscle.getJointSet().get('back')
    
#Create pin joint for back
#Set parent and child locations
locationInParent = origLumbar.get_frames(0).get_translation()
orientationInParent = origLumbar.get_frames(0).get_orientation()
locationInChild = origLumbar.get_frames(1).get_translation()
orientationInChild = origLumbar.get_frames(1).get_orientation()
#Get bodies for joint
pelvisBody = model2D.getBodySet().get('pelvis')
torsoBody = model2D.getBodySet().get('torso')
#Create joint
lumbarJoint = osim.PinJoint('back', pelvisBody, locationInParent,
                            orientationInParent, torsoBody, locationInChild,
                            orientationInChild)
#Update additional joint parameters
lumbarJoint.getCoordinate().setName('lumbar_extension') #set coordinate name
lumbarJoint.getCoordinate().setRangeMin(origLumbar.get_coordinates(0).getRangeMin()) #get min from old joint
lumbarJoint.getCoordinate().setRangeMax(origLumbar.get_coordinates(0).getRangeMax()) #get max from old joint
lumbarJoint.getCoordinate().setDefaultValue(0) #set default value to zero
lumbarJoint.getCoordinate().setDefaultSpeedValue(0) #set default speed to zero
lumbarJoint.getCoordinate().set_clamped(True) #set joint to clamped
lumbarJoint.getCoordinate().set_locked(False) #set locked to false
lumbarJoint.getCoordinate().set_prescribed(False) #set prescribed to false

#Remove existing back joint from model
model2D.getJointSet().remove(model2D.getJointSet().get('back'))

#Add new back joint
model2D.addJoint(lumbarJoint)

#Remove existing lumbar coordinate actuators
#Get force names
forceNames = list()
for ff in range(model2D.getForceSet().getSize()):
    forceNames.append(model2D.getForceSet().get(ff).getName())
#Find indices of actuators to remove
indRemove = list()
indRemove.append(forceNames.index('tau_lumbar_bend'))
indRemove.append(forceNames.index('tau_lumbar_rot'))
#Sort descending to avoid changing index values
indRemove.sort(reverse = True)
#Remove forces
for ff in range(0,len(indRemove)):
    model2D.getForceSet().remove(indRemove[ff])
    
#Get original pelvis from scaled model to work from
origPelvis = scaledModelMuscle.getJointSet().get('ground_pelvis')

#Create planar joint for ground-pelvis
pelvisJoint = osim.PlanarJoint('ground_pelvis',
                                model2D.getGround(), osim.Vec3(0,0,0), osim.Vec3(0,0,0),
                                model2D.getBodySet().get('pelvis'), osim.Vec3(0,0,0), osim.Vec3(0,0,0))
#Update additional joint parameters
#Pelvis tilt
pelvisJoint.getCoordinate(0).setName('pelvis_tilt') #set coordinate name
pelvisJoint.getCoordinate(0).setRangeMin(origPelvis.get_coordinates(0).getRangeMin()) #get min from old joint
pelvisJoint.getCoordinate(0).setRangeMax(origPelvis.get_coordinates(0).getRangeMax()) #get max from old joint
pelvisJoint.getCoordinate(0).setDefaultValue(0) #set default value to zero
pelvisJoint.getCoordinate(0).setDefaultSpeedValue(0) #set default speed to zero
pelvisJoint.getCoordinate(0).set_clamped(True) #set joint to clamped
pelvisJoint.getCoordinate(0).set_locked(False) #set locked to false
pelvisJoint.getCoordinate(0).set_prescribed(False) #set prescribed to false
#Pelvis tx
pelvisJoint.getCoordinate(1).setName('pelvis_tx') #set coordinate name
pelvisJoint.getCoordinate(1).setRangeMin(origPelvis.get_coordinates(3).getRangeMin()) #get min from old joint
pelvisJoint.getCoordinate(1).setRangeMax(origPelvis.get_coordinates(3).getRangeMax()) #get max from old joint
pelvisJoint.getCoordinate(1).setDefaultValue(0) #set default value to zero
pelvisJoint.getCoordinate(1).setDefaultSpeedValue(0) #set default speed to zero
pelvisJoint.getCoordinate(1).set_clamped(True) #set joint to clamped
pelvisJoint.getCoordinate(1).set_locked(False) #set locked to false
pelvisJoint.getCoordinate(1).set_prescribed(False) #set prescribed to false
#Pelvis ty
pelvisJoint.getCoordinate(2).setName('pelvis_ty') #set coordinate name
pelvisJoint.getCoordinate(2).setRangeMin(origPelvis.get_coordinates(4).getRangeMin()) #get min from old joint
pelvisJoint.getCoordinate(2).setRangeMax(origPelvis.get_coordinates(4).getRangeMax()) #get max from old joint
pelvisJoint.getCoordinate(2).setDefaultValue(0.95) #set default value to zero
pelvisJoint.getCoordinate(2).setDefaultSpeedValue(0) #set default speed to zero
pelvisJoint.getCoordinate(2).set_clamped(True) #set joint to clamped
pelvisJoint.getCoordinate(2).set_locked(False) #set locked to false
pelvisJoint.getCoordinate(2).set_prescribed(False) #set prescribed to false

#Remove existing pelvis joint from model
model2D.getJointSet().remove(model2D.getJointSet().get('ground_pelvis'))

#Add new pelvis joint
model2D.addJoint(pelvisJoint)

#Update ankle joints orientation to planar
model2D.getJointSet().get('ankle_r').get_frames(0).set_orientation(osim.Vec3(0,0,0))
model2D.getJointSet().get('ankle_l').get_frames(0).set_orientation(osim.Vec3(0,0,0))
model2D.getJointSet().get('ankle_r').get_frames(1).set_orientation(osim.Vec3(0,0,0))
model2D.getJointSet().get('ankle_l').get_frames(1).set_orientation(osim.Vec3(0,0,0))

#Update mtp joints orientation to planar
model2D.getJointSet().get('mtp_r').get_frames(0).set_orientation(osim.Vec3(0,math.radians(180),0))
model2D.getJointSet().get('mtp_l').get_frames(0).set_orientation(osim.Vec3(0,math.radians(180),0))
model2D.getJointSet().get('mtp_r').get_frames(1).set_orientation(osim.Vec3(0,math.radians(180),0))
model2D.getJointSet().get('mtp_l').get_frames(1).set_orientation(osim.Vec3(0,math.radians(180),0))

#Update subtalar joints orientation to planar
model2D.getJointSet().get('subtalar_r').get_frames(0).set_orientation(osim.Vec3(0,math.radians(180),0))
model2D.getJointSet().get('subtalar_l').get_frames(0).set_orientation(osim.Vec3(0,math.radians(180),0))
model2D.getJointSet().get('subtalar_r').get_frames(1).set_orientation(osim.Vec3(0,math.radians(180),0))
model2D.getJointSet().get('subtalar_l').get_frames(1).set_orientation(osim.Vec3(0,math.radians(180),0))

#Finalize model connections
model2D.finalizeConnections()

#Dump into a model processor to weld the subtalar and mtp joints
#Create the processor
editProcessor = osim.ModelProcessor(model2D)
#Create the string array of joints to weld
weldJoints = osim.StdVectorString()
weldJoints.append('subtalar_r')
weldJoints.append('subtalar_l')
weldJoints.append('mtp_r')
weldJoints.append('mtp_l')
#Append model operator for welding
editProcessor.append(osim.ModOpReplaceJointsWithWelds(weldJoints))
#Append model operator to convert muscles
editProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
#Process model output
model2D_final = editProcessor.process()

#Remove the marker set as it's no use in this 2D context
model2D_final.updMarkerSet().clearAndDestroy()

#Print 2D model output
model2D_final.printToXML('scaledModelMuscle2D_complex.osim')

# %% Set simulation parameters

#Identify time parameters for simulation
#This will be based off a half gait cycle of the right limb, that being from
#right foot heel strike to left foot heel strike
[startTime,endTime] = osimHelper.getHalfGaitCycle('..\\Data\\refGRF.mot')
#Or a full gait cycle
# [startTime,endTime] = osimHelper.getFullGaitCycle('refGRF.mot')

#Add the virtual torso, pelvis and hip joint markers to the .trc file
osimHelper.addVirtualMarkersDynamic(staticTRC = '..\\Data\\staticVirtualMarkers.trc',
                                    dynamicTRC = '..\\Data\\maxSprint.trc',
                                    outputTRC = '..\\Data\\maxSprintVirtualMarkers.trc')

# %% Inverse kinematics

#Navigate to IK directory
os.chdir('..\\IK')

#Initialise IK tool
ikTool = osim.InverseKinematicsTool()

#Set model
ikTool.setModel(scaledModelMuscle)

#Set task set
ikTool.set_IKTaskSet(osim.IKTaskSet('ikTasks.xml'))

#Set marker file
ikTool.set_marker_file('..\\Data\\maxSprintVirtualMarkers.trc')

#Set times
ikTool.setStartTime(startTime)
ikTool.setEndTime(endTime)

#Set output filename
ikTool.set_output_motion_file('ikResults.mot')

#Run IK
ikTool.run()

# %% Residual reduction algorithm

#Navigate to RRA directory
os.chdir('..\\RRA')

#Initialise list for RRA tool
rraTool = []
runRRA = []

#Loop through and perform rra three times
#Each time we'll adjust the mass specified by the RRA tool, and then re-use
#this adjusted model in the next iteration. On the first iteration, we'll
#use the IK motion data - but on subsequent iterations we'll use the
#adjusted rra kinematics. Similarly, on the first iteration we'll set the
#model to our scaled model - but then use the adjusted version on
#subsequent iterations.
for rr in range(1): ##### normally would do multiple iterations, but RRA bugging out with loop?
    
    #Make directory for current iteration
    #Check if directory exists (errors if making when already there)
    if not os.path.isdir('rra'+str(rr+1)):
        os.mkdir('rra'+str(rr+1))
    os.chdir('rra'+str(rr+1))
    
    #Create tool in list
    rraTool.append(osim.RRATool())

    #Set tool to replace model force set
    rraTool[rr].setReplaceForceSet(True)
    
    #Set actuators file
    forceSetFiles = osim.ArrayStr()
    forceSetFiles.set(0,'..\\rraActuators.xml')
    rraTool[rr].setForceSetFiles(forceSetFiles)
    
    #Set tracking tasks file
    rraTool[rr].setTaskSetFileName('..\\rraTasks.xml')
    
    #Set a low pass filter frequency on the kinematics data
    rraTool[rr].setLowpassCutoffFrequency(12)
    
    #Set to adjust the COM to reduce residuals
    rraTool[rr].setAdjustCOMToReduceResiduals(True)
    
    #Set the torso body COM to be the one that gets adjusted
    rraTool[rr].setAdjustedCOMBody('torso')
    
    #Set external loads file
    rraTool[rr].setExternalLoadsFileName('..\\..\\Data\\refGRF.xml')

    #Set tool name based on iteration
    rraTool[rr].setName('rra'+str(rr+1))
    
    #Set desired kinematics file
    if rr == 0:

        #Use IK data
        rraTool[rr].setDesiredKinematicsFileName('..\\..\\IK\\ikResults.mot')

        #Set initial and final time
        #You need to use the IK first and last times here as I don't think the tool
        #likes if the IK data doesn't have the times you put in
        rraTool[rr].setStartTime(osim.Storage('..\\..\\IK\\ikResults.mot').getFirstTime())
        rraTool[rr].setFinalTime(osim.Storage('..\\..\\IK\\ikResults.mot').getLastTime())

        #Set to original scaled model
        rraTool[rr].setModelFilename('..\\..\\Scaling\\scaledModelMuscle.osim')
        
    else:
        
        #Use previous rra data
        rraTool[rr].setDesiredKinematicsFileName('..\\rra'+str(rr)+'\\rra1_Kinematics_q.sto')
        
        #Set initial and final time
        rraTool[rr].setStartTime(osim.Storage('..\\rra'+str(rr)+'\\rra1_Kinematics_q.sto').getFirstTime())
        rraTool[rr].setFinalTime(osim.Storage('..\\rra'+str(rr)+'\\rra1_Kinematics_q.sto').getLastTime())
        
        #Set to last iteration rra model
        rraTool[rr].setModelFilename('..\\rra'+str(rr)+'\\rraAdjustedModel_'+str(rr)+'.osim')
        
    #Set output model file
    rraTool[rr].setOutputModelFileName('rraAdjustedModel_'+str(rr+1)+'.osim')

    #Print out the tool
    rraTool[rr].printToXML('setupRRA'+str(rr+1)+'.xml')

    #Run RRA
    #Reloading the setup tool seems to help with running the tool smoothly
    runRRA.append(osim.RRATool('setupRRA'+str(rr+1)+'.xml'))
    runRRA[rr].run()
    
    #Get the suggested mass change
    #Read back through the file lines until getting to the total mass change
    #Open text log and get lines
    txtLog = open('..\\..\\opensim.log')
    txtLine = txtLog.readlines()
    #Loop through lines in reverse to identify last occurring mass change
    for line in reversed(txtLine):
        if line.startswith('*  Total mass change:'):
            #Get the mass change recommendation
            massChange = float(line.split(sep = ': ')[1])
            break
    #Set suggested new mass
    if rr == 0:
        currMass = osimHelper.getMassOfModel('..\\..\\Scaling\\scaledModelMuscle.osim') 
        
    #Set new mass and ratio
    newMass = currMass + massChange
    massScaleFac = newMass / currMass
        
    #Apply mass change to rra model
    rraModel = osim.Model('rraAdjustedModel_'+str(rr+1)+'.osim')
    allBodies = rraModel.getBodySet()
    for ii in range(allBodies.getSize()):
        currBodyMass = allBodies.get(ii).getMass()
        newBodyMass = currBodyMass * massScaleFac
        allBodies.get(ii).setMass(newBodyMass)
    
    #Print over model
    rraModel.printToXML('rraAdjustedModel_'+str(rr+1)+'.osim')
    
    #Return to RRA directory
    os.chdir('..')
    
#Get the final number of rra iterations
finalRRA = rr+1

#Get the filepath for the final RRA kinematics
rraKinematics = os.getcwd()+'\\rra'+str(finalRRA)+'\\rra'+str(finalRRA)+'_Kinematics_q.sto'

#Get the filepath to the final RRA model
rraModel = os.getcwd()+'\\rra'+str(finalRRA)+'\\rraAdjustedModel_'+str(finalRRA)+'.osim'

#Convert kinematic results to a states file and store in simulation directory
osimHelper.kinematicsToStates(kinematicsFileName = rraKinematics,
                              osimModelFileName = rraModel,
                              outputFileName = '..\\..\\Simulations\\trackingSims\\refQ.sto',
                              inDegrees = True, outDegrees = False)

#Convert states to 2D model format
osimHelper.statesTo2D(statesFileName = '..\\..\\Simulations\\trackingSims\\refQ.sto',
                      outputFileName = '..\\..\\Simulations\\trackingSims\\refQ_2D.sto',
                      renameWalkerKnee = True)

#Copy final rra model to simulation directory
shutil.copy(rraModel,'..\\..\\Simulations\\trackingSims\\gaitModel3D.osim')

#Adjust mass of 2D model coinciding with RRA changes
#Get model masses and ratios
rraMass = osimHelper.getMassOfModel(rraModel)
mass2D = osimHelper.getMassOfModel(scaledModelMuscle2D)
massRatio2D = rraMass / mass2D
#Adjust body masses
allBodies = scaledModelMuscle2D.getBodySet()
for ii in range(allBodies.getSize()):
    currBodyMass = allBodies.get(ii).getMass()
    newBodyMass = currBodyMass * massRatio2D
    allBodies.get(ii).setMass(newBodyMass)
#Re-print adjusted mass model to simulation directory
scaledModelMuscle2D.printToXML('..\\..\\Simulations\\trackingSims\\gaitModel2D.osim')

#Repeat for complex 2D model
mass2D = osimHelper.getMassOfModel(model2D_final)
massRatio2D = rraMass / mass2D
#Adjust body masses
allBodies = model2D_final.getBodySet()
for ii in range(allBodies.getSize()):
    currBodyMass = allBodies.get(ii).getMass()
    newBodyMass = currBodyMass * massRatio2D
    allBodies.get(ii).setMass(newBodyMass)
#Re-print adjusted mass model to simulation directory
model2D_final.printToXML('..\\..\\Simulations\\trackingSims\\gaitModel2D_complex.osim')

# %% Do some final data copies for subsequent simulations

#Copy geometry folder across to simulations directory
os.chdir('..\\..\\Simulations\\trackingSims')
if not os.path.isdir('Geometry'):
    os.mkdir('Geometry')
copy_tree('..\\..\\ExpData\\Scaling\\Geometry',
          os.getcwd()+'\\Geometry')

#Copy GRF data files to simulations directory
os.chdir('..\\..\\ExpData\\Data')
shutil.copy('refGRF_2D.mot', '..\\..\\Simulations\\trackingSims')
shutil.copy('refGRF_2D.xml', '..\\..\\Simulations\\trackingSims')
shutil.copy('refGRF.mot', '..\\..\\Simulations\\trackingSims')
shutil.copy('refGRF.xml', '..\\..\\Simulations\\trackingSims')

# %% Finish up
print('----- Processing of experimental data complete -----')

# %% OLD CODE THAT CAN EVENTUALLY GO BELOW...

# %% Extract scaled model parameters for building external functions

# ##### TODO: check appropriateness of this against the TrackSim_1/2 variants for tracking...

# # This extracts model values for a 2D model. It's important to note that these
# # extracted parameters are related to the models used in the Falisse et al.
# # algorithmic differentiation pipelines external functions.
# #
# # These external functions contain construction of a model that is effectively
# # in OpenSim 3.3 format --- which differs to the 4.x versions used for extracting
# # experimental data here. Hopefully it works anyway...
# #
# # The end of this section prints the outputs to a text file that theoretically
# # can be copied into an external function template

# #Set a dictionary to store the body values
# bodyVals = {'bodyName': [], 'bodyMass': [], 'massCentre': [], 'inertia': []}

# #Set a dictionary to store the joint values
# jointVals = {'jointName': [], 'locInParent': []}

# #Set a list for the bodies to extract
# extractBodies = ['pelvis', 'femur_l', 'femur_r', 'tibia_l', 'tibia_r',
#                  'talus_l', 'talus_r', 'calcn_l', 'calcn_r',
#                  'toes_l', 'toes_r', 'torso']

# #Set a list for the joints to extract
# extractJoints = ['ground_pelvis', 'hip_l', 'hip_r',
#                  'walker_knee_l', 'walker_knee_r',
#                  'ankle_l', 'ankle_r', 'subtalar_l', 'subtalar_r',
#                  'mtp_l', 'mtp_r', 'back']
# nameJoints = ['ground_pelvis', 'hip_l', 'hip_r',
#               'knee_l', 'knee_r',
#               'ankle_l', 'ankle_r', 'subtalar_l', 'subtalar_r',
#               'mtp_l', 'mtp_r', 'back']

# #Get the body parameters
# #Loop through bodies
# for bb in range(len(extractBodies)):
    
#     #Set body name in dictionary
#     bodyVals['bodyName'].append(extractBodies[bb])
    
#     #Get body mass
#     bodyVals['bodyMass'].append(model2D_final.getBodySet().get(extractBodies[bb]).getMass())
    
#     #Get mass centre
#     bodyVals['massCentre'].append(np.array([model2D_final.getBodySet().get(extractBodies[bb]).getMassCenter().get(0),
#                                             model2D_final.getBodySet().get(extractBodies[bb]).getMassCenter().get(1),
#                                             model2D_final.getBodySet().get(extractBodies[bb]).getMassCenter().get(2)]))
    
#     #Get inertia
#     bodyVals['inertia'].append(np.array([model2D_final.getBodySet().get(extractBodies[bb]).getInertia().getMoments().get(0),
#                                          model2D_final.getBodySet().get(extractBodies[bb]).getInertia().getMoments().get(1),
#                                          model2D_final.getBodySet().get(extractBodies[bb]).getInertia().getMoments().get(2),
#                                          model2D_final.getBodySet().get(extractBodies[bb]).getInertia().getProducts().get(0),
#                                          model2D_final.getBodySet().get(extractBodies[bb]).getInertia().getProducts().get(1),
#                                          model2D_final.getBodySet().get(extractBodies[bb]).getInertia().getProducts().get(2)]))

# #Get the joint parameters
# #Loop through bodies
# for jj in range(len(extractJoints)):
    
#     #Set joint name
#     jointVals['jointName'].append(nameJoints[jj])
    
#     #Get location in parent
#     jointVals['locInParent'].append(np.array([model2D_final.getJointSet().get(extractJoints[jj]).get_frames(0).get_translation().get(0),
#                                               model2D_final.getJointSet().get(extractJoints[jj]).get_frames(0).get_translation().get(1),
#                                               model2D_final.getJointSet().get(extractJoints[jj]).get_frames(0).get_translation().get(2)]))

# #Write body set specifications to text file
# #Open file to write body set to
# bodySetExtFile = open('bodySet_externalFunction_2D.txt', 'w')
# #Loop through bodies
# for bb in range(len(extractBodies)):
    
#     #Create current text string for line
#     txtString = extractBodies[bb]+' = new OpenSim::Body("'+extractBodies[bb]+'" , '+\
#         str(bodyVals['bodyMass'][bb])+', Vec3('+str(bodyVals['massCentre'][bb][0])+\
#             ', '+str(bodyVals['massCentre'][bb][1])+', '+str(bodyVals['massCentre'][bb][2])+\
#                 '), Inertia('+str(bodyVals['inertia'][bb][0])+', '+str(bodyVals['inertia'][bb][1])+\
#                     ', '+str(bodyVals['inertia'][bb][2])+', '+str(bodyVals['inertia'][bb][3])+\
#                         ', '+str(bodyVals['inertia'][bb][4])+', '+str(bodyVals['inertia'][bb][5])+\
#                             '));'
        
#     #Write the line to file
#     bodySetExtFile.write('%s\n' % txtString)
    
# #Close file
# bodySetExtFile.close()

# #Write joint set specifications to text file
# jj = 0
# #Open file to write joint set to
# jointSetExtFile = open('jointSet_externalFunction_2D.txt', 'w')

# #These are relatively manual don't work well with a loop
# #Text string for ground to pelvis
# txtString = 'ground_pelvis = new PlanarJoint("ground_pelvis'+\
#     '", model->getGround(), Vec3(0), Vec3(0), *pelvis'+\
#         ', Vec3(0), Vec3(0));'
# jointSetExtFile.write('%s\n' % txtString)        
# jj = jj + 1 #increment on joint value index

# #Text string for custom hip joints
# txtString = 'hip_l = new CustomJoint("hip_l", *pelvis, Vec3('+\
#     str(jointVals['locInParent'][jj][0])+', '+str(jointVals['locInParent'][jj][1])+\
#         ', '+str(jointVals['locInParent'][jj][2])+'), Vec3(0), *femur_l,'+\
#             ' Vec3(0), Vec3(0), st_hip_l);'
# jointSetExtFile.write('%s\n' % txtString)            
# jj = jj + 1 #increment on joint value index
# txtString = 'hip_r = new CustomJoint("hip_r", *pelvis, Vec3('+\
#     str(jointVals['locInParent'][jj][0])+', '+str(jointVals['locInParent'][jj][1])+\
#         ', '+str(jointVals['locInParent'][jj][2])+'), Vec3(0), *femur_r,'+\
#             ' Vec3(0), Vec3(0), st_hip_r);'
# jointSetExtFile.write('%s\n' % txtString)            
# jj = jj + 1 #increment on joint value index
                
# #Text string for custom knee joints
# txtString = 'knee_l = new CustomJoint("knee_l", *femur_l, Vec3('+\
#     str(jointVals['locInParent'][jj][0])+', '+str(jointVals['locInParent'][jj][1])+\
#         ', '+str(jointVals['locInParent'][jj][2])+'), Vec3(0), *tibia_l,'+\
#             ' Vec3(0), Vec3(0), st_knee_l);'
# jointSetExtFile.write('%s\n' % txtString)            
# jj = jj + 1 #increment on joint value index
# txtString = 'knee_r = new CustomJoint("knee_l", *femur_r, Vec3('+\
#     str(jointVals['locInParent'][jj][0])+', '+str(jointVals['locInParent'][jj][1])+\
#         ', '+str(jointVals['locInParent'][jj][2])+'), Vec3(0), *tibia_r,'+\
#             ' Vec3(0), Vec3(0), st_knee_r);'
# jointSetExtFile.write('%s\n' % txtString)            
# jj = jj + 1 #increment on joint value index 

# #Text string for pin ankle joints            
# txtString = 'ankle_l = new PinJoint("ankle_l", *tibia_l, Vec3('+\
#     str(jointVals['locInParent'][jj][0])+', '+str(jointVals['locInParent'][jj][1])+\
#         ', '+str(jointVals['locInParent'][jj][2])+'), Vec3(0), *talus_l,'+\
#             ' Vec3(0), Vec3(0));'
# jointSetExtFile.write('%s\n' % txtString)            
# jj = jj + 1 #increment on joint value index
# txtString = 'ankle_r = new PinJoint("ankle_r", *tibia_r, Vec3('+\
#     str(jointVals['locInParent'][jj][0])+', '+str(jointVals['locInParent'][jj][1])+\
#         ', '+str(jointVals['locInParent'][jj][2])+'), Vec3(0), *talus_r,'+\
#             ' Vec3(0), Vec3(0));'
# jointSetExtFile.write('%s\n' % txtString)            
# jj = jj + 1 #increment on joint value index

# #Text string for weld subtalar joints
# txtString = 'subtalar_l = new WeldJoint("subtalar_l", *talus_l, Vec3('+\
#     str(jointVals['locInParent'][jj][0])+', '+str(jointVals['locInParent'][jj][1])+\
#         ', '+str(jointVals['locInParent'][jj][2])+'), Vec3(0), *calcn_l,'+\
#             ' Vec3(0), Vec3(0));'
# jointSetExtFile.write('%s\n' % txtString)            
# jj = jj + 1 #increment on joint value index
# txtString = 'subtalar_r = new WeldJoint("subtalar_r", *talus_r, Vec3('+\
#     str(jointVals['locInParent'][jj][0])+', '+str(jointVals['locInParent'][jj][1])+\
#         ', '+str(jointVals['locInParent'][jj][2])+'), Vec3(0), *calcn_r,'+\
#             ' Vec3(0), Vec3(0));'
# jointSetExtFile.write('%s\n' % txtString)            
# jj = jj + 1 #increment on joint value index

# #Text string for weld mtp joints
# txtString = 'mtp_l = new WeldJoint("mtp_l", *calcn_l, Vec3('+\
#     str(jointVals['locInParent'][jj][0])+', '+str(jointVals['locInParent'][jj][1])+\
#         ', '+str(jointVals['locInParent'][jj][2])+'), Vec3(0), *toes_l,'+\
#             ' Vec3(0), Vec3(0));'
# jointSetExtFile.write('%s\n' % txtString)            
# jj = jj + 1 #increment on joint value index
# txtString = 'mtp_r = new WeldJoint("mtp_r", *calcn_r, Vec3('+\
#     str(jointVals['locInParent'][jj][0])+', '+str(jointVals['locInParent'][jj][1])+\
#         ', '+str(jointVals['locInParent'][jj][2])+'), Vec3(0), *toes_r,'+\
#             ' Vec3(0), Vec3(0));'
# jointSetExtFile.write('%s\n' % txtString)            
# jj = jj + 1 #increment on joint value index

# #text string for back pin joint
# txtString = 'back = new PinJoint("back", *pelvis, Vec3('+\
#     str(jointVals['locInParent'][jj][0])+', '+str(jointVals['locInParent'][jj][1])+\
#         ', '+str(jointVals['locInParent'][jj][2])+'), Vec3(0), *torso,'+\
#             ' Vec3(0), Vec3(0));'
# jointSetExtFile.write('%s\n' % txtString)            

# #Close file
# jointSetExtFile.close()

# #Convert Falisse et al. 2D contact parameter sizes and locations to present
# #model based on variations in body size from scaling. This relates to spheres
# #placed at the heel and the toe in the 2D model

# #Set the original parameters from the PredSim_2D model
# radiusSphere_heel = 0.035
# radiusSphere_front = 0.015
# locSphere_heel_l = np.array([0.031307527581931796, 0.010435842527310599, 0])
# locSphere_front_l = np.array([0.1774093229642802, -0.015653763790965898, -0.005217921263655299])
# locSphere_heel_r = np.array([0.031307527581931796, 0.010435842527310599, 0])
# locSphere_front_r = np.array([0.1774093229642802, -0.015653763790965898, 0.005217921263655299])

# #Our generic model is a relatively similar size (not so much mass) to the 2D
# #model from Falisse et al., so we'll scale these parameters based on the original
# #scaling factors. Sphere size and placement doesn't matter so much as optimising
# #kinematics when trying to track GRFs --- so the placement here doesn't need to
# #be super refined

# #Use the same scale factors as before to scale the size and location of the 2D spheres
# #Left heel location
# locSphereNew_heel_l = np.array([locSphere_heel_l[0] * scaleFactors['calcn_l'][0],
#                                 locSphere_heel_l[1] * scaleFactors['calcn_l'][1],
#                                 locSphere_heel_l[2] * scaleFactors['calcn_l'][2]])
# #Right heel location
# locSphereNew_heel_r = np.array([locSphere_heel_r[0] * scaleFactors['calcn_r'][0],
#                                 locSphere_heel_r[1] * scaleFactors['calcn_r'][1],
#                                 locSphere_heel_r[2] * scaleFactors['calcn_r'][2]])
# #Left front location
# locSphereNew_front_l = np.array([locSphere_front_l[0] * scaleFactors['toes_l'][0],
#                                  locSphere_front_l[1] * scaleFactors['toes_l'][1],
#                                  locSphere_front_l[2] * scaleFactors['toes_l'][2]])
# #Right front location
# locSphereNew_front_r = np.array([locSphere_front_r[0] * scaleFactors['toes_r'][0],
#                                  locSphere_front_r[1] * scaleFactors['toes_r'][1],
#                                  locSphere_front_r[2] * scaleFactors['toes_r'][2]])
# #Heel size
# radiusSphereNew_heel = radiusSphere_heel * ((heelSphereScale_l + heelSphereScale_r) / 2)
# #Front size
# radiusSphereNew_front = radiusSphere_front * ((toesSphereScale_l + toesSphereScale_r) / 2)

# #Write this contact sphere parameters section to a text file
# contactParamsExtFile = open('contactParameters_externalFunction_2D.txt', 'w')
# contactParamsExtFile.write('osim_double_adouble radiusSphere_heel = '+str(radiusSphereNew_heel)+';\n')
# contactParamsExtFile.write('osim_double_adouble radiusSphere_front = '+str(radiusSphereNew_front)+';\n')
# contactParamsExtFile.write('osim_double_adouble stiffness_heel = 3067776;\n')          
# contactParamsExtFile.write('osim_double_adouble stiffness_front = 3067776;\n')
# contactParamsExtFile.write('osim_double_adouble dissipation = 2.0;\n')
# contactParamsExtFile.write('osim_double_adouble staticFriction = 0.8;\n')
# contactParamsExtFile.write('osim_double_adouble dynamicFriction = 0.8;\n')
# contactParamsExtFile.write('osim_double_adouble viscousFriction = 0.5;\n')
# contactParamsExtFile.write('osim_double_adouble transitionVelocity = 0.2;\n')
# contactParamsExtFile.write('Vec3 halfSpaceLocation(0);\n')
# contactParamsExtFile.write('Vec3 halfSpaceOrientation(0, 0, -0.5*SimTK::Pi);\n')
# contactParamsExtFile.write('Vec3 locSphere_heel_l = Vec3('+\
#                                str(locSphereNew_heel_l[0])+', '+str(locSphereNew_heel_l[1])+\
#                                    ', '+str(locSphereNew_heel_l[2])+');\n')
# contactParamsExtFile.write('Vec3 locSphere_front_l = Vec3('+\
#                                str(locSphereNew_front_l[0])+', '+str(locSphereNew_front_l[1])+\
#                                    ', '+str(locSphereNew_front_l[2])+');\n')
# contactParamsExtFile.write('Vec3 locSphere_heel_r = Vec3('+\
#                                str(locSphereNew_heel_r[0])+', '+str(locSphereNew_heel_r[1])+\
#                                    ', '+str(locSphereNew_heel_r[2])+');\n')
# contactParamsExtFile.write('Vec3 locSphere_front_r = Vec3('+\
#                                str(locSphereNew_front_r[0])+', '+str(locSphereNew_front_r[1])+\
#                                    ', '+str(locSphereNew_front_r[2])+');\n')
# #Close file
# contactParamsExtFile.close()

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

# %% Run a torque driven tracking problem to get consistent kinematics with contact
# At this point we convert to a 2D problem.
# Following this step we should have 2D kinematics that match the experimental 
# 2D GRFs coming from contact sphere estimation

#Set model and parameters
#Create a model processor
modelProcessor = osim.ModelProcessor(osim.Model('model2D.osim'))
#Delete muscles
modelProcessor.append(osim.ModOpRemoveMuscles())

#Get model from the processor
torqueModel = modelProcessor.process()

#Add reserves for torque actuation
#Set list of coordinates to actuate
actuateList = ['hip_flexion_r', 'hip_flexion_l',
               'knee_angle_r', 'knee_angle_l',
               'ankle_angle_r', 'ankle_angle_l',
               'mtp_angle_r', 'mtp_angle_l']
               #'lumbar_extension'] # not needed - already there
#Add actuators
for aa in range(0,len(actuateList)):
    osimHelper.addCoordinateActuator(osimModel = torqueModel,
                                     coordName = actuateList[aa],
                                     optForce = 10)

#Finalise model connections
torqueModel.finalizeConnections()

#Define the motion tracking problem
track = osim.MocoTrack()
track.setName('torqueDriven_GRFtracking_2D')

#Set model
track.setModel(osim.ModelProcessor(torqueModel))

#Set kinematics
tableProcessor = osim.TableProcessor('refQ.sto')

#Set states reference details
track.setStatesReference(tableProcessor) # Apply the target data to the tracking problem
track.set_states_global_tracking_weight(1) # Default tracking weight (is changed below)
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
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tx/value',  (1/(1*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_ty/value',  (1/(2*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tilt/value',(1/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/value',           (1/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/value',     (1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/walker_knee_r/knee_angle_r/value',     (1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/value',   (1/(1*0.0574))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/mtp_r/mtp_angle_r/value',       (1/(2*0.0873))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/value',     (1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/walker_knee_l/knee_angle_l/value',     (1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/value',   (1/(1*0.0574))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/mtp_l/mtp_angle_l/value',       (1/(2*0.0873))**2))
#Joint speeds
w = 0.001 #Scale the generalized speed tracking errors by this constant
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tx/speed',  w*(0/(1*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_ty/speed',  w*(0/(2*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tilt/speed',w*(0/(1*0.0585))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/speed',           w*(0/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/speed',     w*(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/walker_knee_r/knee_angle_r/speed',     w*(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/speed',   w*(1/(1*0.0574))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/mtp_r/mtp_angle_r/speed',       w*(1/(1*0.0873))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/speed',     w*(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/walker_knee_l/knee_angle_l/speed',     w*(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/speed',   w*(1/(1*0.0574))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/mtp_l/mtp_angle_l/speed',       w*(1/(1*0.0873))**2))
# #Add to tracking problem
track.set_states_weight_set(stateWeights)

#Define the Moco study and problem
study = track.initialize()
problem = study.updProblem()

#Define the periodicity goal
#Keeping in mind that the simulation is of a half gait cycle, similar to that
#proposed in Umberger's Moco examples
periodicityGoal = osim.MocoPeriodicityGoal('symmetryGoal')
problem.addGoal(periodicityGoal)
torqueModel.initSystem()

#Set symmetric coordinate values
#Set symmetry pairs
symPairs = [['hip_flexion_r','hip_flexion_l'],
            ['knee_angle_r','knee_angle_l'],
            ['ankle_angle_r','ankle_angle_l'],
            ['mtp_angle_r','mtp_angle_l']]
for ii in range(0,len(symPairs)):
    #Set the jointsets depending on current pair
    if 'hip' in symPairs[ii][0]:
        jointName = ['/jointset/hip_r/','/jointset/hip_l/']
    elif 'knee' in symPairs[ii][0]:
        jointName = ['/jointset/walker_knee_r/','/jointset/walker_knee_l/']
    elif 'ankle' in symPairs[ii][0]:
        jointName = ['/jointset/ankle_r/','/jointset/ankle_l/']
    elif 'mtp' in symPairs[ii][0]:
        jointName = ['/jointset/mtp_r/','/jointset/mtp_l/']
    
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

#Singular periodic kinematics
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_ty/value'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_tilt/value'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/back/lumbar_extension/value'))
        
#Symmetric coordinate actuator controls
periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair('/forceset/tau_lumbar_ext'))
for aa in range(0,len(actuateList)):
    periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair('/forceset/tau_'+actuateList[aa]))

#Regularization term on MocoTrack problem (minimize squared muscle excitations)
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
effort.setWeight(0.01)

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
problem.setStateInfo('/jointset/mtp_r/mtp_angle_r/value', [math.radians(-60), math.radians(60)])
problem.setStateInfo('/jointset/mtp_l/mtp_angle_l/value', [math.radians(-60), math.radians(60)])

#Add contact tracking goal
#Create GRF tracking goal
contactTracking = osim.MocoContactTrackingGoal('contact',1)
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

#Define the solver and set its options
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
#solver.set_multibody_dynamics_mode('implicit')
#solver.set_minimize_implicit_multibody_accelerations(true)
#solver.set_implicit_multibody_accelerations_weight(0.00001)
solver.set_optim_max_iterations(1000)
solver.set_num_mesh_intervals(100)
solver.set_optim_constraint_tolerance(1e-2) ####TODO: might need to be lower
solver.set_optim_convergence_tolerance(1e-2) ####TODO: might need to be lower
# solver.set_minimize_implicit_auxiliary_derivatives(True)
# solver.set_implicit_auxiliary_derivatives_weight(0.001)
solver.resetProblem(problem)

#Solve!
grfTorqueSolution = study.solve()

#Calculate predicted GRFs from 2D tracking simulation
#Add contact elements to extract from
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
#Create forces table
externalForcesTableFlat = osim.createExternalLoadsTableForGait(torqueModel,
                                                               grfTorqueSolution,
                                                               contact_r,
                                                               contact_l)
#Write table to file
osim.writeTableToFile(externalForcesTableFlat,'predictGRF_2D.sto')

##### TODO: compare kinematics and GRFs in figure

# %% Torque driven predictive at higher speeds

#Define the study
study = osim.MocoStudy()
study.setName('torqueDriven_prediction')

#Get problem
problem = study.updProblem()

#Set the model
problem.setModelProcessor(osim.ModelProcessor(torqueModel))

#Define the periodicity goal
#Keeping in mind that the simulation is of a half gait cycle, similar to that
#proposed in Umberger's Moco examples
periodicityGoal = osim.MocoPeriodicityGoal('symmetryGoal')
problem.addGoal(periodicityGoal)
torqueModel.initSystem()

#Set symmetric coordinate values
#Set symmetry pairs
symPairs = [['hip_flexion_r','hip_flexion_l'],
            ['knee_angle_r','knee_angle_l'],
            ['ankle_angle_r','ankle_angle_l'],
            ['mtp_angle_r','mtp_angle_l']]
for ii in range(0,len(symPairs)):
    #Set the jointsets depending on current pair
    if 'hip' in symPairs[ii][0]:
        jointName = ['/jointset/hip_r/','/jointset/hip_l/']
    elif 'knee' in symPairs[ii][0]:
        jointName = ['/jointset/walker_knee_r/','/jointset/walker_knee_l/']
    elif 'ankle' in symPairs[ii][0]:
        jointName = ['/jointset/ankle_r/','/jointset/ankle_l/']
    elif 'mtp' in symPairs[ii][0]:
        jointName = ['/jointset/mtp_r/','/jointset/mtp_l/']
    
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

#Singular periodic kinematics
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_ty/value'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_tilt/value'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/back/lumbar_extension/value'))
        
#Symmetric coordinate actuator controls
periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair('/forceset/tau_lumbar_ext'))
for aa in range(0,len(actuateList)):
    periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair('/forceset/tau_'+actuateList[aa]))

#Prescribed average gait speed
#Note this is upped from trial which is ~6.6
speedGoal = osim.MocoAverageSpeedGoal('speed')
problem.addGoal(speedGoal)
speedGoal.set_desired_average_speed(8.0)

#Set for effort over distance
effortGoal = osim.MocoControlGoal('effort', 10)
problem.addGoal(effortGoal)
effortGoal.setExponent(2)
effortGoal.setDivideByDisplacement(True)

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
problem.setStateInfo('/jointset/mtp_r/mtp_angle_r/value', [math.radians(-60), math.radians(60)])
problem.setStateInfo('/jointset/mtp_l/mtp_angle_l/value', [math.radians(-60), math.radians(60)])

#Define the solver and set its options
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
#solver.set_multibody_dynamics_mode('implicit')
#solver.set_minimize_implicit_multibody_accelerations(true)
#solver.set_implicit_multibody_accelerations_weight(0.00001)
solver.set_optim_max_iterations(1000)
solver.set_num_mesh_intervals(100)
solver.set_optim_constraint_tolerance(1e-2) ####TODO: might need to be lower
solver.set_optim_convergence_tolerance(1e-2) ####TODO: might need to be lower
# solver.set_minimize_implicit_auxiliary_derivatives(True)
# solver.set_implicit_auxiliary_derivatives_weight(0.001)
solver.resetProblem(problem)
solver.setGuess(grfTorqueSolution) #Use tracking solution as initial guess

#Solve!
torquePredictionSolution = study.solve()

#### Torque prediction does what it's meant to, but is
#### not that realistic

#### Try muscle driven tracking sim of torque driven kinematics...

##### TODO: ensure unlocking speeds for joint coordinates and + actuator forces


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


    








# %%
