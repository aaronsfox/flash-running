# -*- coding: utf-8 -*-
"""

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
import opensim as osim
import osimHelper
import numpy as np
from scipy.signal import butter, filtfilt
from distutils.dir_util import copy_tree
import shutil

# %% Set-up

#Set-up initial logger for processing experimental data
osim.Logger.removeFileSink()
osim.Logger.addFileSink('prepExpDataLog.log')

### TODO: add any set-up code here

# #Set matplotlib parameters
# from matplotlib import rcParams
# # rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = 'Arial'
# rcParams['font.weight'] = 'bold'
# rcParams['axes.labelsize'] = 12
# rcParams['axes.titlesize'] = 16
# rcParams['axes.linewidth'] = 1.5
# rcParams['axes.labelweight'] = 'bold'
# rcParams['legend.fontsize'] = 10
# rcParams['xtick.major.width'] = 1.5
# rcParams['ytick.major.width'] = 1.5
# rcParams['legend.framealpha'] = 0.0
# rcParams['savefig.dpi'] = 300
# rcParams['savefig.format'] = 'pdf'

# %% Process the static trial

#Set static file
staticFile = 'Data\\static.c3d'

#Construct opensim 3d object
c3dFile = osim.C3DFileAdapter()
c3dFile.setLocationForForceExpression(osim.C3DFileAdapter.ForceLocation_CenterOfPressure)

#Read in the static trial
staticC3D = c3dFile.read(staticFile)

#Get markers table
staticMarkers = c3dFile.getMarkersTable(staticC3D)

#Rotate marker data
#Create the two rotations needed
markerRot1 = osim.Rotation(np.deg2rad(-90), osim.Vec3(0,0,1))
markerRot2 = osim.Rotation(np.deg2rad(-90), osim.Vec3(1,0,0))
#Rotate the data
for iRow in range(staticMarkers.getNumRows()):
    #Apply the two rotations
    staticMarkers.setRowAtIndex(iRow, markerRot1.multiply(staticMarkers.getRowAtIndex(iRow)))
    staticMarkers.setRowAtIndex(iRow, markerRot2.multiply(staticMarkers.getRowAtIndex(iRow)))
    
#Write static markers to TRC file
osim.TRCFileAdapter().write(staticMarkers, f'{staticFile[0:-4]}.trc')
                            
#Add the virtual markers and create a new .trc file to use in scaling
osimHelper.addVirtualMarkersStatic(f'{staticFile[0:-4]}.trc',
                                   f'{staticFile[0:-4]}_virtualMarkers.trc')

#Get force data to calculate mass
#Participant static trial is done on force plate 4
staticForces = c3dFile.getForcesTable(staticC3D).flatten()
massKg = np.mean(staticForces.getDependentColumn('f4_3').to_numpy()) / 9.81

# %% Process the experimental trial

#Set the dynamic file name
dynamicFile = 'Data\\maxSprint.c3d'
    
#Construct opensim 3d object
c3dFile = osim.C3DFileAdapter()
c3dFile.setLocationForForceExpression(osim.C3DFileAdapter.ForceLocation_CenterOfPressure)

#Read in the c3d file
dynamicC3D = c3dFile.read(dynamicFile)

#Get markers table
dynamicMarkers = c3dFile.getMarkersTable(dynamicC3D)

#Rotate the data
#Use the same two rotations as earlier
for iRow in range(dynamicMarkers.getNumRows()):
    #Apply the two rotations
    dynamicMarkers.setRowAtIndex(iRow, markerRot1.multiply(dynamicMarkers.getRowAtIndex(iRow)))
    dynamicMarkers.setRowAtIndex(iRow, markerRot2.multiply(dynamicMarkers.getRowAtIndex(iRow)))
    
#Write markers to TRC file
osim.TRCFileAdapter().write(dynamicMarkers, f'{dynamicFile[0:-4]}.trc')

#Add the virtual markers and create a new .trc file to use in IK
osimHelper.addVirtualMarkersDynamic(staticTRC = f'{staticFile[0:-4]}_virtualMarkers.trc',
                                    dynamicTRC = f'{dynamicFile[0:-4]}.trc',
                                    outputTRC = f'{dynamicFile[0:-4]}_virtualMarkers.trc')

#Extract the GRF data

#Get forces table
dynamicForces = c3dFile.getForcesTable(dynamicC3D)

#Rotate forces data
#Use the same rotations as earlier
for iRow in range(dynamicForces.getNumRows()):
    dynamicForces.setRowAtIndex(iRow, markerRot1.multiply(dynamicForces.getRowAtIndex(iRow)))
    dynamicForces.setRowAtIndex(iRow, markerRot2.multiply(dynamicForces.getRowAtIndex(iRow)))
    
#Flatten forces data
forcesFlat = dynamicForces.flatten()

#Convert to numpy array
#Pre-allocate numpy array based on data size
dataArray = np.zeros((forcesFlat.getNumRows(),
                      forcesFlat.getNumColumns()))
#Extract data
for forceInd in range(forcesFlat.getNumColumns()):
    dataArray[:,forceInd] = forcesFlat.getDependentColumn(forcesFlat.getColumnLabels()[forceInd]).to_numpy()
    
#Replace nan's for COP and moment data with zeros
np.nan_to_num(dataArray, copy = False, nan = 0.0)

#Convert force point data from mm to m
for forceName in list(forcesFlat.getColumnLabels()):
    if forceName.startswith('p') or forceName.startswith('m'):
        #Get force index
        forceInd = list(forcesFlat.getColumnLabels()).index(forceName)
        #Convert to m units in data array
        dataArray[:,forceInd] = dataArray[:,forceInd] / 1000

#Filter force data
    
#Get the sampling rate
fs = float(dynamicForces.getTableMetaDataAsString('DataRate'))

#Define low-pass Butterworth filter
filtFreq = 50
nyq = 0.5 * fs
normCutoff = filtFreq / nyq
b, a = butter(4, normCutoff, btype = 'low', analog = False)

#Apply lowpass filter to data
for forceName in list(forcesFlat.getColumnLabels()):
    #Get force index
    forceInd = list(forcesFlat.getColumnLabels()).index(forceName)
    #Apply filter
    dataArray[:,forceInd] = filtfilt(b, a, dataArray[:,forceInd])
    
#Build the new time series table
forcesStorage = osim.Storage()

#Get the time data
time = forcesFlat.getIndependentColumn()

#Create maps to replace text from force labels with
#Force plate and type identifiers
forceType = {}
for ii in range(1,9):
    forceType[f'f{ii}'] = f'ground_force_{ii}_v'
    forceType[f'p{ii}'] = f'ground_force_{ii}_p'
    forceType[f'm{ii}'] = f'ground_force_{ii}_m'
#Axis identifiers
forceAxis = {'1': 'x',
             '2': 'y',
             '3': 'z'}

#Set labels in table
newLabels = osim.ArrayStr()
newLabels.append('time')
for forceLabel in forcesFlat.getColumnLabels():
    #Split the label to get parts
    labelSplit = forceLabel.split('_')
    #Create new label
    forceLabel = f'{forceType[labelSplit[0]]}{forceAxis[labelSplit[1]]}'
    #Append to labels vector
    newLabels.append(forceLabel)
forcesStorage.setColumnLabels(newLabels)

#Add data
for iRow in range(dataArray.shape[0]):
    row = osim.ArrayDouble()
    for iCol in range(dataArray.shape[1]):
        row.append(dataArray[iRow,iCol])
    #Add data to storage
    forcesStorage.append(time[iRow], row)

#Set name for storage object
forcesStorage.setName('sprint_grf')

#Write to file
forcesStorage.printResult(forcesStorage, 'refGRF', 'Data', 0.001, '.mot')

#Create the external loads .xml file
#Set the different force plates to the varying foot contacts
#have which contacts
#Right foot = fp7,fp3 
#Left foot = fp5,fp1
forceXML = osim.ExternalLoads()

#Create and append the right GRF external forces
#FP7
rightGRF1 = osim.ExternalForce()
rightGRF1.setName('RightGRF1')
rightGRF1.setAppliedToBodyName('calcn_r')
rightGRF1.setForceExpressedInBodyName('ground')
rightGRF1.setPointExpressedInBodyName('ground')
rightGRF1.setForceIdentifier('ground_force_7_v')
rightGRF1.setPointIdentifier('ground_force_7_p')
rightGRF1.setTorqueIdentifier('ground_force_7_m')
forceXML.cloneAndAppend(rightGRF1)
#FP3
rightGRF2 = osim.ExternalForce()
rightGRF2.setName('RightGRF2')
rightGRF2.setAppliedToBodyName('calcn_r')
rightGRF2.setForceExpressedInBodyName('ground')
rightGRF2.setPointExpressedInBodyName('ground')
rightGRF2.setForceIdentifier('ground_force_3_v')
rightGRF2.setPointIdentifier('ground_force_3_p')
rightGRF2.setTorqueIdentifier('ground_force_3_m')
forceXML.cloneAndAppend(rightGRF2)

#Create and append the left GRF external forces
#FP5
leftGRF1 = osim.ExternalForce()
leftGRF1.setName('LeftGRF1')
leftGRF1.setAppliedToBodyName('calcn_l')
leftGRF1.setForceExpressedInBodyName('ground')
leftGRF1.setPointExpressedInBodyName('ground')
leftGRF1.setForceIdentifier('ground_force_5_v')
leftGRF1.setPointIdentifier('ground_force_5_p')
leftGRF1.setTorqueIdentifier('ground_force_5_m')
forceXML.cloneAndAppend(leftGRF1)
#FP1
leftGRF2 = osim.ExternalForce()
leftGRF2.setName('LeftGRF2')
leftGRF2.setAppliedToBodyName('calcn_l')
leftGRF2.setForceExpressedInBodyName('ground')
leftGRF2.setPointExpressedInBodyName('ground')
leftGRF2.setForceIdentifier('ground_force_1_v')
leftGRF2.setPointIdentifier('ground_force_1_p')
leftGRF2.setTorqueIdentifier('ground_force_1_m')
forceXML.cloneAndAppend(leftGRF2)

#Set GRF datafile
forceXML.setDataFileName('refGRF.mot')

#Write to file
forceXML.printToXML('Data\\refGRF.xml')

# %% Model scaling (3D Model)

#Navigate to scaling directory
os.chdir('Scaling')

#Set-up logger for scaling
osim.Logger.removeFileSink()
osim.Logger.addFileSink('scalingLog.log')

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

#Place in scale tool
scaleTool.getMarkerPlacer().setMarkerFileName('..\\Data\\static_virtualMarkers.trc')
scaleTool.getModelScaler().setMarkerFileName('..\\Data\\static_virtualMarkers.trc')

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

#### NOTE: this 3D model is simply for IK purposes, so no muscle-based or contact
#### sphere changes need to be made for now

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
scaleTool2D.getMarkerPlacer().setMarkerFileName('..\\Data\\static_virtualMarkers.trc')
scaleTool2D.getModelScaler().setMarkerFileName('..\\Data\\static_virtualMarkers.trc')

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

#Load model to edit
scaledModel2D = osim.Model('scaledModelAdjusted2D.osim')

#The locations and size of the contact sphere parameters go unchanged with
#standard model scaling, so these need to be edited to ensure they are in
#an appropriate place. This can be done based on the scale set parameters
#for the representative bodies.

#Load in the 2D scale set
scaleSet2D = osim.ScaleSet('scaleSet2D.xml')

#Grab average of the scale factors for the heel and toes
feetScaleFactors = {}
for segmentInd in range(scaleSet2D.getSize()):
    #Check if segment is one we need scale factors for
    if scaleSet2D.get(segmentInd).getSegmentName() in ['calcn_l', 'calcn_r', 'toes_l', 'toes_r']:
        feetScaleFactors[scaleSet2D.get(segmentInd).getSegmentName()] = scaleSet2D.get(segmentInd).getScaleFactors().to_numpy()

#Scale the radii to each of their respective factors
#While accessing the spheres, also adjust their position based on the scale
#factor for the respective axes
#Create a list of the sheres to loop through and edit
sphereList2D = ['heel_r','midfoot_r', 'toe_r',
                'heel_r','midfoot_l', 'toe_l']
#Loop through sphere list
for sphere in sphereList2D:
    #Get the current sphere
    currSphere = osim.ContactSphere.safeDownCast(scaledModel2D.getContactGeometrySet().get(sphere))
    #Set the current scaling factor based on sphere name
    if '_r' in sphere:
        if 'toe' in sphere or 'midfoot' in sphere:
            scaleFac = feetScaleFactors['toes_r'].mean()
            scaleBod = 'toes_r'
        else:
            scaleFac = feetScaleFactors['calcn_r'].mean()
            scaleBod = 'calcn_r'
    elif '_l' in sphere:
        if 'toe' in sphere or 'midfoot' in sphere:
            scaleFac = feetScaleFactors['toes_l'].mean()
            scaleBod = 'toes_l'
        else:
            scaleFac = feetScaleFactors['calcn_l'].mean()
            scaleBod = 'calcn_l'
    #Rescale the radius
    #Option to scale size here
    # currSphere.setRadius(currSphere.getRadius() * scaleFac)
    #Set new location
    newLoc = [currSphere.getLocation().get(ii) * feetScaleFactors[scaleBod][ii] for ii in [0,1,2]]
    currSphere.setLocation(osim.Vec3(newLoc[0],newLoc[1],newLoc[2]))

#Ensure all spheres lie on a flat plane with the foot in a neutral position (i.e. default pose)
#This will done by finding minimum base of sphere and ensuring each aligns with that
#on the Y-axis
#Loop through sphere list to identify minimum Y-value
minY = 999 #default unnecessarily high starting value
for sphere in sphereList2D:
    #Get the current sphere
    currSphere = osim.ContactSphere.safeDownCast(scaledModel2D.getContactGeometrySet().get(sphere))
    #Check the y-value minimum location taking centre position and radius
    currY = currSphere.getLocation().get(1) - currSphere.getRadius()
    #Replace min-Y if appropriate
    if currY < minY:
        minY = currSphere.getLocation().get(1) - currSphere.getRadius()
#Re-loop through spheres and adjust Y-location if it is higher than the min
for sphere in sphereList2D:
    #Get the current sphere
    currSphere = osim.ContactSphere.safeDownCast(scaledModel2D.getContactGeometrySet().get(sphere))
    #Check min location and adjust if necessary
    #Need to do absolutes as Y-value is negative
    if abs((currSphere.getLocation().get(1) - currSphere.getRadius())) < abs(minY):
        #Drop the y-value on the sphere
        #Calculate new y-loc
        diffY = abs(minY) - abs((currSphere.getLocation().get(1) - currSphere.getRadius()))
        newY = currSphere.getLocation().get(1) - diffY
        currSphere.setLocation(osim.Vec3(currSphere.getLocation().get(0),
                                         newY,
                                         currSphere.getLocation().get(2)))
    
#Reset model name
scaledModel2D.setName('gaitModel2D')

#Remove the marker set as it's no use in this 2D context
scaledModel2D.updMarkerSet().clearAndDestroy()

#Print out updated model to then adjust muscle forces with function
scaledModel2D.finalizeConnections()
scaledModel2D.printToXML('gaitModel2D.osim')

#Scale muscle strength based on linear function presented in Handsfield
#et al. (2014). This uses some convenience functions that are packaged
#with the Rajagopal et al. (2016) gait model. Note that the height of
#the generic model is 1.700; the current participant height is 1.79
osimHelper.scaleOptimalForceSubjectSpecific(genericModelFileName = 'gait9dof18musc_Ong_et_al_Moco.osim',
                                            scaledModelFileName = 'gaitModel2D.osim',
                                            genericHeight = 1.700, scaledHeight = 1.79,
                                            outputModelFileName = 'gaitModel2D.osim')

#Navigate back to experimental directory
os.chdir('..')

# %% Inverse kinematics

#Navigate to IK directory
os.chdir('IK')

#Set-up logger for IK
osim.Logger.removeFileSink()
osim.Logger.addFileSink('ikLog.log')

#Get start and end times for a half gait cycle
[startTime, endTime] = osimHelper.getGaitTimings('..\\Data\\refGRF.mot',
                                                 '..\\Data\\refGRF.xml',
                                                 'RightGRF1', 'LeftGRF1')

#Initialise IK tool
ikTool = osim.InverseKinematicsTool()

#Set model
ikTool.set_model_file('..\\Scaling\\scaledModelAdjusted.osim')

#Set task set
ikTool.set_IKTaskSet(osim.IKTaskSet('ikTasks.xml'))

#Set marker file
ikTool.set_marker_file('..\\Data\\maxSprint_virtualMarkers.trc')

#Set times
ikTool.setStartTime(startTime)
ikTool.setEndTime(endTime)

#Set output filename
ikTool.set_output_motion_file('ikResults.mot')

#Save setup file
ikTool.printToXML('setupIK.xml')

#Run IK
ikTool.run()

#Rename marker error file
#Check if this file needs to be deleted first
if os.path.isfile('ikMarkerErrors.sto'):
    os.remove('ikMarkerErrors.sto')
os.rename('_ik_marker_errors.sto', 'ikMarkerErrors.sto')

#Navigate back to experimental directory
os.chdir('..')

# %% Residual reduction algorithm

#Navigate to RRA directory
os.chdir('RRA')

#Set-up logger for RRA
osim.Logger.removeFileSink()
osim.Logger.addFileSink('rraLog.log')

#Create rra tool
rraTool = osim.RRATool()

#Set tool to replace model force set
rraTool.setReplaceForceSet(True)

#Set actuators file
forceSetFiles = osim.ArrayStr()
forceSetFiles.set(0,'rraActuators.xml')
rraTool.setForceSetFiles(forceSetFiles)

#Set tracking tasks file
rraTool.setTaskSetFileName('rraTasks.xml')

#Set a low pass filter frequency on the kinematics data
rraTool.setLowpassCutoffFrequency(12)

#Set to adjust the COM to reduce residuals
rraTool.setAdjustCOMToReduceResiduals(True)

#Set the torso body COM to be the one that gets adjusted
rraTool.setAdjustedCOMBody('torso')

#Set external loads file
rraTool.setExternalLoadsFileName('..\\Data\\refGRF.xml')

#Set tool name based on iteration
rraTool.setName('rra')

#Set desired kinematics file
rraTool.setDesiredKinematicsFileName('..\\IK\\ikResults.mot')

#Set initial and final time
#You need to use the IK first and last times here as I don't think the tool
#likes if the IK data doesn't have the times you put in
rraTool.setStartTime(osim.Storage('..\\IK\\ikResults.mot').getFirstTime())
rraTool.setFinalTime(osim.Storage('..\\IK\\ikResults.mot').getLastTime())

#Set to original scaled model
rraTool.setModelFilename('..\\Scaling\\scaledModelAdjusted.osim')
    
#Set output model file
rraTool.setOutputModelFileName('rraAdjustedModel.osim')

#Print out the tool
rraTool.printToXML('setupRRA.xml')

#Run RRA
#Reloading the setup tool seems to help with running the tool smoothly
runRRA = osim.RRATool('setupRRA.xml')
runRRA.run()

#Get the suggested mass change

#Open log and get lines
rraLog = open('rraLog.log')
rraLines = rraLog.readlines()

#Loop through lines in reverse to identify last occurring mass change
for line in reversed(rraLines):
    if line.startswith('*  Total mass change:'):
        #Get the mass change recommendation
        massChange = float(line.split(sep = ': ')[1])
        break
    
#Get original mass of model
currMass = osimHelper.getMassOfModel('..\\Scaling\\scaledModelAdjusted.osim')
    
#Set new mass and ratio
newMass = currMass + massChange
massScaleFac = newMass / currMass
    
#Apply mass change to rra model
#Adjust body masses
rraModel = osim.Model('rraAdjustedModel.osim')
allBodies = rraModel.getBodySet()
for ii in range(allBodies.getSize()):
    currBodyMass = allBodies.get(ii).getMass()
    newBodyMass = currBodyMass * massScaleFac
    allBodies.get(ii).setMass(newBodyMass)
#Finalise connections
rraModel.finalizeConnections()
#Print over model
rraModel.printToXML('rraAdjustedModel.osim')

#Apply mass change to 2D gait model
#Adjust body masses
gaitModel = osim.Model('..\\Scaling\\gaitModel2D.osim')
allBodies = gaitModel.getBodySet()
for ii in range(allBodies.getSize()):
    currBodyMass = allBodies.get(ii).getMass()
    newBodyMass = currBodyMass * massScaleFac
    allBodies.get(ii).setMass(newBodyMass)
#Finalise connections
gaitModel.finalizeConnections()
#Print new adjusted model here
gaitModel.printToXML('gaitModel2D.osim')

#Convert RRA kinematic results to a states file
osimHelper.kinematicsToStates(kinematicsFileName = 'rra_Kinematics_q.sto',
                              osimModelFileName = 'rraAdjustedModel.osim',
                              outputFileName = 'refQ.sto',
                              inDegrees = True, outDegrees = False)

#Convert states to 2D model format
osimHelper.statesTo2D(statesFileName = 'refQ.sto',
                      outputFileName = 'refQ_2D.sto',
                      renameWalkerKnee = True)

#Navigate back to experimental directory
os.chdir('..')

# %% Copy relevant files across to the tracking sim folder

#Navigate to tracking sims directory
os.chdir('..\\Simulations\\trackingSim')

#Copy geometry folder across to simulations directory
if not os.path.isdir('Geometry'):
    os.mkdir('Geometry')
copy_tree('..\\..\\ExpData\\Scaling\\Geometry',
          os.getcwd()+'\\Geometry')

#Copy GRF data files to simulations directory
os.chdir('..\\..\\ExpData\\Data')
shutil.copy('refGRF.mot', '..\\..\\Simulations\\trackingSim')
shutil.copy('refGRF.xml', '..\\..\\Simulations\\trackingSim')

#Copy kinematic data and model from RRA
os.chdir('..\\RRA')
shutil.copy('refQ_2D.sto', '..\\..\\Simulations\\trackingSim')
shutil.copy('gaitModel2D.osim', '..\\..\\Simulations\\trackingSim')

#Change back to experimental directory
os.chdir('..')

# %% Finish up
print('----- Processing of experimental data complete -----')

# %% ----- End of processExpData.py -----