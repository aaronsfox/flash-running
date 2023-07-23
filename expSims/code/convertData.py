# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 20:25:31 2020

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This script converts the data provided by Dorn et al. (2012) into 
    OpenSim format for futher processing.
    
"""

# %% Import packages

import opensim as osim
import numpy as np
from scipy.signal import butter, filtfilt

# %% Set-up

#Set-up initial logger for processing experimental data
osim.Logger.removeFileSink()
osim.Logger.addFileSink('convertDataLog.log')

# %% Convert c3d data to OpenSim formats

# %% Static trial

#Set static file
staticFile = '..\\data\\JA1Static05.c3d'

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
    
#Remove black markers with '*' notation
#Loop through marker labels
markerLabels = staticMarkers.getColumnLabels()
for mm in range(len(markerLabels)):
    #Check if current marker starts with *
    if markerLabels[mm].startswith('*'):
        #Remove it
        staticMarkers.removeColumn(markerLabels[mm])
    
#Write static markers to TRC file
osim.TRCFileAdapter().write(staticMarkers, '..\\data\\static.trc')

# %% Dynamic trial

#Set the dynamic file name
dynamicFile = '..\\data\\JA1Gait35_9ms.c3d'
    
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
    
#Remove black markers with '*' notation
#Loop through marker labels
markerLabels = dynamicMarkers.getColumnLabels()
for mm in range(len(markerLabels)):
    #Check if current marker starts with *
    if markerLabels[mm].startswith('*'):
        #Remove it
        dynamicMarkers.removeColumn(markerLabels[mm])
    
#Write markers to TRC file
osim.TRCFileAdapter().write(dynamicMarkers, '..\\data\\sprint.trc')

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
forcesStorage.printResult(forcesStorage, 'sprint_grf', '..\\data', 0.001, '.mot')

#Create the external loads .xml file
#Set the different force plates to the varying foot contacts
#have which contacts
#Right foot = fp2,fp7
#Left foot = fp4,fp8
forceXML = osim.ExternalLoads()

#Create and append the right GRF external forces
#FP2
rightGRF1 = osim.ExternalForce()
rightGRF1.setName('RightGRF1')
rightGRF1.setAppliedToBodyName('calcn_r')
rightGRF1.setForceExpressedInBodyName('ground')
rightGRF1.setPointExpressedInBodyName('ground')
rightGRF1.setForceIdentifier('ground_force_2_v')
rightGRF1.setPointIdentifier('ground_force_2_p')
rightGRF1.setTorqueIdentifier('ground_force_2_m')
forceXML.cloneAndAppend(rightGRF1)
#FP7
rightGRF2 = osim.ExternalForce()
rightGRF2.setName('RightGRF2')
rightGRF2.setAppliedToBodyName('calcn_r')
rightGRF2.setForceExpressedInBodyName('ground')
rightGRF2.setPointExpressedInBodyName('ground')
rightGRF2.setForceIdentifier('ground_force_7_v')
rightGRF2.setPointIdentifier('ground_force_7_p')
rightGRF2.setTorqueIdentifier('ground_force_7_m')
forceXML.cloneAndAppend(rightGRF2)

#Create and append the left GRF external forces
#FP4
leftGRF1 = osim.ExternalForce()
leftGRF1.setName('LeftGRF1')
leftGRF1.setAppliedToBodyName('calcn_l')
leftGRF1.setForceExpressedInBodyName('ground')
leftGRF1.setPointExpressedInBodyName('ground')
leftGRF1.setForceIdentifier('ground_force_4_v')
leftGRF1.setPointIdentifier('ground_force_4_p')
leftGRF1.setTorqueIdentifier('ground_force_4_m')
forceXML.cloneAndAppend(leftGRF1)
#FP8
leftGRF2 = osim.ExternalForce()
leftGRF2.setName('LeftGRF2')
leftGRF2.setAppliedToBodyName('calcn_l')
leftGRF2.setForceExpressedInBodyName('ground')
leftGRF2.setPointExpressedInBodyName('ground')
leftGRF2.setForceIdentifier('ground_force_8_v')
leftGRF2.setPointIdentifier('ground_force_8_p')
leftGRF2.setTorqueIdentifier('ground_force_8_m')
forceXML.cloneAndAppend(leftGRF2)

#Set GRF datafile
forceXML.setDataFileName('sprint_grf.mot')

#Write to file
forceXML.printToXML('..\\data\\sprint_grf.xml')

# %% ----- End of convertData.py ----- %% #