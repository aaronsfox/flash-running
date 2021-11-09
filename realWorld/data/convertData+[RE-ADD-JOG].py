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
import btk
import numpy as np
import os
import pandas as pd
from scipy.signal import butter, lfilter
from scipy.spatial.transform import Rotation as R

# %% Convert c3d data to OpenSim formats

# %% Static trial

#Convert static c3d to .trc

#Set up data adapters
staticC3DFile = osim.C3DFileAdapter()
staticTRCFile = osim.TRCFileAdapter()

#Get markers
staticData = staticC3DFile.read('JA1Static05.c3d')
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
    
#Remove black markers with '*' notation
#Loop through marker labels
markerLabels = staticMarkers.getColumnLabels()
for mm in range(len(markerLabels)):
    #Check if current marker starts with *
    if markerLabels[mm].startswith('*'):
        #Remove it
        staticMarkers.removeColumn(markerLabels[mm])
    
#Write static data to file
staticTRCFile.write(staticMarkers, 'static.trc')

# %% Dynamic trials

# %% Sprint trial

#Convert dynamic c3d to .trc

#Set up data adapters
dynamicC3DFile = osim.C3DFileAdapter()
dynamicTRCFile = osim.TRCFileAdapter()

#Get markers
dynamicData = dynamicC3DFile.read('JA1Gait35.c3d')
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

#Remove black markers with '*' notation
#Loop through marker labels
markerLabels = dynamicMarkers.getColumnLabels()
for mm in range(len(markerLabels)):
    #Check if current marker starts with *
    if markerLabels[mm].startswith('*'):
        #Remove it
        dynamicMarkers.removeColumn(markerLabels[mm])    

#Write static data to file
dynamicTRCFile.write(dynamicMarkers, 'sprint.trc')

#Convert c3d force data to .mot and .xml files

#Load in the c3d data via btk
dynamicC3D = btk.btkAcquisitionFileReader()
dynamicC3D.SetFilename('JA1Gait35.c3d')
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
fs = dynamicAcq.GetAnalogFrequency()
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
#Right foot = fp2,fp7
#Left foot = fp4,fp8
#Right forces
ground_force_r_vx = df_forces['Fx2'].values + df_forces['Fx7'].values
ground_force_r_vy = df_forces['Fy2'].values + df_forces['Fy7'].values
ground_force_r_vz = df_forces['Fz2'].values + df_forces['Fz7'].values
#Left forces
ground_force_l_vx = df_forces['Fx4'].values + df_forces['Fx8'].values
ground_force_l_vy = df_forces['Fy4'].values + df_forces['Fy8'].values
ground_force_l_vz = df_forces['Fz4'].values + df_forces['Fz8'].values
#Right torques
ground_torque_r_x = df_forces['Mx2'].values + df_forces['Mx7'].values
ground_torque_r_y = df_forces['My2'].values + df_forces['My7'].values
ground_torque_r_z = df_forces['Mz2'].values + df_forces['Mz7'].values
#Left torques
ground_torque_l_x = df_forces['Mx4'].values + df_forces['Mx8'].values
ground_torque_l_y = df_forces['My4'].values + df_forces['My8'].values
ground_torque_l_z = df_forces['Mz4'].values + df_forces['Mz8'].values
#Right position
ground_force_r_px = df_forces['Px2'].values + df_forces['Px7'].values
ground_force_r_py = df_forces['Py2'].values + df_forces['Py7'].values
ground_force_r_pz = df_forces['Pz2'].values + df_forces['Pz7'].values
#Left position
ground_force_l_px = df_forces['Px4'].values + df_forces['Px8'].values
ground_force_l_py = df_forces['Py4'].values + df_forces['Py8'].values
ground_force_l_pz = df_forces['Pz4'].values + df_forces['Pz8'].values

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
forcesStorage.setName('grf')
#Print to file
forcesStorage.printResult(forcesStorage, 'sprint_grf', os.getcwd(), 0.001, '.mot')

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
forceXML.setDataFileName('sprint_grf.mot')
#Set filtering for kinematics
forceXML.setLowpassCutoffFrequencyForLoadKinematics(12)
#Write to file
forceXML.printToXML('sprint_grf.xml')

# %% ----- End of convertData.py ----- %% #