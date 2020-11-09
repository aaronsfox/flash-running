# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 21:12:27 2020

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This script process the data provided by Dorn et al. (2012) through:
        - Inverse kinematics
    
"""

# %% Import packages

import opensim as osim

# %% Inverse kinematics

#Load model
osimModel = osim.Model('JA1_SCALED_Osim40.osim')

#Lock the coordinates that we don't need data from for the Falisse pipeline
lockCoordinates = ['mtp_angle_r', 'mtp_angle_l',
                   'pro_sup_r', 'pro_sup_l',
                   'wrist_flex_r', 'wrist_flex_l',
                   'wrist_dev_r', 'wrist_dev_l']
for cc in range(len(lockCoordinates)):
    osimModel.getCoordinateSet().get(lockCoordinates[cc]).set_locked(True)
osimModel.finalizeConnections()

#Initialise IK tool
ikTool = osim.InverseKinematicsTool()
ikTool.setName('sprint')

#Set model
ikTool.setModel(osimModel)

#Set marker file
ikTool.set_marker_file('sprint.trc')

#Create and set task set based off .trc markers
#TODO: make this less generic
#Create task set
ikTaskSet = osim.IKTaskSet()
#Get marker labels
trcData = osim.TimeSeriesTableVec3('sprint.trc')
markerLabels = trcData.getColumnLabels()
#Loop through and set a task value
#TODO: change from default value
for mm in range(len(markerLabels)):
    #Create blank marker task
    markerTask = osim.IKMarkerTask()
    #Set marker name
    markerTask.setName(markerLabels[mm])
    #Set marker weight
    markerTask.setWeight(1)
    #Clone and append to task set
    ikTaskSet.cloneAndAppend(markerTask)
    
#Set task set
ikTool.set_IKTaskSet(ikTaskSet)

#Set times
#Manually identified times from experimental data signifying first to second
#right foot strike
ikTool.setStartTime(0.359)
ikTool.setEndTime(0.819)

#Set output filename
ikTool.set_output_motion_file('sprint_ik.mot')

#Run IK
ikTool.run()

# %% Inverse dynamics

#Initialise ID tool
idTool = osim.InverseDynamicsTool()
idTool.setName('sprint')

#Set model
idTool.setModel(osimModel)

#Set external loads file
idTool.setExternalLoadsFileName('grf.xml')

#Set coordinates file
idTool.setCoordinatesFileName('sprint_ik.mot')

#Set filter of 12Hz for coordinates (matches some sprinting studies)
idTool.setLowpassCutoffFrequency(12)

#Set times
#Start and end of IK
idTool.setStartTime(osim.Storage('sprint_ik.mot').getFirstTime())
idTool.setEndTime(osim.Storage('sprint_ik.mot').getLastTime())

#Set output file
idTool.setOutputGenForceFileName('sprint_id.sto')

#Run ID
idTool.run()

# %% ----- End of processData.py ----- %% #