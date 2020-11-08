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

#Initialise IK tool
ikTool = osim.InverseKinematicsTool()
ikTool.setName('sprint_ik')

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
#Start and end time of .trc file
ikTool.setStartTime(trcData.getIndependentColumn()[0])
ikTool.setEndTime(trcData.getIndependentColumn()[-1])

#Set output filename
ikTool.set_output_motion_file('sprint_ik.mot')

#Run IK
ikTool.run()

###### Things seem to be working up to here...

