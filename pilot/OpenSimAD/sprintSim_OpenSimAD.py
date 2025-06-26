# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This script uses similar concepts to the test of OpenSimAD function but 
    simplifies the process in using functions like the original examples.


"""

# %% Import packages

import os
import sys
import logging
import opensim as osim
import numpy as np
import glob
import shutil

# %% Import OpenSimAD tools

#Set base directory to current dir
baseDir = os.getcwd()

#Append directory to path
sys.path.append(baseDir)

#Set OpenSim AD utilities directory and append to path
opensimADDir = os.path.join(baseDir, 'supplementary', 'UtilsDynamicSimulations', 'OpenSimAD')
sys.path.append(opensimADDir)

#Set general utilities directory and append to path
utilsDir = os.path.join(baseDir, 'supplementary', 'UtilsGeneral')
sys.path.append(utilsDir)

#Import the various tools from the OpenSim AD utilities
from utilsOpenSimAD import processInputsOpenSimAD, plotResultsOpenSimAD
from mainOpenSimAD import run_tracking
# from utilsAuthentication import get_token

# %% Set-up

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

# Insert the name of the trial you want to simulate.
trial_name = 'coordinates'

# Insert the type of activity you want to simulate. We have pre-defined settings
# for different activities (more details above). Visit 
# ./UtilsDynamicSimulations/OpenSimAD/settingsOpenSimAD.py to see all available
# activities and their settings. If your activity is not in the list, select
# 'other' to use generic settings or set your own settings.
motion_type = 'running_torque_driven' ###note this was changed

# Insert the time interval you want to simulate. It is recommended to simulate
# trials shorter than 2s (more details above). Set to [] to simulate full trial.
# Note the time window is extracted from some earlier simulation tests
time_window = [0.352, 0.582]

# Insert the speed of the treadmill in m/s. A positive value indicates that the
# subject is moving forward. You should ignore this parameter or set it to 0 if
# the trial was not measured on a treadmill.
treadmill_speed = 0

# Insert a string to "name" you case.
case = 'torque_driven'

# %% Scale the JA1 model to the LaiUhlrich model

#Set the XYZ marker pairs for each body in the model
scaleMarkerPairs = {
    'pelvis': [   [['RASI', 'SACR'], ['LASI', 'SACR']],   [['RASI', 'LASI']],   [['RASI', 'LASI']]   ],
    'femur_r': [   [['RLEPI', 'RMEPI']],   [['RASI', 'RLEPI']],   [['RLEPI', 'RMEPI']]   ],
    'tibia_r': [   [['RLMAL', 'RMMAL']],   [['RLEPI', 'RLMAL'], ['RMEPI', 'RMMAL']],   [['RLEPI', 'RMEPI']]   ],
    'patella_r': [   [['RLEPI', 'RMEPI']],   [['RLEPI', 'RMEPI']],   [['RLEPI', 'RMEPI']]   ],
    'talus_r': [   [['RMMAL', 'RLMAL']],   [['RMMAL', 'RLMAL']],   [['RMMAL', 'RLMAL']],   ],
    'calcn_r': [   [['RHEEL', 'RLMAL']],   [['RHEEL', 'RLMAL']],   [['RHEEL', 'RLMAL']],   ],
    'toes_r': [   [['RHEEL', 'RTOE']],   [['RHEEL', 'RTOE']],   [['RP1MT', 'RP5MT']],   ],
    'femur_l': [   [['LLEPI', 'LMEPI']],   [['LASI', 'LLEPI']],   [['LLEPI', 'LMEPI']]   ],
    'tibia_l': [   [['LLMAL', 'LMMAL']],   [['LLEPI', 'LLMAL'], ['LMEPI', 'LMMAL']],   [['LLEPI', 'LMEPI']]   ],
    'patella_l': [   [['LLEPI', 'LMEPI']],   [['LLEPI', 'LMEPI']],   [['LLEPI', 'LMEPI']]   ],
    'talus_l': [   [['LMMAL', 'LLMAL']],   [['LMMAL', 'LLMAL']],   [['LMMAL', 'LLMAL']],   ],
    'calcn_l': [   [['LHEEL', 'LLMAL']],   [['LHEEL', 'LLMAL']],   [['LHEEL', 'LLMAL']],   ],
    'toes_l': [   [['LHEEL', 'LTOE']],   [['LHEEL', 'LTOE']],   [['LP1MT', 'LP5MT']],   ],
    'torso': [   [['MAN', 'C7']],   [['RSH', 'RASI'],['LSH', 'LASI']],   [['RSH', 'LSH']],   ],
    'humerus_r': [   [['RSH', 'RELB']],   [['RSH', 'RELB']],   [['RSH', 'RELB']],   ],
    'ulna_r': [   [['RELB', 'RWR']],   [['RELB', 'RWR']],   [['RELB', 'RWR']],   ],
    'radius_r': [   [['RELB', 'RWR']],   [['RELB', 'RWR']],   [['RELB', 'RWR']],   ],
    'humerus_l': [   [['LSH', 'LELB']],   [['LSH', 'LELB']],   [['LSH', 'LELB']],   ],
    'ulna_l': [   [['LELB', 'LWR']],   [['LELB', 'LWR']],   [['LELB', 'LWR']],   ],
    'radius_l': [   [['LELB', 'LWR']],   [['LELB', 'LWR']],   [['LELB', 'LWR']],   ],
    }

#Set model path
pathModelFolder = os.getcwd()+'\\data\\Model\\'

#Set up logging.
logPath = os.path.join(pathModelFolder,'scaling.log')
if os.path.exists(logPath):
    os.remove(logPath)
    
#Remove all handlers associated with the root logger object.
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
logging.shutdown()
logging.basicConfig(filename = logPath,format='%(message)s',
                    level = logging.INFO)
osim.Logger.setLevelString('error')

#Load in the two models
laiUhlrich = osim.Model(pathModelFolder+'LaiUhlrich2022.osim')
laiUhlrichState = laiUhlrich.initSystem()
jaScaled = osim.Model(pathModelFolder+'JA1_SCALED_Osim40_Muscles.osim')
jaScaledState = jaScaled.initSystem()

#Calculate the scaling factors between marker pairs
scaleFactors = {}

#Loop through bodies
for body in scaleMarkerPairs.keys():
    
    #Add empty array in scale factors dictionary
    scaleFactors[body] = np.zeros((3,))
    
    #Loop through three axes
    for axInd in range(3):
        
        #Calculate distances for each model for the current body and axis
        #JA model
        jaModelDist = np.array([np.linalg.norm(jaScaled.updMarkerSet().get(scaleMarkerPairs[body][axInd][markerPairInd][0]).getLocationInGround(jaScaledState).to_numpy() - \
                                               jaScaled.updMarkerSet().get(scaleMarkerPairs[body][axInd][markerPairInd][1]).getLocationInGround(jaScaledState).to_numpy()) for \
                                markerPairInd in range(len(scaleMarkerPairs[body][axInd]))])
        #LaiUhlrich model
        laiModelDist = np.array([np.linalg.norm(laiUhlrich.updMarkerSet().get(scaleMarkerPairs[body][axInd][markerPairInd][0]).getLocationInGround(laiUhlrichState).to_numpy() - \
                                                laiUhlrich.updMarkerSet().get(scaleMarkerPairs[body][axInd][markerPairInd][1]).getLocationInGround(laiUhlrichState).to_numpy()) for \
                                 markerPairInd in range(len(scaleMarkerPairs[body][axInd]))])
            
        #Calculate scale factor average and append to scaling dictionary
        scaleFactors[body][axInd] = np.mean(jaModelDist / laiModelDist)
        
#Create scale tool
scaleTool = osim.ScaleTool()

#Set model mass
jaMass = np.sum([jaScaled.updBodySet().get(bodyInd).getMass() for bodyInd in range(jaScaled.updBodySet().getSize())])
scaleTool.setSubjectMass(jaMass)

#Set unscaled model to LaiUlrich
scaleTool.getGenericModelMaker().setModelFileName('data\Model\\LaiUhlrich2022.osim')

#Set scaling order to use manual measurements
scalingOrder = osim.ArrayStr()
scalingOrder.append('manualScale')
scaleTool.getModelScaler().setScalingOrder(scalingOrder)

#Set measurement scales in scale set
for body in scaleFactors.keys():
    
    #Create scale
    scale = osim.Scale()
    
    #Set parameters
    scale.setScaleFactors(osim.Vec3(scaleFactors[body]))
    scale.setSegmentName(body)
    
    #Append to scale set in tool
    scaleTool.getModelScaler().getScaleSet().cloneAndAppend(scale)

#Set output model file
scaleTool.getModelScaler().setOutputModelFileName('data\Model\\LaiUhlrich2022_JA1_SCALED.osim')

#Set marker placer to false
scaleTool.getMarkerPlacer().setApply(False)

#Run scale tool
scaleTool.run()

#Read in model to make some updates
scaledModel = osim.Model('data\Model\\LaiUhlrich2022_JA1_SCALED.osim')

#Delete markerset as this isn't needed
scaledModel.updMarkerSet().clearAndDestroy()

#Update name
scaledModel.setName('LaiUhlrich2022_JA1_SCALED')

#Re-print to file
scaledModel.printToXML('data\Model\\LaiUhlrich2022_JA1_SCALED.osim')

# %% Set data

#Adjust kinematic datafile to invert knee angle to match with new model
#Convert to joint angle headers .mot file too

#Read in data
coordinatesData = osim.TimeSeriesTable(os.path.join('data', 'Kinematics', trial_name+'.sto'))
    
#Get column headers
colHeaders = coordinatesData.getColumnLabels()

#Rename to remove the joint coordinate state details
newColHeaders = []
for col in colHeaders:
    if col.endswith('/value'):
        #Split the string to the 3rd component which is the joint angle
        newColHeaders.append(col.split('/')[3])
    else:
        newColHeaders.append(col)
        
#Set new column headers in coordinates data
coordinatesData.setColumnLabels(newColHeaders)

#Invert knee angle data
knee_angle_r_inverted = coordinatesData.getDependentColumn('knee_angle_r').to_numpy() * -1
knee_angle_l_inverted = coordinatesData.getDependentColumn('knee_angle_l').to_numpy() * -1

#Set in time series table by extracting vector and then row vector
for ii in range(len(coordinatesData.getIndependentColumn())):
    #Right knee angle at index
    coordinatesData.getDependentColumn('knee_angle_r').row(ii).setTo(knee_angle_r_inverted[ii])
    #Left knee angle at index
    coordinatesData.getDependentColumn('knee_angle_l').row(ii).setTo(knee_angle_l_inverted[ii])
    
#Write to .mot file
osim.STOFileAdapter().write(coordinatesData, os.path.join('data', 'Kinematics', trial_name+'.mot'))

# %% Process inputs for OpenSim AD

"""

The following section of code calls the processInputsOpenSimAD function which
is included the OpenSimAD utilities.
This function:
    > Adjusts wrapping surfaces
    > Adds foot ground contacts
    > Generates external functions for OpenSim AD

"""

#Create data processing and session folder

#Global processing folder
os.makedirs('processed', exist_ok = True)

#Session folder within processed directory
os.makedirs('processed\\testSession', exist_ok = True)

#Copy any model files to a model directory within session
os.makedirs('processed\\testSession\\OpenSimData\\Model', exist_ok = True)
for osimFile in glob.glob('data\\Model\\*.osim'):
    shutil.copy2(osimFile, osimFile.replace('data\\Model\\',
                                            'processed\\testSession\\OpenSimData\\Model\\'))
    
#Copy trial motion files across to session directory
os.makedirs('processed\\testSession\\OpenSimData\\Kinematics', exist_ok = True)
shutil.copy2('data\\Kinematics\\' + trial_name + '.mot',
             'processed\\testSession\\OpenSimData\\Kinematics\\' + trial_name + '.mot')

#Process inputs for OpenSim AD
settings = processInputsOpenSimAD(baseDir, 'processed', 'testSession',
                                  trial_name, motion_type,
                                  OpenSimModel = 'LaiUhlrich2022', downloadData = False,
                                  time_window = time_window,
                                  treadmill_speed = treadmill_speed,
                                  massKg = jaMass, heightM = 1.75,  ### TODO: get accurate height...
                                  unscaledModelFile = 'LaiUhlrich2022.osim',
                                  scaledModelFile = 'LaiUhlrich2022_JA1_SCALED.osim')

# Adjust settings for this example.
# Set the model to be torque-driven.
settings['torque_driven_model'] = True

# Adjust the weights of the objective function and remove the default
# muscle-related weigths. The objective function contains terms for tracking
# coordinate values (positionTrackingTerm), speeds (velocityTrackingTerm), and
# accelerations (accelerationTrackingTerm), as well as terms for minimizing
# excitations of the ideal torque actuators at the arms (armExcitationTerm),
# lumbar (lumbarExcitationTerm), and lower-extremity (coordinateExcitationTerm)
# joints. The objective function also contains a regularization term that
# minimizes the coordinate accelerations (jointAccelerationTerm).
settings['weights'] = {
    'positionTrackingTerm': 10,
    'velocityTrackingTerm': 10,
    'accelerationTrackingTerm': 50,
    'armExcitationTerm': 0.001,
    'lumbarExcitationTerm': 0.001,
    'coordinateExcitationTerm': 1,
    'jointAccelerationTerm': 0.001,}

# Add periodic constraints to the problem. This will constrain initial and
# final states of the problem to be the same. This is useful for obtaining
# faster convergence. Please note that the walking trial we selected might not
# be perfectly periodic. We here add periodic constraints to show how to do it.
# We here add periodic constraints for the coordinate values (coordinateValues)
# and coordinate speeds (coordinateSpeeds) of the lower-extremity joints
# (lowerLimbJoints). We also add periodic constraints for the activations of the
# ideal torque actuators at the lower-extremity (lowerLimbJointActivations) and
# lumbar (lumbarJointActivations) joints. 
settings['periodicConstraints'] = {
    'coordinateValues': ['lowerLimbJoints'],
    'coordinateSpeeds': ['lowerLimbJoints'],
    'lowerLimbJointActivations': ['all'],
    'lumbarJointActivations': ['all']}

# Filter the data to be tracked.
settings['filter_Qs_toTrack'] = True
settings['cutoff_freq_Qs'] = 12
settings['filter_Qds_toTrack'] = True
settings['cutoff_freq_Qds'] = 12
settings['filter_Qdds_toTrack'] = True
settings['cutoff_freq_Qdds'] = 12
settings['splineQds'] = True

# We set the mesh density to 50. We recommend using a mesh density of 100 by
# default, but we here use a lower value to reduce the computation time.
settings['meshDensity'] = 50

# %% Run the dynamic simulation

### TODO:
    ##### > Consider other appropriate settings to include...
    ##### Seems to work - but problem may be infeasible for some reason?
    
    ##### Checks?
        
        ##### Taking out acceleration terms - which were contributing highly to 
        ##### cost function - didn't seem to work...
        
        ##### Infeasibility is generally due to something being bounded too strictly.
        ##### The infesible results point solution is the model floating and the arms
        ##### flying around - it is possible that something to do with the settings
        ##### means the model can't achieve the desired positions and velocities...
            
        #### Turn off the yCalcnToes constraint?
            ##### Still seems to end up with infeasible solution
        
        
        #### Joint velocity bounds?

#Run the torque tracking problem
run_tracking(baseDir, 'processed', 'testSession', settings, case = case,
             model_full_name = 'LaiUhlrich2022_JA1_SCALED_adjusted_contacts')





















# %%% ----- End of sprintSim_OpenSimAD.py ----- %% #