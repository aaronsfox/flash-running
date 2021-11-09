# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: 
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    TODO:
        
        - ADD runOpts DICTIONARY TO IDENTIFY WHETHER TO RUN CERTAIN STEPS
        - ADD DESCRIPTION OF STEPS (SEE BELOW)
        
        - Model creation = adapt the 3D model to convert to 2D for later simulation use
        - Model creation = scale Ong et al. model based on static trial to map muscles [models folder]
        - Run IK with 3D model [IK folder]
        - Run RRA with 3D model to get better representative kinematics [RRA folder]
        - Create 2D version of GRF and kinematic data (states version) [data folder]
        - Run torque tracking simulation to better match kinematics w/ GRF contact [simulations folder]
        - Run muscle driven tracking simulation to get muscle function guess [simulations folder]
    
"""

# %% Import packages

import opensim

# %%% Define functions

#Get full gait cycle for right leg
def getFullGaitCycleRight(grfFile = None):
    
    # Convenience function for getting the first full gait cycle based on GRF
    # data. Note that this function is only applicable to the way the current
    # data is structured, but could easily be edited (e.g. for getting a full
    # gait cycle, or right to right foot strike etc.)
    #
    # This function is currently only applicable to going from right contact
    # to another right contact. This could easily be edited with an input 
    # variable asking for limb.
    #
    # Input:    grfFile - .mot file containing GRF time history
    
    #Check for input
    if grfFile is None:
        raise ValueError('A GRF file in .mot format is required')
        
    #Load the GRF data
    grfTable = osim.TimeSeriesTable(grfFile)
    
    #Convert to more easily readable object for easier manipulation
    #Get number of column labels and data rows
    nLabels = grfTable.getNumColumns()
    nRows = grfTable.getNumRows()
    #Pre-allocate numpy array based on data size
    grfArray = np.zeros((nRows,nLabels))
    #Loop through labels and rows and get data
    for iLabel in range(0,nLabels):
        #Get current column label
        currLabel = grfTable.getColumnLabel(iLabel)
        #Store column index value if one of the vertical GRF data
        if currLabel == 'ground_force_r_vy':
            rightGrfCol = int(iLabel)
        elif  currLabel == 'ground_force_l_vy':
            leftGrfCol = int(iLabel)
        for iRow in range(0,nRows):
            grfArray[iRow,iLabel] = grfTable.getDependentColumn(currLabel).getElt(0,iRow)
    
    #Create a logical for where the right & left foot is in contact with the ground
    #based on 10N threshold
    #Identify columns for right and left vertical GRF
    rightGrfOn = grfArray[:,rightGrfCol] > 10
    leftGrfOn = grfArray[:,leftGrfCol] > 10
    
    #Identify the index where right and left GRF starts
    #Identify where change in boolean value is present
    rightGrfDiff = np.where(rightGrfOn[:-1] != rightGrfOn[1:])[0]
    leftGrfDiff = np.where(leftGrfOn[:-1] != leftGrfOn[1:])[0]
    #Identify where the index one above the change is true for foot strikes
    rightGrfOnInd = list()
    for gg in range(0,len(rightGrfDiff)):
        if rightGrfOn[rightGrfDiff[gg]+1]:
            rightGrfOnInd.append(rightGrfDiff[gg]+1)
    leftGrfOnInd = list()
    for gg in range(0,len(leftGrfDiff)):
        if leftGrfOn[leftGrfDiff[gg]+1]:
            leftGrfOnInd.append(leftGrfDiff[gg]+1)
    
    #Find the first right strike and then get the next. With sprint trials there are generally
    #only two within an overground trial
    rightStartInd = rightGrfOnInd[0] - 1
    rightStopInd = rightGrfOnInd[1] - 1
    
    #Identify the times corresponding to these foot strikes
    startTime = list(grfTable.getIndependentColumn())[rightStartInd]
    endTime = list(grfTable.getIndependentColumn())[rightStopInd]
    
    #Print outputs
    print('Start Time: '+str(startTime))
    print('End Time: '+str(endTime))
    
    return startTime,endTime

# %%% Set-up

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

# %% Create Models

#This step adapts the base 3D model to create a 2D model that includes lower limb
#muscles from the Ong et al. model included in the 'models' directory

#Load the base 3D model
base3DModel = osim.Model('models\\base3DModel.osim')

# %% Create 2D version of model

#Set a copy of this model to edit
model2D = osim.Model('models\\base3DModel.osim')

#Convert hip joints to 2D

#Get hip joints from scaled model to work from
origHip_r = base3DModel.getJointSet().get('hip_r')
origHip_l = base3DModel.getJointSet().get('hip_l')
    
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

#Update ankle to pin joint

#Get hip joints from scaled model to work from
origAnkle_r = base3DModel.getJointSet().get('ankle_r')
origAnkle_l = base3DModel.getJointSet().get('ankle_l')
    
#Create pin joint for right ankle
#Set parent and child locations
locationInParent = origAnkle_r.get_frames(0).get_translation()
orientationInParent = origAnkle_r.get_frames(0).get_orientation()
locationInChild = origAnkle_r.get_frames(1).get_translation()
orientationInChild = origAnkle_r.get_frames(1).get_orientation()
#Get bodies for joint
tibiaBody = model2D.getBodySet().get('tibia_r')
talusBody = model2D.getBodySet().get('talus_r')
#Create joint
ankleJoint_r = osim.PinJoint('ankle_r', tibiaBody, locationInParent,
                            orientationInParent, talusBody, locationInChild,
                            orientationInChild)
#Update additional joint parameters
ankleJoint_r.getCoordinate().setName('ankle_angle_r') #set coordinate name
ankleJoint_r.getCoordinate().setRangeMin(origHip_r.get_coordinates(0).getRangeMin()) #get min from old joint
ankleJoint_r.getCoordinate().setRangeMax(origHip_r.get_coordinates(0).getRangeMax()) #get max from old joint
ankleJoint_r.getCoordinate().setDefaultValue(0) #set default value to zero
ankleJoint_r.getCoordinate().setDefaultSpeedValue(0) #set default speed to zero
ankleJoint_r.getCoordinate().set_clamped(True) #set joint to clamped
ankleJoint_r.getCoordinate().set_locked(False) #set locked to false
ankleJoint_r.getCoordinate().set_prescribed(False) #set prescribed to false

#Create pin joint for left ankle
#Set parent and child locations
locationInParent = origAnkle_l.get_frames(0).get_translation()
orientationInParent = origAnkle_l.get_frames(0).get_orientation()
locationInChild = origAnkle_l.get_frames(1).get_translation()
orientationInChild = origAnkle_l.get_frames(1).get_orientation()
#Get bodies for joint
tibiaBody = model2D.getBodySet().get('tibia_l')
talusBody = model2D.getBodySet().get('talus_l')
#Create joint
ankleJoint_l = osim.PinJoint('ankle_l', tibiaBody, locationInParent,
                            orientationInParent, talusBody, locationInChild,
                            orientationInChild)
#Update additional joint parameters
ankleJoint_l.getCoordinate().setName('ankle_angle_l') #set coordinate name
ankleJoint_l.getCoordinate().setRangeMin(origHip_r.get_coordinates(0).getRangeMin()) #get min from old joint
ankleJoint_l.getCoordinate().setRangeMax(origHip_r.get_coordinates(0).getRangeMax()) #get max from old joint
ankleJoint_l.getCoordinate().setDefaultValue(0) #set default value to zero
ankleJoint_l.getCoordinate().setDefaultSpeedValue(0) #set default speed to zero
ankleJoint_l.getCoordinate().set_clamped(True) #set joint to clamped
ankleJoint_l.getCoordinate().set_locked(False) #set locked to false
ankleJoint_l.getCoordinate().set_prescribed(False) #set prescribed to false

#Remove existing hip joints from model
model2D.getJointSet().remove(model2D.getJointSet().get('ankle_r'))
model2D.getJointSet().remove(model2D.getJointSet().get('ankle_l'))

#Add new hip joints
model2D.addJoint(ankleJoint_r)
model2D.addJoint(ankleJoint_l)

#Convert lumbar joint to 2D

#Get back joint from scaled model to work from
origLumbar = base3DModel.getJointSet().get('back')
    
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

#Update pelvis model to 2D

#Get original pelvis from scaled model to work from
origPelvis = base3DModel.getJointSet().get('ground_pelvis')

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

#Convert ankle and foot joints to planar variants

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

#Convert arm joints to 2D

#Get acromial joints from scaled model to work from
origArm_r = base3DModel.getJointSet().get('acromial_r')
origArm_l = base3DModel.getJointSet().get('acromial_l')
    
#Create pin joint for right arm
#Set parent and child locations
locationInParent = origArm_r.get_frames(0).get_translation()
orientationInParent = origArm_r.get_frames(0).get_orientation()
locationInChild = origArm_r.get_frames(1).get_translation()
orientationInChild = origArm_r.get_frames(1).get_orientation()
#Get bodies for joint
torsoBody = model2D.getBodySet().get('torso')
humerusBody = model2D.getBodySet().get('humerus_r')
#Create joint
armJoint_r = osim.PinJoint('acromial_r', torsoBody, locationInParent,
                            orientationInParent, humerusBody, locationInChild,
                            orientationInChild)
#Update additional joint parameters
armJoint_r.getCoordinate().setName('arm_flex_r') #set coordinate name
armJoint_r.getCoordinate().setRangeMin(origArm_r.get_coordinates(0).getRangeMin()) #get min from old joint
armJoint_r.getCoordinate().setRangeMax(origArm_r.get_coordinates(0).getRangeMax()) #get max from old joint
armJoint_r.getCoordinate().setDefaultValue(0) #set default value to zero
armJoint_r.getCoordinate().setDefaultSpeedValue(0) #set default speed to zero
armJoint_r.getCoordinate().set_clamped(True) #set joint to clamped
armJoint_r.getCoordinate().set_locked(False) #set locked to false
armJoint_r.getCoordinate().set_prescribed(False) #set prescribed to false

#Create pin joint for left arm
#Set parent and child locations
locationInParent = origArm_l.get_frames(0).get_translation()
orientationInParent = origArm_l.get_frames(0).get_orientation()
locationInChild = origArm_l.get_frames(1).get_translation()
orientationInChild = origArm_l.get_frames(1).get_orientation()
#Get bodies for joint
torsoBody = model2D.getBodySet().get('torso')
humerusBody = model2D.getBodySet().get('humerus_l')
#Create joint
armJoint_l = osim.PinJoint('acromial_l', torsoBody, locationInParent,
                            orientationInParent, humerusBody, locationInChild,
                            orientationInChild)
#Update additional joint parameters
armJoint_l.getCoordinate().setName('arm_flex_l') #set coordinate name
armJoint_l.getCoordinate().setRangeMin(origArm_r.get_coordinates(0).getRangeMin()) #get min from old joint
armJoint_l.getCoordinate().setRangeMax(origArm_r.get_coordinates(0).getRangeMax()) #get max from old joint
armJoint_l.getCoordinate().setDefaultValue(0) #set default value to zero
armJoint_l.getCoordinate().setDefaultSpeedValue(0) #set default speed to zero
armJoint_l.getCoordinate().set_clamped(True) #set joint to clamped
armJoint_l.getCoordinate().set_locked(False) #set locked to false
armJoint_l.getCoordinate().set_prescribed(False) #set prescribed to false

#Remove existing hip joints from model
model2D.getJointSet().remove(model2D.getJointSet().get('acromial_r'))
model2D.getJointSet().remove(model2D.getJointSet().get('acromial_l'))

#Add new hip joints
model2D.addJoint(armJoint_r)
model2D.addJoint(armJoint_l)

#Update radioulnar orientation to better align with 2D model
#Supinate for level forearms
#This still results in the hand passing through the legs, but is just a visualisation issue
model2D.getJointSet().get('radioulnar_r').get_frames(1).set_orientation(osim.Vec3(0,math.radians(-90),0))
model2D.getJointSet().get('radioulnar_l').get_frames(1).set_orientation(osim.Vec3(0,math.radians(90),0))

#Finalize model connections
model2D.finalizeConnections()

#Dump into a model processor to weld the foot and hand joints
#Create the processor
editProcessor = osim.ModelProcessor(model2D)
#Create the string array of joints to weld
weldList = ['subtalar_r', 'subtalar_l', 'mtp_r', 'mtp_l',
            'radioulnar_r', 'radioulnar_l', 'radius_hand_r', 'radius_hand_l']
weldJoints = osim.StdVectorString()
for joint in weldList:
    weldJoints.append(joint)
#Append model operator for welding
editProcessor.append(osim.ModOpReplaceJointsWithWelds(weldJoints))
#Process model output
final2DModel = editProcessor.process()

#Remove the marker set as it's no use in this 2D context
final2DModel.updMarkerSet().clearAndDestroy()

#Reset name
final2DModel.setName('base2DModel')

#Print 2D model output
final2DModel.printToXML('models\\base2DModel.osim')

# %% Append muscle set to 2D model

#### TODO: scaling procedures

# %% Inverse Kinematics

#This step uses the .trc data from the sprint trial combined with the base 3D
#model to extract the kinematics for the trial

#Load the model to use for IK
ikModel = osim.Model('models\\base3DModel.osim')

#Initialise IK tool
ikTool = osim.InverseKinematicsTool()
ikTool.setName('sprint')

#Set model
ikTool.setModel(ikModel)

#Set marker file
ikTool.set_marker_file('data\\sprint.trc')

#Create and set task set based off .trc markers

#Set-up a dictionary with various marker task weights
#Note that this is somewhat generic and could be edited
ikTaskWeights = {'C7': 10, 'RSH': 10, 'LSH': 10, 'MAN': 10, 'T7': 10, 'LARM': 2.5,
                 'LELB': 5, 'LFOREARM': 10, 'LWR': 10, 'RARM': 2.5, 'RELB': 5, 'RFOREARM': 10,
                 'RWR': 10, 'RASI': 5, 'LASI': 5, 'SACR': 10, 'P1': 7.5, 'P2': 7.5,
                 'P3': 7.5, 'LTHLP': 5, 'LTHAP': 5, 'LTHAD': 5, 'LTHLD': 5, 'LLEPI': 10, 'LPAT': 10,
                 'LTIAP': 5, 'LTIAD': 5, 'LTILAT': 5, 'LLMAL': 10, 'LHEEL': 5, 'LMFS': 2.5,
                 'LMFL': 2.5, 'LP1MT': 5, 'LTOE': 2.5, 'LP5MT': 5, 'RTHLD': 5, 'RTHAP': 5,
                 'RTHLP': 5, 'RTHAD': 5, 'RLEPI': 10, 'RPAT': 10, 'RTIAP': 5, 'RTIAD': 5, 'RTILAT': 5,
                 'RLMAL': 10, 'RHEEL': 5, 'RMFS': 2.5, 'RMFL': 2.5, 'RP1MT': 5, 'RTOE': 2.5,
                 'RP5MT': 5}
#Create task set
ikTaskSet = osim.IKTaskSet()
#Loop through markers and append to task set
for marker in list(ikTaskWeights.keys()):
    #Create blank marker task
    markerTask = osim.IKMarkerTask()
    #Set marker name
    markerTask.setName(marker)
    #Set marker weight
    markerTask.setWeight(ikTaskWeights[marker])
    #Clone and append to task set
    ikTaskSet.cloneAndAppend(markerTask)
    
#Set task set
ikTool.set_IKTaskSet(ikTaskSet)

#Set times
#Identify right foot strike to right foot strike for analysis from GRF data
#Use the pre-defined function to do so
[startTime,endTime] = getFullGaitCycleRight('data\\sprint_grf.mot')

#Set times
ikTool.setStartTime(startTime)
ikTool.setEndTime(endTime)

#Set output filename
#Create folder if it isn't present
if not os.path.isdir('ik'):
    os.mkdir('ik')
#Set the output file
ikTool.set_output_motion_file('ik\\sprint_ik.mot')

#Print to file
ikTool.printToXML('ik\\sprint_ik_setup.xml')

#Run IK
ikTool.run()




# %% Test ID for opensim working...


#Initialise ID tool
idTool = osim.InverseDynamicsTool()
idTool.setName('sprint')

#Set model
idTool.setModel(osim.Model('models\\base3DModel.osim'))

#Set external loads file
idTool.setExternalLoadsFileName('data\\sprint_grf.xml')

#Set coordinates file
idTool.setCoordinatesFileName('ik\\sprint_ik.mot')

#Set filter of 12Hz for coordinates (matches some sprinting studies)
idTool.setLowpassCutoffFrequency(12)

#Set times
#Start and end of IK
idTool.setStartTime(osim.Storage('ik\\sprint_ik.mot').getFirstTime())
idTool.setEndTime(osim.Storage('ik\\sprint_ik.mot').getLastTime())

#Set output file
idTool.setOutputGenForceFileName('sprint_id.sto')

#Run ID
idTool.run()

# %%% ----- End of X.py -----