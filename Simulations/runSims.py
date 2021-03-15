# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 15:23:54 2021

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    TODO: add notes...tracking sim of experimental data first...
    
"""

# %% Import packages

import os
import opensim as osim
import numpy as np
import pandas as pd

# %% Define functions

#readSTO
def readSTO(fileName):
    
    #Read storage
    sto = osim.Storage(fileName)
        
    #Get column list
    labels = osimArrayToList(sto.getColumnLabels())

    #Get time
    time = osim.ArrayDouble()
    sto.getTimeColumn(time)  
    time = osimArrayToList(time)
    
    #Get data
    data = []
    for i in range(sto.getSize()):
        temp = osimArrayToList(sto.getStateVector(i).getData())
        temp.insert(0, time[i])
        data.append(temp)

    df = pd.DataFrame(data, columns=labels)
    df.index = df.time
    
    return df

#osimArrayToList
def osimArrayToList(array):
    temp = []
    for i in range(array.getSize()):
        temp.append(array.get(i))

    return temp

#addReserve
def addReserve(model, coord, optimal_force, max_control):
    actu = osim.ActivationCoordinateActuator()
    actu.set_coordinate(coord)
    if coord.startswith('pelvis'):
        prefix = 'residual_'
    elif (coord.startswith('lumbar') or
          coord.startswith('arm') or
          coord.startswith('elbow') or
          coord.startswith('pro_sup')):
        prefix = 'torque_'
    else:
        prefix = 'reserve_'
    actu.setName(prefix + coord)
    actu.setOptimalForce(optimal_force)
    actu.setMinControl(-max_control)
    actu.setMaxControl(max_control)
    model.addForce(actu)

# %% Notes

# Tendon compliance in a couple of muscles seems problematic
# Inverse solution was troubled and had NaN's with tendon compliance and no implicit
# dynamics --- and also didn't look like converging with both of these on
# Taking these out of inverse solution started producing better iteration results
# It will probably converge to a solution
# Could this be due to the tendon slack length (or associted tendon) parameters
# not being ideal for a sprinting motion with such higher degrees of motion?

# %% 2D basic tracking sim of experimental data (torque driven)

##### TODO: could package each sim into function?

# This step is designed to generate relevant controls for the 2D model that
# appropriately track the sprint kinematics and GRF from the experimental data
# using torque actuators

#Set simulation parameters

#Timing and mesh data
initialTime = osim.Storage('refQ_2D.sto').getFirstTime()
finalTime = osim.Storage('refQ_2D.sto').getLastTime()
duration = finalTime - initialTime
meshNo = 50 #set as per mesh number for half gait cycle in Falisse et al. 2019, J R Soc Interface
meshInterval = duration / meshNo

#Weights
trackingWeight = 10
effortWeight = 0.1
grfWeight = 1
speedWeight = 1
symmetryWeight = 1

#Create the model processor for the tracking problem

#Load the model and create processor
osimModel = osim.Model('gaitModel2D_doubleStrengthVel.osim')
# osimModel = osim.Model('gaitModel2D_complex_doubleStrength.osim')
modelProcessor = osim.ModelProcessor(osimModel)

#Add processing step to remove muscles
modelProcessor.append(osim.ModOpRemoveMuscles())

#Process model for additional parameters
torqueModel = modelProcessor.process()
torqueModel.initSystem()

#Add torque actuators to support model
#Set list to add reserves to
reserveList = ['pelvis_tilt',
                'hip_flexion_l', 'hip_flexion_r',
                'knee_angle_l','knee_angle_r',
                'ankle_angle_l', 'ankle_angle_r']
#Set optimal force and max torque
optimalForce = 100
maxTorque = np.inf
#Add reserves
for tt in range(len(reserveList)):
    addReserve(torqueModel, reserveList[tt], optimalForce, maxTorque)

#Finalise connections
torqueModel.finalizeConnections()

#Get the track model as a processor object
torqueModelProcessor = osim.ModelProcessor(torqueModel)

#Set-up tracking problem

#Process current model
trackModel = torqueModelProcessor.process()
trackModel.initSystem()
      
#Construct the tracking object and set basic parameters
track = osim.MocoTrack()
track.setName('sprintTracking_torqueDriven')
track.setModel(torqueModelProcessor)

#Set kinematic data and parameters
tableProcessor = osim.TableProcessor('refQ_2D.sto')
# tableProcessor.append(osim.TabOpLowPassFilter(12)) ### data already filtered via RRA
tableProcessor.append(osim.TabOpUseAbsoluteStateNames())
track.setStatesReference(tableProcessor)
track.set_states_global_tracking_weight(trackingWeight)
#Set some coordinates not to be tracked to avoid poor motion
#Similar process done in Dembia et al. Moco paper
stateWeights = osim.MocoWeightSet()
weightList = list()
weightList.append(('/jointset/ground_pelvis/pelvis_ty', 0))
weightList.append(('/jointset/ground_pelvis/pelvis_tilt', 0.1))
for weight in weightList:
    stateWeights.cloneAndAppend(osim.MocoWeight(weight[0] + '/value', weight[1]))
    stateWeights.cloneAndAppend(osim.MocoWeight(weight[0] + '/speed', weight[1]))        
track.set_states_weight_set(stateWeights)
#Set tracked states to be used in guess
track.set_apply_tracked_states_to_guess(True)
#Set unused state references to be allowed in case of file errors
track.set_allow_unused_references(True)
#Set position derivatives to be used as speeds
track.set_track_reference_position_derivatives(True)

#Set control weights
track.set_control_effort_weight(effortWeight)

#Set timing parameters
track.set_initial_time(initialTime)
track.set_final_time(finalTime)
track.set_mesh_interval(meshInterval)

#Customise the base tracking problem with relevant goals
study = track.initialize()
problem = study.updProblem()

#Adjust time bounds to allow for subtle fluctations in finish time
problem.setTimeBounds(initialTime, [finalTime - 0.1, finalTime + 0.1])

#Update the control effort goal to a cost of transport type cost
effort = osim.MocoControlGoal().safeDownCast(problem.updGoal('control_effort'))
effort.setDivideByDisplacement(True)

#Set to not highly penalise the lumbar actuator
effort.setWeightForControl('/forceset/lumbarAct', 0.001)
# effort.setWeightForControl('/forceset/tau_lumbar_ext', 0.001)

#Set an average speed goal based on sprinting data
df_kinematics = readSTO('refQ_2D.sto')
startPos = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'].to_numpy()[0]
endPos = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'].to_numpy()[-1]
sprintSpeed = (endPos - startPos) / (finalTime - initialTime)
speedGoal = osim.MocoAverageSpeedGoal('speed')
speedGoal.set_desired_average_speed(sprintSpeed)
speedGoal.setWeight(speedWeight)
problem.addGoal(speedGoal)

#Set symmetry constraint
#Create symmetry goal
symmetryGoal = osim.MocoPeriodicityGoal('symmetryGoal')
#Set symmetric coordinate values (except for pelvis_tx) and speeds
for cc in range(trackModel.getCoordinateSet().getSize()):
    #Get current coordinate name
    coordName = trackModel.getCoordinateSet().get(cc).getName()
    #Set symmetry parameters based on coordinate
    if coordName.endswith('_r'):
        #Joint angle value
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/value'))
        #Joint speed
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/speed'))
    elif coordName.endswith('_l'):
        #Joint angle value
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/value'))
        #Joint speed
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/speed'))
    elif coordName.endswith('_extension') or coordName.endswith('_tilt') or coordName.endswith('_ty'):
        #Joint angle value
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value'))
        #Joint speed
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed'))
#Add a symmetry pair for pelvis_tx speed
symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
    '/jointset/ground_pelvis/pelvis_tx/speed'))
#Add symmetry goal
symmetryGoal.setWeight(symmetryWeight)
problem.addGoal(symmetryGoal)

#Add contact tracking goal
#Set force names
forceNamesRightFoot = ['forceset/contactHeel_r',
                        'forceset/contactToe1_r',
                        'forceset/contactToe2_r']
forceNamesLeftFoot = ['forceset/contactHeel_l',
                      'forceset/contactToe1_l',
                      'forceset/contactToe2_l']
# forceNamesRightFoot = ['forceset/contactHeel_r',
#                        'forceset/contactMH1_r',
#                        'forceset/contactMH3_r',
#                        'forceset/contactMH5_r',
#                        'forceset/contactHallux_r',
#                        'forceset/contactOtherToes_r']
# forceNamesLeftFoot = ['forceset/contactHeel_l',
#                        'forceset/contactMH1_l',
#                        'forceset/contactMH3_l',
#                        'forceset/contactMH5_l',
#                        'forceset/contactHallux_l',
#                        'forceset/contactOtherToes_l']
#Create contact tracking goal
contactGoal = osim.MocoContactTrackingGoal('contactGoal', grfWeight)
#Set external loads
contactGoal.setExternalLoadsFile('refGRF_2D.xml')
#Set force name groups
forceNames_r = osim.StdVectorString()
forceNames_l = osim.StdVectorString()
for ff in range(len(forceNamesRightFoot)):
    forceNames_r.append(forceNamesRightFoot[ff])
    forceNames_l.append(forceNamesLeftFoot[ff])
#Create and add tracking groups
trackRightGRF = osim.MocoContactTrackingGoalGroup(forceNames_r, 'RightGRF')
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r');
contactGoal.addContactGroup(trackRightGRF)
trackLeftGRF = osim.MocoContactTrackingGoalGroup(forceNames_l, 'LeftGRF')
trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l')
contactGoal.addContactGroup(trackLeftGRF)
#Set parameters
contactGoal.setProjection('plane')
contactGoal.setProjectionVector(osim.Vec3(0, 0, 1))
#Add contact tracking goal
problem.addGoal(contactGoal)

#Add state bounds
#This also sets initial bounds on lumbar and pelvis coordinates that reflect
#'normal' running motions - based off the experimental data
problem.setStateInfo('/jointset/back/lumbar_extension/value', [np.deg2rad(-30),0],
                     [np.deg2rad(-10),0])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tilt/value', [-10*np.pi/180, 0*np.pi/180],
                     [np.deg2rad(-5),0])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tx/value', [0, 3])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_ty/value', [0.8, 1.25])
problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [-25*np.pi/180, 75*np.pi/180])
problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [-25*np.pi/180, 75*np.pi/180])
problem.setStateInfo('/jointset/knee_l/knee_angle_l/value', [-140*np.pi/180, 0])
problem.setStateInfo('/jointset/knee_r/knee_angle_r/value', [-140*np.pi/180, 0])
# problem.setStateInfo('/jointset/walker_knee_l/knee_angle_l/value', [0, 140*np.pi/180])
# problem.setStateInfo('/jointset/walker_knee_r/knee_angle_r/value', [0, 140*np.pi/180])
problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [-20*np.pi/180, 30*np.pi/180])
problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [-20*np.pi/180, 30*np.pi/180])

#Configure the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.resetProblem(problem)
solver.set_optim_constraint_tolerance(1e-2)
solver.set_optim_convergence_tolerance(1e-2)

#Solve
solution = study.solve()

#Option to visualise
# study.visualize(solution)

#Create a full gait cycle trajectory from the periodic solution.
addPatterns = [".*pelvis_tx/value"]
fullTraj = osim.createPeriodicTrajectory(solution, addPatterns)
fullTraj.write('sprintTracking_torqueDriven_solution_fullTrajectory.sto')

#Compute ground reaction forces generated by contact sphere from the 
#full gait cycle trajectory
externalLoads = osim.createExternalLoadsTableForGait(trackModel,
                                                     fullTraj,
                                                     forceNames_r,
                                                     forceNames_l)
osim.STOFileAdapter.write(externalLoads,'trackedGRF_torqueDriven_2D.mot')

##### torque solution is pretty good with these weights and symmetry goal included...

# %% 2D inverse sim of experimental  (nope...)

# This step obtains a solution that is potentially useful as an initial guess to
# the tracking problem

#Set simulation parameters

#Timing and mesh data
# initialTime = osim.Storage('refQ_2D.sto').getFirstTime()
# finalTime = osim.Storage('refQ_2D.sto').getLastTime()
initialTime = osim.Storage('sprintTracking_torqueDriven_solution.sto').getFirstTime()
finalTime = osim.Storage('sprintTracking_torqueDriven_solution.sto').getLastTime()
duration = finalTime - initialTime
meshNo = 200 #set as per mesh number for half gait cycle in Falisse et al. 2019, J R Soc Interface
meshInterval = duration / meshNo

#Model parameters
complexModel = False

#Muscle parameters
passiveForces = False
implicitTendonCompliance = False

#Create the model processor for the inverse problem

#### TODO: conditional for 2D vs. 3D

#Load the model and create processor
if complexModel:
    osimModel = osim.Model('gaitModel2D_complex_doubleStrengthVel.osim')
else:
    osimModel = osim.Model('gaitModel2D_doubleStrengthVel.osim')
# osimModel = osim.Model('gaitModel3D_doubleStrength.osim')

#Remove contact elements from model for inverse problem
forceSet = osimModel.getForceSet()
numContacts = 0
for ii in np.arange(forceSet.getSize()):
    forceName = forceSet.get(int(ii)).getName()
    if 'contact' in forceName:
        numContacts += 1
print('Removing contact force elements from model...')
contactsRemoved = 0
while contactsRemoved < numContacts:
    for ii in np.arange(forceSet.getSize()):
        forceName = forceSet.get(int(ii)).getName()
        if 'contact' in forceName: 
            print('  --> removed', forceSet.get(int(ii)).getName())
            forceSet.remove(int(ii))
            contactsRemoved += 1
            break
print('\n')

#Add reserve actuators to support model
#Set list to add reserves to
reserveList = ['pelvis_tilt', 'pelvis_tx', 'pelvis_ty',
                'hip_flexion_l', 'hip_flexion_r',
                'knee_angle_l','knee_angle_r',
                'ankle_angle_l', 'ankle_angle_r']
# reserveList = ['pelvis_tilt', 'pelvis_list', 'pelvis_rotation',
#                'pelvis_tx', 'pelvis_ty', 'pelvis_tz',
#                'hip_flexion_l', 'hip_adduction_l', 'hip_rotation_l',
#                'hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r',
#                'knee_angle_l','knee_angle_r',
#                'ankle_angle_l', 'ankle_angle_r']
#Set optimal force and max torque
optimalForce = 1
maxTorque = np.inf
#Add reserves
for tt in range(len(reserveList)):
    addReserve(osimModel, reserveList[tt], optimalForce, maxTorque)
    
#Add to model processor
modelProcessor = osim.ModelProcessor(osimModel)
    
#Disable tendon compliance (will be re-enabled for certain muscles later)
modelProcessor.append(osim.ModOpIgnoreTendonCompliance())

#Adjust fiber damping as per Dembia et al. Moco paper
modelProcessor.append(osim.ModOpFiberDampingDGF(0.01))

#Turn off passive forces if appropriate
if not passiveForces:
    modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())

#Add external loads to processor
modelProcessor.append(osim.ModOpAddExternalLoads('refGRF_2D.xml'))
# modelProcessor.append(osim.ModOpAddExternalLoads('refGRF.xml'))
    
#Process model for additional parameters
inverseModel = modelProcessor.process()
inverseModel.initSystem()

# #Enable tendon compliance in the gastrocnemius and soleus
# muscles = inverseModel.updMuscles()
# for mm in np.arange(muscles.getSize()):
#     currMusc = osim.DeGrooteFregly2016Muscle.safeDownCast(muscles.get(int(mm)))
#     muscName = currMusc.getName()
#     #Enable tendon compliance dynamics in the plantarflexors
#     if ('gastroc' in muscName) or ('soleus' in muscName):
#         currMusc.set_ignore_tendon_compliance(False)

# #Lock the subtalar and mtp joints for 3D model
# inverseModel.getCoordinateSet().get('subtalar_angle_r').set_locked(True)
# inverseModel.getCoordinateSet().get('subtalar_angle_l').set_locked(True)
# inverseModel.getCoordinateSet().get('mtp_angle_r').set_locked(True)
# inverseModel.getCoordinateSet().get('mtp_angle_l').set_locked(True)
# inverseModel.finalizeConnections()

#Get the inverse model as a processor object
inverseModelProcessor = osim.ModelProcessor(inverseModel)

#Set model to use implicit tendon compliance if appropriate
if implicitTendonCompliance:
    inverseModelProcessor.append(osim.ModOpUseImplicitTendonComplianceDynamicsDGF())

# #Replace welds for 3D model
# inverseModelProcessor.append(osim.ModOpReplaceJointsWithWelds())
#     # ['subtalar_angle_r', 'mtp_angle_r',
#     #  'subtalar_angle_l', 'mtp_angle_l']))

#Set up inverse problem
inverse = osim.MocoInverse()
inverse.setName('inverseSprintTracking')
inverse.setModel(inverseModelProcessor)

#Set kinematic data and parameters
# tableProcessor = osim.TableProcessor('refQ_2D.sto')
tableProcessor = osim.TableProcessor('sprintTracking_torqueDriven_solution.sto') # try using torque driven tracked kinematics?
# tableProcessor = osim.TableProcessor('refQ.sto')
# tableProcessor.append(osim.TabOpLowPassFilter(12)) ### data already filtered via RRA
tableProcessor.append(osim.TabOpUseAbsoluteStateNames())
inverse.setKinematics(tableProcessor)
inverse.set_kinematics_allow_extra_columns(True)
inverse.set_initial_time(initialTime)
inverse.set_final_time(finalTime)
inverse.set_mesh_interval(meshInterval)

#Solve inverse problem
inverseSolution = inverse.solve()
# inverseSolution.getMocoSolution().unseal()
inverseSolution.getMocoSolution().write('inverseSolution.sto')

#Remove default Moco solution file
os.remove('MocoStudy_solution.sto')

#### got bored and stopped inverse simulation with no tendon dynamics after 920
#### iterations --- it was borderline converging anyway...

#### inverse solution gets better with increased mesh number (i.e. 200) and
#### with the increased contraction velocity (probably needed for better muscle
#### functioning) --- stopped at approximately 550 iters
####
#### still some pretty high activations --- but these models may work better in
#### a tracking sim...

#### inverse with complex model and torque driven tracking kinematics would probably
#### converge eventually, but it had already done ~1500 iterations and therefore
#### isn't tracking to a good solution 
####
#### activations are 'shark finning' which might suggest too big of a mesh interval
####
#### contraction velocity and muscle parameters are perhaps an issue here too...

# %% 2D tracking sim of experimental data (muscle driven)

##### TODO: could package each sim into function?

# This step is designed to generate relevant controls for the 2D model that
# appropriately track the sprint kinematics and GRF from the experimental data

#Set simulation parameters

#Timing and mesh data
# initialTime = osim.Storage('refQ_2D.sto').getFirstTime()
# finalTime = osim.Storage('refQ_2D.sto').getLastTime()
initialTime = osim.Storage('sprintTracking_torqueDriven_solution.sto').getFirstTime()
finalTime = osim.Storage('sprintTracking_torqueDriven_solution.sto').getLastTime()
duration = finalTime - initialTime
meshNo = 100 #set as per mesh number for half gait cycle in Falisse et al. 2019, J R Soc Interface
meshInterval = duration / meshNo

#Model parameters
complexModel = False

#Muscle parameters
passiveForces = False
implicitTendonCompliance = True

### TODO: conditional for tendons

#Weights (more generic weights used than Moco Paper example...)
trackingWeight = 10
effortWeight = 0.1
grfWeight = 1
speedWeight = 1
symmetryWeight = 1

#Create the model processor for the tracking problem

#Load the model and create processor
if complexModel:
    osimModel = osim.Model('gaitModel2D_complex_doubleStrengthVel.osim')
else:
    osimModel = osim.Model('gaitModel2D_doubleStrengthVel.osim')    
modelProcessor = osim.ModelProcessor(osimModel)

# #Convert muscles to DeGrooteFregley variant
# modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
##### NOTE: don't need to do this as already done earlier
##### Is there still the wrapping object issue with this though???

#Disable tendon compliance (will be re-enabled for certain muscles later)
modelProcessor.append(osim.ModOpIgnoreTendonCompliance())

#Adjust fiber damping as per Dembia et al. Moco paper
modelProcessor.append(osim.ModOpFiberDampingDGF(0.01))

#Turn off passive forces if appropriate
if not passiveForces:
    modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
    
#Modify active force width
#This appears most important for "de-stiffening" muscles to allow appropriate movement
modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))

# #Scale isometric forces (even more)
# modelProcessor.append(osim.ModOpScaleMaxIsometricForce(10))
    
#Process model for additional parameters
trackModel = modelProcessor.process()
trackModel.initSystem()

#Enable tendon compliance in the gastrocnemius and soleus
muscles = trackModel.updMuscles()
for mm in np.arange(muscles.getSize()):
    currMusc = osim.DeGrooteFregly2016Muscle.safeDownCast(muscles.get(int(mm)))
    muscName = currMusc.getName()
    #Enable tendon compliance dynamics in the plantarflexors
    if ('gastroc' in muscName) or ('soleus' in muscName):
        currMusc.set_ignore_tendon_compliance(False)
        
#Get the track model as a processor object
trackModelProcessor = osim.ModelProcessor(trackModel)

#Set model to use implicit tendon compliance if appropriate
if implicitTendonCompliance:
    trackModelProcessor.append(osim.ModOpUseImplicitTendonComplianceDynamicsDGF())

#Set-up tracking problem

#Process current model
trackModel = trackModelProcessor.process()
trackModel.initSystem()

# #Count the number of Force objects in the model. This is used to normalise
# #control effort
# numForces = 0
# for ff in range(trackModel.getForceSet().getSize()):
#     if (trackModel.getForceSet().get(ff).getConcreteClassName().endswith('Muscle') or
#         trackModel.getForceSet().get(ff).getConcreteClassName().endswith('Actuator')):
#         numForces += 1
        
#Construct the tracking object and set basic parameters
track = osim.MocoTrack()
track.setName('sprintTracking_muscleDriven')
track.setModel(trackModelProcessor)

#Set kinematic data and parameters
# tableProcessor = osim.TableProcessor('refQ_2D.sto')
tableProcessor = osim.TableProcessor('sprintTracking_torqueDriven_solution.sto') # try using torque driven tracked kinematics?
# tableProcessor.append(osim.TabOpLowPassFilter(12)) ### data already filtered via RRA
tableProcessor.append(osim.TabOpUseAbsoluteStateNames())
track.setStatesReference(tableProcessor)
track.set_states_global_tracking_weight(trackingWeight)
# track.set_states_global_tracking_weight(trackingWeight / (2 * trackModel.getNumCoordinates()))
##### Remove this as already done in torque tracking???
# #Set some coordinates to have low tracking weights
# #Similar process done in Dembia et al. Moco paper
# stateWeights = osim.MocoWeightSet()
# weightList = list()
# weightList.append(('/jointset/ground_pelvis/pelvis_ty', 0.25))
# weightList.append(('/jointset/ground_pelvis/pelvis_tilt', 0.1))
# for weight in weightList:
#     stateWeights.cloneAndAppend(osim.MocoWeight(weight[0] + '/value', weight[1]))
#     stateWeights.cloneAndAppend(osim.MocoWeight(weight[0] + '/speed', weight[1]))        
# track.set_states_weight_set(stateWeights)
#Set tracked states to be used in guess
track.set_apply_tracked_states_to_guess(True)
#Set unused state references to be allowed in case of file errors
track.set_allow_unused_references(True)
#Set position derivatives to be used as speeds
track.set_track_reference_position_derivatives(True)

#Set control weights
track.set_control_effort_weight(effortWeight)
# track.set_control_effort_weight(effortWeight / numForces)

#Set timing parameters
track.set_initial_time(initialTime)
track.set_final_time(finalTime)
track.set_mesh_interval(meshInterval)

#Customise the base tracking problem with relevant goals
study = track.initialize()
problem = study.updProblem()

#Adjust time bounds to allow for subtle fluctations in finish time
problem.setTimeBounds(initialTime, [finalTime - 0.1, finalTime + 0.1])

#Update the control effort goal to a cost of transport type cost
effort = osim.MocoControlGoal().safeDownCast(problem.updGoal('control_effort'))
effort.setDivideByDisplacement(True)

#Set to not highly penalise the lumbar actuator
if complexModel:
    effort.setWeightForControl('/forceset/tau_lumbar_ext', 0.001)
else:
    effort.setWeightForControl('/forceset/lumbarAct', 0.001)

#Set an average speed goal based on sprinting data
df_kinematics = readSTO('refQ_2D.sto')
startPos = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'].to_numpy()[0]
endPos = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'].to_numpy()[-1]
sprintSpeed = (endPos - startPos) / (finalTime - initialTime)
speedGoal = osim.MocoAverageSpeedGoal('speed')
speedGoal.set_desired_average_speed(sprintSpeed)
speedGoal.setWeight(speedWeight)
problem.addGoal(speedGoal)

#Set symmetry constraint
#Create symmetry goal
symmetryGoal = osim.MocoPeriodicityGoal('symmetryGoal')
#Set symmetric coordinate values (except for pelvis_tx) and speeds
for cc in range(trackModel.getCoordinateSet().getSize()):
    #Get current coordinate name
    coordName = trackModel.getCoordinateSet().get(cc).getName()
    #Set symmetry parameters based on coordinate
    if coordName.endswith('_r'):
        #Joint angle value
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/value'))
        #Joint speed
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/speed'))
    elif coordName.endswith('_l'):
        #Joint angle value
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/value'))
        #Joint speed
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/speed'))
    elif coordName.endswith('_extension') or coordName.endswith('_tilt') or coordName.endswith('_ty'):
        #Joint angle value
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value'))
        #Joint speed
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            trackModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed'))
#Add a symmetry pair for pelvis_tx speed
symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
    '/jointset/ground_pelvis/pelvis_tx/speed'))
#Set symmetric activations
for ff in range(trackModel.getForceSet().getSize()):
    #Check for muscle or actuator
    if (trackModel.getForceSet().get(ff).getConcreteClassName().endswith('Muscle') or
        trackModel.getForceSet().get(ff).getConcreteClassName().endswith('Actuator')):
        #Get force name
        forceName = trackModel.getForceSet().get(ff).getName()
        #Set symmetry parameters based on coordinate
        if forceName.endswith('_r'):
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getForceSet().get(ff).getAbsolutePathString()+'/activation',
                trackModel.getForceSet().get(ff).getAbsolutePathString().replace('_r','_l')+'/activation'))
        elif forceName.endswith('_l'):
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                trackModel.getForceSet().get(ff).getAbsolutePathString()+'/activation',
                trackModel.getForceSet().get(ff).getAbsolutePathString().replace('_l','_r')+'/activation'))
        # elif forceName == 'lumbarAct':
        #     symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
        #         trackModel.getForceSet().get(ff).getAbsolutePathString()))
#Add symmetry goal
symmetryGoal.setWeight(symmetryWeight)
problem.addGoal(symmetryGoal)

#Add contact tracking goal
#Set force names
if complexModel:
    forceNamesRightFoot = ['forceset/contactHeel_r',
                           'forceset/contactMH1_r',
                           'forceset/contactMH3_r',
                           'forceset/contactMH5_r',
                           'forceset/contactHallux_r',
                           'forceset/contactOtherToes_r']
    forceNamesLeftFoot = ['forceset/contactHeel_l',
                           'forceset/contactMH1_l',
                           'forceset/contactMH3_l',
                           'forceset/contactMH5_l',
                           'forceset/contactHallux_l',
                           'forceset/contactOtherToes_l']
else:    
    forceNamesRightFoot = ['forceset/contactHeel_r',
                           'forceset/contactToe1_r',
                           'forceset/contactToe2_r']
    forceNamesLeftFoot = ['forceset/contactHeel_l',
                          'forceset/contactToe1_l',
                          'forceset/contactToe2_l']
#Create contact tracking goal
contactGoal = osim.MocoContactTrackingGoal('contactGoal', grfWeight)
#Set external loads
contactGoal.setExternalLoadsFile('refGRF_2D.xml')
#Set force name groups
forceNames_r = osim.StdVectorString()
forceNames_l = osim.StdVectorString()
for ff in range(len(forceNamesRightFoot)):
    forceNames_r.append(forceNamesRightFoot[ff])
    forceNames_l.append(forceNamesLeftFoot[ff])
#Create and add tracking groups
trackRightGRF = osim.MocoContactTrackingGoalGroup(forceNames_r, 'RightGRF')
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r');
contactGoal.addContactGroup(trackRightGRF)
trackLeftGRF = osim.MocoContactTrackingGoalGroup(forceNames_l, 'LeftGRF')
trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l')
contactGoal.addContactGroup(trackLeftGRF)
#Set parameters
contactGoal.setProjection('plane')
contactGoal.setProjectionVector(osim.Vec3(0, 0, 1))
#Add contact tracking goal
problem.addGoal(contactGoal)

#Add state bounds
#Set the initial values for the lumbar and pelvis coordinates that
#reflect 'normal' running motions - bounds based off experimental data
##### TODO: could be more particular about this about bounding initial values to a guess?
problem.setStateInfo('/jointset/back/lumbar_extension/value', [np.deg2rad(-30),0],
                     [np.deg2rad(-10),0])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tilt/value', [-10*np.pi/180, 0*np.pi/180],
                     [np.deg2rad(-5),0])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tx/value', [0, 3])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_ty/value', [0.8, 1.25],
                     [0.95, 1.05]) #likely reasonable starting value
problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [-25*np.pi/180, 75*np.pi/180])
problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [-25*np.pi/180, 75*np.pi/180])
if complexModel:
    problem.setStateInfo('/jointset/walker_knee_l/knee_angle_l/value', [0, 140*np.pi/180])
    problem.setStateInfo('/jointset/walker_knee_r/knee_angle_r/value', [0, 140*np.pi/180])
else:
    problem.setStateInfo('/jointset/knee_l/knee_angle_l/value', [-140*np.pi/180, 0])
    problem.setStateInfo('/jointset/knee_r/knee_angle_r/value', [-140*np.pi/180, 0])
problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [-20*np.pi/180, 30*np.pi/180])
problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [-20*np.pi/180, 30*np.pi/180])

#Set reasonable bounds for tendon force
problem.setStateInfoPattern('/forceset/.*/normalized_tendon_force', [0, 1.8], [], []);

#Configure the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_optim_constraint_tolerance(1e-2)
solver.set_optim_convergence_tolerance(1e-2)
# solver.set_multibody_dynamics_mode('implicit')
# solver.set_minimize_implicit_multibody_accelerations(True)
# solver.set_implicit_multibody_accelerations_weight(1e-3 / trackModel.getNumCoordinates())
# solver.set_implicit_multibody_acceleration_bounds(osim.MocoBounds(-500, 500))
# solver.set_implicit_multibody_acceleration_bounds(osim.MocoBounds(-200, 200))
# solver.set_minimize_implicit_auxiliary_derivatives(True)
# solver.set_implicit_auxiliary_derivatives_weight(0.00001)
solver.resetProblem(problem)

##### acceleration derivatives end up in problem but not guess if using some 
##### of the above settings...

##### seems a guess needs to be created in some way for solver to progress?
##### didn't seem to be the case with some elements commented out...
#solver.setGuess(solver.createGuess())
#does this keep tracked states in guess though?

#Set the normalized tendon forces if not loading initial guess from file
#### TODO: conditional for tendons
guess = solver.getGuess()
numRows = guess.getNumTimes()
stateNames = trackModel.getStateVariableNames()
for ii in range(trackModel.getNumStateVariables()):
    currState = stateNames.get(ii)
    if 'normalized_tendon_force' in currState:
        guess.setState(currState, np.linspace(0.2,0.2,numRows))

#Set guess
solver.setGuess(guess)

#Solve
solution = study.solve()

# solution.unseal()

#### some highly active muscles - muscle strength?
#### feet projecting into ground - contact tracking/pelvis_ty?
#### didn't set anay bounds to help
#### inverse tracking first?

#### seems to go OK when pelvis ty tracking included a bit more - perhaps still
#### weighted a little high --- balance out with other state weights perhaps?

#### perhaps use torque driven solution that tracked well as initial kinematics
#### guess??? --- or actually just track these kinematics with muscle-driven!!!!

#### muscles are definitely not all good --- increase contraction velocity, might
#### need to still increase strength too? 

#### contraction velocity seems to have helped, but there's still an issue with
#### the interplay of the hamstrings, quads and plantarflexors in the stance
#### limb that is preventing adequate motion. could be:
    #### - tendon slack length?
    #### - plantarflexors in 2D model having altered moment arm at knee?
        #### > i think the gastroc is extending the knee at toe-off...
    #### - tendon forces (off-shoot of TSL)?
    
#### including passive forces definitely an issue...

#### scaling the active force width combined with isometric force x 10 pushed
#### through many more iterations (> 1000) before I stopped it --- earlier convergence
#### without this was much simpler. It probably would have converged eventually
#### but it's a seemingly infeasible concept to take that long...
    ##### this did however reduce the activation levels below saturation points
    ##### unclear I suppose if this was due to active force width or strength (or both?)
    
#### upped strength x 10 with no active force width scaling
#### also increased mesh no from 50 to 100...
    ##### obviously double mesh is slower
    ##### stopped this at about 250 iters (~30 mins) to have a look --- optimisation seemed to be going OK though
        ##### toe-pattern was back though --- so active force width may be the fix for this...?
        ##### activations weren't saturating - but r gastroc was still close
        ##### 100 mesh smooths out activations
        
#### scaled active force width to 1.5 with no increase in force
#### retained mesh no of 100...
    ##### stopped this after about 50 iters and 5 mins --- but was going well
    ##### looking pretty good - muscle forces not saturated (high though) even
    ##### without added force
    #####
    ##### this eventually solves in 1 hour 16 mins & ~600 iters
    ##### looks pretty solid - could perhaps speed convergence up by constraining
    ##### initial bound of pelvix_tx to exp. data --- this moved around early
    ##### in iterations...
    

#Option to visualise
study.visualize(solution)

#Create a full gait cycle trajectory from the periodic solution.
addPatterns = [".*pelvis_tx/value"]
fullTraj = osim.createPeriodicTrajectory(solution, addPatterns)
fullTraj.write('sprintTracking_muscleDriven_solution_fullTrajectory.sto')

#Compute ground reaction forces generated by contact sphere from the 
#full gait cycle trajectory
externalLoads = osim.createExternalLoadsTableForGait(trackModel,
                                                     fullTraj,
                                                     forceNames_r,
                                                     forceNames_l)
osim.STOFileAdapter.write(externalLoads,'trackedGRF_muscleDriven_2D.mot')

##### solution is getting better..pelvis_ty symmetry doesn't seem to be there 
##### and that is seemingly the biggest problem??? the tracked GRF is not even
##### present as the spheres are penetrated into the ground to begin with...?
##### weights seem slightly better scaled now...I think
##### only took 15 mins to converge...so not overly complex
##### muscles still don't seem strong enough as signals are quite saturated
##### still not using tendon compliance...

##### to get an actual OK solution:
    #### - double strength and velocity for muscles
    #### - tendon dynamics in implicit mode for plantarflexors
    #### - default solver parameters (not implicit or derivatives stuff)
    
##### in the muscle driven solution there is more plantarflexion though (i.e.
##### too much tension by the plantarflexors, weak dorsiflexors) and more trunk
##### movement...

# %% 2D predictive at higher speed (torque driven)

# This section attempts to generate a predictive simulation of a half gait cycle
# at a slightly faster speed than the tracking simulation using torque actuators

#Set simulation parameters

#Timing and mesh data
initialTime = osim.Storage('refQ_2D.sto').getFirstTime()
finalTime = osim.Storage('refQ_2D.sto').getLastTime()
duration = finalTime - initialTime
meshNo = 50 #set as per mesh number for half gait cycle in Falisse et al. 2019, J R Soc Interface
meshInterval = duration / meshNo

#Weights
effortWeight = 0.1
speedWeight = 1
symmetryWeight = 1

#Create the model processor for the tracking problem

#Load the model and create processor
osimModel = osim.Model('gaitModel2D_doubleStrength.osim')
modelProcessor = osim.ModelProcessor(osimModel)

#Add processing step to remove muscles
modelProcessor.append(osim.ModOpRemoveMuscles())

#Process model for additional parameters
torqueModel = modelProcessor.process()
torqueModel.initSystem()

#Add torque actuators to support model
#Set list to add reserves to
reserveList = ['pelvis_tilt',
                'hip_flexion_l', 'hip_flexion_r',
                'knee_angle_l','knee_angle_r',
                'ankle_angle_l', 'ankle_angle_r']
#Set optimal force and max torque
optimalForce = 100
maxTorque = np.inf
#Add reserves
for tt in range(len(reserveList)):
    addReserve(torqueModel, reserveList[tt], optimalForce, maxTorque)

#Finalise connections
torqueModel.finalizeConnections()

#Get the track model as a processor object
torqueModelProcessor = osim.ModelProcessor(torqueModel)

#Set-up the predictive problem

#Process current model
studyModel = torqueModelProcessor.process()
studyModel.initSystem()
      
#Construct the study object and set basic parameters
study = osim.MocoStudy()
study.setName('sprintPrediction_7ms_torqueDriven')

#Update problem
problem = study.updProblem()

#Set model
problem.setModelCopy(studyModel)

#Set the time bounds (note for higher speed could be more refined)
problem.setTimeBounds(initialTime, [finalTime - 0.2, finalTime])

#Add and set the effort goal to the problem
problem.addGoal(osim.MocoControlGoal('effort'))
effort = osim.MocoControlGoal().safeDownCast(problem.updGoal('effort'))
effort.setDivideByDisplacement(True)
effort.setWeight(effortWeight)
effort.setWeightForControl('/forceset/lumbarAct', 0.001)

#Set the speed goal
#Tracking speed was ~6.5, so going for 7m.s here
problem.addGoal(osim.MocoAverageSpeedGoal('speed'))
speedGoal = osim.MocoAverageSpeedGoal().safeDownCast(problem.updGoal('speed'))
speedGoal.set_desired_average_speed(7.0)
speedGoal.setWeight(speedWeight)

#Set the symmetry goal
problem.addGoal(osim.MocoPeriodicityGoal('symmetry'))
symmetryGoal = osim.MocoPeriodicityGoal().safeDownCast(problem.updGoal('symmetry'))
#Set symmetric coordinate values (except for pelvis_tx) and speeds
for cc in range(studyModel.getCoordinateSet().getSize()):
    #Get current coordinate name
    coordName = studyModel.getCoordinateSet().get(cc).getName()
    #Set symmetry parameters based on coordinate
    if coordName.endswith('_r'):
        #Joint angle value
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/value'))
        #Joint speed
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/speed'))
    elif coordName.endswith('_l'):
        #Joint angle value
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/value'))
        #Joint speed
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/speed'))
    elif coordName.endswith('_extension') or coordName.endswith('_tilt'):
        #Joint angle value
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value'))
        #Joint speed
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed'))
#Add a symmetry pair for pelvis_tx speed
symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
    '/jointset/ground_pelvis/pelvis_tx/speed'))
#Add symmetry goal
symmetryGoal.setWeight(symmetryWeight)

#Add state bounds
#This also sets initial bounds on all coordinates that reflect
#'normal' running motions - based off the experimental data
#These bounds constrain the starting position fairly strongly in a similar fashion
#to what was observed in experimental data...
problem.setStateInfo('/jointset/back/lumbar_extension/value', [np.deg2rad(-30),0],
                     [np.deg2rad(-10),0])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tilt/value', [-10*np.pi/180, 0*np.pi/180],
                     [np.deg2rad(-5),0])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tx/value', [0, 3],
                     1.061) #starting at x = 1.061 from exp. data
problem.setStateInfo('/jointset/ground_pelvis/pelvis_ty/value', [0.8, 1.25],
                     [1, 1.1])
problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [-25*np.pi/180, 75*np.pi/180],
                     [np.deg2rad(-20),np.deg2rad(-10)])
problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [-25*np.pi/180, 75*np.pi/180],
                     [np.deg2rad(30),np.deg2rad(40)])
problem.setStateInfo('/jointset/knee_l/knee_angle_l/value', [-140*np.pi/180, 0],
                     [np.deg2rad(-85),np.deg2rad(-75)])
problem.setStateInfo('/jointset/knee_r/knee_angle_r/value', [-140*np.pi/180, 0],
                     [np.deg2rad(-25),np.deg2rad(-20)])
problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [-20*np.pi/180, 30*np.pi/180],
                     [np.deg2rad(-15),np.deg2rad(-5)])
problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [-20*np.pi/180, 30*np.pi/180],
                     [np.deg2rad(0),np.deg2rad(5)])

#Configure the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.resetProblem(problem)
solver.set_optim_constraint_tolerance(1e-2)
solver.set_optim_convergence_tolerance(1e-2)
solver.set_num_mesh_intervals(meshNo)

#Set guess from tracking simulation
solver.setGuessFile('sprintTracking_torqueDriven_solution.sto')

#Solve
solution = study.solve()

#Option to visualise
study.visualize(solution)

#Create a full gait cycle trajectory from the periodic solution.
addPatterns = [".*pelvis_tx/value"]
fullTraj = osim.createPeriodicTrajectory(solution, addPatterns)
fullTraj.write('sprintPrediction_7ms_torqueDriven_solution_fullTrajectory.sto')

#Compute ground reaction forces generated by contact sphere from the 
#full gait cycle trajectory
#Set force names
forceNamesRightFoot = ['forceset/contactHeel_r',
                       'forceset/contactToe1_r',
                       'forceset/contactToe2_r']
forceNamesLeftFoot = ['forceset/contactHeel_l',
                      'forceset/contactToe1_l',
                      'forceset/contactToe2_l']
#Set force name groups
forceNames_r = osim.StdVectorString()
forceNames_l = osim.StdVectorString()
for ff in range(len(forceNamesRightFoot)):
    forceNames_r.append(forceNamesRightFoot[ff])
    forceNames_l.append(forceNamesLeftFoot[ff])
#Create loads
externalLoads = osim.createExternalLoadsTableForGait(studyModel,
                                                     fullTraj,
                                                     forceNames_r,
                                                     forceNames_l)
osim.STOFileAdapter.write(externalLoads,'predictGRF_7ms_torqueDriven_2D.mot')

#### solution looks a little odd, and the GRFs are almost nothing, which doesn't
#### necessarily suggest a very good or realistic solution...
#### something to do with the torque actuators being able to cheat...

# %% 2D predictive at higher speed (muscle driven)

# This section attempts to generate a predictive simulation of a half gait cycle
# at a slightly faster speed than the tracking simulation using muscle actuators

#Set simulation parameters

#Timing and mesh data
# initialTime = osim.Storage('refQ_2D.sto').getFirstTime()
# finalTime = osim.Storage('refQ_2D.sto').getLastTime()
initialTime = osim.Storage('sprintTracking_torqueDriven_solution.sto').getFirstTime()
finalTime = osim.Storage('sprintTracking_torqueDriven_solution.sto').getLastTime()
duration = finalTime - initialTime
meshNo = 100 #set as per mesh number for half gait cycle in Falisse et al. 2019, J R Soc Interface
meshInterval = duration / meshNo

#Weights
effortWeight = 0.1
speedWeight = 1
symmetryWeight = 1

#Model parameters
complexModel = False

#Muscle parameters
passiveForces = False
implicitTendonCompliance = True

### TODO: conditional for tendons

#Create the model processor for the tracking problem

#Load the model and create processor
if complexModel:
    osimModel = osim.Model('gaitModel2D_complex_doubleStrengthVel.osim')
else:
    osimModel = osim.Model('gaitModel2D_doubleStrengthVel.osim')    
modelProcessor = osim.ModelProcessor(osimModel)

# #Convert muscles to DeGrooteFregley variant
# modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
##### NOTE: don't need to do this as already done earlier
##### Is there still the wrapping object issue with this though???

#Disable tendon compliance (will be re-enabled for certain muscles later)
modelProcessor.append(osim.ModOpIgnoreTendonCompliance())

#Adjust fiber damping as per Dembia et al. Moco paper
modelProcessor.append(osim.ModOpFiberDampingDGF(0.01))

#Turn off passive forces if appropriate
if not passiveForces:
    modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
    
#Modify active force width
#This appears most important for "de-stiffening" muscles to allow appropriate movement
modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))

#Scale isometric forces (even more)
#### Infinite strength vs. infinite activation level?
modelProcessor.append(osim.ModOpScaleMaxIsometricForce(10))
    
#Process model for additional parameters
studyModel = modelProcessor.process()
studyModel.initSystem()

#Enable tendon compliance in the gastrocnemius and soleus
muscles = studyModel.updMuscles()
for mm in np.arange(muscles.getSize()):
    currMusc = osim.DeGrooteFregly2016Muscle.safeDownCast(muscles.get(int(mm)))
    muscName = currMusc.getName()
    #Enable tendon compliance dynamics in the plantarflexors
    if ('gastroc' in muscName) or ('soleus' in muscName):
        currMusc.set_ignore_tendon_compliance(False)
        
#Upscale contraction velocity for super speed
#Also turn off activation dynamics
for mm in range(studyModel.getMuscles().getSize()):
    #Get the current muscle
    currMusc = studyModel.getMuscles().get(mm)
    #Get current max contraction velocity
    currVel = currMusc.getMaxContractionVelocity()
    #Set new max contraction velocity
    currMusc.set_max_contraction_velocity(currVel*10)
    # #Turn off activation dynamics
    # currMusc.set_ignore_activation_dynamics(True)
    #Turning off seems to change states - just set to zero?
    #Can't set to zero --- set really low...
    osim.DeGrooteFregly2016Muscle.safeDownCast(currMusc).set_activation_time_constant(1e-10)
    osim.DeGrooteFregly2016Muscle.safeDownCast(currMusc).set_deactivation_time_constant(1e-10)

#Finalise model connections
studyModel.finalizeConnections()    
        
#Get the track model as a processor object
studyModelProcessor = osim.ModelProcessor(studyModel)

#Set model to use implicit tendon compliance if appropriate
if implicitTendonCompliance:
    studyModelProcessor.append(osim.ModOpUseImplicitTendonComplianceDynamicsDGF())
 
#Construct the study object and set basic parameters
study = osim.MocoStudy()
study.setName('sprintPrediction_3xSpeed_muscleDriven')

#Update problem
problem = study.updProblem()

#Set model
studyModel = studyModelProcessor.process()
studyModel.initSystem()
problem.setModelCopy(studyModel)

#Set the speed goal
#Tracking speed was ~6.5, so going for 7m.s here
#7m.s was ok - trying 9m.s
#9m.s was ok - really pushing to 15m.s

#Get the original speed
df_kinematics = readSTO('refQ_2D.sto')
startPos = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'].to_numpy()[0]
endPos = df_kinematics['/jointset/ground_pelvis/pelvis_tx/value'].to_numpy()[-1]
sprintSpeed = (endPos - startPos) / (finalTime - initialTime)
#Add a speed goal to go at three times as fast
problem.addGoal(osim.MocoAverageSpeedGoal('speed'))
speedGoal = osim.MocoAverageSpeedGoal().safeDownCast(problem.updGoal('speed'))
speedGoal.set_desired_average_speed(sprintSpeed * 3)
speedGoal.setWeight(speedWeight)

#Set the time bounds (note for higher speed could be more refined)
# problem.setTimeBounds(initialTime, [finalTime - 0.2, finalTime])
# problem.setTimeBounds(initialTime, [initialTime + 0.1, finalTime])
#Scale based on speed times three with 10% lee-way either side
problem.setTimeBounds(initialTime, [initialTime + (duration / 3) - (duration / 3 * 0.10),
                                    initialTime + (duration / 3) + (duration / 3 * 0.10)])

#Add and set the effort goal to the problem
problem.addGoal(osim.MocoControlGoal('effort'))
effort = osim.MocoControlGoal().safeDownCast(problem.updGoal('effort'))
effort.setDivideByDisplacement(True)
effort.setWeight(effortWeight)
effort.setWeightForControl('/forceset/lumbarAct', 0.001)

#Set the symmetry goal
problem.addGoal(osim.MocoPeriodicityGoal('symmetry'))
symmetryGoal = osim.MocoPeriodicityGoal().safeDownCast(problem.updGoal('symmetry'))
#Set symmetric coordinate values (except for pelvis_tx) and speeds
for cc in range(studyModel.getCoordinateSet().getSize()):
    #Get current coordinate name
    coordName = studyModel.getCoordinateSet().get(cc).getName()
    #Set symmetry parameters based on coordinate
    if coordName.endswith('_r'):
        #Joint angle value
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/value'))
        #Joint speed
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_r','_l')+'/speed'))
    elif coordName.endswith('_l'):
        #Joint angle value
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value',
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/value'))
        #Joint speed
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed',
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString().replace('_l','_r')+'/speed'))
    elif coordName.endswith('_extension') or coordName.endswith('_tilt') or coordName.endswith('_ty'):
        #Joint angle value
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/value'))
        #Joint speed
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
            studyModel.getCoordinateSet().get(cc).getAbsolutePathString()+'/speed'))
#Add a symmetry pair for pelvis_tx speed
symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
    '/jointset/ground_pelvis/pelvis_tx/speed'))
#Set symmetric activations
for ff in range(studyModel.getForceSet().getSize()):
    #Check for muscle or actuator
    if (studyModel.getForceSet().get(ff).getConcreteClassName().endswith('Muscle') or
        studyModel.getForceSet().get(ff).getConcreteClassName().endswith('Actuator')):
        #Get force name
        forceName = studyModel.getForceSet().get(ff).getName()
        #Set symmetry parameters based on coordinate
        if forceName.endswith('_r'):
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                studyModel.getForceSet().get(ff).getAbsolutePathString()+'/activation',
                studyModel.getForceSet().get(ff).getAbsolutePathString().replace('_r','_l')+'/activation'))
        elif forceName.endswith('_l'):
            symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
                studyModel.getForceSet().get(ff).getAbsolutePathString()+'/activation',
                studyModel.getForceSet().get(ff).getAbsolutePathString().replace('_l','_r')+'/activation'))
        # elif forceName == 'lumbarAct':
        #     symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
        #         trackModel.getForceSet().get(ff).getAbsolutePathString()))
#Add symmetry goal
symmetryGoal.setWeight(symmetryWeight)

#Add state bounds
#This also sets initial bounds on all coordinates that reflect
#'normal' running motions - based off the experimental data
#These bounds constrain the starting position fairly strongly in a similar fashion
#to what was observed in experimental data...
#opened up for faster running (15m.s) --- are these appropriate, normal kinematics for sprinting?
problem.setStateInfo('/jointset/back/lumbar_extension/value', [np.deg2rad(-30),0])#,
                     # [np.deg2rad(-10),0])
                     # [np.deg2rad(-20),0])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tilt/value', [-10*np.pi/180, 0*np.pi/180],
                     # [np.deg2rad(-5),0])
                     [np.deg2rad(-10),0])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tx/value', [0, 10],
                     1.061) #starting at x = 1.061 from exp. data
problem.setStateInfo('/jointset/ground_pelvis/pelvis_ty/value', [0.8, 1.25],
                     [1, 1.1])
problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [-35*np.pi/180, 75*np.pi/180],
                     # [np.deg2rad(-20),np.deg2rad(-10)])
                     [np.deg2rad(-35),np.deg2rad(-5)])
problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [-25*np.pi/180, 75*np.pi/180],
                     # [np.deg2rad(30),np.deg2rad(40)])
                     [np.deg2rad(20),np.deg2rad(50)])
problem.setStateInfo('/jointset/knee_l/knee_angle_l/value', [-140*np.pi/180, 0],
                     # [np.deg2rad(-85),np.deg2rad(-75)])
                     [np.deg2rad(-100),np.deg2rad(-70)])
problem.setStateInfo('/jointset/knee_r/knee_angle_r/value', [-140*np.pi/180, 0],
                     # [np.deg2rad(-25),np.deg2rad(-20)])
                     [np.deg2rad(-40),np.deg2rad(-10)])
problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [-25*np.pi/180, 30*np.pi/180],
                     # [np.deg2rad(-15),np.deg2rad(-5)])
                     [np.deg2rad(-25),np.deg2rad(0)])
problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [-20*np.pi/180, 30*np.pi/180],
                     # [np.deg2rad(0),np.deg2rad(5)])
                     [np.deg2rad(-5),np.deg2rad(15)])

#Configure the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.resetProblem(problem)
solver.set_optim_constraint_tolerance(1e-2)
solver.set_optim_convergence_tolerance(1e-2)
solver.set_num_mesh_intervals(meshNo)

#Set guess from tracking simulation
solver.setGuessFile('sprintTracking_muscleDriven_solution.sto')
# solver.setGuessFile('sprintPrediction_9ms_muscleDriven_solution.sto') #piggy back off predictions

#Solve
solution = study.solve()

#### predictive sim converged in 24 mins and 70 iters (7m.s)
#### predictive sim covered in 34 mins and 114 iters (9m.s)
#### running at 15m.s was struggling at an hour - stopped to check
#### solution looks pretty good...

#Option to visualise
study.visualize(solution)

#Create a full gait cycle trajectory from the periodic solution.
addPatterns = [".*pelvis_tx/value"]
fullTraj = osim.createPeriodicTrajectory(solution, addPatterns)
fullTraj.write('sprintPrediction_3xSpeed_muscleDriven_solution_fullTrajectory.sto')

#Compute ground reaction forces generated by contact sphere from the 
#full gait cycle trajectory
#Set force names
forceNamesRightFoot = ['forceset/contactHeel_r',
                       'forceset/contactToe1_r',
                       'forceset/contactToe2_r']
forceNamesLeftFoot = ['forceset/contactHeel_l',
                      'forceset/contactToe1_l',
                      'forceset/contactToe2_l']
#Set force name groups
forceNames_r = osim.StdVectorString()
forceNames_l = osim.StdVectorString()
for ff in range(len(forceNamesRightFoot)):
    forceNames_r.append(forceNamesRightFoot[ff])
    forceNames_l.append(forceNamesLeftFoot[ff])
#Compute external loads
externalLoads = osim.createExternalLoadsTableForGait(studyModel,
                                                     fullTraj,
                                                     forceNames_r,
                                                     forceNames_l)
osim.STOFileAdapter.write(externalLoads,'predictGRF_3xSpeed_muscleDriven_2D.mot')

#### predicted forces look a little off, but not too bad...
#### there is too much of a double peak - probably due to lack of spheres in mid-foot
#### large braking forces - perhaps due to too large heel sphere --- maybe just move and don't scale?

#### adjust final time to allow for faster running...
#### perhaps set initial and final bounds for pelvis_tx to speed goal? with a little flexibility...
    ##### don't know if there is a pelvis_tx/speed state?
#### minimise lumbar motion? a bit of flexion and extension throughout? or keep posture upright?

#### 3 x speed can be achieved with really high contractile velocities and strength
#### muscle activations are extremely noisy due to being switched on and off all the time
    #### due to activation/deactivation dynamics being so low
        #### suggest best bet IS to turn of activation dynamics and resolve the guess/problem...
#### mini-flight phase during stance from GRF data probably due to no mid-foot spheres...