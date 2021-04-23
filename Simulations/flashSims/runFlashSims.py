# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 16:30:31 2021

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
"""

# %% SAMPLE from original predictive sim tests

import opensim as osim
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

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

# %% 2D predictive at higher speed (muscle driven)

# This section attempts to generate a predictive simulation of a half gait cycle
# at a slightly faster speed than the tracking simulation using muscle actuators

#Set simulation parameters

#Timing and mesh data
# initialTime = osim.Storage('refQ_2D.sto').getFirstTime()
# finalTime = osim.Storage('refQ_2D.sto').getLastTime()
initialTime = osim.Storage('..\\trackingSims\\sprintTracking_torqueDriven_solution.sto').getFirstTime()
finalTime = osim.Storage('..\\trackingSims\\sprintTracking_torqueDriven_solution.sto').getLastTime()
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
tendonDynamics = True
activationDynamics = False
contractVelScale = 10 #upscale for sprinting
maxForceScale = 10 #upscale for sprinting 


#Create the model processor for the tracking problem

#Load the model and create processor
if complexModel:
    osimModel = osim.Model('..\\trackingSims\\gaitModel2D_complex.osim')
else:
    osimModel = osim.Model('..\\trackingSims\\gaitModel2D.osim')
    
#Update contractile velocity maximum to assist with sprinting sim
for mm in range(osimModel.getMuscles().getSize()):
    #Get the current muscle
    currMusc = osimModel.getMuscles().get(mm)
    #Get current max contraction velocity
    currVel = currMusc.getMaxContractionVelocity()
    #Set new max contraction velocity
    currMusc.set_max_contraction_velocity(currVel*contractVelScale)
    
#Finalise connections
osimModel.finalizeConnections()
    
#Put the model in a processor
modelProcessor = osim.ModelProcessor(osimModel)

#Scale the maximum isometric force of muscles to assist with sprinting sim
modelProcessor.append(osim.ModOpScaleMaxIsometricForce(maxForceScale))

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
    
#Process model for additional parameters
studyModel = modelProcessor.process()
studyModel.initSystem()

#Enable tendon compliance in the gastrocnemius and soleus
#Check activation dynamics
muscles = studyModel.updMuscles()
for mm in np.arange(muscles.getSize()):
    #Get current muscle
    currMusc = osim.DeGrooteFregly2016Muscle.safeDownCast(muscles.get(int(mm)))
    #Switch off activation dynamics if necessary
    if not activationDynamics:
        currMusc.set_ignore_activation_dynamics(True)    
    #Get muscle name
    muscName = currMusc.getName()
    #Enable tendon compliance dynamics in the plantarflexors
    if ('gastroc' in muscName) or ('soleus' in muscName):
        currMusc.set_ignore_tendon_compliance(False)

#Finalise connections
studyModel.finalizeConnections()

#Get the model as a processor object
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
problem.setModelAsCopy(studyModel)

#Set the speed goal
#Tracking speed was ~6.5, so going for 7m.s here
#7m.s was ok - trying 9m.s
#9m.s was ok - really pushing to 15m.s

#Get the original speed
df_kinematics = readSTO('..\\trackingSims\\refQ_2D.sto')
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
##### TODO: such low weight here might be why you get lumbar fluctuation in predictive???

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
#### TODO: without activation dynamics, might need something different here...
# #Set symmetric activations
# for ff in range(studyModel.getForceSet().getSize()):
#     #Check for muscle or actuator
#     if (studyModel.getForceSet().get(ff).getConcreteClassName().endswith('Muscle') or
#         studyModel.getForceSet().get(ff).getConcreteClassName().endswith('Actuator')):
#         #Get force name
#         forceName = studyModel.getForceSet().get(ff).getName()
#         #Set symmetry parameters based on coordinate
#         #Check if using activation dynamics within here
#         if forceName.endswith('_r'):
#             if activationDynamics:
#                 symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
#                     studyModel.getForceSet().get(ff).getAbsolutePathString()+'/activation',
#                     studyModel.getForceSet().get(ff).getAbsolutePathString().replace('_r','_l')+'/activation'))
#             else:
#                 symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
#                     studyModel.getForceSet().get(ff).getAbsolutePathString(),
#                     studyModel.getForceSet().get(ff).getAbsolutePathString().replace('_r','_l')))                    
#         elif forceName.endswith('_l'):
#             if activationDynamics:
#                 symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
#                     studyModel.getForceSet().get(ff).getAbsolutePathString()+'/activation',
#                     studyModel.getForceSet().get(ff).getAbsolutePathString().replace('_l','_r')+'/activation'))
#             else:
#                 symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
#                     studyModel.getForceSet().get(ff).getAbsolutePathString(),
#                     studyModel.getForceSet().get(ff).getAbsolutePathString().replace('_l','_r')))
#         # elif forceName == 'lumbarAct':
#         #     symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(
#         #         trackModel.getForceSet().get(ff).getAbsolutePathString()))
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
#Different if activation dynamics taken out
guess = solver.createGuess()
origGuess = osim.MocoTrajectory('..\\trackingSims\\sprintTracking_muscleDriven_solution.sto')

#Match up guess nodes
#Times theoretically already match up?
if guess.getNumTimes() != origGuess.getNumTimes():
    guess.setNumTimes(origGuess.getNumTimes())
    
#Fill the guess with data from the old guess where column headers match
#States
for gg in range(len(guess.getStateNames())):
    #Get current state name
    currState = guess.getStateNames()[gg]
    #Fill the value from the original guess
    guess.setState(currState,origGuess.getStateMat(currState))
#Controls
for gg in range(len(guess.getControlNames())):
    #Get current control name
    currControl = guess.getControlNames()[gg]
    #Fill the value from the original guess
    guess.setControl(currControl,origGuess.getControlMat(currControl))
#Derivatives
for gg in range(len(guess.getDerivativeNames())):
    #Get current derivative name
    currDerivative = guess.getDerivativeNames()[gg]
    #Fill the value from the original guess
    guess.setDerivative(currDerivative,origGuess.getDerivativesTrajectoryMat()[:,gg])
    
#Reset time due to it getting filled with nan's
guess.setTime(origGuess.getTime())

#Set guess
solver.setGuess(guess)

# solver.setGuessFile('sprintTracking_muscleDriven_solution.sto')
# solver.setGuessFile('sprintPrediction_9ms_muscleDriven_solution.sto') #piggy back off predictions

#Solve
solution = study.solve()

#### predictive sim converged in 24 mins and 70 iters (7m.s)
#### predictive sim covered in 34 mins and 114 iters (9m.s)
#### running at 15m.s was struggling at an hour - stopped to check
#### solution looks pretty good...

#### above at 3 x speed solves pretty quickly in 10 mins
#### quick little strides so will be interesting to see further speed progression

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
		##### just figure out how far tx is going to travel with that speed and put the bounds on that
#### minimise lumbar motion? a bit of flexion and extension throughout? or keep posture upright?

#### 3 x speed can be achieved with really high contractile velocities and strength
#### muscle activations are extremely noisy due to being switched on and off all the time
    #### due to activation/deactivation dynamics being so low
        #### suggest best bet IS to turn of activation dynamics and resolve the guess/problem...
#### mini-flight phase during stance from GRF data probably due to no mid-foot spheres...