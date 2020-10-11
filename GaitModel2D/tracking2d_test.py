# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 22:40:53 2020

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    Test of 2D simulation with new model
    
    TODO: balance weights of goals better; scaling?
    
    NOTES:
        - Restoration fails > maybe passive forces with hip/knee flexion,
        inadequate muscle strength, active force width, lack of loads, lack of
        reserves, metabolics?
        
        - Without passive forces, tendon compliance, and wider active force width
        the problem solves quite quickly (7 minutes). The movement is not tracked
        well though, it's a weird bounding sort of thing. Higher state tracking
        weight is probably necessary...
        
        - Adding individualised state weights to above produce a good kinematic
        result and also reduced the metabolic cost -- took ~15 mins to solve
        
        - Add GRF tracking -- results in slightly longer sim, but converges in
        ~26 mins with ~770 iterations. COntact tracking makes up about 9 of the
        objective function cost, while state tracking is at 20. There's some weird
        stuff going on with the trunk angle, seemingly because of the need to
        be periodic...could try removing this...GRF looks good, but not smooth
        though given the mesh interval used (i.e. 25)
        
        - Took out symmetric pelvis motion to try and fix trunk issue, symmetric
        trunk activation remained though; also went up to a 50 mesh solution and
        used the previous solution as a guess
        
        - Weird kink in joining full stride might actually be lack of periodicity
        goal on trunk -- seemed to struggle with when the trunk periodicity was added? 
        This corresponded to using the ignore tendon compliance model operator
        but I wouldn't think that would make a difference -- it definitely seems
        like it's the periodic trunk stuff. Solution might be to just sim a whole
        gait cycle, although that could prove time costly...
        
        - Test out having tendon compliance on, can add a model processor to 
        ignore this if wanting to shift back with model -- struggles, could be
        to do with the metabolics, test without -- still struggles...
        
        - Trunk issue are potentially stemming from fact that runk is lurching 
        forward in half simulations -- test out if this happens with full cycle
        
        - Contraction velocity already needs to be increased...
        
        - Kinemtics aren't really that periodic in gait cycle -- maybe save these
        for predictive simulations
        
        - It makes sense that the GRF can't quite match up with the new model,
        withou RRA the kinematics might be off and the mass/inertia of models is
        different
        
        - Maybe try full blown predictive sim - or could try a predictive sim
        with some basic kinematics (i.e. hip, knee, ankle, trunk -- not pelvis)
        prescribed in starting guess...
    
"""

import opensim as osim
import math
# import numpy as np



# %% Make shift RRA tracking simulation
#    Align kinematics with torque driven simulation w/ GRF tracking

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName('torque_driven_tracking')

#Construct a ModelProcessor and add it to the tool.
modelProcessor = osim.ModelProcessor('testModel.osim')
#Add CoordinateActuators to the model degrees-of-freedom. This
modelProcessor.append(osim.ModOpAddReserves(1000))
#Ste model processor
track.setModel(modelProcessor)

#Setup of MocoTrack
tabProcessor = osim.TableProcessor('refQ_converted.sto')
tabProcessor.append(osim.TabOpLowPassFilter(12))
track.setStatesReference(tabProcessor)
track.set_states_global_tracking_weight(10)
track.set_allow_unused_references(True)
track.set_track_reference_position_derivatives(True)
track.set_apply_tracked_states_to_guess(True)
track.set_initial_time(osim.Storage('refQ_converted.sto').getFirstTime())
track.set_final_time(osim.Storage('refQ_converted.sto').getLastTime())

#Set individual state weights
#These are based off Ross Miller's UMocoD project for now
stateWeights = osim.MocoWeightSet()
#elevated pelvis_tx weight to try and account for weirdness? (100 > 200)
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tx/value',(1/(1*0.2000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_ty/value',(1/(2*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tilt/value',(1/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/value',(1/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/value',(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_r/knee_angle_r/value',(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/value',(1/(1*0.0574))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/value',(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_l/knee_angle_l/value',(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/value',(1/(1*0.0574))**2))
w = 0.001 #Scale the generalized speed tracking errors by this constant
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tx/speed',w*(0/(1*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_ty/speed',w*(0/(2*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tilt/speed',w*(0/(1*0.0585))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/speed',w*(0/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/speed',w*(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_r/knee_angle_r/speed',w*(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/speed',w*(1/(1*0.0574))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/speed',w*(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_l/knee_angle_l/speed',w*(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/speed',w*(1/(1*0.0574))**2))
track.set_states_weight_set(stateWeights)

#Initialise study and problem
study = track.initialize()
problem = study.updProblem()

#Goals

#Effort
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
effort.setWeight(1)

#Put a large weight on the pelvis CoordinateActuators, which act as the
#residual, or 'hand-of-god', forces which we would like to keep as small
#as possible.
model = modelProcessor.process()
model.initSystem()
forceSet = model.getForceSet()
for ii in range(0,forceSet.getSize()):
   forcePath = forceSet.get(ii).getAbsolutePathString()
   if 'pelvis' in forcePath:
       effort.setWeightForControl(forcePath,10)

#GRF contact tracking goal
contactTracking = osim.MocoContactTrackingGoal('contact',1)
contactTracking.setExternalLoadsFile('refGRF_2D.xml')
#Right foot
forceNamesRightFoot = osim.StdVectorString()
forceNamesRightFoot.append('contactHeel_r')
forceNamesRightFoot.append('contactToe1_r')
forceNamesRightFoot.append('contactToe2_r')
trackRightGRF = osim.MocoContactTrackingGoalGroup(forceNamesRightFoot,'RightGRF')
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r')
contactTracking.addContactGroup(trackRightGRF)
#Left foot
forceNamesLeftFoot = osim.StdVectorString()
forceNamesLeftFoot.append('contactHeel_l')
forceNamesLeftFoot.append('contactToe1_l')
forceNamesLeftFoot.append('contactToe2_l')
trackLeftGRF = osim.MocoContactTrackingGoalGroup(forceNamesLeftFoot,'LeftGRF')
trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l')
contactTracking.addContactGroup(trackLeftGRF)
#Properties
contactTracking.setProjection('plane')
contactTracking.setProjectionVector(osim.Vec3(0, 0, 1))
problem.addGoal(contactTracking)

#Bounds
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tilt/value', [-20*math.pi/180, 5*math.pi/180])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tx/value', [-5, 5])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_ty/value', [0.75, 1.25])
problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [-40*math.pi/180, 80*math.pi/180])
problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [-40*math.pi/180, 80*math.pi/180])
problem.setStateInfo('/jointset/knee_l/knee_angle_l/value', [-140*math.pi/180, 0])
problem.setStateInfo('/jointset/knee_r/knee_angle_r/value', [-140*math.pi/180, 0])
problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [-30*math.pi/180, 30*math.pi/180])
problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [-30*math.pi/180, 30*math.pi/180])
problem.setStateInfo('/jointset/back/lumbar_extension/value', [-30, 0*math.pi/180])

#Configure the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_num_mesh_intervals(50)
solver.set_verbosity(2)
solver.set_optim_solver('ipopt')
solver.set_optim_convergence_tolerance(1e-4)
solver.set_optim_constraint_tolerance(1e-4)
solver.set_optim_max_iterations(1000)
# solver.setGuessFile('MocoStudy_solution.sto')

#Solve
solution = study.solve()

# %%

#Setup tracking object
track = osim.MocoTrack()

#Get model
baseModel = osim.Model('gait9dof18musc_Ong_et_al_Moco.osim')

#Setup metabolics
metabolics = osim.Bhargava2004Metabolics()
metabolics.setName('metabolics')
metabolics.set_use_smoothing(True)
metabolics.addMuscle('hamstrings_r',baseModel.getComponent('hamstrings_r'))
metabolics.addMuscle('hamstrings_l',baseModel.getComponent('hamstrings_l'))
metabolics.addMuscle('bifemsh_r',baseModel.getComponent('bifemsh_r'))
metabolics.addMuscle('bifemsh_l',baseModel.getComponent('bifemsh_l'))
metabolics.addMuscle('glut_max_r',baseModel.getComponent('glut_max_r'))
metabolics.addMuscle('glut_max_l',baseModel.getComponent('glut_max_l'))
metabolics.addMuscle('iliopsoas_r',baseModel.getComponent('iliopsoas_r'))
metabolics.addMuscle('iliopsoas_l',baseModel.getComponent('iliopsoas_l'))
metabolics.addMuscle('rect_fem_r',baseModel.getComponent('rect_fem_r'))
metabolics.addMuscle('rect_fem_l',baseModel.getComponent('rect_fem_l'))
metabolics.addMuscle('vasti_r',baseModel.getComponent('vasti_r'))
metabolics.addMuscle('vasti_l',baseModel.getComponent('vasti_l'))
metabolics.addMuscle('gastroc_r',baseModel.getComponent('gastroc_r'))
metabolics.addMuscle('gastroc_l',baseModel.getComponent('gastroc_l'))
metabolics.addMuscle('soleus_r',baseModel.getComponent('soleus_r'))
metabolics.addMuscle('soleus_l',baseModel.getComponent('soleus_l'))
metabolics.addMuscle('tib_ant_r',baseModel.getComponent('tib_ant_r'))
metabolics.addMuscle('tib_ant_l',baseModel.getComponent('tib_ant_l'))
baseModel.addComponent(metabolics)
# baseModel.finalizeConnections()

#Turn off contact forces to make it standard simulation
baseModel.getComponent('contactHeel_r').set_appliesForce(False)
baseModel.getComponent('contactToe1_r').set_appliesForce(False)
baseModel.getComponent('contactToe2_r').set_appliesForce(False)
baseModel.getComponent('contactHeel_l').set_appliesForce(False)
baseModel.getComponent('contactToe1_l').set_appliesForce(False)
baseModel.getComponent('contactToe2_l').set_appliesForce(False)
baseModel.finalizeConnections()

#Process model
modelProcessor = osim.ModelProcessor(baseModel)
#Ignore passive forces
modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
#Scale active force width
modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
#Ignore tendon compliance
modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
#Add ground reaction external loads in lieu of a ground-contact model.
modelProcessor.append(osim.ModOpAddExternalLoads('refGRF_2D.xml'))
# #Add reserves
# modelProcessor.append(osim.ModOpAddReserves(1))

#Setup of MocoTrack
track.setModel(modelProcessor)
tabProcessor = osim.TableProcessor('refQ_converted.sto')
tabProcessor.append(osim.TabOpLowPassFilter(12))
track.setStatesReference(tabProcessor)
track.set_states_global_tracking_weight(10)
track.set_allow_unused_references(True)
track.set_track_reference_position_derivatives(True)
track.set_apply_tracked_states_to_guess(True)
track.set_initial_time(osim.Storage('refQ_converted.sto').getFirstTime())
track.set_final_time(osim.Storage('refQ_converted.sto').getLastTime())

#Set individual state weights
#These are based off Ross Miller's UMocoD project for now
stateWeights = osim.MocoWeightSet()
#elevated pelvis_tx weight to try and account for weirdness? (100 > 200)
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tx/value',(1/(1*0.2000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_ty/value',(1/(2*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tilt/value',(1/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/value',(1/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/value',(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_r/knee_angle_r/value',(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/value',(1/(1*0.0574))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/value',(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_l/knee_angle_l/value',(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/value',(1/(1*0.0574))**2))
w = 0.001 #Scale the generalized speed tracking errors by this constant
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tx/speed',w*(0/(1*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_ty/speed',w*(0/(2*0.1000))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tilt/speed',w*(0/(1*0.0585))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/speed',w*(0/(1*0.1745))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/speed',w*(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_r/knee_angle_r/speed',w*(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/speed',w*(1/(1*0.0574))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/speed',w*(1/(1*0.0647))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_l/knee_angle_l/speed',w*(1/(1*0.0889))**2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/speed',w*(1/(1*0.0574))**2))
track.set_states_weight_set(stateWeights)

#Initialise study and problem
study = track.initialize()
problem = study.updProblem()

#Goals

# #Symmetry
# periodicityGoal = osim.MocoPeriodicityGoal('symmetryGoal')
# problem.addGoal(periodicityGoal)
model = modelProcessor.process()
model.initSystem()

#Set symmetric coordinate values

#Half gait cycle

# #Set symmetry pairs
# symPairs = [['hip_flexion_r','hip_flexion_l'],
#             ['knee_angle_r','knee_angle_l'],
#             ['ankle_angle_r','ankle_angle_l']]
# for ii in range(0,len(symPairs)):
#     #Set the jointsets depending on current pair
#     if 'hip' in symPairs[ii][0]:
#         jointName = ['/jointset/hip_r/','/jointset/hip_l/']
#     elif 'knee' in symPairs[ii][0]:
#         jointName = ['/jointset/knee_r/','/jointset/knee_l/']
#     elif 'ankle' in symPairs[ii][0]:
#         jointName = ['/jointset/ankle_r/','/jointset/ankle_l/']
    
#     #Set the pair for coordinate value
#     pair = osim.MocoPeriodicityGoalPair(jointName[0]+symPairs[ii][0]+'/value',
#                                         jointName[1]+symPairs[ii][1]+'/value')
#     #Add to the goal
#     periodicityGoal.addStatePair(pair)
    
#     #Set the pair for coordinate speed
#     pair = osim.MocoPeriodicityGoalPair(jointName[0]+symPairs[ii][0]+'/speed',
#                                         jointName[1]+symPairs[ii][1]+'/speed')
#     #Add to the goal
#     periodicityGoal.addStatePair(pair)

# #Symmetric pelvis joint coordinates
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_ty/value'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_ty/speed'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_tilt/value'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_tilt/speed'))
# ##### Seems to throw out solution for some reason?
# # periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/back/lumbar_extension/value'))
# # periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/back/lumbar_extension/speed'))
    
# #Full gait cycle

# #Set symmetry across relevant states from start to end
# ##### Efficiency of this could be improved
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_ty/value'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_ty/speed'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_tilt/value'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_tilt/speed'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/back/lumbar_extension/value'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/back/lumbar_extension/speed'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/hip_r/hip_flexion_r/value'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/hip_r/hip_flexion_r/speed'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/knee_r/knee_angle_r/value'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/knee_r/knee_angle_r/speed'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ankle_r/ankle_angle_r/value'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ankle_r/ankle_angle_r/speed'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/hip_l/hip_flexion_l/value'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/hip_l/hip_flexion_l/speed'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/knee_l/knee_angle_l/value'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/knee_l/knee_angle_l/speed'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ankle_l/ankle_angle_l/value'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ankle_l/ankle_angle_l/speed'))

# #Set symmetric muscle activations

# #Half gait cycle

# # #Set symmetry pairs
# # symPairs = [['hamstrings_r','hamstrings_l'],
# #             ['bifemsh_r','bifemsh_l'],
# #             ['glut_max_r','glut_max_l'],
# #             ['iliopsoas_r','iliopsoas_l'],
# #             ['rect_fem_r','rect_fem_l'],
# #             ['vasti_r','vasti_l'],
# #             ['gastroc_r','gastroc_l'],
# #             ['soleus_r','soleus_l'],
# #             ['tib_ant_r','tib_ant_l']]
# # for ii in range(0,len(symPairs)):
    
# #     #Set the pair for coordinate value
# #     pair = osim.MocoPeriodicityGoalPair('/'+symPairs[ii][0]+'/activation',
# #                                         '/'+symPairs[ii][1]+'/activation')
# #     #Add to the goal
# #     periodicityGoal.addStatePair(pair)

# # ##### Turn off if not having matched periodic lumbar kinematics???
# # # #Symmetric coordinate actuator controls
# # # periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair('/lumbarAct'))

# #Full gait cycle

# #Set symmetry across relevant states from start to end
# ##### Efficiency of this could be improved
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/hamstrings_r/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/bifemsh_r/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/glut_max_r/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/iliopsoas_r/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/rect_fem_r/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/vasti_r/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/gastroc_r/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/soleus_r/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/tib_ant_r/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/hamstrings_l/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/bifemsh_l/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/glut_max_l/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/iliopsoas_l/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/rect_fem_l/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/vasti_l/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/gastroc_l/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/soleus_l/activation'))
# periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/tib_ant_l/activation'))
# periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair('/lumbarAct'))

#Effort
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
effort.setWeight(0.1)

#Metabolics
metGoal = osim.MocoOutputGoal('met',0.1)
metGoal.setOutputPath('/metabolics|total_metabolic_rate')
metGoal.setDivideByDisplacement(True)
metGoal.setDivideByMass(True)
problem.addGoal(metGoal)

# #GRF contact tracking goal
# contactTracking = osim.MocoContactTrackingGoal('contact',1)
# contactTracking.setExternalLoadsFile('refGRF_2D.xml')
# #Right foot
# forceNamesRightFoot = osim.StdVectorString()
# forceNamesRightFoot.append('contactHeel_r')
# forceNamesRightFoot.append('contactToe1_r')
# forceNamesRightFoot.append('contactToe2_r')
# trackRightGRF = osim.MocoContactTrackingGoalGroup(forceNamesRightFoot,'RightGRF')
# trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r')
# contactTracking.addContactGroup(trackRightGRF)
# #Left foot
# forceNamesLeftFoot = osim.StdVectorString()
# forceNamesLeftFoot.append('contactHeel_l')
# forceNamesLeftFoot.append('contactToe1_l')
# forceNamesLeftFoot.append('contactToe2_l')
# trackLeftGRF = osim.MocoContactTrackingGoalGroup(forceNamesLeftFoot,'LeftGRF')
# trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l')
# contactTracking.addContactGroup(trackLeftGRF)
# #Properties
# contactTracking.setProjection('plane')
# contactTracking.setProjectionVector(osim.Vec3(0, 0, 1))
# problem.addGoal(contactTracking)

#Bounds
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tilt/value', [-20*math.pi/180, 0*math.pi/180])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tx/value', [-5, 5])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_ty/value', [0.75, 1.25])
problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [-40*math.pi/180, 80*math.pi/180])
problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [-40*math.pi/180, 80*math.pi/180])
problem.setStateInfo('/jointset/knee_l/knee_angle_l/value', [-140*math.pi/180, 0])
problem.setStateInfo('/jointset/knee_r/knee_angle_r/value', [-140*math.pi/180, 0])
problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [-30*math.pi/180, 30*math.pi/180])
problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [-30*math.pi/180, 30*math.pi/180])
problem.setStateInfo('/jointset/back/lumbar_extension/value', [-30, 0*math.pi/180])

#Configure the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_num_mesh_intervals(50)
solver.set_verbosity(2)
solver.set_optim_solver('ipopt')
solver.set_optim_convergence_tolerance(1e-4)
solver.set_optim_constraint_tolerance(1e-4)
solver.set_optim_max_iterations(1000)
# solver.setGuessFile('MocoStudy_solution.sto')

# #Set the normalized tendon forces if not loading initial guess from file
# #This is for when tendon compliance is used
# guess = solver.getGuess()
# numRows = guess.getNumTimes();
# stateNames = model.getStateVariableNames()
# for ii in range(0,model.getNumStateVariables()):
#     currentStateName = stateNames.getitem(ii)
#     if 'normalized_tendon_force' in currentStateName:
#         guess.setState(currentStateName, np.linspace(0.2,0.2,numRows))

#Solve
solution = study.solve()

# #Create full stride
# full = osim.createPeriodicTrajectory(solution)

#Print metabolics output
print(str(round(solution.getObjectiveTerm('met'),3))+' J kg-1 m-1')

#Visualise
# study.visualize(full)
study.visualize(solution)

#Write tracked GRF to file
contact_r = osim.StdVectorString()
contact_l = osim.StdVectorString()
contact_r.append('contactHeel_r')
contact_r.append('contactToe1_r')
contact_r.append('contactToe2_r')
contact_l.append('contactHeel_l')
contact_l.append('contactToe1_l')
contact_l.append('contactToe2_l')
externalForcesTableFlat = osim.createExternalLoadsTableForGait(model,
                                                               solution,
                                                               contact_r,
                                                               contact_l)
osim.writeTableToFile(externalForcesTableFlat,
                      'solutionGRF_fullStride.sto')

# %% Predictive simulation

# Uses a tracking simulation as a guess...

#Setup study
study = osim.MocoStudy()
study.setName('gaitPrediction')
problem = study.updProblem()

#Get model
baseModel = osim.Model('gait9dof18musc_Ong_et_al_Moco.osim')

#Setup metabolics
metabolics = osim.Bhargava2004Metabolics()
metabolics.setName('metabolics')
metabolics.set_use_smoothing(True)
metabolics.addMuscle('hamstrings_r',baseModel.getComponent('hamstrings_r'))
metabolics.addMuscle('hamstrings_l',baseModel.getComponent('hamstrings_l'))
metabolics.addMuscle('bifemsh_r',baseModel.getComponent('bifemsh_r'))
metabolics.addMuscle('bifemsh_l',baseModel.getComponent('bifemsh_l'))
metabolics.addMuscle('glut_max_r',baseModel.getComponent('glut_max_r'))
metabolics.addMuscle('glut_max_l',baseModel.getComponent('glut_max_l'))
metabolics.addMuscle('iliopsoas_r',baseModel.getComponent('iliopsoas_r'))
metabolics.addMuscle('iliopsoas_l',baseModel.getComponent('iliopsoas_l'))
metabolics.addMuscle('rect_fem_r',baseModel.getComponent('rect_fem_r'))
metabolics.addMuscle('rect_fem_l',baseModel.getComponent('rect_fem_l'))
metabolics.addMuscle('vasti_r',baseModel.getComponent('vasti_r'))
metabolics.addMuscle('vasti_l',baseModel.getComponent('vasti_l'))
metabolics.addMuscle('gastroc_r',baseModel.getComponent('gastroc_r'))
metabolics.addMuscle('gastroc_l',baseModel.getComponent('gastroc_l'))
metabolics.addMuscle('soleus_r',baseModel.getComponent('soleus_r'))
metabolics.addMuscle('soleus_l',baseModel.getComponent('soleus_l'))
metabolics.addMuscle('tib_ant_r',baseModel.getComponent('tib_ant_r'))
metabolics.addMuscle('tib_ant_l',baseModel.getComponent('tib_ant_l'))
baseModel.addComponent(metabolics)
baseModel.finalizeConnections()

#Process model
modelProcessor = osim.ModelProcessor(baseModel)
#Ignore passive forces
modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
#Scale active force width
modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
#Ignore tendon compliance
modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
#Set model processor
problem.setModelProcessor(modelProcessor)

#Goals

#Symmetry goal - full gait cycle
periodicityGoal = osim.MocoPeriodicityGoal('symmetryGoal')
problem.addGoal(periodicityGoal)
model = modelProcessor.process()
model.initSystem()

#Set symmetry across relevant states from start to end
##### Efficiency of this could be improved
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_ty/value'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_ty/speed'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_tilt/value'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_tilt/speed'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/back/lumbar_extension/value'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/back/lumbar_extension/speed'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/hip_r/hip_flexion_r/value'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/hip_r/hip_flexion_r/speed'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/knee_r/knee_angle_r/value'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/knee_r/knee_angle_r/speed'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ankle_r/ankle_angle_r/value'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ankle_r/ankle_angle_r/speed'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/hip_l/hip_flexion_l/value'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/hip_l/hip_flexion_l/speed'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/knee_l/knee_angle_l/value'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/knee_l/knee_angle_l/speed'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ankle_l/ankle_angle_l/value'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/jointset/ankle_l/ankle_angle_l/speed'))

#Set symmetry across relevant states from start to end
##### Efficiency of this could be improved
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/hamstrings_r/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/bifemsh_r/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/glut_max_r/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/iliopsoas_r/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/rect_fem_r/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/vasti_r/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/gastroc_r/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/soleus_r/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/tib_ant_r/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/hamstrings_l/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/bifemsh_l/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/glut_max_l/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/iliopsoas_l/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/rect_fem_l/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/vasti_l/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/gastroc_l/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/soleus_l/activation'))
periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair('/tib_ant_l/activation'))
periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair('/lumbarAct'))

#Metabolics
metGoal = osim.MocoOutputGoal('met',1)
metGoal.setOutputPath('/metabolics|total_metabolic_rate')
metGoal.setDivideByDisplacement(True)
metGoal.setDivideByMass(True)
problem.addGoal(metGoal)

#Prescribed average gait speed
speedGoal = osim.MocoAverageSpeedGoal('speed')
problem.addGoal(speedGoal)
speedGoal.set_desired_average_speed(6.0)

#Effort over distance
effortGoal = osim.MocoControlGoal('effort', 10)
problem.addGoal(effortGoal)
effortGoal.setExponent(2)
effortGoal.setDivideByDisplacement(True)
       
#Bounds
problem.setTimeBounds(0, [0.1, 0.5])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tilt/value', [-20*math.pi/180, 5*math.pi/180])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_tx/value', [0, 5])
problem.setStateInfo('/jointset/ground_pelvis/pelvis_ty/value', [0.75, 1.25])
problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [-40*math.pi/180, 80*math.pi/180])
problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [-40*math.pi/180, 80*math.pi/180])
problem.setStateInfo('/jointset/knee_l/knee_angle_l/value', [-140*math.pi/180, 0])
problem.setStateInfo('/jointset/knee_r/knee_angle_r/value', [-140*math.pi/180, 0])
problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [-30*math.pi/180, 30*math.pi/180])
problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [-30*math.pi/180, 30*math.pi/180])
problem.setStateInfo('/jointset/back/lumbar_extension/value', [-30, 0*math.pi/180])

#Configure the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_num_mesh_intervals(50)
solver.set_verbosity(2)
solver.set_optim_solver('ipopt')
solver.set_optim_convergence_tolerance(1e-4)
solver.set_optim_constraint_tolerance(1e-4)
solver.set_optim_max_iterations(1000)
solver.setGuessFile('MocoStudy_solution.sto') #fix up

#Solve
solutionPred = study.solve()

#Visualise
study.visualize(solutionPred)
