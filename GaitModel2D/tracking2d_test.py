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

# %% Import packages

import opensim as osim
import math
import xml.etree.ElementTree as et
# import numpy as np

# %% Scale the 2D model to match the experimental data

#### TODO: adjust model properties to be original Ong et al. size

#Spheres are maybe quite large in the first place?

#The models essentially have the same geometry, so the same scale set can be applied

#Set up the scale tool
scaleTool = osim.ScaleTool()

#Set participant mass (taken from model) TODO: automate
scaleTool.setSubjectMass(69.36)

#Set generic model file
genModelFileName = 'gait9dof18musc_Ong_et_al_Moco.osim'
scaleTool.getGenericModelMaker().setModelFileName(genModelFileName)

#Set options
scaleTool.getModelScaler().setPreserveMassDist(True)
scaleOrder = osim.ArrayStr(); scaleOrder.set(0,'manualScale')
scaleTool.getModelScaler().setScalingOrder(scaleOrder)

#Set the previous scale set for scaling
scaleTool.getModelScaler().setScaleSetFile('..\\ExpData\\scaleSet.xml')

#Set output files
scaleTool.getModelScaler().setOutputModelFileName('expData_2Dmodel.osim')

#Set marker placer to false
scaleTool.getMarkerPlacer().setApply(False)

#Run scale tool
scaleTool.run()

#### fails trying to adjust markers?

#The locations and size of the contact sphere parameters go unchanged with
#standard model scaling, so these need to be edited to ensure they are in
#an appropriate place. This can be done based on the scale set parameters
#for the representative bodies.

#Load in the scale set, parsed from the XML tree
xmlTree = et.parse('..\\ExpData\\scaleSet.xml')
xmlRoot = xmlTree.getroot()

#Load in the scaled model
scaledModel = osim.Model('expData_2Dmodel.osim')

#Create a dictionary with segments and scale factors to access for calculations
scaleFactors = {}

#Loop through segments and place  in dictionary
for segment in range(0,len(xmlRoot.findall('./ScaleSet/objects/Scale/segment'))):
    #Get current segment name
    currSegment = xmlRoot.findall('./ScaleSet/objects/Scale/segment')[segment].text
    #Get current scale name and parse to 0 x 3 array
    currScale = xmlRoot.findall('./ScaleSet/objects/Scale/scales')[segment].text
    currScale = str.split(currScale)
    scaleFactors[currSegment] = [float(currScale[0]),float(currScale[1]),float(currScale[2])]
    
#Get the 3D scale factors for the relevant bodies and average to scale the
#sphere radii. Note that these scale factors for each foot will be the same
#for heel and toes as it looks like same scale factors are applied.
heelSphereScale_r = sum(scaleFactors['calcn_r']) / 3
heelSphereScale_l = sum(scaleFactors['calcn_l']) / 3
toesSphereScale_r = sum(scaleFactors['toes_r']) / 3
toesSphereScale_l = sum(scaleFactors['toes_l']) / 3

#Scale the radii to each of their respective factors
#While accessing the spheres, also adjust their position based on the scale
#factor for the respective axes
#Create a list of the sheres to loop through and edit
sphereList = ['heel_r','toe1_r','toe2_r','heel_l','toe1_l','toe2_l']
#Loop through sphere list
for ss in range(0,len(sphereList)):
    #Get the current sphere
    currSphere = osim.ContactSphere.safeDownCast(scaledModel.getContactGeometrySet().get(sphereList[ss]))
    #Set the current scaling factor based on sphere name
    if '_r' in sphereList[ss]:
        if 'hallux' in sphereList[ss] or 'othertoes' in sphereList[ss]:
            scaleFac = toesSphereScale_r
            scaleBod = 'toes_r'
        else:
            scaleFac = heelSphereScale_r
            scaleBod = 'calcn_r'
    elif '_l' in sphereList[ss]:
        if 'hallux' in sphereList[ss] or 'othertoes' in sphereList[ss]:
            scaleFac = toesSphereScale_l
            scaleBod = 'toes_l'
        else:
            scaleFac = heelSphereScale_l
            scaleBod = 'calcn_l'
    #Rescale the radius
    currSphere.setRadius(currSphere.getRadius() * scaleFac)
    #Set new location
    newLoc = []
    newLoc.append(currSphere.getLocation().get(0) * scaleFactors[scaleBod][0])
    newLoc.append(currSphere.getLocation().get(1) * scaleFactors[scaleBod][1])
    newLoc.append(currSphere.getLocation().get(2) * scaleFactors[scaleBod][2])
    currSphere.setLocation(osim.Vec3(newLoc[0],newLoc[1],newLoc[2]))    
    
#Reset model name
scaledModel.setName('scaledModel_expData')

#Finalise connections and update scaled model
scaledModel.finalizeConnections()
scaledModel.printToXML('expData_2Dmodel.osim')

# %% Make shift RRA tracking simulation

# This section uses a torque driven simulation to adjust the kinematics so that
# the GRFs from the contact spheres match the experimental data

##### TODO: could try and incorporate mass and CoM parameters in this?

#Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName('torque_driven_rra_tracking')

##### Model processor with muscles removed, DON'T add actuators
##### Process model
##### Use a pre-specified actuator set (like RRAactiators) to add torque actuators to model
##### Set this model with the new actuators in the problem
##### Use a state and GRF tracking problem w/out pelvis actuators

#Construct a ModelProcessor and add it to the tool.
modelProcessor = osim.ModelProcessor('expData_2Dmodel.osim')
#Remove the muscles
modelProcessor.append(osim.ModOpRemoveMuscles())
# #Add CoordinateActuators to the model degrees-of-freedom. This
# modelProcessor.append(osim.ModOpAddReserves(500))

#Process to get model object
osimModel = modelProcessor.process()

#Append torque actuators to model
#Load in predefined force set
forceSet = osim.ForceSet('torque_actuators.xml')
#Loop through and append forces to model
for ff in range(0,forceSet.getSize()):
    osimModel.updForceSet().cloneAndAppend(osim.CoordinateActuator.safeDownCast(forceSet.get(ff)))

#Finalize model connections
osimModel.finalizeConnections()    

#Set model processor
track.setModel(osim.ModelProcessor(osimModel))

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
#Set each weight in a dictionary
weights = {'/jointset/ground_pelvis/pelvis_tx/value': 250,
           '/jointset/ground_pelvis/pelvis_ty/value': 50,
           '/jointset/ground_pelvis/pelvis_tilt/value': 25,
           '/jointset/back/lumbar_extension/value': 150,
           '/jointset/hip_r/hip_flexion_r/value': 100,
           '/jointset/knee_r/knee_angle_r/value': 250,
           '/jointset/ankle_r/ankle_angle_r/value': 25,
           '/jointset/hip_l/hip_flexion_l/value': 150,
           '/jointset/knee_l/knee_angle_l/value': 250,
           '/jointset/ankle_l/ankle_angle_l/value': 25}
#Initialise state weights object
stateWeights = osim.MocoWeightSet()
#Loop through and set weights
for kk in range(0,len(list(weights.keys()))):
    stateWeights.cloneAndAppend(osim.MocoWeight(list(weights.keys())[kk],weights[list(weights.keys())[kk]]))
#Set weights in tracking object
track.set_states_weight_set(stateWeights)

#Initialise study and problem
study = track.initialize()
problem = study.updProblem()

#Goals

#Effort
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
effort.setWeight(0.01)

# #Put a large weight on the pelvis CoordinateActuators, which act as the
# #residual, or 'hand-of-god', forces which we would like to keep as small
# #as possible.
# forceSet = osimModel.getForceSet()
# for ii in range(0,forceSet.getSize()):
#     forcePath = forceSet.get(ii).getAbsolutePathString()
#     if 'pelvis' in forcePath:
#         effort.setWeightForControl(forcePath,10)

#GRF contact tracking goal
contactTracking = osim.MocoContactTrackingGoal('contact',1)
contactTracking.setExternalLoadsFile('refGRF_2D.xml')
#Right foot
forceNamesRightFoot = osim.StdVectorString()
forceNamesRightFoot.append('/forceset/contactHeel_r')
forceNamesRightFoot.append('/forceset/contactToe1_r')
forceNamesRightFoot.append('/forceset/contactToe2_r')
trackRightGRF = osim.MocoContactTrackingGoalGroup(forceNamesRightFoot,'RightGRF')
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r')
contactTracking.addContactGroup(trackRightGRF)
#Left foot
forceNamesLeftFoot = osim.StdVectorString()
forceNamesLeftFoot.append('/forceset/contactHeel_l')
forceNamesLeftFoot.append('/forceset/contactToe1_l')
forceNamesLeftFoot.append('/forceset/contactToe2_l')
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
problem.setStateInfo('/jointset/back/lumbar_extension/value', [-30*math.pi/180, 0*math.pi/180])

#Configure the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_num_mesh_intervals(100)
solver.set_verbosity(2)
solver.set_optim_solver('ipopt')
solver.set_optim_convergence_tolerance(1e-4)
solver.set_optim_constraint_tolerance(1e-4)
solver.set_optim_max_iterations(2000)
# solver.setGuessFile('torque_driven_rra_tracking_solution.sto')

#Solve
solution = study.solve()

#Visualise solution
study.visualize(solution)

#Write tracked GRF to file
contact_r = osim.StdVectorString()
contact_l = osim.StdVectorString()
contact_r.append('/forceset/contactHeel_r')
contact_r.append('/forceset/contactToe1_r')
contact_r.append('/forceset/contactToe2_r')
contact_l.append('/forceset/contactHeel_l')
contact_l.append('/forceset/contactToe1_l')
contact_l.append('/forceset/contactToe2_l')
externalForcesTableFlat = osim.createExternalLoadsTableForGait(osimModel,
                                                               solution,
                                                               contact_r,
                                                               contact_l)
osim.writeTableToFile(externalForcesTableFlat,
                      'torque_driven_rra_tracking_GRF.sto')


####Solution looks good but doesn't solve in 1000 iterations
####Needed ~1770 iterations - GRFs are mostly good, although still seem to drop
#### off quite quickly rather than a smooth tail off at toe-off

#### Could filter ground reaction forces a little more - AP is noisy?

#Get kinematic states of the solution
outputTable = solution.exportToStatesTable()
fileAdapter = osim.STOFileAdapter()
fileAdapter.write(outputTable,'torque_drive_rra_tracking_kinematics.sto')

# %% [no good] Test out inverse simulation with RRA tracking kinematics

#Setup inverse object
inverse = osim.MocoInverse()
inverse.setName('expInverse')

#Load in the model and turn off the contact forces
osimModel = osim.Model('expData_2Dmodel.osim')
osimModel.updForceSet().get('contactHeel_r').set_appliesForce(False)
osimModel.updForceSet().get('contactToe1_r').set_appliesForce(False)
osimModel.updForceSet().get('contactToe2_r').set_appliesForce(False)
osimModel.updForceSet().get('contactHeel_l').set_appliesForce(False)
osimModel.updForceSet().get('contactToe1_l').set_appliesForce(False)
osimModel.updForceSet().get('contactToe2_l').set_appliesForce(False)
osimModel.finalizeConnections()

#Construct model processor to set on the tool
modelProcessor = osim.ModelProcessor(osimModel)
#Add external loads
modelProcessor.append(osim.ModOpAddExternalLoads('refGRF_2D.xml'))
#Ignore passive forces
modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
#Scale active force width
modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
#Ignore tendon compliance
modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
#Set on model
inverse.setModel(modelProcessor)

#Set the kinematics
tabProcessor = osim.TableProcessor('torque_drive_rra_tracking_kinematics.sto')
tabProcessor.append(osim.TabOpLowPassFilter(12))
inverse.setKinematics(tabProcessor)
inverse.set_initial_time(osim.Storage('torque_drive_rra_tracking_kinematics.sto').getFirstTime())
inverse.set_final_time(osim.Storage('torque_drive_rra_tracking_kinematics.sto').getLastTime())
inverse.set_mesh_interval(0.01)
inverse.set_kinematics_allow_extra_columns(True)

#Solve
inverseSolution = inverse.solve()
inverseSolution.getMocoSolution().write('inverse_solution.sto')

##### Doesn't work well

# %% [not tracking great] Basic tracking simulation of experimental data

#Setup tracking object
track = osim.MocoTrack()
track.setName('expDataTracking')

#Load in the model and turn off the contact forces
osimModel = osim.Model('expData_2Dmodel.osim')
osimModel.updForceSet().get('contactHeel_r').set_appliesForce(False)
osimModel.updForceSet().get('contactToe1_r').set_appliesForce(False)
osimModel.updForceSet().get('contactToe2_r').set_appliesForce(False)
osimModel.updForceSet().get('contactHeel_l').set_appliesForce(False)
osimModel.updForceSet().get('contactToe1_l').set_appliesForce(False)
osimModel.updForceSet().get('contactToe2_l').set_appliesForce(False)

#Setup metabolics
metabolics = osim.Bhargava2004Metabolics()
metabolics.setName('metabolics')
metabolics.set_use_smoothing(True)
#Get muscle list
muscleNames = list()
for mm in range(osimModel.getMuscles().getSize()):
    muscleNames.append(osimModel.getMuscles().get(mm).getName())
#Add muscles to metabolics
for mm in range(len(muscleNames)):
    metabolics.addMuscle(muscleNames[mm],
                         osim.DeGrooteFregly2016Muscle.safeDownCast(osimModel.getForceSet().get(muscleNames[mm])))
#Add to model
osimModel.addComponent(metabolics)

#Finalise model edits
osimModel.finalizeConnections()

#Process model
modelProcessor = osim.ModelProcessor(osimModel)
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
tabProcessor = osim.TableProcessor('torque_driven_rra_tracking_solution.sto')
tabProcessor.append(osim.TabOpLowPassFilter(12))
track.setStatesReference(tabProcessor)
track.set_states_global_tracking_weight(1.0)
track.set_allow_unused_references(True)
track.set_track_reference_position_derivatives(True)
track.set_apply_tracked_states_to_guess(True)
track.set_initial_time(osim.Storage('torque_driven_rra_tracking_solution.sto').getFirstTime())
track.set_final_time(osim.Storage('torque_driven_rra_tracking_solution.sto').getLastTime())

#Set individual state weights
#Set each weight in a dictionary
# weights = {'/jointset/ground_pelvis/pelvis_tx/value': 250,
#            '/jointset/ground_pelvis/pelvis_ty/value': 50,
#            '/jointset/ground_pelvis/pelvis_tilt/value': 25,
#            '/jointset/back/lumbar_extension/value': 150,
#            '/jointset/hip_r/hip_flexion_r/value': 100,
#            '/jointset/knee_r/knee_angle_r/value': 250,
#            '/jointset/ankle_r/ankle_angle_r/value': 25,
#            '/jointset/hip_l/hip_flexion_l/value': 150,
#            '/jointset/knee_l/knee_angle_l/value': 250,
#            '/jointset/ankle_l/ankle_angle_l/value': 25}
w = 0.001
weights = {'/jointset/ground_pelvis/pelvis_tx/value': (1/(1*0.1000))**2,
           '/jointset/ground_pelvis/pelvis_ty/value': (1/(2*0.1000))**2,
           '/jointset/ground_pelvis/pelvis_tilt/value': (1/(1*0.1745))**2,
           '/jointset/back/lumbar_extension/value': (1/(1*0.1745))**2,
           '/jointset/hip_r/hip_flexion_r/value': (1/(1*0.0647))**2,
           '/jointset/knee_r/knee_angle_r/value': (1/(1*0.0889))**2,
           '/jointset/ankle_r/ankle_angle_r/value': (1/(1*0.0574))**2,
           '/jointset/hip_l/hip_flexion_l/value': (1/(1*0.0647))**2,
           '/jointset/knee_l/knee_angle_l/value': (1/(1*0.0889))**2,
           '/jointset/ankle_l/ankle_angle_l/value': (1/(1*0.0574))**2,
           '/jointset/ground_pelvis/pelvis_tx/speed': w*(1/(1*0.1000))**2,
           '/jointset/ground_pelvis/pelvis_ty/speed': w*(1/(2*0.1000))**2,
           '/jointset/ground_pelvis/pelvis_tilt/speed': w*(1/(1*0.0585))**2,
           '/jointset/back/lumbar_extension/speed': w*(1/(1*0.1745))**2,
           '/jointset/hip_r/hip_flexion_r/speed': w*(1/(1*0.0647))**2,
           '/jointset/knee_r/knee_angle_r/speed': w*(1/(1*0.0889))**2,
           '/jointset/ankle_r/ankle_angle_r/speed': w*(1/(1*0.0574))**2,
           '/jointset/hip_l/hip_flexion_l/speed': w*(1/(1*0.0647))**2,
           '/jointset/knee_l/knee_angle_l/speed': w*(1/(1*0.0889))**2,
           '/jointset/ankle_l/ankle_angle_l/speed': w*(1/(1*0.0574))**2}
#Initialise state weights object
stateWeights = osim.MocoWeightSet()
#Loop through and set weights
for kk in range(0,len(list(weights.keys()))):
    stateWeights.cloneAndAppend(osim.MocoWeight(list(weights.keys())[kk],weights[list(weights.keys())[kk]]))
#Set weights in tracking object
track.set_states_weight_set(stateWeights)

#Initialise study and problem
study = track.initialize()
problem = study.updProblem()

#Goals

#Define the periodicity goal
#Keeping in mind that the simulation is of a full gait cycle, similar to that
#proposed in Umberger's Moco examples
periodicityGoal = osim.MocoPeriodicityGoal('symmetryGoal')
problem.addGoal(periodicityGoal)

#Set symmetric coordinate values (except for pelvis tx)
jointLabels = ['/jointset/ground_pelvis/pelvis_ty/value',
           '/jointset/ground_pelvis/pelvis_tilt/value',
           '/jointset/back/lumbar_extension/value',
           '/jointset/hip_r/hip_flexion_r/value',
           '/jointset/knee_r/knee_angle_r/value',
           '/jointset/ankle_r/ankle_angle_r/value',
           '/jointset/hip_l/hip_flexion_l/value',
           '/jointset/knee_l/knee_angle_l/value',
           '/jointset/ankle_l/ankle_angle_l/value',
           '/jointset/ground_pelvis/pelvis_tx/speed',
           '/jointset/ground_pelvis/pelvis_ty/speed',
           '/jointset/ground_pelvis/pelvis_tilt/speed',
           '/jointset/back/lumbar_extension/speed',
           '/jointset/hip_r/hip_flexion_r/speed',
           '/jointset/knee_r/knee_angle_r/speed',
           '/jointset/ankle_r/ankle_angle_r/speed',
           '/jointset/hip_l/hip_flexion_l/speed',
           '/jointset/knee_l/knee_angle_l/speed',
           '/jointset/ankle_l/ankle_angle_l/speed']

#Loop through joints and set as own periodicity pair
for jj in range(len(jointLabels)):

    #Set the pair for coordinate value
    pair = osim.MocoPeriodicityGoalPair(jointLabels[jj])
    
    #Add to the goal
    periodicityGoal.addStatePair(pair)
    
#Set periodic controls
#Muscles
muscleLabels = list()
for mm in range(osimModel.getMuscles().getSize()):
    muscleLabels.append(osimModel.getMuscles().get(mm).getAbsolutePathString())

#Loop through muscles and set as own periodicity pair
for mm in range(len(muscleLabels)):

    #Set the pair for coordinate value
    pair = osim.MocoPeriodicityGoalPair(muscleLabels[mm]+'/activation')
    
    #Add to the goal
    periodicityGoal.addStatePair(pair)
    
#Symmetric coordinate actuator control
periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair('/forceset/lumbarAct'))

# #Muscle excitations
# #Loop through muscles and set as own periodicity pair
# for mm in range(len(muscleLabels)):

#     #Set the pair for coordinate value
#     pair = osim.MocoPeriodicityGoalPair(muscleLabels[mm])
    
#     #Add to the goal
#     periodicityGoal.addStatePair(pair)
    
#Average speed goal
#Based on average speed of experimental data along x-axis
speedGoal = osim.MocoAverageSpeedGoal('speed')
problem.addGoal(speedGoal)
#Calculate average speed
#Load in kinematic storage
kinSto = osim.Storage('torque_driven_rra_tracking_solution.sto')
#Get state index for pelvis tx
txInd = kinSto.getStateIndex('/jointset/ground_pelvis/pelvis_tx/value')
#Get first data vector and extract pelvis tx value
startPos = kinSto.getStateVector(0).getData().get(txInd)
#Get end data vector and extract pelvis tx value
endPos = kinSto.getStateVector(kinSto.getSize()-1).getData().get(txInd)
#Calculate velocity
runSpeed = (endPos - startPos) / (kinSto.getLastTime() - kinSto.getFirstTime())
#Set speed in goal
speedGoal.set_desired_average_speed(runSpeed)

#Effort
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
effort.setWeight(0.1)

#Metabolics
metGoal = osim.MocoOutputGoal('met',0.1)
metGoal.setOutputPath('/metabolics|total_metabolic_rate')
metGoal.setDivideByDisplacement(True)
metGoal.setDivideByMass(True)
problem.addGoal(metGoal)

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
problem.setStateInfo('/jointset/back/lumbar_extension/value', [-30*math.pi/180, 0*math.pi/180])

#Configure the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_num_mesh_intervals(100)
solver.set_verbosity(2)
solver.set_optim_solver('ipopt')
solver.set_optim_convergence_tolerance(1e-2)
solver.set_optim_constraint_tolerance(1e-2)
solver.set_optim_max_iterations(1000)
solver.set_minimize_implicit_auxiliary_derivatives(True)
solver.set_implicit_auxiliary_derivatives_weight(0.001)
solver.resetProblem(problem)

#Solve
solution = study.solve()

#Visualise
study.visualize(solution)

# %% Tracking simulation with experimental 2D model

#Setup tracking object
track = osim.MocoTrack()
track.setName('expTracking')

#Get model
baseModel = osim.Model('expData_2Dmodel.osim')

# #Setup metabolics
# metabolics = osim.Bhargava2004Metabolics()
# metabolics.setName('metabolics')
# metabolics.set_use_smoothing(True)
# metabolics.addMuscle('hamstrings_r',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('hamstrings_r')))
# metabolics.addMuscle('hamstrings_l',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('hamstrings_l')))
# metabolics.addMuscle('bifemsh_r',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('bifemsh_r')))
# metabolics.addMuscle('bifemsh_l',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('bifemsh_l')))
# metabolics.addMuscle('glut_max_r',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('glut_max_r')))
# metabolics.addMuscle('glut_max_l',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('glut_max_l')))
# metabolics.addMuscle('iliopsoas_r',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('iliopsoas_r')))
# metabolics.addMuscle('iliopsoas_l',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('iliopsoas_l')))
# metabolics.addMuscle('rect_fem_r',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('rect_fem_r')))
# metabolics.addMuscle('rect_fem_l',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('rect_fem_l')))
# metabolics.addMuscle('vasti_r',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('vasti_r')))
# metabolics.addMuscle('vasti_l',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('vasti_l')))
# metabolics.addMuscle('gastroc_r',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('gastroc_r')))
# metabolics.addMuscle('gastroc_l',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('gastroc_l')))
# metabolics.addMuscle('soleus_r',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('soleus_r')))
# metabolics.addMuscle('soleus_l',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('soleus_l')))
# metabolics.addMuscle('tib_ant_r',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('tib_ant_r')))
# metabolics.addMuscle('tib_ant_l',osim.DeGrooteFregly2016Muscle.safeDownCast(baseModel.getForceSet().get('tib_ant_l')))
# baseModel.addComponent(metabolics)
# baseModel.finalizeConnections()

#Turn off contact forces to make it standard simulation
# baseModel.getComponent('contactHeel_r').set_appliesForce(False)
# baseModel.getComponent('contactToe1_r').set_appliesForce(False)
# baseModel.getComponent('contactToe2_r').set_appliesForce(False)
# baseModel.getComponent('contactHeel_l').set_appliesForce(False)
# baseModel.getComponent('contactToe1_l').set_appliesForce(False)
# baseModel.getComponent('contactToe2_l').set_appliesForce(False)
# baseModel.finalizeConnections()

#Process model
modelProcessor = osim.ModelProcessor(baseModel)
#Ignore passive forces
modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
#Scale active force width
modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
#Ignore tendon compliance
modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
#Add ground reaction external loads in lieu of a ground-contact model.
# modelProcessor.append(osim.ModOpAddExternalLoads('refGRF_2D.xml'))
# #Add reserves
# modelProcessor.append(osim.ModOpAddReserves(1))

#Setup of MocoTrack
track.setModel(modelProcessor)
tabProcessor = osim.TableProcessor('torque_drive_rra_tracking_kinematics.sto')
tabProcessor.append(osim.TabOpLowPassFilter(12))
track.setStatesReference(tabProcessor)
track.set_states_global_tracking_weight(10)
track.set_allow_unused_references(True)
track.set_track_reference_position_derivatives(True)
track.set_apply_tracked_states_to_guess(True)
track.set_initial_time(osim.Storage('torque_drive_rra_tracking_kinematics.sto').getFirstTime())
track.set_final_time(osim.Storage('torque_drive_rra_tracking_kinematics.sto').getLastTime())

#Set individual state weights
#Set each weight in a dictionary
weights = {'/jointset/ground_pelvis/pelvis_tx/value': 250,
            '/jointset/ground_pelvis/pelvis_ty/value': 250,
            '/jointset/ground_pelvis/pelvis_tilt/value': 25,
            '/jointset/back/lumbar_extension/value': 150,
            '/jointset/hip_r/hip_flexion_r/value': 100,
            '/jointset/knee_r/knee_angle_r/value': 250,
            '/jointset/ankle_r/ankle_angle_r/value': 25,
            '/jointset/hip_l/hip_flexion_l/value': 150,
            '/jointset/knee_l/knee_angle_l/value': 250,
            '/jointset/ankle_l/ankle_angle_l/value': 25}
#Initialise state weights object
stateWeights = osim.MocoWeightSet()
#Loop through and set weights
for kk in range(0,len(list(weights.keys()))):
    stateWeights.cloneAndAppend(osim.MocoWeight(list(weights.keys())[kk],weights[list(weights.keys())[kk]]))
#Set weights in tracking object
track.set_states_weight_set(stateWeights)

#Initialise study and problem
study = track.initialize()
problem = study.updProblem()

#Goals

# #Symmetry
# periodicityGoal = osim.MocoPeriodicityGoal('symmetryGoal')
# problem.addGoal(periodicityGoal)
# model = modelProcessor.process()
# model.initSystem()

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

# #Metabolics
# metGoal = osim.MocoOutputGoal('met',0.1)
# metGoal.setOutputPath('/metabolics|total_metabolic_rate')
# metGoal.setDivideByDisplacement(True)
# metGoal.setDivideByMass(True)
# problem.addGoal(metGoal)

#GRF contact tracking goal
contactTracking = osim.MocoContactTrackingGoal('contact',5)
contactTracking.setExternalLoadsFile('refGRF_2D.xml')
#Right foot
forceNamesRightFoot = osim.StdVectorString()
forceNamesRightFoot.append('/forceset/contactHeel_r')
forceNamesRightFoot.append('/forceset/contactToe1_r')
forceNamesRightFoot.append('/forceset/contactToe2_r')
trackRightGRF = osim.MocoContactTrackingGoalGroup(forceNamesRightFoot,'RightGRF')
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r')
contactTracking.addContactGroup(trackRightGRF)
#Left foot
forceNamesLeftFoot = osim.StdVectorString()
forceNamesLeftFoot.append('/forceset/contactHeel_l')
forceNamesLeftFoot.append('/forceset/contactToe1_l')
forceNamesLeftFoot.append('/forceset/contactToe2_l')
trackLeftGRF = osim.MocoContactTrackingGoalGroup(forceNamesLeftFoot,'LeftGRF')
trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l')
contactTracking.addContactGroup(trackLeftGRF)
#Properties
contactTracking.setProjection('plane')
contactTracking.setProjectionVector(osim.Vec3(0, 0, 1))
problem.addGoal(contactTracking)

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
problem.setStateInfo('/jointset/back/lumbar_extension/value', [-30*math.pi/180, 0*math.pi/180])

#Configure the solver
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_num_mesh_intervals(50)
solver.set_verbosity(2)
solver.set_optim_solver('ipopt')
solver.set_optim_convergence_tolerance(1e-4)
solver.set_optim_constraint_tolerance(1e-4)
solver.set_optim_max_iterations(1000)
#Set guess from basic tracking
# solver.setGuessFile('expDataTracking_solution.sto')

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

# #Print metabolics output
# print(str(round(solution.getObjectiveTerm('met'),3))+' J kg-1 m-1')

#Visualise
# study.visualize(full)
study.visualize(solution)

##### There's a clear need to add extra tracking on the trunk, as it otherwise 
##### can just have no torque input given the need to minimise it's actuator 
##### function. Another solution could be to zero weight the torque actuator and
##### and just focus on the trunk mechanics that minimise lower body functions?
##### The bounds also don't seem to be working on that trunk position, given a 
##### successful solve had the lumbar coordinate outside of the bounds (unless
##### the specified value is wrong...?) --- it is, it's not in radians and likely
##### isn't across a number of the problems!

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

# %% Tracking simulation

#Setup tracking object
track = osim.MocoTrack()

#Get model
baseModel = osim.Model('gait9dof18musc_Ong_et_al_Moco.osim')


#### TODO: adapt 2D model - tendon force incr. to 10% for plantarflexors

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
# baseModel.getComponent('contactHeel_r').set_appliesForce(False)
# baseModel.getComponent('contactToe1_r').set_appliesForce(False)
# baseModel.getComponent('contactToe2_r').set_appliesForce(False)
# baseModel.getComponent('contactHeel_l').set_appliesForce(False)
# baseModel.getComponent('contactToe1_l').set_appliesForce(False)
# baseModel.getComponent('contactToe2_l').set_appliesForce(False)
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
# modelProcessor.append(osim.ModOpAddExternalLoads('refGRF_2D.xml'))
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
#Set each weight in a dictionary
weights = {'/jointset/ground_pelvis/pelvis_tx/value': 100,
           '/jointset/ground_pelvis/pelvis_ty/value': 25,
           '/jointset/ground_pelvis/pelvis_tilt/value': 25,
           '/jointset/back/lumbar_extension/value': 150,
           '/jointset/hip_r/hip_flexion_r/value': 250,
           '/jointset/knee_r/knee_angle_r/value': 500,
           '/jointset/ankle_r/ankle_angle_r/value': 75,
           '/jointset/hip_l/hip_flexion_l/value': 250,
           '/jointset/knee_l/knee_angle_l/value': 500,
           '/jointset/ankle_l/ankle_angle_l/value': 75}
#Initialise state weights object
stateWeights = osim.MocoWeightSet()
#Loop through and set weights
for kk in range(0,len(list(weights.keys()))):
    stateWeights.cloneAndAppend(osim.MocoWeight(list(weights.keys())[kk],weights[list(weights.keys())[kk]]))
#Set weights in tracking object
track.set_states_weight_set(stateWeights)

# #These are based off Ross Miller's UMocoD project for now
# stateWeights = osim.MocoWeightSet()
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tx/value',(1/(1*0.2000))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_ty/value',(1/(2*0.1000))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tilt/value',(1/(1*0.1745))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/value',(1/(1*0.1745))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/value',(1/(1*0.0647))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_r/knee_angle_r/value',(1/(1*0.0889))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/value',(1/(1*0.0574))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/value',(1/(1*0.0647))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_l/knee_angle_l/value',(1/(1*0.0889))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/value',(1/(1*0.0574))**2))
# w = 0.001 #Scale the generalized speed tracking errors by this constant
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tx/speed',w*(0/(1*0.1000))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_ty/speed',w*(0/(2*0.1000))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ground_pelvis/pelvis_tilt/speed',w*(0/(1*0.0585))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/speed',w*(0/(1*0.1745))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/speed',w*(1/(1*0.0647))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_r/knee_angle_r/speed',w*(1/(1*0.0889))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/speed',w*(1/(1*0.0574))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/speed',w*(1/(1*0.0647))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_l/knee_angle_l/speed',w*(1/(1*0.0889))**2))
# stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/speed',w*(1/(1*0.0574))**2))
# track.set_states_weight_set(stateWeights)

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
problem.setStateInfo('/jointset/back/lumbar_extension/value', [-30*math.pi/180, 0*math.pi/180])

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


# #Test out locking torso joint
# baseModel.getCoordinateSet().get('lumbar_extension').set_locked(True)

#Finalise model connections
baseModel.finalizeConnections()

#Process model
modelProcessor = osim.ModelProcessor(baseModel)
#Ignore passive forces
modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
#Scale active force width
modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
#Ignore tendon compliance
modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
# #Convert locked joints (now the back) to weld
# modelProcessor.append(osim.ModOpReplaceJointsWithWelds())

### Consider passive force operators if using?

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
# problem.setStateInfo('/jointset/back/lumbar_extension/value', [-30*math.pi/180, 0*math.pi/180])
problem.setStateInfo('/jointset/back/lumbar_extension/value', [0, 0*math.pi/180])

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
