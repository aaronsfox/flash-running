# This example attempts to run some Moco problems in Python all from the one
# folder (similar to the paackaged Moco examples) to see whether the folder 
# navigation is what is causing issues.
#
# See the README.MD file next to this for more info.

import os
import opensim as osim
import osimHelper
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xml.etree.ElementTree as et

# %% Model scaling

#Set up the scale tool
scaleTool = osim.ScaleTool()

#Set participant mass (64.3kg)
scaleTool.setSubjectMass(64.3)

#Set generic model file
genModelFileName = os.getcwd()+'\\pfjLoading_GenericOsimModel_LaiArnold.osim'
scaleTool.getGenericModelMaker().setModelFileName(genModelFileName)

#Set the measurement set
measurementSetObject = osim.OpenSimObject.makeObjectFromFile('scaleMeasurementSet.xml')
measurementSet = osim.MeasurementSet.safeDownCast(measurementSetObject)
scaleTool.getModelScaler().setMeasurementSet(measurementSet)

#Set scale tasks
taskSet = osim.IKTaskSet('scaleTasks.xml')
for k in range(0,taskSet.getSize()-1):
    scaleTool.getMarkerPlacer().getIKTaskSet().adoptAndAppend(taskSet.get(k))

#Set marker file

#Add the virtual markers and create a new .trc file to use in scaling
osimHelper.addVirtualMarkersStatic('static.trc','staticVirtualMarkers.trc')

#Place in scale tool
scaleTool.getMarkerPlacer().setMarkerFileName('staticVirtualMarkers.trc')
scaleTool.getModelScaler().setMarkerFileName('staticVirtualMarkers.trc')

#Set options
scaleTool.getModelScaler().setPreserveMassDist(True)
scaleOrder = osim.ArrayStr(); scaleOrder.set(0,'measurements')
scaleTool.getModelScaler().setScalingOrder(scaleOrder)

#Set time ranges
timeRange = osim.ArrayDouble()
timeRange.set(0,0.5); timeRange.set(1,1.5)
scaleTool.getMarkerPlacer().setTimeRange(timeRange)
scaleTool.getModelScaler().setTimeRange(timeRange)

#Set output files
scaleTool.getModelScaler().setOutputModelFileName('scaledModel.osim')
scaleTool.getModelScaler().setOutputScaleFileName('scaleSet.xml')

#Set marker adjuster parameters
scaleTool.getMarkerPlacer().setOutputMotionFileName('static_motion.mot')
scaleTool.getMarkerPlacer().setOutputModelFileName('scaledModelAdjusted.osim')

#Save and run scale tool
scaleTool.printToXML('scaleSetup.xml')
scaleTool.run()

#Load in scaled model
scaledModel = osim.Model('scaledModelAdjusted.osim')

#The locations and size of the contact sphere parameters go unchanged with
#standard model scaling, so these need to be edited to ensure they are in
#an appropriate place. This can be done based on the scale set parameters
#for the representative bodies.

#Load in the scale set, parsed from the XML tree
xmlTree = et.parse('scaleSet.xml')
xmlRoot = xmlTree.getroot()

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
sphereList = ['heel_r','mh1_r','mh3_r','mh5_r','hallux_r','othertoes_r',
              'heel_l','mh1_l','mh3_l','mh5_l','hallux_l','othertoes_l']
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
scaledModel.setName('scaledModel')

#Finalise connections and update scaled model
scaledModel.finalizeConnections()
scaledModel.printToXML('scaledModelAdjusted.osim')

#Scale muscle strength based on linear function presented in Handsfield
#et al. (2014). This uses some convenience functions that are packaged
#with the Rajagopal et al. (2016) gait model. Note that the height of
#the generic model is 1.700; the current participant height is 1.759
osimHelper.scaleOptimalForceSubjectSpecific(genericModelFileName = 'pfjLoading_GenericOsimModel_LaiArnold.osim',
                                            scaledModelFileName = 'scaledModelAdjusted.osim',
                                            genericHeight = 1.700, scaledHeight = 1.759,
                                            outputModelFileName = 'scaledModelMuscle.osim')
#Load in new scaled model
scaledModelMuscle = osim.Model('scaledModelMuscle.osim')

# %% Simulation parameters

#Identify time parameters for dynamic simulations
[startTime,endTime] = osimHelper.getHalfGaitCycle('Jog05_grf.mot')

# %% Replicate the torque driven marker tracking Moco example

##### NOTE: the marker tracking doesn't necessarily offer a great benefit over
##### inverse kinematics, so we may as well use IK for its better speed anyway...

# def torqueDrivenMarkerTracking(startTime = 0.0, endTime = 1.0, meshInterval = 0.05,
#                                scaledModelFile = None, grfXml = None, trcFile = None,
#                                ikTasksXml = None):
    
#     #Check for inputs
#     if scaledModelFile is None:
#         raise ValueError('A .osim model filename is required.')
#     if grfXml is None:
#         raise ValueError('A .xml file of ground reaction forces is required.')
#     if trcFile is None:
#         raise ValueError('A .trc filename is required.')
#     if ikTasksXml is None:
#         raise ValueError('A .xml file for marker task weights is required.')
    
#     #Create and name an instance of the MocoTrack tool.
#     track = osim.MocoTrack()
#     track.setName('torque_driven_marker_tracking')

#     #Construct a ModelProcessor and add it to the tool. ModelProcessors
#     #accept a base model and allow you to easily modify the model by appending
#     #ModelOperators. Operations are performed in the order that they are
#     #appended to the model.
#     #Create the base Model by passing in the model file.
#     modelProcessor = osim.ModelProcessor(scaledModelFile)
#     #Add ground reaction external loads in lieu of a ground-contact model.
#     modelProcessor.append(osim.ModOpAddExternalLoads(grfXml))
#     #Remove all the muscles in the model's ForceSet.
#     modelProcessor.append(osim.ModOpRemoveMuscles())
#     #Add CoordinateActuators to the model degrees-of-freedom.
#     modelProcessor.append(osim.ModOpAddReserves(300))
#     track.setModel(modelProcessor)

#     #Use this convenience function to set the MocoTrack markers reference
#     #directly from a TRC file. By default, the markers data is filtered at
#     #6 Hz and if in millimeters, converted to meters.
    
#     ##### TODO: change this filter frequency...
#     ##### TODO: potentially add virtual markers to dynamic trial...
    
#     track.setMarkersReferenceFromTRC(trcFile)

#     #There is marker data in the 'marker_trajectories.trc' associated with
#     #model markers that no longer exists (i.e. markers on the arms). Set this
#     #flag to avoid an exception from being thrown.
#     track.set_allow_unused_references(True)

#     # Increase the global marker tracking weight, which is the weight
#     # associated with the internal MocoMarkerTrackingCost term.
#     track.set_markers_global_tracking_weight(10)

#     #Set the individual marker tracking weights based on the IK tasks
    
#     #Load in the IK task set
#     ikTaskSet = osim.IKTaskSet(ikTasksXml)
    
#     #Create the marker weight set
#     markerWeights = osim.MocoWeightSet()
    
#     #Loop through task set and add the markers that are set to apply
#     for ii in range(0,ikTaskSet.getSize()-1):
#         if ikTaskSet.get(ii).getApply():
#             #Get marker name
#             markerName = ikTaskSet.get(ii).getName()
#             #Get weight
#             markerWeight = ikTaskSet.get(ii).getWeight()
#             #Set in weight set
#             markerWeights.cloneAndAppend(osim.MocoWeight(markerName,markerWeight))
    
#     #Set marker weights in tracking problem
#     track.set_markers_weight_set(markerWeights)

#     # Initial time, final time, and mesh interval. The number of mesh points
#     # used to discretize the problem is computed internally using these values.
#     track.set_initial_time(startTime)
#     track.set_final_time(endTime)
#     track.set_mesh_interval(meshInterval)

#     # Solve! The boolean argument indicates to visualize the solution.
#     solution = track.solve(True)
    
#Add the virtual torso, pelvis and hip joint markers to the .trc file
osimHelper.addVirtualMarkersDynamic(staticTRC = 'staticVirtualMarkers.trc',
                                    dynamicTRC = 'Jog05.trc',
                                    outputTRC = 'Jog05_virtualMarkers.trc')
    
# #Run the torque driven marker tracking problem
# torqueDrivenMarkerTracking(startTime = startTime,endTime = endTime,meshInterval = 0.05,
#                            scaledModelFile = 'scaledModelAdjusted.osim',
#                            grfXml = 'Jog05_grf.xml', trcFile = 'Jog05_virtualMarkers.trc',
#                            ikTasksXml = 'ikTasks.xml')

##### Marker tracking looks fairly good - compare to IK?
##### Still getting visualiser error though - weird considering it works for Moco examples?
    
# %% Inverse kinematics

##### TODO: could set this in osimHelper function

#Initialise IK tool
ikTool = osim.InverseKinematicsTool()

#Set model
ikTool.setModel(scaledModel)

#Set task set
ikTool.set_IKTaskSet(osim.IKTaskSet('ikTasks.xml'))

#Set marker file
ikTool.set_marker_file('Jog05_virtualMarkers.trc')

#Set times
ikTool.setStartTime(startTime)
ikTool.setEndTime(endTime)

#Set output filename
ikTool.set_output_motion_file('ikResults.mot')

#Run IK
ikTool.run()

# %% Compare marker tracking vs. IK solution

##### NOTE: without marker tracking results, we can just plot IK

#Read in IK results
#Need to skip first 10 rows of file
##### TODO: automate process of identifying 'endheader' string
ik_df = pd.read_table('ikResults.mot',
                      skiprows = list(range(0,10)))

# #Read in marker tracking results
# #Need to skip first 17 rows
# markerTracking_df = pd.read_table('torque_driven_marker_tracking_solution.sto',
#                                   skiprows = list(range(0,17)))
# #Convert marker tracking results to degrees
# for cc in range(0,len(list(markerTracking_df))):
#     #Get current column header
#     currHeader = list(markerTracking_df)[cc]
#     #Check if it contains the jointset and value keywords
#     if 'jointset' in currHeader and 'value' in currHeader:
#         #Check if it is a pelvis translation value not to convert
#         if '_tx' not in currHeader and '_ty' not in currHeader and '_tz' not in currHeader:
#             #Convert column to degrees
#             markerTracking_df[currHeader] = np.rad2deg(markerTracking_df[currHeader])

#Create list of coordinates to plot
coordPlot = ['lumbar_extension','lumbar_bending','lumbar_rotation',
             'pelvis_tx','pelvis_ty','pelvis_tz',
             'pelvis_list','pelvis_tilt','pelvis_rotation',
             'hip_flexion_r','hip_adduction_r','hip_rotation_r',
             'knee_angle_r','ankle_angle_r']
coordPlotName = ['Lumbar Extension','Lumbar Bending','Lumbar Rotation',
                 'Pelvis X Translation','Pelvis Y Translation','Pelvis Z Translation',
                 'Pelvis List','Pelvis Tilt','Pelvis Rotation',
                 'R. Hip Flexion','R. Hip Adduction','R. Hip Rotation',
                 'R. Knee Angle','R. Ankle Angle']

#Plot results

#Set subplot shape
fig, axs = plt.subplots(5,3,figsize=(8,6))
#Set a list to access the axes coordinates in a loop
whichAx = [[0,0],[0,1],[0,2],
           [1,0],[1,1],[1,2],
           [2,0],[2,1],[2,2],
           [3,0],[3,1],[3,2],
           [4,0],[4,1],[4,2]]
#Loop through coordinates
for cc in range(0,len(coordPlot)):
    #Plot IK data
    ik_df.plot.line(x = 'time', y = coordPlot[cc],
                    legend = False, xlim = [startTime,endTime],
                    color = 'black', linestyle = '-',
                    linewidth = 1.5, ax = axs[whichAx[cc][0],whichAx[cc][1]])
    # #Plot marker tracking data
    # #Get the absolute path to the current coordinate
    # markerTracking_df.plot.line(x = 'time', y = scaledModel.updCoordinateSet().get(coordPlot[cc]).getAbsolutePathString()+'/value',
    #                             legend = False, xlim = [startTime,endTime],
    #                             color = 'red', linestyle = '--',
    #                             linewidth = 1.5, ax = axs[whichAx[cc][0],whichAx[cc][1]])
    #Set title
    axs[whichAx[cc][0],whichAx[cc][1]].set_title(coordPlotName[cc])
#Turn off the final extra subplot axes
axs[4,2].axis('off')
#Set tight layout for figure
fig.tight_layout()
#Show figure
plt.show()

##### Marker tracking follows IK pretty well. Slight differences in some areas
##### Some Y and Z rotations are a bit noisy, probably improve with shorter mesh interval
    ##### Shorter mesh interval doesn't necessarily work better, potentially given sample rate of data?
##### Marker tracking is cumbersome and can sometimes have weird results...cost-benefit?

# %% Run an inverse torque driven problem

#First, convert IK results to a states file
osimHelper.kinematicsToStates(kinematicsFileName = 'ikResults.mot',
                              osimModelFileName = 'scaledModelAdjusted.osim',
                              outputFileName = 'ikResults_states.sto',
                              inDegrees = True, outDegrees = False)

#Use MocoInverse to identify the torque controls that replicate the experimental data
##### TODO: could wrap this as a function with relevant inputs
def solveMocoInverse():

    #Construct the MocoInverse tool.
    inverse = osim.MocoInverse()
    inverse.setName('inverseTorqueTracking')

    #Construct a ModelProcessor and set it on the tool. The muscles are removed
    #and reserve actuators added to generate a torque driven solution
    modelProcessor = osim.ModelProcessor('scaledModelMuscle.osim')
    modelProcessor.append(osim.ModOpAddExternalLoads('Jog05_grf.xml'))
    modelProcessor.append(osim.ModOpRemoveMuscles())
    modelProcessor.append(osim.ModOpAddReserves(300))
    inverse.setModel(modelProcessor)
    
    #Construct a TableProcessor of the coordinate data and pass it to the
    #inverse tool. TableProcessors can be used in the same way as
    #ModelProcessors by appending TableOperators to modify the base table.
    #A TableProcessor with no operators, as we have here, simply returns the
    #base table.
    inverse.setKinematics(osim.TableProcessor('ikResults_states.sto'))

    # Initial time, final time, and mesh interval.
    inverse.set_initial_time(osim.Storage('ikResults_states.sto').getFirstTime())
    inverse.set_final_time(osim.Storage('ikResults_states.sto').getLastTime())
    inverse.set_mesh_interval(0.02)

    # By default, Moco gives an error if the kinematics contains extra columns.
    # Here, we tell Moco to allow (and ignore) those extra columns.
    inverse.set_kinematics_allow_extra_columns(True)

    # Solve the problem and write the solution to a Storage file.
    solution = inverse.solve()
    
    #Return the solution as an object that we can use
    return solution

#Run solution
inverseTorqueSolution = solveMocoInverse()

#Write as a MocoSolution
inverseTorqueSolution.getMocoSolution().write('inverseTorqueSolution_fromIK.sto')

#Remove the originally written solution
os.remove('MocoStudy_solution.sto')

# %% Run a torque driven tracking problem, with experimental GRFs

##### TODO: could wrap this as a function with relevant inputs
def torqueDrivenStateTracking():

    #Create and name an instance of the MocoTrack tool.
    track = osim.MocoTrack()
    track.setName('torqueDrivenStateTracking')

    #Construct a ModelProcessor and set it on the tool. Remove the muscles and
    #add reserve actuators
    modelProcessor = osim.ModelProcessor('scaledModelMuscle.osim')
    modelProcessor.append(osim.ModOpAddExternalLoads('Jog05_grf.xml'))
    modelProcessor.append(osim.ModOpRemoveMuscles())
    modelProcessor.append(osim.ModOpAddReserves(300))
    track.setModel(modelProcessor)

    #Construct a TableProcessor of the coordinate data and pass it to the 
    #tracking tool. TableProcessors can be used in the same way as
    #ModelProcessors by appending TableOperators to modify the base table.
    #A TableProcessor with no operators, as we have here, simply returns the
    #base table.
    track.setStatesReference(osim.TableProcessor('ikResults_states.sto'))
    track.set_states_global_tracking_weight(5)

    ##### TODO: add more specific weights for different state coordinates like in RRA

    #This setting allows extra data columns contained in the states
    #reference that don't correspond to model coordinates.
    track.set_allow_unused_references(True)

    # Since there is only coordinate position data the states references, this
    # setting is enabled to fill in the missing coordinate speed data using
    # the derivative of splined position data.
    track.set_track_reference_position_derivatives(True)

    # Initial time, final time, and mesh interval.
    track.set_initial_time(osim.Storage('ikResults_states.sto').getFirstTime())
    track.set_final_time(osim.Storage('ikResults_states.sto').getLastTime())
    track.set_mesh_interval(0.08)

    #Instead of calling solve(), call initialize() to receive a pre-configured
    #MocoStudy object based on the settings above. Use this to customize the
    #problem beyond the MocoTrack interface.
    study = track.initialize()

    #Get a reference to the MocoControlCost that is added to every MocoTrack
    #problem by default.
    problem = study.updProblem()
    effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
    effort.setWeight(1)

    # # Put a large weight on the pelvis CoordinateActuators, which act as the
    # # residual, or 'hand-of-god', forces which we would like to keep as small
    # # as possible.
    # model = modelProcessor.process()
    # model.initSystem()
    # forceSet = model.getForceSet()
    # for ii in range(forceSet.getSize()):
    #     forcePath = forceSet.get(ii).getAbsolutePathString()
    #     if 'pelvis' in str(forcePath):
    #         effort.setWeightForControl(forcePath, 10)
            
    #Get solver and set the mesh interval
    solver = study.initCasADiSolver()
    #50 mesh intervals for half gait cycle recommended, keep smaller for now
    #19 will produce the same as inverse solution above
    # solver.set_num_mesh_intervals(19)
    
    #Set solver parameters
    solver.set_optim_constraint_tolerance(1e-3)
    solver.set_optim_convergence_tolerance(1e-3)
    
    #Set guess from inverse solution using the convenience function to transfer
    #an inverse solution to a solver guess
    osimHelper.inverseSolutionToTrackGuess(guessFile = 'inverseTorqueSolution_fromIK.sto',
                                           mocoSolver = solver)
    
    # Solve and visualize.
    solution = study.solve()

    #Return the solution as an object that we can use
    return solution

#Run solution
stateTrackingTorqueSolution = torqueDrivenStateTracking()

##### TODO: plots to compare IK to torque drive state tracking

##### Torque tracking is still quite time consuming, perhaps because of the
##### infinite space for controls to lie within

# %% Use MocoInverse to solve for muscle driven solution

#Use MocoInverse to identify the muscle controls that replicate the experimental data
##### TODO: could wrap this as a function with relevant inputs
def solveMocoInverseMuscle():

    #Construct the MocoInverse tool.
    inverse = osim.MocoInverse()
    inverse.setName('inverseMuscleTracking')

    #Construct a ModelProcessor and set it on the tool.
    #Currently the coordinate actuators for the pelvis, along with the reserve
    #actuators are fairly generally set, and not weighted at all in cost function.
    #These actuators could be more carefully considered to generate an appropriate
    #muscle driven simulation. For example, if they max and min control to 1 - the 
    #pelvis actuators may not produce enough torque.
    modelProcessor = osim.ModelProcessor('scaledModelMuscle.osim')
    modelProcessor.append(osim.ModOpAddExternalLoads('Jog05_grf.xml'))
    modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
    modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
    # Only valid for DeGrooteFregly2016Muscles.
    # modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
    # Only valid for DeGrooteFregly2016Muscles.
    modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
    modelProcessor.append(osim.ModOpAddReserves(2))
    inverse.setModel(modelProcessor)
    
    #Construct a TableProcessor of the coordinate data and pass it to the
    #inverse tool. TableProcessors can be used in the same way as
    #ModelProcessors by appending TableOperators to modify the base table.
    #A TableProcessor with no operators, as we have here, simply returns the
    #base table.
    inverse.setKinematics(osim.TableProcessor('ikResults_states.sto'))

    #Initial time, final time, and mesh interval.
    inverse.set_initial_time(osim.Storage('ikResults_states.sto').getFirstTime())
    inverse.set_final_time(osim.Storage('ikResults_states.sto').getLastTime())
    inverse.set_mesh_interval(0.02)

    # By default, Moco gives an error if the kinematics contains extra columns.
    # Here, we tell Moco to allow (and ignore) those extra columns.
    inverse.set_kinematics_allow_extra_columns(True)

    # Solve the problem and write the solution to a Storage file.
    solution = inverse.solve()
    
    #Return the solution as an object that we can use
    return solution

#Run solution
inverseMuscleSolution = solveMocoInverseMuscle()

#Write as a MocoSolution
inverseMuscleSolution.getMocoSolution().unseal().write('inverseMuscleSolution_fromIK.sto')

#Remove the originally written solution
os.remove('MocoStudy_solution.sto')

##### Inverse solution not converging well, and has high residuals
##### Perhaps RRA is still appropriate here...

# %% Run a muscle drive tracking problem

##### TODO: muscle driven state tracking is very, very slow (44 minutes for 4 iterations!) to begin with
##### Consider things to speed up:
    ##### - Finite differences?
    ##### - Muscles on one leg?
    ##### - Initial guess from inverse solution?

##### TODO: could wrap this as a function with relevant inputs
def muscleDrivenStateTracking():

    #Create and name an instance of the MocoTrack tool.
    track = osim.MocoTrack()
    track.setName('muscleDrivenStateTracking')

    #Construct a ModelProcessor and set it on the tool.
    #Currently the coordinate actuators for the pelvis, along with the reserve
    #actuators are fairly generally set, and not weighted at all in cost function.
    #These actuators could be more carefully considered to generate an appropriate
    #muscle driven simulation. For example, if they max and min control to 1 - the 
    #pelvis actuators may not produce enough torque.
    modelProcessor = osim.ModelProcessor('scaledModelMuscle.osim')
    modelProcessor.append(osim.ModOpAddExternalLoads('Jog05_grf.xml'))
    modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
    modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
    # Only valid for DeGrooteFregly2016Muscles.
    modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
    # Only valid for DeGrooteFregly2016Muscles.
    modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
    modelProcessor.append(osim.ModOpAddReserves(2))
    track.setModel(modelProcessor)

    #Construct a TableProcessor of the coordinate data and pass it to the 
    #tracking tool. TableProcessors can be used in the same way as
    #ModelProcessors by appending TableOperators to modify the base table.
    #A TableProcessor with no operators, as we have here, simply returns the
    #base table.
    track.setStatesReference(osim.TableProcessor('ikResults_states.sto'))
    track.set_states_global_tracking_weight(5)

    ##### TODO: add more specific weights for different state coordinates like in RRA

    #This setting allows extra data columns contained in the states
    #reference that don't correspond to model coordinates.
    track.set_allow_unused_references(True)

    # Since there is only coordinate position data the states references, this
    # setting is enabled to fill in the missing coordinate speed data using
    # the derivative of splined position data.
    track.set_track_reference_position_derivatives(True)

    # Initial time, final time, and mesh interval.
    track.set_initial_time(osim.Storage('ikResults_states.sto').getFirstTime())
    track.set_final_time(osim.Storage('ikResults_states.sto').getLastTime())
    track.set_mesh_interval(0.08)

    #Instead of calling solve(), call initialize() to receive a pre-configured
    #MocoStudy object based on the settings above. Use this to customize the
    #problem beyond the MocoTrack interface.
    study = track.initialize()

    #Get a reference to the MocoControlCost that is added to every MocoTrack
    #problem by default.
    problem = study.updProblem()
    effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
    effort.setWeight(1)

    # # Put a large weight on the pelvis CoordinateActuators, which act as the
    # # residual, or 'hand-of-god', forces which we would like to keep as small
    # # as possible.
    # model = modelProcessor.process()
    # model.initSystem()
    # forceSet = model.getForceSet()
    # for ii in range(forceSet.getSize()):
    #     forcePath = forceSet.get(ii).getAbsolutePathString()
    #     if 'pelvis' in str(forcePath):
    #         effort.setWeightForControl(forcePath, 10)
            
    #Get solver and set the mesh interval
    solver = study.initCasADiSolver()
    #50 mesh intervals for half gait cycle recommended, keep smaller for now
    #19 will produce the same as inverse solution above
    # solver.set_num_mesh_intervals(19)
    
    #Set solver parameters
    solver.set_optim_constraint_tolerance(1e-3)
    solver.set_optim_convergence_tolerance(1e-3)
    
    # Solve and visualize.
    solution = study.solve()

    #Return the solution as an object that we can use
    return solution

#Run solution
stateTrackingMuscleSolution = muscleDrivenStateTracking()





# %% Replicate the muscle driven state tracking example...but make it torque driven...

#Replicate moco track example
##### TODO: set in function

# Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName("torque_driven_state_tracking")

#Probably need to turn off contact element forceset elements
# trackModel = osim.Model("scaledModelAdjusted.osim")
# for ii in range(0,trackModel.updForceSet().getSize()-1):
#     #Get current force name
#     currForce = trackModel.updForceSet().get(ii).getName()
#     #Turn off if a contact element
#     if 'contact' in currForce:
#         trackModel.updForceSet().get(ii).set_appliesForce(False)
##### TODO: keep the spheres, but these haven't been adjusted yet...

# #Finalize
# trackModel.finalizeConnections()

# Construct a ModelProcessor and set it on the tool. The default
# muscles in the model are replaced with optimization-friendly
# DeGrooteFregly2016Muscles, and adjustments are made to the default muscle
# parameters.
modelProcessor = osim.ModelProcessor("scaledModelAdjusted.osim")
modelProcessor.append(osim.ModOpAddExternalLoads("Jog05_grf.xml"))
# modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
# modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
#Remove all the muscles in the model's ForceSet.
modelProcessor.append(osim.ModOpRemoveMuscles())
# Only valid for DeGrooteFregly2016Muscles.
# modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
# Only valid for DeGrooteFregly2016Muscles.
# modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
#Add reserves
modelProcessor.append(osim.ModOpAddReserves(300))
track.setModel(modelProcessor)

# Construct a TableProcessor of the coordinate data and pass it to the 
# tracking tool. TableProcessors can be used in the same way as
# ModelProcessors by appending TableOperators to modify the base table.
# A TableProcessor with no operators, as we have here, simply returns the
# base table.
track.setStatesReference(osim.TableProcessor("coordinates.sto"))
track.set_states_global_tracking_weight(5)

# This setting allows extra data columns contained in the states
# reference that don't correspond to model coordinates.
track.set_allow_unused_references(True)

# Since there is only coordinate position data the states references, this
# setting is enabled to fill in the missing coordinate speed data using
# the derivative of splined position data.
track.set_track_reference_position_derivatives(True)

# Initial time, final time, and mesh interval.
# The IK data doesn't match the GRF data, so specific times need ti be taken from this
track.set_initial_time(osim.Storage('coordinates.sto').getFirstTime())
track.set_final_time(osim.Storage('coordinates.sto').getLastTime())
# track.set_mesh_interval(0.05)

# Instead of calling solve(), call initialize() to receive a pre-configured
# MocoStudy object based on the settings above. Use this to customize the
# problem beyond the MocoTrack interface.
study = track.initialize()

# Get a reference to the MocoControlCost that is added to every MocoTrack
# problem by default.
problem = study.updProblem()
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal("control_effort"))

# Put a large weight on the pelvis CoordinateActuators, which act as the
# residual, or 'hand-of-god', forces which we would like to keep as small
# as possible.
model = modelProcessor.process()
model.initSystem()
forceSet = model.getForceSet()
for i in range(forceSet.getSize()):
    forcePath = forceSet.get(i).getAbsolutePathString()
    if 'pelvis' in str(forcePath):
        effort.setWeightForControl(forcePath, 10)

#Set more specific weights for each coordinate in states tracking weight
#Get a reference to the already added state tracking goal
stateTracking = osim.MocoStateTrackingGoal.safeDownCast(problem.updGoal("state_tracking"))
#Loop through model coordinates
for cc in range(0,model.updCoordinateSet().getSize()-1):
    #Check if current coordinate is unlocked to proceed
    if not model.updCoordinateSet().get(cc).get_locked():
        #Get the absolute path string for the current coordinate and add the value key
        currState = model.updCoordinateSet().get(cc).getAbsolutePathString()+'/value'
        #Set the weight for the state
        #Set to lower value if for pelvis, higher for others
        if 'pelvis' in currState:
            stateTracking.setWeightForState(currState,1.0)
        else:
            stateTracking.setWeightForState(currState,5.0)

#Add contact tracking goal

#Create the GRF tracking goal
contactTracking = osim.MocoContactTrackingGoal('contact',1.0)

#Set the external loads file
contactTracking.setExternalLoadsFile('Jog05_grf.xml')

#Set the force names for the right foot
forceNamesRightFoot = osim.StdVectorString()
forceNamesRightFoot.append('/forceset/contactHeel_r')
forceNamesRightFoot.append('/forceset/contactMH1_r')
forceNamesRightFoot.append('/forceset/contactMH3_r')
forceNamesRightFoot.append('/forceset/contactMH5_r')
forceNamesRightFoot.append('/forceset/contactOtherToes_r')
forceNamesRightFoot.append('/forceset/contactHallux_r')

#Set the force names for the left foot
forceNamesLeftFoot = osim.StdVectorString()
forceNamesLeftFoot.append('/forceset/contactHeel_l')
forceNamesLeftFoot.append('/forceset/contactMH1_l')
forceNamesLeftFoot.append('/forceset/contactMH3_l')
forceNamesLeftFoot.append('/forceset/contactMH5_l')
forceNamesLeftFoot.append('/forceset/contactOtherToes_l')
forceNamesLeftFoot.append('/forceset/contactHallux_l')

#Add the tracking groups
#Right foot
trackRightGRF = osim.MocoContactTrackingGoalGroup(forceNamesRightFoot,'calcn_r')
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r')
contactTracking.addContactGroup(trackRightGRF)
#Left foot
trackLeftGRF = osim.MocoContactTrackingGoalGroup(forceNamesLeftFoot,'calcn_l')
trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l')
contactTracking.addContactGroup(trackLeftGRF)

#Add to problem
problem.addGoal(contactTracking)
        
#Set mesh intervals and convergence tolerance in solver
solver = study.initCasADiSolver()
solver.set_optim_constraint_tolerance(1e-3)
solver.set_optim_convergence_tolerance(1e-3)
solver.set_num_mesh_intervals(50)
solver.set_optim_max_iterations(1500)
# #Use an average performing guess from past iteration
# solver.setGuessFile('torque_driven_state_tracking_solution.sto')

# Solve and visualize.
solution = study.solve()
# study.visualize(solution)

# problem.printToXML('test.omoco')

### trunk rotation has weird movement in tracked solution???
### adding GRF tracking did not work well...

#Plot against IK to compare

#Read in state tracking results
#Need to skip first 17 rows
stateTracking_df = pd.read_table('torque_driven_state_tracking_solution.sto',
                                 skiprows = list(range(0,17)))
#Convert marker tracking results to degrees
for cc in range(0,len(list(stateTracking_df))):
    #Get current column header
    currHeader = list(stateTracking_df)[cc]
    #Check if it contains the jointset and value keywords
    if 'jointset' in currHeader and 'value' in currHeader:
        #Check if it is a pelvis translation value not to convert
        if '_tx' not in currHeader and '_ty' not in currHeader and '_tz' not in currHeader:
            #Convert column to degrees
            stateTracking_df[currHeader] = np.rad2deg(stateTracking_df[currHeader])

#Create list of coordinates to plot
coordPlot = ['lumbar_extension','lumbar_bending','lumbar_rotation',
             'pelvis_tx','pelvis_ty','pelvis_tz',
             'pelvis_list','pelvis_tilt','pelvis_rotation',
             'hip_flexion_r','hip_adduction_r','hip_rotation_r',
             'knee_angle_r','ankle_angle_r']
coordPlotName = ['Lumbar Extension','Lumbar Bending','Lumbar Rotation',
                 'Pelvis X Translation','Pelvis Y Translation','Pelvis Z Translation',
                 'Pelvis List','Pelvis Tilt','Pelvis Rotation',
                 'R. Hip Flexion','R. Hip Adduction','R. Hip Rotation',
                 'R. Knee Angle','R. Ankle Angle']

#Plot results

#Set subplot shape
fig, axs = plt.subplots(5,3,figsize=(12,10))
#Set a list to access the axes coordinates in a loop
whichAx = [[0,0],[0,1],[0,2],
           [1,0],[1,1],[1,2],
           [2,0],[2,1],[2,2],
           [3,0],[3,1],[3,2],
           [4,0],[4,1],[4,2]]
#Loop through coordinates
for cc in range(0,len(coordPlot)):
    #Plot IK data
    ik_df.plot.line(x = 'time', y = coordPlot[cc],
                    legend = False, xlim = [startTime,endTime],
                    color = 'black', linestyle = '-',
                    linewidth = 1.5, ax = axs[whichAx[cc][0],whichAx[cc][1]])
    #Plot marker tracking data
    #Get the absolute path to the current coordinate
    stateTracking_df.plot.line(x = 'time', y = scaledModel.updCoordinateSet().get(coordPlot[cc]).getAbsolutePathString()+'/value',
                               legend = False, xlim = [startTime,endTime],
                               color = 'red', linestyle = '--',
                               linewidth = 1.5, ax = axs[whichAx[cc][0],whichAx[cc][1]])
    #Set title
    axs[whichAx[cc][0],whichAx[cc][1]].set_title(coordPlotName[cc])
#Turn off the final extra subplot axes
axs[4,2].axis('off')
#Set tight layout for figure
fig.tight_layout()
#Show figure
plt.show()