
%%%% TODO: add starting notes...

%Import OpenSim libraries
import org.opensim.modeling.*

%Add supplementary functions to path
addpath('functions');

%%

% 
% %  Set-up
% 

%Set base directory
baseDir = pwd;

%%%%%% TODO: set runOpts variable; add progressive displays from run
%%%%%% processes

%% Inverse Kinematics

%This step uses the .trc data from the sprint trial combined with the base 3D
%model to extract the kinematics for the trial

%Create folder if it isn't present
if ~exist('ik', 'dir')
    mkdir('ik')
end

%Add file logger
Logger.addFileSink('ik\\output_log.log');

%Load the model to use for IK
ikModel = Model('models\\base3DModel.osim');

%Initialise IK tool
ikTool = InverseKinematicsTool();
ikTool.setName('sprint');
% % % ikTool.setName('jog');

%Set model
ikTool.setModel(ikModel);

%Set marker file
ikTool.set_marker_file('data\\sprint.trc');
% % % ikTool.set_marker_file('data\\jog.trc');

%Create and set task set based off .trc markers

%Set-up a dictionary with various marker task weights
%Note that this is somewhat generic and could be edited
ikTaskWeights = struct('C7', 10, 'RSH', 10, 'LSH', 10, 'MAN', 10, 'T7', 10, 'LARM', 2.5,...
    'LELB', 5, 'LFOREARM', 10, 'LWR', 10, 'RARM', 2.5, 'RELB', 5, 'RFOREARM', 10,...
    'RWR', 10, 'RASI', 5, 'LASI', 5, 'SACR', 10, 'P1', 7.5, 'P2', 7.5,...
    'P3', 7.5, 'LTHLP', 5, 'LTHAP', 5, 'LTHAD', 5, 'LTHLD', 5, 'LLEPI', 10, 'LPAT', 10,...
    'LTIAP', 5, 'LTIAD', 5, 'LTILAT', 5, 'LLMAL', 10, 'LHEEL', 5, 'LMFS', 2.5,...
    'LMFL', 2.5, 'LP1MT', 5, 'LTOE', 2.5, 'LP5MT', 5, 'RTHLD', 5, 'RTHAP', 5,...
    'RTHLP', 5, 'RTHAD', 5, 'RLEPI', 10, 'RPAT', 10, 'RTIAP', 5, 'RTIAD', 5, 'RTILAT', 5,...
    'RLMAL', 10, 'RHEEL', 5, 'RMFS', 2.5, 'RMFL', 2.5, 'RP1MT', 5, 'RTOE', 2.5,...
    'RP5MT', 5);
markers = fieldnames(ikTaskWeights);
%Create task set
ikTaskSet = IKTaskSet();
%Loop through markers and append to task set
for mm = 1:length(markers)
    %Create blank marker task
    markerTask = IKMarkerTask();
    %Set marker name
    markerTask.setName(markers{mm});
    %Set marker weight
    markerTask.setWeight(ikTaskWeights.(markers{mm}));
    %Clone and append to task set
    ikTaskSet.cloneAndAppend(markerTask);
end
    
%Set task set
ikTool.set_IKTaskSet(ikTaskSet);

%Set times
%Identify right foot strike to right foot strike for analysis from GRF data
%Use the pre-defined function to do so
%Full vs. half gait cycle options available
[startTime,endTime] = getHalfGaitCycleRight('data\\sprint_grf.mot');
% % % [startTime,endTime] = getFullGaitCycleRight('data\\sprint_grf.mot');
% % % [startTime,endTime] = getHalfGaitCycleRight('data\\jog_grf.mot');

%Set times
ikTool.setStartTime(startTime)
ikTool.setEndTime(endTime)

%Set output filename
ikTool.set_output_motion_file('ik\\sprint_ik.mot')
% % % ikTool.set_output_motion_file('ik\\jog_ik.mot')

%Set results directory
ikTool.set_results_directory('ik\\');

%Print to file
ikTool.print('ik\\sprint_ik_setup.xml');
% % % ikTool.print('ik\\jog_ik_setup.xml');

%Run IK
ikTool.run()

%Remove file sink
Logger.removeFileSink();

%% Convert IK to states

%Use function
kinematicsToStates('ik\\sprint_ik.mot', 'models\\base3DModel.osim', 'ik\\sprint_ik_states.sto');
% % % kinematicsToStates('ik\\jog_ik.mot', 'models\\base3DModel.osim', 'ik\\jog_ik_states.sto');

%% Residual Reduction

%This step takes the IK results combined with GRF data to run the residual
%reduction algorithm to produce more dynamically consistent kinematics and
%models with respect to the external loads data

%Create folder if it isn't present
if ~exist('rra', 'dir')
    mkdir('rra')
end

%Create actuators set to use in RRA

%Set the base RRA model
rraModel = Model('models\\base3DModel.osim');

%Create lists for point, torque and coordinate actuators
%Point residuals
pointActuators = [{'FX'}, {'FY'}, {'FZ'}];
pointDirections = [Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,1)];
%Torque residuals
torqueActuators = [{'MX'}, {'MY'}, {'MZ'}];
torqueAxes = [Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,1)];
%Coordinate actuators
nAct = 1;
for cc = 0:rraModel.getCoordinateSet().getSize()-1
    %Don't need pelvis actuators
    if ~contains(char(rraModel.getCoordinateSet().get(cc).getName()),'pelvis')
        coordinateActuators{nAct} = char(rraModel.getCoordinateSet().get(cc).getName());
        nAct = nAct + 1;
    end
end

%Set structure for optimal forces
rraOptForce = struct('hip_flexion_r', 300, 'hip_adduction_r', 200, 'hip_rotation_r', 100, ...
    'knee_angle_r', 300, 'ankle_angle_r', 300, 'subtalar_angle_r', 100, 'mtp_angle_r', 100, ...
    'hip_flexion_l', 300, 'hip_adduction_l', 200, 'hip_rotation_l', 100, ...
    'knee_angle_l', 300, 'ankle_angle_l', 300, 'subtalar_angle_l', 100, 'mtp_angle_l', 100, ...
    'lumbar_extension', 300, 'lumbar_bending', 200, 'lumbar_rotation', 100, ...
    'arm_flex_r', 100, 'arm_add_r', 100, 'arm_rot_r', 100, ...
    'arm_flex_l', 100, 'arm_add_l', 100, 'arm_rot_l', 100, ...
    'elbow_flex_r', 100, 'pro_sup_r', 50, 'wrist_flex_r', 50, 'wrist_dev_r', 50, ...
    'elbow_flex_l', 100, 'pro_sup_l', 50, 'wrist_flex_l', 50, 'wrist_dev_l', 50);
    

%Create rra actuators force set
rraActuators = ForceSet();

%%%% TODO: package in function???

%Loop through point actuators to create and append to force set
for pp = 1:length(pointActuators)
    %Create point actuator
    actu = PointActuator();
    %Set name
    actu.setName(pointActuators{pp});
    %Set body name
    actu.set_body('pelvis');
    %Set point as pelvis CoM
    actu.set_point(rraModel.updBodySet().get('pelvis').get_mass_center());
    %Set direction
    actu.set_direction(pointDirections(pp));
    %Set optimal force
    actu.set_optimal_force(5.0);
    %Clone to force set
    rraActuators.cloneAndAppend(actu);
end

%Loop through torque actuators to create and append to force set
for pp = 1:length(torqueActuators)
    %Create torque actuator
    actu = TorqueActuator();
    %Set name
    actu.setName(torqueActuators{pp});
    %Set body names
    actu.set_bodyA('pelvis');
    actu.set_bodyB('ground');
    %Set torque axis
    actu.set_axis(torqueAxes(pp));
    %Set optimal force
    actu.set_optimal_force(2.5);
    %Clone to force set
    rraActuators.cloneAndAppend(actu);
end

%Loop through coordinate actuators to create and append to force set
for pp = 1:length(coordinateActuators)
    %Create coordinate actuator
    actu = CoordinateActuator();
    %Set name
    actu.setName(coordinateActuators{pp});
    %Set coordinate
    actu.set_coordinate(coordinateActuators{pp});
    %Set optimal force
    actu.set_optimal_force(rraOptForce.(coordinateActuators{pp}));
    %Clone to force set
    rraActuators.cloneAndAppend(actu);
end

%Write rra actuators force set to file
rraActuators.print('rra\rraActuators.xml');

%Create the rra task set
rraTasks = CMC_TaskSet();

%Set the weights for the task set in a struct
%Somewhat generic and could be altered
rraTaskWeights = struct('pelvis_tx', 5, 'pelvis_ty', 5, 'pelvis_tz', 5, ...
    'pelvis_tilt', 1000, 'pelvis_list', 500, 'pelvis_rotation', 100, ...
    'hip_flexion_r', 50, 'hip_adduction_r', 25, 'hip_rotation_r', 25, ...
    'knee_angle_r', 50, 'ankle_angle_r', 50, 'subtalar_angle_r', 20, 'mtp_angle_r', 20, ...    
    'hip_flexion_l', 50, 'hip_adduction_l', 25, 'hip_rotation_l', 25, ...
    'knee_angle_l', 50, 'ankle_angle_l', 50, 'subtalar_angle_l', 20, 'mtp_angle_l', 20, ...
    'lumbar_extension', 30, 'lumbar_bending', 30, 'lumbar_rotation', 30, ...
    'arm_flex_r', 25, 'arm_add_r', 25, 'arm_rot_r', 10, ...
    'arm_flex_l', 25, 'arm_add_l', 25, 'arm_rot_l', 10, ...
    'elbow_flex_r', 25, 'pro_sup_r', 10, 'wrist_flex_r', 25, 'wrist_dev_r', 10, ...
    'elbow_flex_l', 25, 'pro_sup_l', 10, 'wrist_flex_l', 25, 'wrist_dev_l', 10);

%Loop through model coordinates to create and append tasks
for cc = 0:rraModel.getCoordinateSet().getSize()-1
    %Create task
    task = CMC_Joint();
    %Set active
    task.setActive(true);
    %Get coordinate name as string
    coord = char(rraModel.getCoordinateSet().get(cc).getName());
    %Set name
    task.setName(coord);
    %Set coordinate
    task.setCoordinateName(coord);
    %Set weight
    task.setWeight(rraTaskWeights.(coord));
    %Set kp and kv to standard values
    task.setKP(100);
    task.setKV(20);
    %Clone and append to task set
    rraTasks.cloneAndAppend(task);    
end

%Write rra tasks set to file
rraTasks.print('rra\rraTasks.xml');

%%%% TODO: loop through RRA iterations???

%Make current directory for RRA if it doesn't yet exist
%%%% TODO: looping...
if ~exist('rra\\iter1', 'dir')
    mkdir('rra\\iter1')
end

%Add file logger
%%%% TODO: looping...
Logger.addFileSink('rra\\iter1\\output_log.log');

%Create RRA tool
%Note that directives to filenames in the tool use filepaths relative to
%where the setup file is saved
rraTool = RRATool();

%Set actuators file
forceSetFiles = ArrayStr(); forceSetFiles.append('..\\rraActuators.xml');
rraTool.setForceSetFiles(forceSetFiles);

%Set tracking tasks file
rraTool.setTaskSetFileName('..\\rraTasks.xml');

%Set a low pass filter frequency on the kinematics data
rraTool.setLowpassCutoffFrequency(12);

%Set to adjust the COM to reduce residuals
rraTool.setAdjustCOMToReduceResiduals(true);

%Set the torso body COM to be the one that gets adjusted
rraTool.setAdjustedCOMBody('torso');

%Set external loads file
rraTool.setExternalLoadsFileName('..\\..\\data\\sprint_grf.xml');

%Set tool name based on iteration
%%%% TODO: looping...
rraTool.setName('rra1');

%Set desired kinematics
%%%% TODO: looping...

%Use IK data
rraTool.setDesiredKinematicsFileName('..\\..\\ik\\sprint_ik.mot');

%Set initial and final time
%You need to use the IK first and last times here as I don't think the tool
%likes if the IK data doesn't have the times you put in
rraTool.setStartTime(Storage('ik\\sprint_ik.mot').getFirstTime());
rraTool.setFinalTime(Storage('ik\\sprint_ik.mot').getLastTime());

%Set model
%%%% TODO: use Model object or filename???
%%%% TODO: looping...
rraTool.setModelFilename('..\\..\\models\\base3DModel.osim');

%Set output model file
%%%% TODO: looping...
rraTool.setOutputModelFileName('rraAdjustedModel_1.osim');

%Set print level
%%%% TODO: can't access this parameter???

%Print out the tool
%%%% TODO: looping...
rraTool.print('rra\\iter1\\setupRRA_1.xml');

%Run RRA
%Navigate to rra directory
cd('rra\\iter1');
%Reloading file to run seems necessary to not have OpenSim crack it...
runRRA = RRATool('setupRRA_1.xml');
runRRA.run();


%%%%%%% Mass adjustments - need verbosity for this???


%Return to base directory
cd(baseDir);





%Compare solution to experimental data
if compareRRA
    
    %Kinematics
    
    %Set coordinate states to work through
    %%%% TODO: add more as this is a 3d problem...
    plotStates = [{'/jointset/ground_pelvis/pelvis_ty/value'};
        {'/jointset/ground_pelvis/pelvis_tx/value'};
        {'/jointset/ground_pelvis/pelvis_tilt/value'}
        {'/jointset/back/lumbar_extension/value'};
        {'/jointset/hip_r/hip_flexion_r/value'};
        {'/jointset/hip_l/hip_flexion_l/value'};
        {'/jointset/knee_r/knee_angle_r/value'};
        {'/jointset/knee_l/knee_angle_l/value'};
        {'/jointset/ankle_r/ankle_angle_r/value'};
        {'/jointset/ankle_l/ankle_angle_l/value'};
        {'/jointset/acromial_r/arm_flex_r/value'};
        {'/jointset/acromial_l/arm_flex_l/value'};
        {'/jointset/elbow_r/elbow_flex_r/value'};
        {'/jointset/elbow_l/elbow_flex_l/value'};
        ];
    
    %Read IK states
    ikTable = TimeSeriesTable('ik\\sprint_ik_states.sto');
    
    %Read RRA states
    %%%% TODO: looping
    rraTable = TimeSeriesTable('rra\\iter1\\rra1_states.sto');
    
    %Set plot
    figure; set(gcf,'units','normalized','position',[0.01 0.01 0.6 0.9]);
    
    %Set IK time vector
    for tt = 0:ikTable.getIndependentColumn().size()-1
        ikTime(tt+1) = ikTable.getIndependentColumn().get(tt).doubleValue() - ...
            ikTable.getIndependentColumn().get(0).doubleValue();
    end
    
    %Set RRA time vector
    for tt = 0:rraTable.getIndependentColumn().size()-1
        rraTime(tt+1) = rraTable.getIndependentColumn().get(tt).doubleValue() - ...
            rraTable.getIndependentColumn().get(0).doubleValue();
    end
    
    %Loop through states
    for nn = 1:length(plotStates)
        %Set subplot
        subplot(4,4,nn); hold on
        %Turn box on
        ax = gca;
        box(ax,'on');
        ax.LineWidth = 1;
        set(ax, 'Layer', 'Top')
        %White background
        set(gcf,'Color','w');
        %Get experimental data from kinematic table
        plot(ikTime, ikTable.getDependentColumn(plotStates{nn}).getAsMat(),...
            'k', 'LineWidth', 1.5);
        %Get RRA data
        %%%%% TODO: looping...
        plot(rraTime, rraTable.getDependentColumn(plotStates{nn}).getAsMat(),...
            'r--', 'LineWidth', 1.5);
        %Add plot labels
        %X-label
        xlabel('Time (s)');
        %Y-label
        if contains(plotStates{nn},'pelvis_t')
            ylabel('Pos. (m)')
        else
            ylabel('Angle (rad)')
        end
        %Title
        splitState = split(plotStates{nn},'/');
        title(regexprep(splitState{4},'_',' '));        
    end
    
end

%% Update Models

%This step adapts the base 3D model to create a 2D model that includes lower limb
%muscles from the Ong et al. model included in the 'models' directory

%Load the base 3D model
base3DModel = Model('models\\base3DModel.osim');

%% Create 2D version of model

%Set a copy of this model to edit
model2D = Model('models\\base3DModel.osim');

%Convert hip joints to 2D

%Get hip joints from scaled model to work from
origHip_r = base3DModel.getJointSet().get('hip_r');
origHip_l = base3DModel.getJointSet().get('hip_l');
    
%Create pin joint for right hip
%Set parent and child locations
locationInParent = origHip_r.get_frames(0).get_translation();
orientationInParent = origHip_r.get_frames(0).get_orientation();
locationInChild = origHip_r.get_frames(1).get_translation();
orientationInChild = origHip_r.get_frames(1).get_orientation();
%Get bodies for joint
pelvisBody = model2D.getBodySet().get('pelvis');
femurBody = model2D.getBodySet().get('femur_r');
%Create joint
hipJoint_r = PinJoint('hip_r', pelvisBody, locationInParent,...
    orientationInParent, femurBody, locationInChild,...
    orientationInChild);
%Update additional joint parameters
hipJoint_r.getCoordinate().setName('hip_flexion_r'); %set coordinate name
hipJoint_r.getCoordinate().setRangeMin(origHip_r.get_coordinates(0).getRangeMin()); %get min from old joint
hipJoint_r.getCoordinate().setRangeMax(origHip_r.get_coordinates(0).getRangeMax()); %get max from old joint
hipJoint_r.getCoordinate().setDefaultValue(0); %set default value to zero
hipJoint_r.getCoordinate().setDefaultSpeedValue(0); %set default speed to zero
hipJoint_r.getCoordinate().set_clamped(true); %set joint to clamped
hipJoint_r.getCoordinate().set_locked(false); %set locked to false
hipJoint_r.getCoordinate().set_prescribed(false); %set prescribed to false

%Create pin joint for left hip
%Set parent and child locations
locationInParent = origHip_l.get_frames(0).get_translation();
orientationInParent = origHip_l.get_frames(0).get_orientation();
locationInChild = origHip_l.get_frames(1).get_translation();
orientationInChild = origHip_l.get_frames(1).get_orientation();
%Get bodies for joint
pelvisBody = model2D.getBodySet().get('pelvis');
femurBody = model2D.getBodySet().get('femur_l');
%Create joint
hipJoint_l = PinJoint('hip_l', pelvisBody, locationInParent,...
    orientationInParent, femurBody, locationInChild,...
    orientationInChild);
%Update additional joint parameters
hipJoint_l.getCoordinate().setName('hip_flexion_l'); %set coordinate name
hipJoint_l.getCoordinate().setRangeMin(origHip_l.get_coordinates(0).getRangeMin()); %get min from old joint
hipJoint_l.getCoordinate().setRangeMax(origHip_l.get_coordinates(0).getRangeMax()); %get max from old joint
hipJoint_l.getCoordinate().setDefaultValue(0); %set default value to zero
hipJoint_l.getCoordinate().setDefaultSpeedValue(0); %set default speed to zero
hipJoint_l.getCoordinate().set_clamped(true); %set joint to clamped
hipJoint_l.getCoordinate().set_locked(false); %set locked to false
hipJoint_l.getCoordinate().set_prescribed(false); %set prescribed to false

%Remove existing hip joints from model
model2D.getJointSet().remove(model2D.getJointSet().get('hip_r'));
model2D.getJointSet().remove(model2D.getJointSet().get('hip_l'));

%Add new hip joints
model2D.addJoint(hipJoint_r);
model2D.addJoint(hipJoint_l);

%Update ankle to pin joint

%Get hip joints from scaled model to work from
origAnkle_r = base3DModel.getJointSet().get('ankle_r');
origAnkle_l = base3DModel.getJointSet().get('ankle_l');
    
%Create pin joint for right ankle
%Set parent and child locations
locationInParent = origAnkle_r.get_frames(0).get_translation();
orientationInParent = origAnkle_r.get_frames(0).get_orientation();
locationInChild = origAnkle_r.get_frames(1).get_translation();
orientationInChild = origAnkle_r.get_frames(1).get_orientation();
%Get bodies for joint
tibiaBody = model2D.getBodySet().get('tibia_r');
talusBody = model2D.getBodySet().get('talus_r');
%Create joint
ankleJoint_r = PinJoint('ankle_r', tibiaBody, locationInParent,...
    orientationInParent, talusBody, locationInChild,...
    orientationInChild);
%Update additional joint parameters
ankleJoint_r.getCoordinate().setName('ankle_angle_r'); %set coordinate name
ankleJoint_r.getCoordinate().setRangeMin(origHip_r.get_coordinates(0).getRangeMin()); %get min from old joint
ankleJoint_r.getCoordinate().setRangeMax(origHip_r.get_coordinates(0).getRangeMax()); %get max from old joint
ankleJoint_r.getCoordinate().setDefaultValue(0); %set default value to zero
ankleJoint_r.getCoordinate().setDefaultSpeedValue(0); %set default speed to zero
ankleJoint_r.getCoordinate().set_clamped(true); %set joint to clamped
ankleJoint_r.getCoordinate().set_locked(false); %set locked to false
ankleJoint_r.getCoordinate().set_prescribed(false); %set prescribed to false

%Create pin joint for left ankle
%Set parent and child locations
locationInParent = origAnkle_l.get_frames(0).get_translation();
orientationInParent = origAnkle_l.get_frames(0).get_orientation();
locationInChild = origAnkle_l.get_frames(1).get_translation();
orientationInChild = origAnkle_l.get_frames(1).get_orientation();
%Get bodies for joint
tibiaBody = model2D.getBodySet().get('tibia_l');
talusBody = model2D.getBodySet().get('talus_l');
%Create joint
ankleJoint_l = PinJoint('ankle_l', tibiaBody, locationInParent,...
    orientationInParent, talusBody, locationInChild,...
    orientationInChild);
%Update additional joint parameters
ankleJoint_l.getCoordinate().setName('ankle_angle_l'); %set coordinate name
ankleJoint_l.getCoordinate().setRangeMin(origHip_r.get_coordinates(0).getRangeMin()); %get min from old joint
ankleJoint_l.getCoordinate().setRangeMax(origHip_r.get_coordinates(0).getRangeMax()); %get max from old joint
ankleJoint_l.getCoordinate().setDefaultValue(0); %set default value to zero
ankleJoint_l.getCoordinate().setDefaultSpeedValue(0); %set default speed to zero
ankleJoint_l.getCoordinate().set_clamped(true); %set joint to clamped
ankleJoint_l.getCoordinate().set_locked(false); %set locked to false
ankleJoint_l.getCoordinate().set_prescribed(false); %set prescribed to false

%Remove existing hip joints from model
model2D.getJointSet().remove(model2D.getJointSet().get('ankle_r'));
model2D.getJointSet().remove(model2D.getJointSet().get('ankle_l'));

%Add new hip joints
model2D.addJoint(ankleJoint_r);
model2D.addJoint(ankleJoint_l);

%Convert lumbar joint to 2D

%Get back joint from scaled model to work from
origLumbar = base3DModel.getJointSet().get('back');
    
%Create pin joint for back
%Set parent and child locations
locationInParent = origLumbar.get_frames(0).get_translation();
orientationInParent = origLumbar.get_frames(0).get_orientation();
locationInChild = origLumbar.get_frames(1).get_translation();
orientationInChild = origLumbar.get_frames(1).get_orientation();
%Get bodies for joint
pelvisBody = model2D.getBodySet().get('pelvis');
torsoBody = model2D.getBodySet().get('torso');
%Create joint
lumbarJoint = PinJoint('back', pelvisBody, locationInParent,...
    orientationInParent, torsoBody, locationInChild,...
    orientationInChild);
%Update additional joint parameters
lumbarJoint.getCoordinate().setName('lumbar_extension'); %set coordinate name
lumbarJoint.getCoordinate().setRangeMin(origLumbar.get_coordinates(0).getRangeMin()); %get min from old joint
lumbarJoint.getCoordinate().setRangeMax(origLumbar.get_coordinates(0).getRangeMax()); %get max from old joint
lumbarJoint.getCoordinate().setDefaultValue(0); %set default value to zero
lumbarJoint.getCoordinate().setDefaultSpeedValue(0); %set default speed to zero
lumbarJoint.getCoordinate().set_clamped(true); %set joint to clamped
lumbarJoint.getCoordinate().set_locked(false);  %set locked to false
lumbarJoint.getCoordinate().set_prescribed(false); %set prescribed to false

%Remove existing back joint from model
model2D.getJointSet().remove(model2D.getJointSet().get('back'));

%Add new back joint
model2D.addJoint(lumbarJoint)

%Update pelvis model to 2D

%Get original pelvis from scaled model to work from
origPelvis = base3DModel.getJointSet().get('ground_pelvis');

%Create planar joint for ground-pelvis
pelvisJoint = PlanarJoint('ground_pelvis',...
    model2D.getGround(), Vec3(0,0,0), Vec3(0,0,0),...
    model2D.getBodySet().get('pelvis'), Vec3(0,0,0), Vec3(0,0,0));
%Update additional joint parameters
%Pelvis tilt
pelvisJoint.get_coordinates(0).setName('pelvis_tilt'); %set coordinate name
pelvisJoint.get_coordinates(0).setRangeMin(origPelvis.get_coordinates(0).getRangeMin()); %get min from old joint
pelvisJoint.get_coordinates(0).setRangeMax(origPelvis.get_coordinates(0).getRangeMax()); %get max from old joint
pelvisJoint.get_coordinates(0).setDefaultValue(0); %set default value to zero
pelvisJoint.get_coordinates(0).setDefaultSpeedValue(0); %set default speed to zero
pelvisJoint.get_coordinates(0).set_clamped(true); %set joint to clamped
pelvisJoint.get_coordinates(0).set_locked(false); %set locked to false
pelvisJoint.get_coordinates(0).set_prescribed(false); %set prescribed to false
%Pelvis tx
pelvisJoint.get_coordinates(1).setName('pelvis_tx'); %set coordinate name
pelvisJoint.get_coordinates(1).setRangeMin(origPelvis.get_coordinates(3).getRangeMin()); %get min from old joint
pelvisJoint.get_coordinates(1).setRangeMax(origPelvis.get_coordinates(3).getRangeMax()); %get max from old joint
pelvisJoint.get_coordinates(1).setDefaultValue(0); %set default value to zero
pelvisJoint.get_coordinates(1).setDefaultSpeedValue(0); %set default speed to zero
pelvisJoint.get_coordinates(1).set_clamped(true); %set joint to clamped
pelvisJoint.get_coordinates(1).set_locked(false); %set locked to false
pelvisJoint.get_coordinates(1).set_prescribed(false); %set prescribed to false
%Pelvis ty
pelvisJoint.get_coordinates(2).setName('pelvis_ty'); %set coordinate name
pelvisJoint.get_coordinates(2).setRangeMin(origPelvis.get_coordinates(4).getRangeMin()); %get min from old joint
pelvisJoint.get_coordinates(2).setRangeMax(origPelvis.get_coordinates(4).getRangeMax()); %get max from old joint
pelvisJoint.get_coordinates(2).setDefaultValue(0.95); %set default value to zero
pelvisJoint.get_coordinates(2).setDefaultSpeedValue(0); %set default speed to zero
pelvisJoint.get_coordinates(2).set_clamped(true); %set joint to clamped
pelvisJoint.get_coordinates(2).set_locked(false); %set locked to false
pelvisJoint.get_coordinates(2).set_prescribed(false); %set prescribed to false

%Remove existing pelvis joint from model
model2D.getJointSet().remove(model2D.getJointSet().get('ground_pelvis'));

%Add new pelvis joint
model2D.addJoint(pelvisJoint)

%Convert ankle and foot joints to planar variants

%Update ankle joints orientation to planar
model2D.getJointSet().get('ankle_r').get_frames(0).set_orientation(Vec3(0,0,0));
model2D.getJointSet().get('ankle_l').get_frames(0).set_orientation(Vec3(0,0,0));
model2D.getJointSet().get('ankle_r').get_frames(1).set_orientation(Vec3(0,0,0));
model2D.getJointSet().get('ankle_l').get_frames(1).set_orientation(Vec3(0,0,0));

%Update mtp joints orientation to planar
model2D.getJointSet().get('mtp_r').get_frames(0).set_orientation(Vec3(0,deg2rad(180),0));
model2D.getJointSet().get('mtp_l').get_frames(0).set_orientation(Vec3(0,deg2rad(180),0));
model2D.getJointSet().get('mtp_r').get_frames(1).set_orientation(Vec3(0,deg2rad(180),0));
model2D.getJointSet().get('mtp_l').get_frames(1).set_orientation(Vec3(0,deg2rad(180),0));

%Update subtalar joints orientation to planar
model2D.getJointSet().get('subtalar_r').get_frames(0).set_orientation(Vec3(0,deg2rad(180),0));
model2D.getJointSet().get('subtalar_l').get_frames(0).set_orientation(Vec3(0,deg2rad(180),0));
model2D.getJointSet().get('subtalar_r').get_frames(1).set_orientation(Vec3(0,deg2rad(180),0));
model2D.getJointSet().get('subtalar_l').get_frames(1).set_orientation(Vec3(0,deg2rad(180),0));

%Convert arm joints to 2D

%Get acromial joints from scaled model to work from
origArm_r = base3DModel.getJointSet().get('acromial_r');
origArm_l = base3DModel.getJointSet().get('acromial_l');
    
%Create pin joint for right arm
%Set parent and child locations
locationInParent = origArm_r.get_frames(0).get_translation();
orientationInParent = origArm_r.get_frames(0).get_orientation();
locationInChild = origArm_r.get_frames(1).get_translation();
orientationInChild = origArm_r.get_frames(1).get_orientation();
%Get bodies for joint
torsoBody = model2D.getBodySet().get('torso');
humerusBody = model2D.getBodySet().get('humerus_r');
%Create joint
armJoint_r = PinJoint('acromial_r', torsoBody, locationInParent,...
    orientationInParent, humerusBody, locationInChild,...
    orientationInChild);
%Update additional joint parameters
armJoint_r.getCoordinate().setName('arm_flex_r'); %set coordinate name
armJoint_r.getCoordinate().setRangeMin(origArm_r.get_coordinates(0).getRangeMin()); %get min from old joint
armJoint_r.getCoordinate().setRangeMax(origArm_r.get_coordinates(0).getRangeMax()); %get max from old joint
armJoint_r.getCoordinate().setDefaultValue(0); %set default value to zero
armJoint_r.getCoordinate().setDefaultSpeedValue(0); %set default speed to zero
armJoint_r.getCoordinate().set_clamped(true); %set joint to clamped
armJoint_r.getCoordinate().set_locked(false); %set locked to false
armJoint_r.getCoordinate().set_prescribed(false); %set prescribed to false

%Create pin joint for left arm
%Set parent and child locations
locationInParent = origArm_l.get_frames(0).get_translation();
orientationInParent = origArm_l.get_frames(0).get_orientation();
locationInChild = origArm_l.get_frames(1).get_translation();
orientationInChild = origArm_l.get_frames(1).get_orientation();
%Get bodies for joint
torsoBody = model2D.getBodySet().get('torso');
humerusBody = model2D.getBodySet().get('humerus_l');
%Create joint
armJoint_l = PinJoint('acromial_l', torsoBody, locationInParent,...
    orientationInParent, humerusBody, locationInChild,...
    orientationInChild);
%Update additional joint parameters
armJoint_l.getCoordinate().setName('arm_flex_l'); %set coordinate name
armJoint_l.getCoordinate().setRangeMin(origArm_r.get_coordinates(0).getRangeMin()); %get min from old joint
armJoint_l.getCoordinate().setRangeMax(origArm_r.get_coordinates(0).getRangeMax()); %get max from old joint
armJoint_l.getCoordinate().setDefaultValue(0); %set default value to zero
armJoint_l.getCoordinate().setDefaultSpeedValue(0); %set default speed to zero
armJoint_l.getCoordinate().set_clamped(true); %set joint to clamped
armJoint_l.getCoordinate().set_locked(false); %set locked to false
armJoint_l.getCoordinate().set_prescribed(false); %set prescribed to false

%Remove existing hip joints from model
model2D.getJointSet().remove(model2D.getJointSet().get('acromial_r'));
model2D.getJointSet().remove(model2D.getJointSet().get('acromial_l'));

%Add new hip joints
model2D.addJoint(armJoint_r);
model2D.addJoint(armJoint_l);

%Update radioulnar orientation to better align with 2D model
%Supinate for level forearms
%This still results in the hand passing through the legs, but is just a visualisation issue
model2D.getJointSet().get('radioulnar_r').get_frames(1).set_orientation(Vec3(0,deg2rad(-90),0));
model2D.getJointSet().get('radioulnar_l').get_frames(1).set_orientation(Vec3(0,deg2rad(90),0));

%Finalize model connections
model2D.finalizeConnections();

%Dump into a model processor to weld the foot and hand joints
%Create the processor
editProcessor = ModelProcessor(model2D);
%Create the string array of joints to weld
weldList = [{'subtalar_r'}, {'subtalar_l'}, {'mtp_r'}, {'mtp_l'},...
    {'radioulnar_r'}, {'radioulnar_l'}, {'radius_hand_r'}, {'radius_hand_l'}];
weldJoints = StdVectorString();
for joint = 1:length(weldList)
    weldJoints.add(char(weldList(joint)));
end
%Append model operator for welding
editProcessor.append(ModOpReplaceJointsWithWelds(weldJoints));
%Process model output
final2DModel = editProcessor.process();

%Remove the marker set as it's no use in this 2D context
final2DModel.updMarkerSet().clearAndDestroy()

%Reset name
final2DModel.setName('base2DModel');

%Print 2D model output
final2DModel.print('models\\base2DModel.osim')

%% Adapt lower limb muscle set to 2D model

%This process scales the Ong et al model to the static trial of the
%experimental data, so that the lower limb muscle set can be appropriately
%adapted to the experimental model

%Set the parameters required to calculate for the manual scaling
%Note that this ignores the arms, as they don't need to be scaled to get
%appropriate anatomy for the lower limb muscles
scaleParameters = struct('torso', [], ...
    'pelvis', [], ...
    'femur_r', [], ...
    'femur_l', [], ...
    'tibia_r', [], ...
    'tibia_l', [], ...
    'talus_r', [], ...
    'talus_l', [], ...
    'calcn_r', [], ...
    'calcn_l', [], ...
    'toes_r', [], ...
    'toes_l', []);

%Set a variable to grab the XYZ label from
axesLabels = [{'X'}; {'Y'}; {'Z'}];

%Set the corresponding marker pairs to be used for the XYZ directions on
%each body
markerPairs = struct('torso', {{'C7','MAN','T7','MAN'},{'RASI','RSH','LASI','LSH'}, {'RSH','LSH'}}, ...
    'pelvis', {{'SACR','RASI','SACR','LASI'},{'RASI','LASI'},{'RASI','LASI'}}, ...
    'femur_r', {{'RLEPI','RMEPI'},{'RASI','RLEPI','RASI','RMEPI'},{'RLEPI','RMEPI'}}, ...
    'femur_l', {{'LLEPI','LMEPI'},{'LASI','LLEPI','LASI','LMEPI'},{'LLEPI','LMEPI'}}, ...
    'tibia_r', {{'RLMAL','RMMAL'},{'RLEPI','RLMAL','RMEPI','RMMAL'},{'RLMAL','RMMAL'}}, ...
    'tibia_l', {{'LLMAL','LMMAL'},{'LLEPI','LLMAL','LMEPI','LMMAL'},{'LLMAL','LMMAL'}}, ...
    'talus_r', {{'RLMAL','RMMAL'},{'RLMAL','RMMAL'},{'RLMAL','RMMAL'}}, ...
    'talus_l', {{'LLMAL','LMMAL'},{'LLMAL','LMMAL'},{'LLMAL','LMMAL'}}, ...
    'calcn_r', {{'RHEEL','RLMAL','RHEEL','RMMAL'},{'RHEEL','RLMAL','RHEEL','RMMAL'},{'RHEEL','RLMAL','RHEEL','RMMAL'}}, ...
    'calcn_l', {{'LHEEL','LLMAL','LHEEL','LMMAL'},{'LHEEL','LLMAL','LHEEL','LMMAL'},{'LHEEL','LLMAL','LHEEL','LMMAL'}}, ...
    'toes_r', {{'RLMAL','RTOE','RMMAL','RTOE'},{'RLMAL','RP5MT','RMMAL','RP1MT'},{'RP1MT','RP5MT'}}, ...
    'toes_l', {{'RLMAL','RTOE','RMMAL','RTOE'},{'RLMAL','RP5MT','RMMAL','RP1MT'},{'RP1MT','RP5MT'}});

%Create scale tool
%Note that this is really just for kinematic scaling, so many parameters
%related to mass and marker movement are ignored in this scaling process
scaleTool = ScaleTool();

%Set the generic model filename
scaleTool.getGenericModelMaker().setModelFileName('models\\gait9dof18musc_Ong_et_al_Moco.osim');

%Loop through the bodies and axes to create measurements
bodies = fieldnames(scaleParameters);
for bb = 1:length(bodies)
    
    %Loop through XYZ components
    for cc = 1:length(axesLabels)
        
        %Create the measurement
        currMeasurement = Measurement();
        
        %Set the name
        currMeasurement.setName([bodies{bb},'_',axesLabels{cc}]);
        
        %Get the current set of marker pairs
        currPairs = markerPairs(cc).(bodies{bb});
        
        %Get number of pairs
        nPairs = length(currPairs) / 2;
        
        %Loop through pairs
        for pp = 1:nPairs            
            %Get the two current markers
            marker1 = currPairs{pp*1}; marker2 = currPairs{pp*1+1};            
            %Create the marker pair object
            currMarkerPair = MarkerPair(marker1,marker2);            
            %Append to marker pair set in measurement
            currMeasurement.getMarkerPairSet().cloneAndAppend(currMarkerPair);           
        end
        
        %Add the body scale for the measurement
        currBodyScale = BodyScale();
        currBodyScale.setName(bodies{bb});
        axisNames = ArrayStr(); axisNames.append(axesLabels{cc});
        currBodyScale.setAxisNames(axisNames);
        currMeasurement.getBodyScaleSet().cloneAndAppend(currBodyScale);
        
        %Add the current measurement to the measurement set
        scaleTool.getModelScaler().getMeasurementSet().cloneAndAppend(currMeasurement);
        
    end

end

%Set marker file from the static trial
scaleTool.getModelScaler().setMarkerFileName('data\\static.trc');

%Set scaling order
scaleOrder = ArrayStr(); scaleOrder.append('measurements'); scaleOrder.append('manualScale');
scaleTool.getModelScaler().setScalingOrder(scaleOrder);

%Set marker placer as false
scaleTool.getMarkerPlacer().setApply(false);

%Create the model from the scale tool
scaledModel = scaleTool.createModel();

%Load in a copy of the base 2D model to add the muscles to
model2D_muscles = Model('models\\base2DModel.osim');

%Get the muscles from the scaled model
muscleSet = scaledModel.getMuscles();

%Loop through the muscles and append to the base 2D model
for mm = 0:muscleSet.getSize()-1
    model2D_muscles.updForceSet().cloneAndAppend(Force.safeDownCast(muscleSet.get(mm)));
end

%Update model name
model2D_muscles.setName('muscles2DModel');

%Finalise model connections
model2D_muscles.finalizeConnections();

%Print to file
model2D_muscles.print('models\\muscles2DModel.osim');



%%

%% Test out torque tracking with 2D

%% Convert the GRF data to 2D for the tracking sims

%Get original GRF data
forcesTable = TimeSeriesTable('data\\sprint_grf.mot');

%Get number of column labels and data rows
nLabels = forcesTable.getNumColumns();
nRows = forcesTable.getNumRows();

%Loop through labels and rows to manipulate data
for iLabel = 0:nLabels-1
    %Get current column label
    currLabel = forcesTable.getColumnLabel(iLabel);
    %Check for z-axis data
    if contains(char(currLabel), '_vz') || contains(char(currLabel), '_pz') || contains(char(currLabel), '_z')
        %Set rows to zero (i.e. no data)
        for iRow = 0:nRows-1
            forcesTable.getDependentColumn(currLabel).set(iRow,0);
        end
    end
end

%Write to new file
STOFileAdapter.write(forcesTable,'data\\sprint_grf_2D.mot');

%Adapt external loads file to match 2D GRF variant

%Read external loads file in
extLoads = ExternalLoads('data\\sprint_grf.xml', true);

%Update the datafile
extLoads.setDataFileName('sprint_grf_2D.mot');

%Write to file
extLoads.print('data\\sprint_grf_2D.xml');

%% Convert IK to states

%Use function
kinematicsToStates('ik\\sprint_ik.mot', 'models\\base3DModel.osim', 'ik\\sprint_ik_states.sto');
% % % kinematicsToStates('ik\\jog_ik.mot', 'models\\base3DModel.osim', 'ik\\jog_ik_states.sto');


%% Torque driven tracking simulation


%%%%%% TODO: Highly penalise pelvis torques!!!

%This step is designed to generate relevant controls for the 2D model that
%appropriately track the sprint kinematics and GRF from the experimental data
%using torque actuators

%Set simulation parameters
    
%Timing and mesh data
initialTime = Storage('ik\\sprint_ik_states.sto').getFirstTime();
finalTime = Storage('ik\\sprint_ik_states.sto').getLastTime();
duration = finalTime - initialTime;
meshNo = 50; %%% mesh number based on some work by Falisse et al. recommending 50 for half gait cycle 
meshInterval = duration / meshNo;

%Weights
trackingWeight = 10;
effortWeight = 0.1;
grfWeight = 2.5;
speedWeight = 1;
symmetryWeight = 1;

%General parameters
visualiseSolution = false;
compareToExp = false;

%Create the model processor for the tracking problem

%Load the model
torqueModel = Model('models\\base2DModel.osim');
torqueModel.initSystem();

%Add torque actuators to support model
%Set list to add reserves to based on coordinate list
%Set optimal force and max torque
optimalForce = 100;
maxTorque = Inf;
%Add torque actuators
for coordinate = 0:torqueModel.getCoordinateSet().getSize()-1
    %Get coordinate name
    coord = char(torqueModel.getCoordinateSet().get(coordinate));
    %Add actuator
    addReserve(torqueModel, coord, optimalForce, maxTorque);
end

%Finalise connections
torqueModel.finalizeConnections();

%Get the track model as a processor object
%Not really necessary but here in case some extra processing steps need to
%be inserted
torqueModelProcessor = ModelProcessor(torqueModel);

%Set-up tracking problem

%Process current model
trackModel = torqueModelProcessor.process();
trackModel.initSystem();

%Construct the tracking object and set basic parameters
track = MocoTrack();
track.setName('sprintTracking_torqueDriven');
track.setModel(torqueModelProcessor);

%Set kinematic data and parameters
tableProcessor = TableProcessor('ik\\sprint_ik_states.sto');
tableProcessor.append(TabOpLowPassFilter(12));
tableProcessor.append(TabOpUseAbsoluteStateNames());
track.setStatesReference(tableProcessor);
track.set_states_global_tracking_weight(trackingWeight);
%Set some coordinates not to be tracked or low weight so that kinematics
%can be better adjusted to reflect GRF tracking problem. Weights here
%reflect confidence in data from inverse kinematics solution
%Similar process done in Dembia et al. Moco paper
stateWeights = MocoWeightSet();
stateList = [{'/jointset/ground_pelvis/pelvis_ty'};
    {'/jointset/ground_pelvis/pelvis_tx'};
    {'/jointset/ground_pelvis/pelvis_tilt'}
    {'/jointset/back/lumbar_extension'};
    {'/jointset/hip_r/hip_flexion_r'};
    {'/jointset/hip_l/hip_flexion_l'};
    {'/jointset/knee_r/knee_angle_r'};
    {'/jointset/knee_l/knee_angle_l'};
    {'/jointset/ankle_r/ankle_angle_r'};
    {'/jointset/ankle_l/ankle_angle_l'};
    {'/jointset/acromial_r/arm_flex_r'};
    {'/jointset/acromial_l/arm_flex_l'};
    {'/jointset/elbow_r/elbow_flex_r'};
    {'/jointset/elbow_l/elbow_flex_l'};
    ];
weightList = [0; 1; 0.1; 0.5; 1; 1; 1; 1; 0.75; 0.75; 0.5; 0.5; 0.75; 0.75];
for ww = 1:length(stateList)
    stateWeights.cloneAndAppend(MocoWeight([stateList{ww},+'/value'], weightList(ww)));
    stateWeights.cloneAndAppend(MocoWeight([stateList{ww},+'/speed'], weightList(ww)));
end
track.set_states_weight_set(stateWeights);
%Set tracked states to be used in guess
track.set_apply_tracked_states_to_guess(true);
%Set unused state references to be allowed in case of file errors
%This is relevant given use of 3D on 2D model
track.set_allow_unused_references(true);
%Set position derivatives to be used as speeds
track.set_track_reference_position_derivatives(true);

%Set control weights
track.set_control_effort_weight(effortWeight);

%Set timing parameters
track.set_initial_time(initialTime);
track.set_final_time(finalTime);
track.set_mesh_interval(meshInterval);

%Customise the base tracking problem with relevant goals
study = track.initialize();
problem = study.updProblem();

%Adjust time bounds to allow for subtle fluctations in finish time
problem.setTimeBounds(initialTime, [finalTime - 0.1, finalTime + 0.1]);

%Update the control effort goal to a cost of transport type cost
effort = MocoControlGoal().safeDownCast(problem.updGoal('control_effort'));
effort.setDivideByDisplacement(true);

%Set to not highly penalise the lumbar actuator
%Seems necessary for a good solution to match GRF tracking
effort.setWeightForControl('/forceset/reserve_lumbar_extension', 0.001)

%Set an average speed goal based on sprinting data
%Use the pelvis_tx displacement to calculate speed
kinematicsTable = TimeSeriesTable('ik\\sprint_ik_states.sto');
startPos = kinematicsTable.getDependentColumn('/jointset/ground_pelvis/pelvis_tx/value').get(0);
endPos = kinematicsTable.getDependentColumn('/jointset/ground_pelvis/pelvis_tx/value').get(kinematicsTable.getNumRows()-1);
sprintSpeed = (endPos - startPos) / (finalTime - initialTime);
%Create and add speed goal
speedGoal = MocoAverageSpeedGoal('speed');
speedGoal.set_desired_average_speed(sprintSpeed);
speedGoal.setWeight(speedWeight);
problem.addGoal(speedGoal);

%Set symmetry constraint
%Create symmetry goal
symmetryGoal = MocoPeriodicityGoal('symmetryGoal');
%Set symmetric coordinate values (except for pelvis_tx) and speeds
for cc = 0:trackModel.getCoordinateSet().getSize()-1
    %Get current coordinate name
    coordName = char(trackModel.getCoordinateSet().get(cc).getName());
    %Set symmetry parameters based on coordinate
    if endsWith(coordName,'_r')
        %Set current state name
        currStateName = char(trackModel.getCoordinateSet().get(cc).getAbsolutePathString());
        %Joint angle value
        symmetryGoal.addStatePair(MocoPeriodicityGoalPair(...
            [currStateName,'/value'], regexprep([currStateName,'/value'],'_r','_l')));
        %Joint speed
        symmetryGoal.addStatePair(MocoPeriodicityGoalPair(...
            [currStateName,'/speed'], regexprep([currStateName,'/speed'],'_r','_l')));
    elseif endsWith(coordName,'_l')
        %Set current state name
        currStateName = char(trackModel.getCoordinateSet().get(cc).getAbsolutePathString());
        %Joint angle value
        symmetryGoal.addStatePair(MocoPeriodicityGoalPair(...
            [currStateName,'/value'], regexprep([currStateName,'/value'],'_l','_r')));
        %Joint speed
        symmetryGoal.addStatePair(MocoPeriodicityGoalPair(...
            [currStateName,'/speed'], regexprep([currStateName,'/speed'],'_l','_r')));
    elseif endsWith(coordName,'_extension') || ...
            endsWith(coordName,'_tilt') || ...
            endsWith(coordName,'_ty')
        %Set current state name
        currStateName = char(trackModel.getCoordinateSet().get(cc).getAbsolutePathString());
        %Joint angle value
        symmetryGoal.addStatePair(MocoPeriodicityGoalPair([currStateName,'/value']));
        %Joint speed
        symmetryGoal.addStatePair(MocoPeriodicityGoalPair([currStateName,'/speed']));
    end
end
%Add a symmetry pair for pelvis_tx speed
symmetryGoal.addStatePair(MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_tx/speed'));
%Add symmetry goal
symmetryGoal.setWeight(symmetryWeight);
problem.addGoal(symmetryGoal);

%Add contact tracking goal
%Set force names
forceNamesRightFoot = [{'forceset/contactHeel_r'};
    {'forceset/contactMidfoot_r'};
    {'forceset/contactToe_r'}];
forceNamesLeftFoot = [{'forceset/contactHeel_l'};
    {'forceset/contactMidfoot_l'};
    {'forceset/contactToe_l'}];

%Create contact tracking goal
contactGoal = MocoContactTrackingGoal('contactGoal', grfWeight);
%Set external loads
contactGoal.setExternalLoadsFile('data\\sprint_grf_2D.xml');
%Set force name groups
forceNames_r = StdVectorString();
forceNames_l = StdVectorString();
for ff = 1:length(forceNamesRightFoot)
    forceNames_r.add(forceNamesRightFoot(ff));
    forceNames_l.add(forceNamesLeftFoot(ff));
end
%Create and add tracking groups
trackRightGRF = MocoContactTrackingGoalGroup(forceNames_r, 'RightGRF');
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r');
contactGoal.addContactGroup(trackRightGRF);
trackLeftGRF = MocoContactTrackingGoalGroup(forceNames_l, 'LeftGRF');
trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l');
contactGoal.addContactGroup(trackLeftGRF);
%Set parameters
contactGoal.setProjection('plane');
contactGoal.setProjectionVector(Vec3(0, 0, 1));
%Add contact tracking goal
problem.addGoal(contactGoal);

%Add state bounds
%This process sets both bounds on the entire problem, as well as initial
%bounds based on the experimental data. Here we allow slight variation
%around both of these by extracting from the experimental data and allowing
%an amount of variation.
%We apply a special case to the pelvis _ty initial state based on some
%exploration that aligns the starting position with the floor
%We also apply a special case for the pelvis_tx to start at the same point
%as the experimental data to simplify the problem
for cc = 0:trackModel.getCoordinateSet().getSize()-1
    %Set current state name
    currStateName = [char(trackModel.getCoordinateSet().get(cc).getAbsolutePathString()),'/value'];
    %Get the current data from the kinematics table
    currStateData = kinematicsTable.getDependentColumn(currStateName).getAsMat();
    %Set overall limits by allowing 10% of the range either way
    currStateMax = max(currStateData);
    currStateMin = min(currStateData);
    currStateRange = abs(currStateMax - currStateMin);
    %Check to see if the proposed max and mins exceed typical
    %coordinate allowed values and reset if so
    %Max
    if (currStateMax + (currStateRange * 0.1)) > trackModel.getCoordinateSet().get(cc).getRangeMax()
        %Use the coordinate max
        maxLimit = trackModel.getCoordinateSet().get(cc).getRangeMax();
    else
        %Use the calculate
        maxLimit = currStateMax + (currStateRange * 0.1);
    end
    %Min
    if (currStateMin - (currStateRange * 0.1)) < trackModel.getCoordinateSet().get(cc).getRangeMin()
        %Use the coordinate max
        minLimit = trackModel.getCoordinateSet().get(cc).getRangeMin();
    else
        %Use the calculate
        minLimit = currStateMin - (currStateRange * 0.1);
    end
    %Set initial limits
    %Check if this is special pelvis_ty state
    if contains(currStateName,'pelvis_ty')
        %Set initial bounds to be within pre-specific points
        initialBounds = MocoInitialBounds(0.94,0.97);
        %Reset min and max limits from those calculated
        minLimit = 0.93; maxLimit = 0.99;
    elseif contains(currStateName,'pelvis_tx')
        %Set initial state to start at the same spot
        initialBounds = MocoInitialBounds(currStateData(1));
    else
        %Set initial state to be within 20% either way of starting point
        initialLimits = [currStateData(1) - (currStateData(1)*0.2),...
            currStateData(1) + (currStateData(1)*0.2)];
        %Need a check in place for negative coordinates as sometimes the
        %order of these values is the wrong way around
        if initialLimits(2) < initialLimits(1)
            minInitial = initialLimits(2); maxInitial = initialLimits(1);
        else
            minInitial = initialLimits(1); maxInitial = initialLimits(2);
        end
        %Check if the minimum & maximum initial bounds exceeds the already
        %approved ranges
        if minInitial < minLimit
            minInitial = minLimit;
        end
        if maxInitial > maxLimit
            maxInitial = maxLimit;
        end
        %Set the initial bounds object
        initialBounds = MocoInitialBounds(minInitial,maxInitial);
    end
    
    %Set the state bounds in the problem
    problem.setStateInfo(currStateName, ...
        MocoBounds(minLimit,maxLimit),...
        initialBounds);    
end

%Configure the solver
solver = MocoCasADiSolver.safeDownCast(study.updSolver());
solver.resetProblem(problem);
solver.set_optim_constraint_tolerance(1e-2) %%% probably a bit high, but doing so for speedy solution
solver.set_optim_convergence_tolerance(1e-2) %%% probably a bit high, but doing so for speedy solution

%Solve
%Solution with 1e-2 parameters takes ~1000 iterations and ~44 mins with 8 threads 
torqueTrackingSolution = study.solve();

%Option to visualise
if visualiseSolution
    study.visualize(torqueTrackingSolution);
end

%Remove the tracked states file as this is no longer needed
delete(dir('*_tracked_states.sto').name);

%Extract ground reaction forces from tracking simulation
%Set the contact points on each foot
contact_r = StdVectorString();
contact_l = StdVectorString();
contact_r.add('/forceset/contactHeel_r');
contact_r.add('/forceset/contactMidfoot_r');
contact_r.add('/forceset/contactToe_r');
contact_l.add('/forceset/contactHeel_l');
contact_l.add('/forceset/contactMidfoot_l');
contact_l.add('/forceset/contactToe_l');
%Create the table of forces
trackingGrfTable = opensimMoco.createExternalLoadsTableForGait(torqueModel,torqueTrackingSolution,contact_r,contact_l);

%Create folder to store if it isn't present
if ~exist('trackingSims', 'dir')
    mkdir('trackingSims')
end
if ~exist('trackingSims\\torqueDrivenTracking','dir')
    mkdir('trackingSims\\torqueDrivenTracking')
end

%Write solution files 
torqueTrackingSolution.write('trackingSims\\torqueDrivenTracking\\sprint_torqueDriven_tracking.sto');
STOFileAdapter.write(trackingGrfTable, 'trackingSims\\torqueDrivenTracking\\sprint_torqueDriven_tracking_grf.sto');

%Compare solution to experimental data
if compareToExp
    
    %Kinematics
    
    %Set coordinate states to work through
    plotStates = [{'/jointset/ground_pelvis/pelvis_ty/value'};
        {'/jointset/ground_pelvis/pelvis_tx/value'};
        {'/jointset/ground_pelvis/pelvis_tilt/value'}
        {'/jointset/back/lumbar_extension/value'};
        {'/jointset/hip_r/hip_flexion_r/value'};
        {'/jointset/hip_l/hip_flexion_l/value'};
        {'/jointset/knee_r/knee_angle_r/value'};
        {'/jointset/knee_l/knee_angle_l/value'};
        {'/jointset/ankle_r/ankle_angle_r/value'};
        {'/jointset/ankle_l/ankle_angle_l/value'};
        {'/jointset/acromial_r/arm_flex_r/value'};
        {'/jointset/acromial_l/arm_flex_l/value'};
        {'/jointset/elbow_r/elbow_flex_r/value'};
        {'/jointset/elbow_l/elbow_flex_l/value'};
        ];
    
    %Convert solution to states table
    solutionTable = torqueTrackingSolution.exportToStatesTable();
    
    %Set plot
    figure; set(gcf,'units','normalized','position',[0.01 0.01 0.6 0.9]);
    
    %Set experimental kinematics time vector
    for tt = 0:kinematicsTable.getIndependentColumn().size()-1
        expTime(tt+1) = kinematicsTable.getIndependentColumn().get(tt).doubleValue() - ...
            kinematicsTable.getIndependentColumn().get(0).doubleValue();
    end
    
    %Set solution time vector
    for tt = 0:solutionTable.getIndependentColumn().size()-1
        solTime(tt+1) = solutionTable.getIndependentColumn().get(tt).doubleValue() - ...
            solutionTable.getIndependentColumn().get(0).doubleValue();
    end
    
    %Loop through states
    for nn = 1:length(plotStates)
        %Set subplot
        subplot(4,4,nn); hold on
        %Turn box on
        ax = gca;
        box(ax,'on');
        ax.LineWidth = 1;
        set(ax, 'Layer', 'Top')
        %White background
        set(gcf,'Color','w');
        %Get experimental data from kinematic table
        plot(expTime, kinematicsTable.getDependentColumn(plotStates{nn}).getAsMat(),...
            'k', 'LineWidth', 1.5);
        %Get tracking solution data
        plot(solTime, solutionTable.getDependentColumn(plotStates{nn}).getAsMat(),...
            'r--', 'LineWidth', 1.5);
        %Add plot labels
        %X-label
        xlabel('Time (s)');
        %Y-label
        if contains(plotStates{nn},'pelvis_t')
            ylabel('Pos. (m)')
        else
            ylabel('Angle (rad)')
        end
        %Title
        splitState = split(plotStates{nn},'/');
        title(regexprep(splitState{4},'_',' '));        
    end
    
    %Save figure
    %Note this uses export fig rather than built in Matlab functions
    export_fig trackingSims\\torqueDrivenTracking\\exp-vs-tracked-kinematics.png -m1
    
    %Close figure
    close(gcf);
    
    %GRFs
    
    %Get experimental GRFs
    grfsTable = TimeSeriesTable('data\\sprint_grf_2D.mot');
    
    %Set grfs to plot states to work through
    plotGrfs = [{'ground_force_r_vx'};
        {'ground_force_r_vy'};
        ];
    
    %Set plot
    figure; set(gcf,'units','normalized','position',[0.1 0.1 0.4 0.3]);
    
    %Set experimental grfs time vector
    %Get the full time
    for tt = 0:grfsTable.getIndependentColumn().size()-1
        grfTime(tt+1) = grfsTable.getIndependentColumn().get(tt).doubleValue();
    end
    %Identify the indices that relate to the start and end times from kinematics
    grfStartInd = find(grfTime > startTime,1) - 1;
    grfEndInd = find(grfTime > finalTime,1) - 1;
    %Reset time to match-up with indices and zero
    grfTime = grfTime(grfStartInd:grfEndInd) - grfTime(grfStartInd);
    
    %Set solution time vector
    for tt = 0:trackingGrfTable.getIndependentColumn().size()-1
        solTime(tt+1) = trackingGrfTable.getIndependentColumn().get(tt).doubleValue() - ...
            trackingGrfTable.getIndependentColumn().get(0).doubleValue();
    end
        
    %Loop through states
    for nn = 1:length(plotGrfs)
        %Set subplot
        subplot(1,2,nn); hold on
        %Turn box on
        ax = gca;
        box(ax,'on');
        ax.LineWidth = 1;
        set(ax, 'Layer', 'Top')
        %White background
        set(gcf,'Color','w');
        %Get experimental data from kinematic table
        grfData = grfsTable.getDependentColumn(plotGrfs{nn}).getAsMat();
        plot(grfTime, grfData(grfStartInd:grfEndInd),...
            'k', 'LineWidth', 1.5);
        %Get tracking solution data
        plot(solTime, trackingGrfTable.getDependentColumn(plotGrfs{nn}).getAsMat(),...
            'r--', 'LineWidth', 1.5);
        %Add plot labels
        %X-label
        xlabel('Time (s)');
        %Y-label
        ylabel('Force (N)');
        %Title
        if strcmp(plotGrfs{nn}(end-1:end),'vx')
            title('Ant. / Post. GRFs')
        else
            title('Vertical GRFs')
        end
    end
    
    %Save figure
    %Note this uses export fig rather than built in Matlab functions
    export_fig trackingSims\\torqueDrivenTracking\\exp-vs-tracked-grfs.png -m1
    
    %Close figure
    close(gcf);

end

%% Muscle driven inverse simulation

%This is a muscle driven simulation of the half gait cycle, using the
%kinematics from the previous torque tracking simulation. As an inverse
%simulation, the kinematics are not tracked but held consistent - with the
%muscle forces required to drive the motion estimated.

%%%% TODO: consider more complex metabolics cost function???

%%%% TODO: model processing seems a little repetitive and could be cleaner

%Set simulation parameters

%Timing and mesh data
initialTime = Storage('trackingSims\\torqueDrivenTracking\\sprint_torqueDriven_tracking.sto').getFirstTime();
finalTime = Storage('trackingSims\\torqueDrivenTracking\\sprint_torqueDriven_tracking.sto').getLastTime();
% % % initialTime = Storage('ik\\jog_ik.mot').getFirstTime();
% % % finalTime = Storage('ik\\jog_ik.mot').getLastTime();
duration = finalTime - initialTime;
meshNo = 100; %50; %%% mesh number based on some work by Falisse et al. recommending 50 for half gait cycle 
meshInterval = duration / meshNo;

%General parameters
visualiseSolution = false;
compareToExp = false;

%Muscle parameters
passiveForces = false;
implicitTendonCompliance = true;
tendonDynamics = true;
contractVelScale = 5; %%%2.5; %upscale contractile velocity for sprinting
maxForceScale = 5; %%%2; %upscale force for sprinting

%Create the model processor for the tracking problem

%Load the model
muscleModel = Model('models\\muscles2DModel.osim');

%Update contractile velocity maximum to assist with sprinting sim
for mm = 0:muscleModel.getMuscles().getSize()-1
    %Get the current muscle
    currMusc = muscleModel.getMuscles().get(mm);
    %Get current max contraction velocity
    currVel = currMusc.getMaxContractionVelocity();
    %Set new max contraction velocity
    currMusc.set_max_contraction_velocity(currVel*contractVelScale);
end

%Finalise connections
muscleModel.finalizeConnections();

%Put the model in a processor
muscleModelProcessor = ModelProcessor(muscleModel);

%Scale the maximum isometric force of muscles to assist with sprinting sim
muscleModelProcessor.append(ModOpScaleMaxIsometricForce(maxForceScale));

%Disable tendon compliance (will be re-enabled for certain muscles later)
muscleModelProcessor.append(ModOpIgnoreTendonCompliance());

%Adjust fiber damping as per Dembia et al. Moco paper
muscleModelProcessor.append(ModOpFiberDampingDGF(0.01));

%Turn off passive forces if appropriate
if ~passiveForces
    muscleModelProcessor.append(ModOpIgnorePassiveFiberForcesDGF());
end

%Modify active force width
%This appears most important for "de-stiffening" muscles to allow muscle
%to function well over sprinting joint ranges
muscleModelProcessor.append(ModOpScaleActiveFiberForceCurveWidthDGF(1.5));

%Process model for additional parameters
inverseModel = muscleModelProcessor.process();
inverseModel.initSystem();

%Enable tendon compliance in the gastrocnemius and soleus
if tendonDynamics
    muscles = inverseModel.updMuscles();
    for mm = 0:muscles.getSize()-1
        currMusc = DeGrooteFregly2016Muscle.safeDownCast(muscles.get(mm));
        muscName = char(currMusc.getName());
        %Enable tendon compliance dynamics in the plantarflexors
        if contains(muscName,'gastroc') || contains(muscName,'soleus')
            currMusc.set_ignore_tendon_compliance(false)
        end
    end
end

%Get the track model as a processor object for extra operations
inverseModelProcessor = ModelProcessor(inverseModel);

%Set model to use implicit tendon compliance if appropriate
if implicitTendonCompliance
    inverseModelProcessor.append(ModOpUseImplicitTendonComplianceDynamicsDGF());
end

%Process and add torque actuators to un-actuated coordinates
inverseModel = inverseModelProcessor.process();
%Set optimal force and max torque
%Split this into ones that are torque drive coordinates vs. reserves
optimalForce = 100;
optimalForce_res = 1;
maxTorque = Inf;
%Add torque actuators
for coordinate = 0:inverseModel.getCoordinateSet().getSize()-1
    %Get coordinate name
    coord = char(inverseModel.getCoordinateSet().get(coordinate));
    %Check if a torque drive or reserve
    if strcmp(coord,'lumbar_extension') || ...
            strcmp(coord,'arm_flex_r') || ...
            strcmp(coord,'arm_flex_l') || ...
            strcmp(coord,'elbow_flex_r') || ...
            strcmp(coord,'elbow_flex_l')
            %Add actuator
            addReserve(inverseModel, coord, optimalForce, maxTorque);
    else
        %Add actuator
        addReserve(inverseModel, coord, optimalForce_res, maxTorque);
    end
end


%%%%%% Turn off contact force geometry
for ff = 0:inverseModel.updForceSet().getSize()-1
    %Check for contact force
    if contains(char(inverseModel.updForceSet().get(ff).getName()),'contact')
        %Turn off force
        inverseModel.updForceSet().get(ff).set_appliesForce(false)
    end
end

%Finalise connections
inverseModel.finalizeConnections();

%Construct the inverse object and set basic parameters
inverse = MocoInverse();
inverse.setName('sprintInverse_muscleDriven');
% % % inverse.setName('jogInverse_muscleDriven');

%%%%% Set model with external loads
inverseProcessor = ModelProcessor(inverseModel);
inverseProcessor.append(ModOpAddExternalLoads('data\\sprint_grf_2D.xml'));
% % % inverseProcessor.append(ModOpAddExternalLoads('data\\jog_gr_2Df.xml'));

%Set model processor
inverse.setModel(inverseProcessor);

%Initial time, final time, and mesh interval.
inverse.set_initial_time(initialTime);
inverse.set_final_time(finalTime);
inverse.set_mesh_interval(meshInterval);

%Set kinematic data and parameters
%Get the states table from the torque driven solution
torqueDrivenStates = MocoTrajectory(...
    'trackingSims\\torqueDrivenTracking\\sprint_torqueDriven_tracking.sto').exportToStatesTable();
%Get the column labels to loop through
statesColumns = torqueDrivenStates.getColumnLabels();
%Loop through and remove the column labels that aren't kinematic values
for cc = 0:statesColumns.size()-1
    if ~endsWith(char(statesColumns.get(cc)),'/value')
        %Remove non kinematic value column
        torqueDrivenStates.removeColumn(char(statesColumns.get(cc)));        
    end    
end
%Set table processor
tableProcessor = TableProcessor(torqueDrivenStates);
% % % tableProcessor = TableProcessor('ik\\jog_ik_states.sto');
% % % tableProcessor.append(TabOpLowPassFilter(12)); %pretty clean from past tracking sim
tableProcessor.append(TabOpUseAbsoluteStateNames());
inverse.setKinematics(tableProcessor);

%Set unused state references to be allowed in case of any issues
inverse.set_kinematics_allow_extra_columns(true);

%Solve the inverse
inverseMuscleSolution = inverse.solve();


%%%%% There's clearly an issue with solving this problem - it's either the
%%%%% mesh interval being too small, the muscle parameters having an issue,
%%%%% or the torques/reserves not being strong enough...

%%%%% inf_du (dual infeasibility) seems to be the biggest issue...whatever
%%%%% that relates to is probably causing the problem...
    %%%% same problem with jog...
    %%%% even very early part of sprint trial generates high infeasibility
    %%%% with increased strength and contractile velocity...

%%%% Test making much smaller mesh interval, mesh interval from 50 to 200...
%%%% Wasn't setting properly, so only had 25 nodes!
    %%%% Node doesn't make a huge difference
    
%%%% Ankle and hip flexion reserves high, pelvis tx reserves high
%%%% Strength as an issue? Residuals...?

%%%% Contact force geometry could be messing around with this, rather than
%%%% the external loads --- as this is different to typical MocoInverse
%%%% solutions...
    %%%% this works better (tried without tendon dynamics)
    %%%% still some decent saturation of muscle signals
    %%%% reserves for residuals are pretty high --- maybe do need that RRA step... 
        %%%% this would defeat purpose of tracking kinematics
    %%%% knee and hip flexion reserves are probably the highest
        %%%% fix tendon slack length --- shouldn't matter given 
    %%%% with 100 mesh intervals still kinda shark fin like


% % % inverseMuscleSolution.getMocoSolution().unseal().write('test.sto');







%% Muscle driven tracking simulation

%This is a muscle driven simulation of the half gait cycle. Considering we
%took the time to run the torque driven simulation to generate kinematics
%that track the GRFs with the contact spheres, we will use these data to
%drive this simulation.

%%%% TODO: consider more complex metabolics cost function???

%%%% TODO: model processing seems a little repetitive and could be cleaner

%Set simulation parameters
    
%Timing and mesh data
initialTime = Storage('trackingSims\\torqueDrivenTracking\\sprint_torqueDriven_tracking.sto').getFirstTime();
finalTime = Storage('trackingSims\\torqueDrivenTracking\\sprint_torqueDriven_tracking.sto').getLastTime();
duration = finalTime - initialTime;
meshNo = 50; %%% mesh number based on some work by Falisse et al. recommending 50 for half gait cycle 
meshInterval = duration / meshNo;

%Weights
trackingWeight = 10;
effortWeight = 0.1;
grfWeight = 2.5;
speedWeight = 1;
symmetryWeight = 1;

%General parameters
visualiseSolution = false;
compareToExp = false;

%Muscle parameters
passiveForces = false; %%%true;
implicitTendonCompliance = true; %%%false;
tendonDynamics = true;
contractVelScale = 5; %upscale contractile velocity for sprinting
maxForceScale = 10; %upscale force for sprinting

%Create the model processor for the tracking problem

%Load the model
muscleModel = Model('models\\muscles2DModel.osim');

%Update contractile velocity maximum to assist with sprinting sim
for mm = 0:muscleModel.getMuscles().getSize()-1
    %Get the current muscle
    currMusc = muscleModel.getMuscles().get(mm);
    %Get current max contraction velocity
    currVel = currMusc.getMaxContractionVelocity();
    %Set new max contraction velocity
    currMusc.set_max_contraction_velocity(currVel*contractVelScale);
end

%Finalise connections
muscleModel.finalizeConnections();

%Put the model in a processor
muscleModelProcessor = ModelProcessor(muscleModel);

%Scale the maximum isometric force of muscles to assist with sprinting sim
muscleModelProcessor.append(ModOpScaleMaxIsometricForce(maxForceScale));

%Disable tendon compliance (will be re-enabled for certain muscles later)
muscleModelProcessor.append(ModOpIgnoreTendonCompliance());

%Adjust fiber damping as per Dembia et al. Moco paper
muscleModelProcessor.append(ModOpFiberDampingDGF(0.01));

%Turn off passive forces if appropriate
if ~passiveForces
    muscleModelProcessor.append(ModOpIgnorePassiveFiberForcesDGF());
end

%Modify active force width
%This appears most important for "de-stiffening" muscles to allow muscle
%to function well over sprinting joint ranges
muscleModelProcessor.append(ModOpScaleActiveFiberForceCurveWidthDGF(1.5));
% % % muscleModelProcessor.append(ModOpScaleActiveFiberForceCurveWidthDGF(3.0));

%Process model for additional parameters
trackModel = muscleModelProcessor.process();
trackModel.initSystem();

%Enable tendon compliance in the gastrocnemius and soleus
if tendonDynamics
    muscles = trackModel.updMuscles();
    for mm = 0:muscles.getSize()-1
        currMusc = DeGrooteFregly2016Muscle.safeDownCast(muscles.get(mm));
        muscName = char(currMusc.getName());
        %Enable tendon compliance dynamics in the plantarflexors
        if contains(muscName,'gastroc') || contains(muscName,'soleus')
            currMusc.set_ignore_tendon_compliance(false)
        end
    end
end

%Get the track model as a processor object for extra operations
trackModelProcessor = ModelProcessor(trackModel);

%Set model to use implicit tendon compliance if appropriate
if implicitTendonCompliance
    trackModelProcessor.append(ModOpUseImplicitTendonComplianceDynamicsDGF());
end

%Process model for additional parameters
trackModel = trackModelProcessor.process();
trackModel.initSystem();

% % % %Add torque actuators to support model
% % % %Set list to add reserves to
% % % reserveTorqueList = [{'lumbar_extension'};
% % %     {'arm_flex_r'};
% % %     {'arm_flex_l'};
% % %     {'elbow_flex_r'};
% % %     {'elbow_flex_l'};
% % %     {'pelvis_tx'};
% % %     {'pelvis_ty'};
% % %     {'pelvis_tilt'}];
% % % %Set optimal force and max torque
% % % optimalForce = 100; %for torque actuated coordinates
% % % optimalForce_res = 1; %for residual coordinates
% % % maxTorque = Inf;
% % % %Add torque actuators
% % % for cc = 1:length(reserveTorqueList)
% % %     %Check which type of actuator to apply
% % %     if contains(reserveTorqueList{cc},'pelvis')
% % %         %Add actuator with low optimal force
% % %         addReserve(trackModel, reserveTorqueList{cc}, optimalForce_res, maxTorque);
% % %     else
% % %         %Add actuator with high optimal force
% % %         addReserve(trackModel, reserveTorqueList{cc}, optimalForce, maxTorque);
% % %     end    
% % % end


%%%%%% Possible that reserve actuators are needed across the board...
%%%%%% This will also help see where problems are...

%Add torque actuators to support model
%Set list to add reserves to based on coordinate list
%Set optimal force and max torque
optimalForce = 100;
optimalForce_res = 1;
maxTorque = Inf;
%Add torque actuators
for coordinate = 0:trackModel.getCoordinateSet().getSize()-1
    %Get coordinate name
    coord = char(trackModel.getCoordinateSet().get(coordinate));
    %Check which level of optimal force to add
    if contains(coord, 'lumbar') || ...
            contains(coord, 'arm') || ...
            contains(coord, 'elbow')
        %Add higher force actuator
        addReserve(trackModel, coord, optimalForce, maxTorque);
    else
        %Add lower force actuator
        addReserve(trackModel, coord, optimalForce_res, maxTorque);
    end
end

%Finalise connections
trackModel.finalizeConnections();

%Construct the tracking object and set basic parameters
track = MocoTrack();
track.setName('sprintTracking_muscleDriven');
track.setModel(ModelProcessor(trackModel));

%Set kinematic data and parameters
%Get the states table from the torque driven solution
torqueDrivenStates = MocoTrajectory(...
    'trackingSims\\torqueDrivenTracking\\sprint_torqueDriven_tracking.sto').exportToStatesTable();
%Get the column labels to loop through
statesColumns = torqueDrivenStates.getColumnLabels();
%Loop through and remove the column labels that aren't kinematic values
for cc = 0:statesColumns.size()-1
    if ~endsWith(char(statesColumns.get(cc)),'/value')
        %Remove non kinematic value column
        torqueDrivenStates.removeColumn(char(statesColumns.get(cc)));        
    end    
end
%Set table processor
tableProcessor = TableProcessor(torqueDrivenStates);
% % % tableProcessor.append(TabOpLowPassFilter(12)); %pretty clean from past tracking sim
tableProcessor.append(TabOpUseAbsoluteStateNames());
track.setStatesReference(tableProcessor);
%The global tracking states weight is enough here, as the torque driven
%tracking problem gives us confidence that these kinematics align with the
%GRF tracking problem - therefore individual state weights aren't used
track.set_states_global_tracking_weight(trackingWeight);
%Set tracked states to be used in guess
track.set_apply_tracked_states_to_guess(true);
%Set unused state references to be allowed in case of any issues
track.set_allow_unused_references(true);
%Set position derivatives to be used as speeds
track.set_track_reference_position_derivatives(true);

%Set control weights
track.set_control_effort_weight(effortWeight);

%Set timing parameters
track.set_initial_time(initialTime);
track.set_final_time(finalTime);
track.set_mesh_interval(meshInterval);

%Customise the base tracking problem with relevant goals
study = track.initialize();
problem = study.updProblem();

%Set solid time bounds from the torque driven simulation
problem.setTimeBounds(initialTime, finalTime);

%Update the control effort goal to a cost of transport type cost
effort = MocoControlGoal().safeDownCast(problem.updGoal('control_effort'));
effort.setDivideByDisplacement(true);

% % % %Set relevant penalties on torque actuators
% % % %High for residual, low for the reserve torque actuators
% % % for ff = 0:trackModel.updForceSet().getSize()-1
% % %     %Check for residual force
% % %     if contains(char(trackModel.updForceSet().get(ff).getName()),'residual')
% % %         %Add high weight to tracking
% % %         effort.setWeightForControl(char(trackModel.updForceSet().get(ff).getAbsolutePathString()),...
% % %             10)
% % %     %Check for reserve
% % %     elseif contains(char(trackModel.updForceSet().get(ff).getName()),'reserve')
% % %         %Add low weight to tracking
% % %         effort.setWeightForControl(char(trackModel.updForceSet().get(ff).getAbsolutePathString()),...
% % %             0.001)
% % %     end
% % % end

%Set relevant penalties on torque actuators
%High for residual, low for the reserve torque actuators
for ff = 0:trackModel.updForceSet().getSize()-1
    %Check for residual force
    if contains(char(trackModel.updForceSet().get(ff).getName()),'residual')
        %Add high weight to tracking
        effort.setWeightForControl(char(trackModel.updForceSet().get(ff).getAbsolutePathString()),...
            10)
    %Check for upper body torque
    elseif contains(char(trackModel.updForceSet().get(ff).getName()),'lumbar') || ...
            contains(char(trackModel.updForceSet().get(ff).getName()),'arm') || ...
            contains(char(trackModel.updForceSet().get(ff).getName()),'elbow')
        %Add low weight to tracking
        effort.setWeightForControl(char(trackModel.updForceSet().get(ff).getAbsolutePathString()),...
            0.001)
    %Check for other reserves
    elseif contains(char(trackModel.updForceSet().get(ff).getName()),'reserve')
        %Add moderate weight to tracking
        effort.setWeightForControl(char(trackModel.updForceSet().get(ff).getAbsolutePathString()),...
            0.1)
    end
end

%Set an average speed goal based on sprinting data
%Use the pelvis_tx displacement to calculate speed
startPos = torqueDrivenStates.getDependentColumn('/jointset/ground_pelvis/pelvis_tx/value').get(0);
endPos = torqueDrivenStates.getDependentColumn('/jointset/ground_pelvis/pelvis_tx/value').get(torqueDrivenStates.getNumRows()-1);
sprintSpeed = (endPos - startPos) / (finalTime - initialTime);
%Create and add speed goal
speedGoal = MocoAverageSpeedGoal('speed');
speedGoal.set_desired_average_speed(sprintSpeed);
speedGoal.setWeight(speedWeight);
problem.addGoal(speedGoal);

%Set symmetry constraint
%Create symmetry goal
symmetryGoal = MocoPeriodicityGoal('symmetryGoal');
%Set symmetric coordinate values (except for pelvis_tx) and speeds
for cc = 0:trackModel.getCoordinateSet().getSize()-1
    %Get current coordinate name
    coordName = char(trackModel.getCoordinateSet().get(cc).getName());
    %Set symmetry parameters based on coordinate
    if endsWith(coordName,'_r')
        %Set current state name
        currStateName = char(trackModel.getCoordinateSet().get(cc).getAbsolutePathString());
        %Joint angle value
        symmetryGoal.addStatePair(MocoPeriodicityGoalPair(...
            [currStateName,'/value'], regexprep([currStateName,'/value'],'_r','_l')));
        %Joint speed
        symmetryGoal.addStatePair(MocoPeriodicityGoalPair(...
            [currStateName,'/speed'], regexprep([currStateName,'/speed'],'_r','_l')));
    elseif endsWith(coordName,'_l')
        %Set current state name
        currStateName = char(trackModel.getCoordinateSet().get(cc).getAbsolutePathString());
        %Joint angle value
        symmetryGoal.addStatePair(MocoPeriodicityGoalPair(...
            [currStateName,'/value'], regexprep([currStateName,'/value'],'_l','_r')));
        %Joint speed
        symmetryGoal.addStatePair(MocoPeriodicityGoalPair(...
            [currStateName,'/speed'], regexprep([currStateName,'/speed'],'_l','_r')));
    elseif endsWith(coordName,'_extension') || ...
            endsWith(coordName,'_tilt') || ...
            endsWith(coordName,'_ty')
        %Set current state name
        currStateName = char(trackModel.getCoordinateSet().get(cc).getAbsolutePathString());
        %Joint angle value
        symmetryGoal.addStatePair(MocoPeriodicityGoalPair([currStateName,'/value']));
        %Joint speed
        symmetryGoal.addStatePair(MocoPeriodicityGoalPair([currStateName,'/speed']));
    end
end
%Add a symmetry pair for pelvis_tx speed
symmetryGoal.addStatePair(MocoPeriodicityGoalPair('/jointset/ground_pelvis/pelvis_tx/speed'));
%Add symmetry goal
symmetryGoal.setWeight(symmetryWeight);
problem.addGoal(symmetryGoal);

%%%%% TODO: can add muscle state symmetry goals...

%Add contact tracking goal
%Set force names
forceNamesRightFoot = [{'forceset/contactHeel_r'};
    {'forceset/contactMidfoot_r'};
    {'forceset/contactToe_r'}];
forceNamesLeftFoot = [{'forceset/contactHeel_l'};
    {'forceset/contactMidfoot_l'};
    {'forceset/contactToe_l'}];

%Create contact tracking goal
contactGoal = MocoContactTrackingGoal('contactGoal', grfWeight);
%Set external loads
contactGoal.setExternalLoadsFile('data\\sprint_grf_2D.xml');
%Set force name groups
forceNames_r = StdVectorString();
forceNames_l = StdVectorString();
for ff = 1:length(forceNamesRightFoot)
    forceNames_r.add(forceNamesRightFoot(ff));
    forceNames_l.add(forceNamesLeftFoot(ff));
end
%Create and add tracking groups
trackRightGRF = MocoContactTrackingGoalGroup(forceNames_r, 'RightGRF');
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r');
contactGoal.addContactGroup(trackRightGRF);
trackLeftGRF = MocoContactTrackingGoalGroup(forceNames_l, 'LeftGRF');
trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l');
contactGoal.addContactGroup(trackLeftGRF);
%Set parameters
contactGoal.setProjection('plane');
contactGoal.setProjectionVector(Vec3(0, 0, 1));
%Add contact tracking goal
problem.addGoal(contactGoal);

%Add state bounds
%This process is a little tighter than the original tracking simulation.
%Here we restrict the states to their exact initial starting position from
%the torque driven simulation, and also set the max and min bounds to that
%from this tracking simulation too.
for cc = 0:trackModel.getCoordinateSet().getSize()-1
    %Set current state name
    currStateName = [char(trackModel.getCoordinateSet().get(cc).getAbsolutePathString()),'/value'];
    %Get the current data from the states table
    currStateData = torqueDrivenStates.getDependentColumn(currStateName).getAsMat();
    %Set the state bounds in the problem
    problem.setStateInfo(currStateName, ...
        MocoBounds(min(currStateData),max(currStateData)),...
        MocoInitialBounds(currStateData(1)));
end

%Set reasonable bounds for tendon force
problem.setStateInfoPattern('/forceset/.*/normalized_tendon_force', [0, 1.8], [], [])

%Configure the solver
solver = MocoCasADiSolver.safeDownCast(study.updSolver());
solver.resetProblem(problem);
solver.set_optim_constraint_tolerance(1e-2) %%% probably a bit high, but doing so for speedy solution
solver.set_optim_convergence_tolerance(1e-2) %%% probably a bit high, but doing so for speedy solution

%Set the normalized tendon forces if not loading initial guess from file
if tendonDynamics

    %Get the current default-ish guess
    guess = solver.getGuess();
    numRows = guess.getNumTimes();

    %Set the tendon force to reasonable values in initial guess to avoid a bad initial guess
    trackModel.initSystem();
    stateNames = trackModel.getStateVariableNames();
    for ii = 0:trackModel.getNumStateVariables()-1
        currState = stateNames.get(ii);
        if contains(char(currState),'normalized_tendon_force')
            guess.setState(char(currState), ones(numRows,1)*0.2);
        end
    end

    %Set guess
    solver.setGuess(guess)
    
end

%Solve
%Solution with 1e-2 parameters takes 505 iterations and 65 mins with 8 threads (no passive forces; no joint reserves)
muscleTrackingSolution = study.solve();

%%%% muscle tracking sim seems to suffer similar issue to inverse -- really
%%%% high inf_du - suggests something infeasible with problem???
    %%%% residual/reserve actuators seem to be the issue here --- it's
    %%%% possible this simulation still needs some assistance
        %%%% adding residual actuators in helps --- probably also need some
        %%%% reserves at joints if the muscles aren't quite right to solve
        %%%% the problem...
        
        %%%% residuals are pretty massive at the beginning, but quite low
        %%%% from that point on...
        
            %%%% definitely it - can solve when these are added, motion
            %%%% still looks a litte weird, perhaps cause reserves are
            %%%% needed - but not too bad overall
                %%%% still done without passive forces...
                %%%% iliopsoas_r fully activated the whole time
                
        %%%% used past guess when tracking with passive forces...
            %%%% passive forces seem to be a real struggle, perhaps
            %%%% coordinate forces (i.e. UMocoD approach) is better? Though
            %%%% this may still wreak some havoc...
            
        %%%% added joint reserves as there still seemed to be problems with
        %%%% the muscles achieving the desired states and hence GRF
        %%%% tracking too...
            %%%% stopped after about 1500 iterations, was progressively
            %%%% getting better but still just taking too many
            %%%% iterations...
        

%Option to visualise
if visualiseSolution
    study.visualize(muscleTrackingSolution);
end

%Remove the tracked states file as this is no longer needed
delete(dir('*_tracked_states.sto').name);

%Extract ground reaction forces from tracking simulation
%Set the contact points on each foot
contact_r = StdVectorString();
contact_l = StdVectorString();
contact_r.add('/forceset/contactHeel_r');
contact_r.add('/forceset/contactMidfoot_r');
contact_r.add('/forceset/contactToe_r');
contact_l.add('/forceset/contactHeel_l');
contact_l.add('/forceset/contactMidfoot_l');
contact_l.add('/forceset/contactToe_l');
%Create the table of forces
trackingGrfTable = opensimMoco.createExternalLoadsTableForGait(trackModel,muscleTrackingSolution,contact_r,contact_l);

%Create folder to store if it isn't present
if ~exist('trackingSims', 'dir')
    mkdir('trackingSims')
end
if ~exist('trackingSims\\muscleDrivenTracking','dir')
    mkdir('trackingSims\\muscleDrivenTracking')
end

%Write solution files 
muscleTrackingSolution.write('trackingSims\\muscleDrivenTracking\\sprint_muscleDriven_tracking.sto');
STOFileAdapter.write(trackingGrfTable, 'trackingSims\\muscleDrivenTracking\\sprint_muscleDriven_tracking_grf.sto');

%%

















