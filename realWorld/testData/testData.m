%% Muscle driven inverse simulation

%%%%% This function crashes MATLAB!!!!

import org.opensim.modeling.*

addpath('..\\functions');

%This is a muscle driven simulation of the half gait cycle, using the
%kinematics from the previous torque tracking simulation. As an inverse
%simulation, the kinematics are not tracked but held consistent - with the
%muscle forces required to drive the motion estimated.

%%%% TODO: consider more complex metabolics cost function???

%%%% TODO: model processing seems a little repetitive and could be cleaner

%Set simulation parameters

%Timing and mesh data
initialTime = Storage('referenceCoordinates.sto').getFirstTime();
finalTime = Storage('referenceCoordinates.sto').getLastTime();
% % % initialTime = Storage('ik\\jog_ik.mot').getFirstTime();
% % % finalTime = Storage('ik\\jog_ik.mot').getLastTime();
duration = finalTime - initialTime;
meshNo = 50;
meshInterval = duration / meshNo;

%General parameters
visualiseSolution = false;
compareToExp = false;

%Muscle parameters
passiveForces = false;
implicitTendonCompliance = true;
tendonDynamics = true;
contractVelScale = 2.5; %upscale contractile velocity for sprinting
maxForceScale = 2; %upscale force for sprinting

%Create the model processor for the tracking problem

%Load the model
muscleModel = Model('..\\models\\muscles2DModel.osim');

%Update contractile velocity maximum to assist with sprinting sim
for mm = 0:muscleModel.getMuscles().getSize()-1
    %Get the current muscle
    currMusc = muscleModel.getMuscles().get(mm);
    %Get current max contraction velocity
    currVel = currMusc.getMaxContractionVelocity();
    %Set new max contraction velocity
    currMusc.set_max_contraction_velocity(currVel*contractVelScale);
end

%Lock arm coordinates
muscleModel.getCoordinateSet().get('arm_flex_r').set_locked(true);
muscleModel.getCoordinateSet().get('arm_flex_l').set_locked(true);
muscleModel.getCoordinateSet().get('elbow_flex_r').set_locked(true);
muscleModel.getCoordinateSet().get('elbow_flex_l').set_locked(true);

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

%Weld joints for locked coordinates
muscleModelProcessor.append(ModOpReplaceJointsWithWelds());

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

%Get the model as a processor object for extra operations
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
    if strcmp(coord,'lumbar_extension')
% % %             strcmp(coord,'arm_flex_r') || ...
% % %             strcmp(coord,'arm_flex_l') || ...
% % %             strcmp(coord,'elbow_flex_r') || ...
% % %             strcmp(coord,'elbow_flex_l')
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
inverse.setName('testData_muscleDriven_noLoads');
% % % inverse.setName('jogInverse_muscleDriven');

%%%%% Set model with external loads
%%%% External loads aren't right, basically running this problem without
%%%% external loads, so the residual forces should be high...

%Set model processor
inverse.setModel(ModelProcessor(inverseModel));

%Initial time, final time, and mesh interval.
inverse.set_initial_time(initialTime);
inverse.set_final_time(finalTime);
inverse.set_mesh_interval(meshInterval);

%Set kinematic data and parameters
%Set table processor
tableProcessor = TableProcessor('referenceCoordinates.sto');
tableProcessor.append(TabOpLowPassFilter(6));
tableProcessor.append(TabOpUseAbsoluteStateNames());
inverse.setKinematics(tableProcessor);

%Set unused state references to be allowed in case of any issues
inverse.set_kinematics_allow_extra_columns(true);

%Solve the inverse
inverseMuscleSolution = inverse.solve();