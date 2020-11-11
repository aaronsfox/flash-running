%% This script performs a muscle-driven predictive simulation of the running.
%  The code here is adapted from the code provided by Falisse et al.
%  (2019). Comments are made throughout for the different aspects of code,
%  but I would defer to the original publication for a thorough
%  explanation. The initial guess for the predictive simulation is based
%  off tracking simulations of data from Dorn et al. (2012).
%
%  This code does have a dependency on CasADi (https://web.casadi.org/),
%  hence this must be installed with MATLAB for the code to work.
%
%  References:
%
%  Dorn et al. (2012). Muscular strategy shift in human running: Dependence
%  of running speed on hip and ankle muscle performance. J Exp Biol, 215:
%  1944-1956.
%
%  Falisse et al. (2019). Rapid predictive simulations with complex
%  musculoskeletal models suggest that diverse healthy and pathological
%  human gaits can emerge from similar control strategies. J R Soc
%  Interface, 16: 20190402.

%%%%% NOTE: first attempt wasn't too bad --- didn't converge and produced
%%%%% more of a skipping pattern. There may be some issues with the
%%%%% bounds/periodicity goals that is causing the motion; or a better
%%%%% initial guess with muscle activations etc. included may lead to a
%%%%% better solution.
%%%%%
%%%%% Elbow joint bounds need to be fixed as these went into
%%%%% hyperextension. Setting elbow angle to 90 degrees may be a valid
%%%%% option just to simplify things...

%% Settings

%%%%% NOTE: some weights and settings changed...

%Import CasADi
import casadi.*

%Import OpenSim
import org.opensim.modeling.*

%Set main path
mainPath = pwd;

%Set weights and options
%%%%% TODO: may not be appropriate for jog/sprinting
v_tgt = 3.5; % average speed
tol_ipopt = 4; % NLP error tolerance: 1*10^(-value)
N = 50; % number of mesh intervals
W.E = 500; % weight metabolic energy rate
W.Ak = 50000; % weight joint accelerations
W.ArmE = 1000000; % weight arm excitations
W.passMom = 1000; % weight passive torques
W.A = 2000; % weight muscle activations
mE = 0; % metabolic energy model identifier Bhargava et al. (2004)
        % 1 Umberger et al. (2003), 2 Umberger (2010), 3 Uchida et al. (2016), 
        % 4 Umberger (2010) treating muscle lengthening heat rate as 
        % Umberger et al. (2003), 5 Umberger (2010) treating negative mechanical 
        % work as Umberger et al. (2003)
exp_E = 2; % power metabolic energy

%Fixed parameter
W.u = 0.001;

%Set file save name
savename = 'test_predictiveSim';

%% Load external functions

%Original notes from Falisse code:

    % The external function performs inverse dynamics through the
    % OpenSim/Simbody C++ API. This external function is compiled as a dll from
    % which we create a Function instance using CasADi in MATLAB. More details
    % about the external function can be found in the documentation.

%%%%% TODO: note that these dll's have been built from CPP files that
%%%%% contain model based information (i.e. body masses, centre of masses,
%%%%% contact sphere locations etc.) --- it's likely that they wouldn't
%%%%% necessarily work the same with a scaled model without being rebuilt
    
%Loading external functions
%%%%% NOTE: currently only set to use generic contact sphere locations
setup.derivatives =  'AD'; % Algorithmic differentiation
if ispc 
    switch setup.derivatives
        case 'AD'     
            cd('..\externalFunctions');
            F = external('F','PredSim.dll'); 
            F1 = external('F','PredSim_pp.dll'); 
    end
elseif ismac
    switch setup.derivatives
        case 'AD'     
            cd('..\externalFunctions');
            F = external('F','PredSim.dylib');
            F1 = external('F','PredSim_pp.dylib'); 
    end
else
    disp('Platform not supported')
end   
cd(mainPath);

%% Indices for external function

%Indices of the elements in the external functions
%External function: F

%First, joint torques. 
jointi.pelvis.tilt  = 1; 
jointi.pelvis.list  = 2; 
jointi.pelvis.rot   = 3; 
jointi.pelvis.tx    = 4;
jointi.pelvis.ty    = 5;
jointi.pelvis.tz    = 6;
jointi.hip_flex.l   = 7;
jointi.hip_add.l    = 8;
jointi.hip_rot.l    = 9;
jointi.hip_flex.r   = 10;
jointi.hip_add.r    = 11;
jointi.hip_rot.r    = 12;
jointi.knee.l       = 13;
jointi.knee.r       = 14;
jointi.ankle.l      = 15;
jointi.ankle.r      = 16;
jointi.subt.l       = 17;
jointi.subt.r       = 18;
jointi.trunk.ext    = 19;
jointi.trunk.ben    = 20;
jointi.trunk.rot    = 21;
jointi.sh_flex.l    = 22;
jointi.sh_add.l     = 23;
jointi.sh_rot.l     = 24;
jointi.sh_flex.r    = 25;
jointi.sh_add.r     = 26;
jointi.sh_rot.r     = 27;
jointi.elb.l        = 28;
jointi.elb.r        = 29;

%Vectors of indices for later use
residualsi          = jointi.pelvis.tilt:jointi.elb.r; % all 
ground_pelvisi      = jointi.pelvis.tilt:jointi.pelvis.tz; % ground-pelvis
trunki              = jointi.trunk.ext:jointi.trunk.rot; % trunk
armsi               = jointi.sh_flex.l:jointi.elb.r; % arms
residuals_noarmsi   = jointi.pelvis.tilt:jointi.trunk.rot; % all but arms

%Number of degrees of freedom for later use
nq.all      = length(residualsi); % all 
nq.abs      = length(ground_pelvisi); % ground-pelvis
nq.trunk    = length(trunki); % trunk
nq.arms     = length(armsi); % arms
nq.leg      = 9; % #joints needed for polynomials

%Second, origins bodies. 
%Calcaneus
calcOr.r    = 30:31;
calcOr.l    = 32:33;
calcOr.all  = [calcOr.r,calcOr.l];
NcalcOr     = length(calcOr.all);

%Femurs
femurOr.r   = 34:35;
femurOr.l   = 36:37;
femurOr.all = [femurOr.r,femurOr.l];
NfemurOr    = length(femurOr.all);

%Hands
handOr.r    = 38:39;
handOr.l    = 40:41;
handOr.all  = [handOr.r,handOr.l];
NhandOr     = length(handOr.all);

%Tibias
tibiaOr.r   = 42:43;
tibiaOr.l   = 44:45;
tibiaOr.all = [tibiaOr.r,tibiaOr.l];
NtibiaOr    = length(tibiaOr.all);

%External function: F1 (post-processing purpose only)

%Ground reaction forces (GRFs)
GRFi.r      = 30:32;
GRFi.l      = 33:35;
GRFi.all    = [GRFi.r,GRFi.l];
NGRF        = length(GRFi.all);

%Origins calcaneus (3D)
calcOrall.r     = 36:38;
calcOrall.l     = 39:41;
calcOrall.all   = [calcOrall.r,calcOrall.l];
NcalcOrall      = length(calcOrall.all);

%% Set model info

%%%%% TODO: investigate capacity of altering this and external function
%%%%% parameters to suit a different model

body_mass = 62;
body_weight = body_mass*9.81;

%% Collocation scheme

%Original notes from Falisse code:

    % We use a pseudospectral direct collocation method, i.e. we use Lagrange
    % polynomials to approximate the state derivatives at the collocation
    % points in each mesh interval. We use d=3 collocation points per mesh
    % interval and Radau collocation points. 
    
%Set path to collocation scheme
pathCollocationScheme = [mainPath,'\..\collocationScheme'];
addpath(genpath(pathCollocationScheme));

%Settings
d = 3; % degree of interpolating polynomial
method = 'radau'; % collocation method

%Create collocation scheme
[tau_root,C,D,B] = CollocationScheme(d,method);

%% Muscle-tendon parameters 

%%%%% TODO: consider checking/updating these in line with Lai et al. to be
%%%%% more suitable for sprinting data

%Muscles from one leg and from the back
muscleNames = {'glut_med1_r','glut_med2_r','glut_med3_r',...
        'glut_min1_r','glut_min2_r','glut_min3_r','semimem_r',...
        'semiten_r','bifemlh_r','bifemsh_r','sar_r','add_long_r',...
        'add_brev_r','add_mag1_r','add_mag2_r','add_mag3_r','tfl_r',...
        'pect_r','grac_r','glut_max1_r','glut_max2_r','glut_max3_r',......
        'iliacus_r','psoas_r','quad_fem_r','gem_r','peri_r',...
        'rect_fem_r','vas_med_r','vas_int_r','vas_lat_r','med_gas_r',...
        'lat_gas_r','soleus_r','tib_post_r','flex_dig_r','flex_hal_r',...
        'tib_ant_r','per_brev_r','per_long_r','per_tert_r','ext_dig_r',...
        'ext_hal_r','ercspn_r','intobl_r','extobl_r','ercspn_l',...
        'intobl_l','extobl_l'};
    
%Set the path to muscle model
pathMuscleModel = [mainPath,'\..\muscleModel'];
addpath(genpath(pathMuscleModel));

%Muscle indices for later use
%(1:end-3), since we do not want to count twice the back muscles
musi = MuscleIndices(muscleNames(1:end-3));

%Total number of muscles
NMuscle = length(muscleNames(1:end-3))*2;

%Muscle-tendon parameters. Row 1: maximal isometric forces; Row 2: optimal
%fiber lengths; Row 3: tendon slack lengths; Row 4: optimal pennation 
%angles; Row 5: maximal contraction velocities

%%%%% TODO: update from the 'subject1' muscle parameters here

load([pathMuscleModel,'\MTparameters_subject1.mat']);
MTparameters_m = [MTparameters(:,musi),MTparameters(:,musi)];

%Adjust the maximal contraction velocity by a factor of 2
MTparameters_m(5,:) = MTparameters_m(5,:) * 2;

%Adjust the maximal isometric force by a factor of 2
MTparameters_m(1,:) = MTparameters_m(1,:) * 2;

%Indices of the muscles actuating the different joints for later use
pathPolynomial = [mainPath,'\..\polynomials'];
addpath(genpath(pathPolynomial));
tl = load([pathPolynomial,'\muscle_spanning_joint_INFO_subject1.mat']);
[~,mai] = MomentArmIndices(muscleNames(1:end-3),...
    tl.muscle_spanning_joint_INFO(1:end-3,:));

%We load some variables for the polynomial approximations
load([pathPolynomial,'/muscle_spanning_joint_INFO_subject1.mat']);
load([pathPolynomial,'/MuscleInfo_subject1.mat']);

%For the polynomials, we want all independent muscles. So we do not need 
%the muscles from both legs, since we assume bilateral symmetry, but want
%all muscles from the back (indices 47:49).
musi_pol = [musi,47,48,49];
NMuscle_pol = NMuscle/2+3;

%By default, the tendon stiffness is 35 and the shift is 0.
aTendon = 35*ones(NMuscle,1);
shift = zeros(NMuscle,1);

%% Metabolic energy model parameters

%Add the path to metabolic energy models and data
pathMetabolicEnergy = [mainPath,'\..\metabolicEnergy'];
addpath(genpath(pathMetabolicEnergy));

%Extract the specific tensions and slow twitch rations.
% (1:end-3), since we do not want to count twice the back muscles
tension = getSpecificTensions(muscleNames(1:end-3)); 
tensions = [tension;tension];
% (1:end-3), since we do not want to count twice the back muscles
pctst = getSlowTwitchRatios(muscleNames(1:end-3)); 
pctsts = [pctst;pctst];

%% CasADi functions

%Create CasADi functions for later use

%Add path to CasADi functions
pathCasADiFunctions = [mainPath,'\..\CasADiFunctions'];
addpath(genpath(pathCasADiFunctions));

%Load the CasADi functions for tracking simulations
CasADiFunctions_all

%% Passive joint torques

%%%%% TODO: consider altering passive moments to make more
%%%%% appropriate/useable with sprinting?

%Extract the parameters for the passive torques of the lower limbs and trunk
pathPassiveMoments = [mainPath,'\..\passiveMoments'];
addpath(genpath(pathPassiveMoments));
PassiveMomentsData

%% Experimental data

%We extract experimental data to set bounds and initial guesses if needed

%%%%% NOTE: currently set to use some data from an OK iteration of a
%%%%% jogging tracking simulation

%Path to experimental tracking simulation data
pathData = [mainPath,'\..\results\trackingSim'];

%Set joints
joints = {'pelvis_tilt','pelvis_list','pelvis_rotation','pelvis_tx',...
    'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
    'hip_rotation_l','hip_flexion_r','hip_adduction_r','hip_rotation_r',...
    'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
    'subtalar_angle_l','subtalar_angle_r','lumbar_extension',...
    'lumbar_bending','lumbar_rotation','arm_flex_l','arm_add_l',...
    'arm_rot_l','arm_flex_r','arm_add_r','arm_rot_r','elbow_flex_l',...
    'elbow_flex_r'};

%Set and add path of the helper functions
pathHelperFunctions = [pwd,'\..\helperFunctions'];
addpath(genpath(pathHelperFunctions));

%Extract joint positions from experimental jogging data
%%%%% TODO: update based on more appropriate file naming
%%%%% TODO: currently uses a 6Hz filter

%Set path to experimental data
pathExpData = [pwd,'\..\..\Dorn2012Data'];

%Set IK filename
% % % fileIK = 'sprint_ik.mot';
fileIK = 'jog_ik.mot';

%Extract joint kinematics
Qs_jog = getIK(pathExpData,fileIK,joints);

%%%%% Subtract 0.15 from the pelvis ty data as was done originally
Qs_jog.pelvis_ty(:,1) = Qs_jog.pelvis_ty - 0.15;
Qs_jog.all(:,6) = Qs_jog.all(:,6) - 0.15;
Qs_jog.allfilt(:,6) = Qs_jog.allfilt(:,6) - 0.15;

%We'll use a previous tracking sim as the initial guess
%%%%% TODO: fix this up with more appropriate naming/processes

%Set path to tracking sim results
pathResults = [pwd,'\..\results\trackingSim'];

%Set filename
%%%%% TODO: update file naming process
trackFileName = 'jog_trackingSim_10-Nov-2020_kinematics.mot';

%Extract joint positions from existing motion (previous results)
Qs_ig = getIK_ig(pathResults,trackFileName,joints);

%When tracking, we save a full gait cycle (2*N) so here we
%only select 1:N to have half a gait cycle
Qs_ig_sel.allfilt = Qs_ig.allfilt(1:N,:);
Qs_ig_sel.time = Qs_ig.time(1:N,:);
Qs_ig_sel.colheaders = Qs_ig.colheaders;  

%% Bounds

%Set path to bounds folder
pathBounds = [mainPath,'\..\bounds'];
addpath(genpath(pathBounds));

%%%%% TODO: bounds seem like they are oriented to the experimental data,
%%%%% but may need some checking to determine appropriateness to sprinting
%%%%% data

%Set bounds with tracking function
[bounds,scaling] = getBounds_all(Qs_jog,NMuscle,nq,jointi,v_tgt);

%% Initial guess

%%%%% TODO: only the motion is included in the initial guess here, and not
%%%%% the activations or any other design vairables --- these could be
%%%%% obtained from tracking results?????

%Add initial guess path
pathIG = [mainPath,'\..\initialGuess'];
addpath(genpath(pathIG));

%Set initial guess to data from the selected motion earlier
time_IC = [Qs_ig_sel.time(1),Qs_ig_sel.time(end)];
guess = getGuess_DI_t(Qs_ig_sel,nq,N,time_IC,NMuscle,jointi,scaling); 

%This allows visualizing the initial guess and the bounds
checkBoundsIG = input('Check bounds? [Y]es or [N]o? ','s');
if strcmpi(checkBoundsIG,'y')    
    %Add the plots path
    pathPlots = [mainPath,'\..\plots'];
    addpath(genpath(pathPlots));    
    %Check bounds of initial guess    
    plot_BoundsVSInitialGuess_all    
end

%% Formulate the NLP

% Start with an empty NLP
w   = {}; % design variables
w0  = []; % initial guess for design variables
lbw = []; % lower bounds for design variables
ubw = []; % upper bounds for design variables
J   = 0;  % initial value of cost function
g   = {}; % constraints
lbg = []; % lower bounds for constraints
ubg = []; % upper bounds for constraints

% Define static parameters
% Final time
tf              = MX.sym('tf',1);
w               = [w {tf}];
lbw             = [lbw; bounds.tf.lower];
ubw             = [ubw; bounds.tf.upper];
w0              = [w0;  guess.tf];

% Define states at first mesh point
% Muscle activations
a0              = MX.sym('a0',NMuscle);
w               = [w {a0}];
lbw             = [lbw; bounds.a.lower'];
ubw             = [ubw; bounds.a.upper'];
w0              = [w0;  guess.a(1,:)'];

% Muscle-tendon forces
FTtilde0        = MX.sym('FTtilde0',NMuscle);
w               = [w {FTtilde0}];
lbw             = [lbw; bounds.FTtilde.lower'];
ubw             = [ubw; bounds.FTtilde.upper'];
w0              = [w0;  guess.FTtilde(1,:)'];

% Qs and Qdots 
X0              = MX.sym('X0',2*nq.all);
w               = [w {X0}];    
lbw             = [lbw; bounds.QsQdots_0.lower'];
ubw             = [ubw; bounds.QsQdots_0.upper'];    
w0              = [w0;  guess.QsQdots(1,:)'];

% Arm activations
a_a0            = MX.sym('a_a0',nq.arms);
w               = [w {a_a0}];
lbw             = [lbw; bounds.a_a.lower'];
ubw             = [ubw; bounds.a_a.upper'];
w0              = [w0;  guess.a_a(1,:)'];    

% We pre-allocate some of the states so that we can provide an
% expression for the distance traveled
for k=0:N
    Xk{k+1,1} = MX.sym(['X_' num2str(k+1)], 2*nq.all);
end     

% "Lift" initial conditions
ak          = a0;
FTtildek    = FTtilde0;
Xk{1,1}     = X0;
a_ak        = a_a0;    

% Provide expression for the distance traveled
pelvis_tx0 = Xk{1,1}(2*jointi.pelvis.tx-1,1).*...
    scaling.QsQdots(2*jointi.pelvis.tx-1); % initial position pelvis_tx     
pelvis_txf = Xk{N+1,1}(2*jointi.pelvis.tx-1,1).*...
    scaling.QsQdots(2*jointi.pelvis.tx-1); % final position pelvis_tx    
dist_trav_tot = pelvis_txf-pelvis_tx0;% distance traveled   

% Time step
h = tf/N; 

% Loop over mesh points
for k=0:N-1
    % Define controls at mesh point (piecewise-constant in interval) 
    % Time derivative of muscle activations (states)
    vAk                 = MX.sym(['vA_' num2str(k)], NMuscle);
    w                   = [w {vAk}];
    lbw                 = [lbw; bounds.vA.lower'];
    ubw                 = [ubw; bounds.vA.upper'];
    w0                  = [w0; guess.vA(k+1,:)'];
    % Time derivative of muscle-tendon forces (states)
    dFTtildek           = MX.sym(['dFTtilde_' num2str(k)], NMuscle);
    w                   = [w {dFTtildek}];
    lbw                 = [lbw; bounds.dFTtilde.lower'];
    ubw                 = [ubw; bounds.dFTtilde.upper'];
    w0                  = [w0; guess.dFTtilde(k+1,:)'];  
    % Time derivative of Qdots (states) 
    Ak                  = MX.sym(['A_' num2str(k)], nq.all);
    w                   = [w {Ak}];
    lbw                 = [lbw; bounds.Qdotdots.lower'];
    ubw                 = [ubw; bounds.Qdotdots.upper'];
    w0                  = [w0; guess.Qdotdots(k+1,:)'];    
    % Arm excitations
    e_ak                = MX.sym(['e_a_' num2str(k)], nq.arms);
    w                   = [w {e_ak}];
    lbw                 = [lbw; bounds.e_a.lower'];
    ubw                 = [ubw; bounds.e_a.upper'];
    w0                  = [w0; guess.e_a(k+1,:)'];      
    % Define states at collocation points    
    % Muscle activations
    akj = {};
    for j=1:d
        akj{j}  = MX.sym(['	a_' num2str(k) '_' num2str(j)], NMuscle);
        w       = {w{:}, akj{j}};
        lbw     = [lbw; bounds.a.lower'];
        ubw     = [ubw; bounds.a.upper'];
        w0      = [w0;  guess.a(k+1,:)'];
    end   
    % Muscle-tendon forces
    FTtildekj = {};
    for j=1:d
        FTtildekj{j} = ...
            MX.sym(['FTtilde_' num2str(k) '_' num2str(j)], NMuscle);
        w            = {w{:}, FTtildekj{j}};
        lbw          = [lbw; bounds.FTtilde.lower'];
        ubw          = [ubw; bounds.FTtilde.upper'];
        w0           = [w0;  guess.FTtilde(k+1,:)'];
    end
    % Qs and Qdots        
    Xkj = {};
    for j=1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], 2*nq.all);
        w      = {w{:}, Xkj{j}};
        lbw    = [lbw; bounds.QsQdots.lower'];
        ubw    = [ubw; bounds.QsQdots.upper'];
        w0     = [w0;  guess.QsQdots(k+1,:)'];
    end   
    % Arm activations
    a_akj = {};
    for j=1:d
        a_akj{j}= MX.sym(['	a_a_' num2str(k) '_' num2str(j)], nq.arms);
        w       = {w{:}, a_akj{j}};
        lbw     = [lbw; bounds.a_a.lower'];
        ubw     = [ubw; bounds.a_a.upper'];
        w0      = [w0;  guess.a_a(k+1,:)'];
    end           
    % Unscale variables for later use
    Xk_nsc          = Xk{k+1,1}.*scaling.QsQdots';
    FTtildek_nsc    = FTtildek.*scaling.FTtilde';
    Ak_nsc          = Ak.*scaling.Qdotdots';
    for j=1:d
        Xkj_nsc{j} = Xkj{j}.*scaling.QsQdots';
        FTtildekj_nsc{j} = FTtildekj{j}.*scaling.FTtilde';
    end                
    % Get muscle-tendon lengths, velocities, and moment arms
    % Left leg
    qin_l = [Xk_nsc(jointi.hip_flex.l*2-1,1),...
        Xk_nsc(jointi.hip_add.l*2-1,1), ...
        Xk_nsc(jointi.hip_rot.l*2-1,1), ...
        Xk_nsc(jointi.knee.l*2-1,1), ...
        Xk_nsc(jointi.ankle.l*2-1,1),...
        Xk_nsc(jointi.subt.l*2-1,1),...
        Xk_nsc(jointi.trunk.ext*2-1,1),...
        Xk_nsc(jointi.trunk.ben*2-1,1),...
        Xk_nsc(jointi.trunk.rot*2-1,1)];  
    qdotin_l = [Xk_nsc(jointi.hip_flex.l*2,1),...
        Xk_nsc(jointi.hip_add.l*2,1),...
        Xk_nsc(jointi.hip_rot.l*2,1),...
        Xk_nsc(jointi.knee.l*2,1),...
        Xk_nsc(jointi.ankle.l*2,1),...
        Xk_nsc(jointi.subt.l*2,1),...
        Xk_nsc(jointi.trunk.ext*2,1),...
        Xk_nsc(jointi.trunk.ben*2,1),...
        Xk_nsc(jointi.trunk.rot*2,1)];  
    [lMTk_l,vMTk_l,MA_l] = f_lMT_vMT_dM(qin_l,qdotin_l);    
    MA.hip_flex.l   =  MA_l(mai(1).mus.l',1);
    MA.hip_add.l    =  MA_l(mai(2).mus.l',2);
    MA.hip_rot.l    =  MA_l(mai(3).mus.l',3);
    MA.knee.l       =  MA_l(mai(4).mus.l',4);
    MA.ankle.l      =  MA_l(mai(5).mus.l',5);  
    MA.subt.l       =  MA_l(mai(6).mus.l',6); 
    % For the back muscles, we want left and right together: left
    % first, right second. In MuscleInfo, we first have the right
    % muscles (44:46) and then the left muscles (47:49). Since the back
    % muscles only depend on back dofs, we do not care if we extract
    % them "from the left or right leg" so here we just picked left.
    MA.trunk_ext    =  MA_l([47:49,mai(7).mus.l]',7);
    MA.trunk_ben    =  MA_l([47:49,mai(8).mus.l]',8);
    MA.trunk_rot    =  MA_l([47:49,mai(9).mus.l]',9);
    % Right leg
    qin_r = [Xk_nsc(jointi.hip_flex.r*2-1,1),...
        Xk_nsc(jointi.hip_add.r*2-1,1),...
        Xk_nsc(jointi.hip_rot.r*2-1,1),...
        Xk_nsc(jointi.knee.r*2-1,1),...
        Xk_nsc(jointi.ankle.r*2-1,1),...
        Xk_nsc(jointi.subt.r*2-1,1),...
        Xk_nsc(jointi.trunk.ext*2-1,1),...
        Xk_nsc(jointi.trunk.ben*2-1,1),...
        Xk_nsc(jointi.trunk.rot*2-1,1)];  
    qdotin_r = [Xk_nsc(jointi.hip_flex.r*2,1),...
        Xk_nsc(jointi.hip_add.r*2,1),...
        Xk_nsc(jointi.hip_rot.r*2,1),...
        Xk_nsc(jointi.knee.r*2,1),...
        Xk_nsc(jointi.ankle.r*2,1),...
        Xk_nsc(jointi.subt.r*2,1),...
        Xk_nsc(jointi.trunk.ext*2,1),...
        Xk_nsc(jointi.trunk.ben*2,1),...
        Xk_nsc(jointi.trunk.rot*2,1)];      
    [lMTk_r,vMTk_r,MA_r] = f_lMT_vMT_dM(qin_r,qdotin_r);
    % Here we take the indices from left since the vector is 1:49
    MA.hip_flex.r   =  MA_r(mai(1).mus.l',1);
    MA.hip_add.r    =  MA_r(mai(2).mus.l',2);
    MA.hip_rot.r    =  MA_r(mai(3).mus.l',3);
    MA.knee.r       =  MA_r(mai(4).mus.l',4);
    MA.ankle.r      =  MA_r(mai(5).mus.l',5);
    MA.subt.r       =  MA_r(mai(6).mus.l',6);
    % Both legs
    % In MuscleInfo, we first have the right back muscles (44:46) and 
    % then the left back muscles (47:49). Here we re-organize so that
    % we have first the left muscles and then the right muscles.
    lMTk_lr = [lMTk_l([1:43,47:49],1);lMTk_r(1:46,1)];
    vMTk_lr = [vMTk_l([1:43,47:49],1);vMTk_r(1:46,1)];   
    % Get muscle-tendon forces and derive Hill-equilibrium
    [Hilldiffk,FTk,Fcek,Fpassk,Fisok,vMmaxk,massMk] = ...
        f_forceEquilibrium_FtildeState_all_tendon(ak,...
            FTtildek.*scaling.FTtilde',dFTtildek.*scaling.dFTtilde,...
            lMTk_lr,vMTk_lr,tensions,aTendon,shift); 
    % Get metabolic energy rate if in the cost function   
    if W.E ~= 0    
        % Get muscle fiber lengths
        [~,lMtildek] = f_FiberLength_TendonForce_tendon(...
            FTtildek.*scaling.FTtilde',lMTk_lr,aTendon,shift); 
        % Get muscle fiber velocities
        [vMk,~] = f_FiberVelocity_TendonForce_tendon(...
            FTtildek.*scaling.FTtilde',dFTtildek.*scaling.dFTtilde,...
            lMTk_lr,vMTk_lr,aTendon,shift);
        % Get metabolic energy rate
        if mE == 0 % Bhargava et al. (2004)
            [e_tot,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(ak,ak,...
                lMtildek,vMk,Fcek,Fpassk,massMk,pctsts,Fisok,...
                MTparameters_m(1,:)',body_mass,10);
        elseif mE == 1 % Umberger et al. (2003)
            % vMtilde defined for this model as vM/lMopt
            vMtildeUmbk = vMk./(MTparameters_m(2,:)');
            [e_tot,~,~,~,~] = fgetMetabolicEnergySmooth2003all(...
                ak,ak,lMtildek,vMtildeUmbk,vMk,Fcek,massMk,...
                pctsts,vMmaxk,Fisok,body_mass,10);
        elseif mE == 2 % Umberger (2010)
            % vMtilde defined for this model as vM/lMopt
            vMtildeUmbk = vMk./(MTparameters_m(2,:)');
            [e_tot,~,~,~,~] = fgetMetabolicEnergySmooth2010all(...
                ak,ak,lMtildek,vMtildeUmbk,vMk,Fcek,massMk,...
                pctsts,vMmaxk,Fisok,body_mass,10);  
        elseif mE == 3 % Uchida et al. (2016)
            % vMtilde defined for this model as vM/lMopt
            vMtildeUmbk = vMk./(MTparameters_m(2,:)');
            [e_tot,~,~,~,~] = fgetMetabolicEnergySmooth2016all(...
                ak,ak,lMtildek,vMtildeUmbk,vMk,Fcek,massMk,...
                pctsts,vMmaxk,Fisok,body_mass,10);  
        elseif mE == 4 % Umberger (2010) treating muscle lengthening 
            % heat rate as Umberger et al. (2003)
            % vMtilde defined for this model as vM/lMopt
            vMtildeUmbk = vMk./(MTparameters_m(2,:)');
            [e_tot,~,~,~,~] = fgetMetabolicEnergySmooth2010all_hl(...
                ak,ak,lMtildek,vMtildeUmbk,vMk,Fcek,massMk,...
                pctsts,vMmaxk,Fisok,body_mass,10); 
        elseif mE == 5 % Umberger (2010) treating negative mechanical 
            % work as Umberger et al. (2003)
            % vMtilde defined for this model as vM/lMopt
        vMtildeUmbk = vMk./(MTparameters_m(2,:)');
        [e_tot,~,~,~,~] = fgetMetabolicEnergySmooth2010all_neg(...
            ak,ak,lMtildek,vMtildeUmbk,vMk,Fcek,massMk,...
            pctsts,vMmaxk,Fisok,body_mass,10); 
        end                
    end        
    % Get passive joint torques
    Tau_passk.hip.flex.l    = f_PassiveMoments(k_pass.hip.flex,...
        theta.pass.hip.flex,Xk_nsc(jointi.hip_flex.l*2-1,1),...
        Xk_nsc(jointi.hip_flex.l*2,1));
    Tau_passk.hip.flex.r    = f_PassiveMoments(k_pass.hip.flex,...
        theta.pass.hip.flex,Xk_nsc(jointi.hip_flex.r*2-1,1),...
        Xk_nsc(jointi.hip_flex.r*2,1));
    Tau_passk.hip.add.l     = f_PassiveMoments(k_pass.hip.add,...
        theta.pass.hip.add,Xk_nsc(jointi.hip_add.l*2-1,1),...
        Xk_nsc(jointi.hip_add.l*2,1));
    Tau_passk.hip.add.r     = f_PassiveMoments(k_pass.hip.add,...
        theta.pass.hip.add,Xk_nsc(jointi.hip_add.r*2-1,1),...
        Xk_nsc(jointi.hip_add.r*2,1));
    Tau_passk.hip.rot.l     = f_PassiveMoments(k_pass.hip.rot,...
        theta.pass.hip.rot,Xk_nsc(jointi.hip_rot.l*2-1,1),...
        Xk_nsc(jointi.hip_rot.l*2,1));
    Tau_passk.hip.rot.r     = f_PassiveMoments(k_pass.hip.rot,...
        theta.pass.hip.rot,Xk_nsc(jointi.hip_rot.r*2-1,1),...
        Xk_nsc(jointi.hip_rot.r*2,1));
    Tau_passk.knee.l        = f_PassiveMoments(k_pass.knee,...
        theta.pass.knee,Xk_nsc(jointi.knee.l*2-1,1),...
        Xk_nsc(jointi.knee.l*2,1));
    Tau_passk.knee.r        = f_PassiveMoments(k_pass.knee,...
        theta.pass.knee,Xk_nsc(jointi.knee.r*2-1,1),...
        Xk_nsc(jointi.knee.r*2,1));
    Tau_passk.ankle.l       = f_PassiveMoments(k_pass.ankle,...
        theta.pass.ankle,Xk_nsc(jointi.ankle.l*2-1,1),...
        Xk_nsc(jointi.ankle.l*2,1));
    Tau_passk.ankle.r       = f_PassiveMoments(k_pass.ankle,...
        theta.pass.ankle,Xk_nsc(jointi.ankle.r*2-1,1),...
        Xk_nsc(jointi.ankle.r*2,1));        
    Tau_passk.subt.l       = f_PassiveMoments(k_pass.subt,...
        theta.pass.subt,Xk_nsc(jointi.subt.l*2-1,1),...
        Xk_nsc(jointi.subt.l*2,1));
    Tau_passk.subt.r       = f_PassiveMoments(k_pass.subt,...
        theta.pass.subt,Xk_nsc(jointi.subt.r*2-1,1),...
        Xk_nsc(jointi.subt.r*2,1));        
    Tau_passk.trunk.ext     = f_PassiveMoments(k_pass.trunk.ext,...
        theta.pass.trunk.ext,Xk_nsc(jointi.trunk.ext*2-1,1),...
        Xk_nsc(jointi.trunk.ext*2,1));
    Tau_passk.trunk.ben     = f_PassiveMoments(k_pass.trunk.ben,...
        theta.pass.trunk.ben,Xk_nsc(jointi.trunk.ben*2-1,1),...
        Xk_nsc(jointi.trunk.ben*2,1));
    Tau_passk.trunk.rot     = f_PassiveMoments(k_pass.trunk.rot,...
        theta.pass.trunk.rot,Xk_nsc(jointi.trunk.rot*2-1,1),...
        Xk_nsc(jointi.trunk.rot*2,1));        
    Tau_passk_all = [Tau_passk.hip.flex.l,Tau_passk.hip.flex.r,...
        Tau_passk.hip.add.l,Tau_passk.hip.add.r,...
        Tau_passk.hip.rot.l,Tau_passk.hip.rot.r,...
        Tau_passk.knee.l,Tau_passk.knee.r,Tau_passk.ankle.l,...
        Tau_passk.ankle.r,Tau_passk.subt.l,Tau_passk.subt.r,...
        Tau_passk.trunk.ext,Tau_passk.trunk.ben,...
        Tau_passk.trunk.rot]';            
    % Loop over collocation points
    Xk_nsc_end          = D(1)*Xk_nsc;
    FTtildek_nsc_end    = D(1)*FTtildek_nsc;
    ak_end              = D(1)*ak;
    a_ak_end            = D(1)*a_ak;
    for j=1:d
        % Expression for the state derivatives at the collocation point
        xp_nsc          = C(1,j+1)*Xk_nsc;
        FTtildep_nsc    = C(1,j+1)*FTtildek_nsc;
        ap              = C(1,j+1)*ak;
        a_ap            = C(1,j+1)*a_ak;
        for r=1:d
            xp_nsc       = xp_nsc + C(r+1,j+1)*Xkj_nsc{r};
            FTtildep_nsc = FTtildep_nsc + C(r+1,j+1)*FTtildekj_nsc{r};
            ap           = ap + C(r+1,j+1)*akj{r};
            a_ap         = a_ap + C(r+1,j+1)*a_akj{r};
        end 
        % Append collocation equations
        % Dynamic constraints are scaled using the same scale
        % factors as was used to scale the states
        % Activation dynamics (implicit formulation)  
        g       = {g{:}, (h*vAk.*scaling.vA - ap)./scaling.a};
        lbg     = [lbg; zeros(NMuscle,1)];
        ubg     = [ubg; zeros(NMuscle,1)]; 
        % Contraction dynamics (implicit formulation)            
        g       = {g{:}, (h*dFTtildek.*scaling.dFTtilde - ...
            FTtildep_nsc)./(scaling.FTtilde')};
        lbg     = [lbg; zeros(NMuscle,1)];
        ubg     = [ubg; zeros(NMuscle,1)];
        % Skeleton dynamics (implicit formulation)       
        xj_nsc  = [...
            Xkj_nsc{j}(2); Ak_nsc(1); Xkj_nsc{j}(4); Ak_nsc(2);...
            Xkj_nsc{j}(6); Ak_nsc(3); Xkj_nsc{j}(8); Ak_nsc(4);...
            Xkj_nsc{j}(10); Ak_nsc(5); Xkj_nsc{j}(12); Ak_nsc(6);...
            Xkj_nsc{j}(14); Ak_nsc(7); Xkj_nsc{j}(16); Ak_nsc(8);...
            Xkj_nsc{j}(18); Ak_nsc(9); Xkj_nsc{j}(20); Ak_nsc(10);...
            Xkj_nsc{j}(22); Ak_nsc(11); Xkj_nsc{j}(24); Ak_nsc(12);...
            Xkj_nsc{j}(26); Ak_nsc(13); Xkj_nsc{j}(28); Ak_nsc(14);...
            Xkj_nsc{j}(30); Ak_nsc(15); Xkj_nsc{j}(32); Ak_nsc(16);...
            Xkj_nsc{j}(34); Ak_nsc(17); Xkj_nsc{j}(36); Ak_nsc(18);...
            Xkj_nsc{j}(38); Ak_nsc(19); Xkj_nsc{j}(40); Ak_nsc(20);...
            Xkj_nsc{j}(42); Ak_nsc(21); Xkj_nsc{j}(44); Ak_nsc(22);...
            Xkj_nsc{j}(46); Ak_nsc(23); Xkj_nsc{j}(48); Ak_nsc(24); ...
            Xkj_nsc{j}(50); Ak_nsc(25); Xkj_nsc{j}(52); Ak_nsc(26); ...
            Xkj_nsc{j}(54); Ak_nsc(27); Xkj_nsc{j}(56); Ak_nsc(28); ...
            Xkj_nsc{j}(58); Ak_nsc(29);];
        g       = {g{:}, (h*xj_nsc - xp_nsc)./(scaling.QsQdots')};
        lbg     = [lbg; zeros(2*nq.all,1)];
        ubg     = [ubg; zeros(2*nq.all,1)];   
        % Arm activation dynamics (explicit formulation)  
        dadt    = f_ArmActivationDynamics(e_ak,a_akj{j});
        g       = {g{:}, (h*dadt - a_ap)./scaling.a_a};
        lbg     = [lbg; zeros(nq.arms,1)];
        ubg     = [ubg; zeros(nq.arms,1)]; 
        % Add contribution to the end state
        Xk_nsc_end       = Xk_nsc_end + D(j+1)*Xkj_nsc{j};
        FTtildek_nsc_end = FTtildek_nsc_end + D(j+1)*FTtildekj_nsc{j};
        ak_end           = ak_end + D(j+1)*akj{j};  
        a_ak_end         = a_ak_end + D(j+1)*a_akj{j};    
        % Add contribution to the quadrature function
        if W.E == 0
        J = J + 1/(dist_trav_tot)*(...
            W.A*B(j+1)      *(f_J92(akj{j}))*h + ...
            W.ArmE*B(j+1)   *(f_J8(e_ak))*h +... 
            W.Ak*B(j+1)     *(f_J21(Ak(residuals_noarmsi,1)))*h + ... 
            W.passMom*B(j+1)*(f_J15(Tau_passk_all))*h + ...
            W.u*B(j+1)      *(f_J92(vAk))*h + ...
            W.u*B(j+1)      *(f_J92(dFTtildek))*h + ...                
            W.u*B(j+1)      *(f_J8(Ak(armsi,1)))*h);
        elseif W.A == 0
        J = J + 1/(dist_trav_tot)*(...
            W.E*B(j+1)      *(f_J92exp(e_tot,exp_E))/body_mass*h + ...                
            W.ArmE*B(j+1)   *(f_J8(e_ak))*h +... 
            W.Ak*B(j+1)     *(f_J21(Ak(residuals_noarmsi,1)))*h + ... 
            W.passMom*B(j+1)*(f_J15(Tau_passk_all))*h + ...
            W.u*B(j+1)      *(f_J92(vAk))*h + ...
            W.u*B(j+1)      *(f_J92(dFTtildek))*h + ...                
            W.u*B(j+1)      *(f_J8(Ak(armsi,1)))*h); 
        elseif W.passMom == 0
            J = J + 1/(dist_trav_tot)*(...
            W.E*B(j+1)      *(f_J92exp(e_tot,exp_E))/body_mass*h + ...                
            W.A*B(j+1)      *(f_J92(akj{j}))*h + ...
            W.ArmE*B(j+1)   *(f_J8(e_ak))*h +... 
            W.Ak*B(j+1)     *(f_J21(Ak(residuals_noarmsi,1)))*h + ... 
            W.u*B(j+1)      *(f_J92(vAk))*h + ...
            W.u*B(j+1)      *(f_J92(dFTtildek))*h + ...                
            W.u*B(j+1)      *(f_J8(Ak(armsi,1)))*h);
        else                
        J = J + 1/(dist_trav_tot)*(...
            W.E*B(j+1)      *(f_J92exp(e_tot,exp_E))/body_mass*h + ...                
            W.A*B(j+1)      *(f_J92(akj{j}))*h + ...
            W.ArmE*B(j+1)   *(f_J8(e_ak))*h +... 
            W.Ak*B(j+1)     *(f_J21(Ak(residuals_noarmsi,1)))*h + ... 
            W.passMom*B(j+1)*(f_J15(Tau_passk_all))*h + ...
            W.u*B(j+1)      *(f_J92(vAk))*h + ...
            W.u*B(j+1)      *(f_J92(dFTtildek))*h + ...                
            W.u*B(j+1)      *(f_J8(Ak(armsi,1)))*h); 
        end
    end            
    % Call external function
    [Tk] = F([Xk_nsc;Ak_nsc]);           
    % Add path constraints
    % Null pelvis residuals
    g               = {g{:},Tk(ground_pelvisi,1)};
    lbg             = [lbg; zeros(nq.abs,1)];
    ubg             = [ubg; zeros(nq.abs,1)];  
    % Muscle-driven joint torques for the lower limbs and the trunk
    % Hip flexion, left
    Ft_hip_flex_l   = FTk(mai(1).mus.l',1);
    T_hip_flex_l    = f_T27(MA.hip_flex.l,Ft_hip_flex_l);
    g               = {g{:},Tk(jointi.hip_flex.l,1)-(T_hip_flex_l + ...
        Tau_passk.hip.flex.l - 0.1*Xk_nsc(jointi.hip_flex.l*2,1))};        
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];    
    % Hip flexion, right
    Ft_hip_flex_r   = FTk(mai(1).mus.r',1);
    T_hip_flex_r    = f_T27(MA.hip_flex.r,Ft_hip_flex_r);
    g               = {g{:},Tk(jointi.hip_flex.r,1)-(T_hip_flex_r + ...
        Tau_passk.hip.flex.r - 0.1*Xk_nsc(jointi.hip_flex.r*2,1))};        
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];   
    % Hip adduction, left
    Ft_hip_add_l    = FTk(mai(2).mus.l',1);
    T_hip_add_l     = f_T27(MA.hip_add.l,Ft_hip_add_l);
    g               = {g{:},Tk(jointi.hip_add.l,1)-(T_hip_add_l + ...
        Tau_passk.hip.add.l - 0.1*Xk_nsc(jointi.hip_add.l*2,1))};
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];    
    % Hip adduction, right
    Ft_hip_add_r    = FTk(mai(2).mus.r',1);
    T_hip_add_r     = f_T27(MA.hip_add.r,Ft_hip_add_r);
    g               = {g{:},Tk(jointi.hip_add.r,1)-(T_hip_add_r + ...
        Tau_passk.hip.add.r - 0.1*Xk_nsc(jointi.hip_add.r*2,1))};
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];  
    % Hip rotation, left
    Ft_hip_rot_l    = FTk(mai(3).mus.l',1);
    T_hip_rot_l     = f_T27(MA.hip_rot.l,Ft_hip_rot_l);
    g               = {g{:},Tk(jointi.hip_rot.l,1)-(T_hip_rot_l + ...
        Tau_passk.hip.rot.l - 0.1*Xk_nsc(jointi.hip_rot.l*2,1))};
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];    
    % Hip rotation, right
    Ft_hip_rot_r    = FTk(mai(3).mus.r',1);
    T_hip_rot_r     = f_T27(MA.hip_rot.r,Ft_hip_rot_r);
    g               = {g{:},Tk(jointi.hip_rot.r,1)-(T_hip_rot_r + ...
        Tau_passk.hip.rot.r - 0.1*Xk_nsc(jointi.hip_rot.r*2,1))};
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];   
    % Knee, left
    Ft_knee_l       = FTk(mai(4).mus.l',1);
    T_knee_l        = f_T13(MA.knee.l,Ft_knee_l);
    g               = {g{:},Tk(jointi.knee.l,1)-(T_knee_l + ...
        Tau_passk.knee.l - 0.1*Xk_nsc(jointi.knee.l*2,1))};
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];    
    % Knee, right
    Ft_knee_r       = FTk(mai(4).mus.r',1);
    T_knee_r        = f_T13(MA.knee.r,Ft_knee_r);
    g               = {g{:},Tk(jointi.knee.r,1)-(T_knee_r + ...
        Tau_passk.knee.r - 0.1*Xk_nsc(jointi.knee.r*2,1))};
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];    
    % Ankle, left
    Ft_ankle_l      = FTk(mai(5).mus.l',1);
    T_ankle_l       = f_T12(MA.ankle.l,Ft_ankle_l);
    g               = {g{:},Tk(jointi.ankle.l,1)-(T_ankle_l + ...
        Tau_passk.ankle.l - 0.1*Xk_nsc(jointi.ankle.l*2,1))};
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];    
    % Ankle, right
    Ft_ankle_r      = FTk(mai(5).mus.r',1);
    T_ankle_r       = f_T12(MA.ankle.r,Ft_ankle_r);
    g               = {g{:},Tk(jointi.ankle.r,1)-(T_ankle_r + ...
        Tau_passk.ankle.r - 0.1*Xk_nsc(jointi.ankle.r*2,1))};
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];
    % Subtalar, left
    Ft_subt_l       = FTk(mai(6).mus.l',1);
    T_subt_l        = f_T12(MA.subt.l,Ft_subt_l);
    g               = {g{:},(Tk(jointi.subt.l,1)-(T_subt_l + ...
        Tau_passk.subt.l - 0.1*Xk_nsc(jointi.subt.l*2,1)))};
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];    
    % Subtalar, right
    Ft_subt_r       = FTk(mai(6).mus.r',1);
    T_subt_r        = f_T12(MA.subt.r,Ft_subt_r);
    g               = {g{:},(Tk(jointi.subt.r,1)-(T_subt_r + ...
        Tau_passk.subt.r - 0.1*Xk_nsc(jointi.subt.r*2,1)))};
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];
    % Lumbar extension
    Ft_trunk_ext    = FTk([mai(7).mus.l,mai(7).mus.r]',1);
    T_trunk_ext     = f_T6(MA.trunk_ext,Ft_trunk_ext);
    g               = {g{:},Tk(jointi.trunk.ext,1)-(T_trunk_ext + ...
        Tau_passk.trunk.ext - 0.1*Xk_nsc(jointi.trunk.ext*2,1))};
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];
    % Lumbar bending
    Ft_trunk_ben    = FTk([mai(8).mus.l,mai(8).mus.r]',1);
    T_trunk_ben     = f_T6(MA.trunk_ben,Ft_trunk_ben);
    g               = {g{:},Tk(jointi.trunk.ben,1)-(T_trunk_ben + ...
        Tau_passk.trunk.ben - 0.1*Xk_nsc(jointi.trunk.ben*2,1))};
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];
    % Lumbar rotation
    Ft_trunk_rot    = FTk([mai(9).mus.l,mai(9).mus.r]',1);
    T_trunk_rot     = f_T6(MA.trunk_rot,Ft_trunk_rot);
    g               = {g{:},Tk(jointi.trunk.rot,1)-(T_trunk_rot + ...
        Tau_passk.trunk.rot - 0.1*Xk_nsc(jointi.trunk.rot*2,1))};
    lbg             = [lbg; 0];
    ubg             = [ubg; 0];
    % Torque-driven joint torques for the arms
    % Arms
    g               = {g{:},Tk(armsi,1)/scaling.ArmTau - a_ak + ...
        0.1/scaling.ArmTau*Xk_nsc(armsi*2,1)};
    lbg             = [lbg; zeros(nq.arms,1)];
    ubg             = [ubg; zeros(nq.arms,1)];
    % Activation dynamics (implicit formulation)
    tact = 0.015;
    tdeact = 0.06;
    act1 = vAk*scaling.vA + ak./(ones(size(ak,1),1)*tdeact);
    act2 = vAk*scaling.vA + ak./(ones(size(ak,1),1)*tact);
    % act1
    g               = {g{:},act1};
    lbg             = [lbg; zeros(NMuscle,1)];
    ubg             = [ubg; inf*ones(NMuscle,1)]; 
    % act2
    g               = {g{:},act2};
    lbg             = [lbg; -inf*ones(NMuscle,1)];
    ubg             = [ubg; ones(NMuscle,1)./(ones(NMuscle,1)*tact)];        
    % Contraction dynamics (implicit formulation)
    g               = {g{:},Hilldiffk};
    lbg             = [lbg; zeros(NMuscle,1)];
    ubg             = [ubg; zeros(NMuscle,1)];  
    % Constraints to prevent parts of the skeleton to penetrate each
    % other.
    % Origins calcaneus (transv plane) at minimum 9 cm from each other.
    g               = {g{:},f_Jnn2(Tk(calcOr.r,1) - Tk(calcOr.l,1))};
    lbg             = [lbg; 0.0081];
    ubg             = [ubg; 4];   
    % Constraint to prevent the arms to penetrate the skeleton       
    % Origins femurs and ipsilateral hands (transv plane) at minimum 
    % 18 cm from each other.
    g               = {g{:},f_Jnn2(Tk(femurOr.r,1) - Tk(handOr.r,1))};
    lbg             = [lbg; 0.0324];
    ubg             = [ubg; 4];     
    g               = {g{:},f_Jnn2(Tk(femurOr.l,1) - Tk(handOr.l,1))};
    lbg             = [lbg; 0.0324];
    ubg             = [ubg; 4]; 
    % Origins tibia (transv plane) at minimum 11 cm from each other.   
    g               = {g{:},f_Jnn2(Tk(tibiaOr.r,1) - Tk(tibiaOr.l,1))};
    lbg             = [lbg; 0.0121];
    ubg             = [ubg; 4]; 
    % We impose periodicity and symmetry and create some helper
    % variables to easily impose these constraints. Here we use these
    % variables when settings initial guesses of states at the end of
    % the current interval.
    % Muscles
    orderMusInv = [NMuscle/2+1:NMuscle,1:NMuscle/2]; 
    % Arms
    orderArmInv = [jointi.sh_flex.r:jointi.sh_rot.r,...
        jointi.sh_flex.l:jointi.sh_rot.l,...
        jointi.elb.r:jointi.elb.r,...
        jointi.elb.l:jointi.elb.l]-jointi.sh_flex.l+1; 
    % Lower limbs and trunk
    % Qs and Qdots are inverted after a half gait cycle BUT 
    % Pelvis: pelvis tilt and pelvis ty should be equal, pelvis
    % list, rot, tz should be opposite and pelvis tx should be equal
    % plus dist traveled.
    % Trunk: lumbar ext should be equal, lumbar bend and lumbar rot
    % should be of opposite.     
    orderQsInv = [jointi.pelvis.tilt:2*jointi.pelvis.tz,...
        2*jointi.hip_flex.r-1:2*jointi.hip_rot.r,...
        2*jointi.hip_flex.l-1:2*jointi.hip_rot.l,...
        2*jointi.knee.r-1:2*jointi.knee.r,...
        2*jointi.knee.l-1:2*jointi.knee.l,...
        2*jointi.ankle.r-1:2*jointi.ankle.r,...
        2*jointi.ankle.l-1:2*jointi.ankle.l,...
        2*jointi.subt.r-1:2*jointi.subt.r,...
        2*jointi.subt.l-1:2*jointi.subt.l,...
        2*jointi.trunk.ext-1:2*jointi.trunk.rot,...
        2*jointi.sh_flex.r-1:2*jointi.sh_rot.r,...
        2*jointi.sh_flex.l-1:2*jointi.sh_rot.l,...
        2*jointi.elb.r-1:2*jointi.elb.r,...
        2*jointi.elb.l-1:2*jointi.elb.l];        
    orderQsOpp = [2*jointi.pelvis.list-1:2*jointi.pelvis.list,...   
        2*jointi.pelvis.rot-1:2*jointi.pelvis.rot,...
        2*jointi.pelvis.tz-1:2*jointi.pelvis.tz,...
        2*jointi.trunk.ben-1:2*jointi.trunk.ben,...
        2*jointi.trunk.rot-1:2*jointi.trunk.rot];
    % New NLP variables for states at end of interval
    if k ~= N-1
        % Muscle activations
        ak              = MX.sym(['a_' num2str(k+1)], NMuscle);
        w               = {w{:}, ak};
        lbw             = [lbw; bounds.a.lower'];
        ubw             = [ubw; bounds.a.upper'];
        w0              = [w0;  guess.a(k+2,:)'];
        % Muscle-tendon forces
        FTtildek        = MX.sym(['FTtilde_' num2str(k+1)], NMuscle);
        w               = {w{:}, FTtildek};
        lbw             = [lbw; bounds.FTtilde.lower'];
        ubw             = [ubw; bounds.FTtilde.upper'];
        w0              = [w0;  guess.FTtilde(k+2,:)'];    
        % Qs and Qdots
        w               = {w{:}, Xk{k+2,1}};
        lbw             = [lbw; bounds.QsQdots.lower'];
        ubw             = [ubw; bounds.QsQdots.upper']; 
        w0              = [w0;  guess.QsQdots(k+2,:)'];
        % Arm activations
        a_ak            = MX.sym(['a_a_' num2str(k+1)], nq.arms);
        w               = {w{:}, a_ak};
        lbw             = [lbw; bounds.a_a.lower'];
        ubw             = [ubw; bounds.a_a.upper'];
        w0              = [w0;  guess.a_a(k+2,:)'];
    else % Periodicty            
        % Muscle activations
        ak              = MX.sym(['a_' num2str(k+1)], NMuscle);
        w               = {w{:}, ak};
        lbw             = [lbw; bounds.a.lower'];
        ubw             = [ubw; bounds.a.upper'];            
        w0              = [w0;  guess.a(1,orderMusInv)'];
        % Muscle-tendon forces
        FTtildek        = MX.sym(['FTtilde_' num2str(k+1)], NMuscle);
        w               = {w{:}, FTtildek};
        lbw             = [lbw; bounds.FTtilde.lower'];
        ubw             = [ubw; bounds.FTtilde.upper'];
        w0              = [w0;  guess.FTtilde(1,orderMusInv)'];    
        % Qs and Qdots
        w               = {w{:}, Xk{k+2,1}};
        lbw             = [lbw; bounds.QsQdots.lower'];
        ubw             = [ubw; bounds.QsQdots.upper'];
        % For "bilateral" joints, we invert right and left
        inv_X           = guess.QsQdots(1,orderQsInv);
        % For other joints, we take the opposite right and left
        inv_X(orderQsOpp) = -guess.QsQdots(1,orderQsOpp);           
        dx = guess.QsQdots(end,2*jointi.pelvis.tx-1) - ...
            guess.QsQdots(end-1,2*jointi.pelvis.tx-1);
        inv_X(2*jointi.pelvis.tx-1) = ...
            guess.QsQdots(end,2*jointi.pelvis.tx-1) + dx;            
        w0                = [w0;  inv_X'];   
        % Arm activations
        a_ak            = MX.sym(['a_a_' num2str(k+1)], nq.arms);
        w               = {w{:}, a_ak};
        lbw             = [lbw; bounds.a_a.lower'];
        ubw             = [ubw; bounds.a_a.upper'];            
        w0              = [w0;  guess.a_a(1,orderArmInv)'];
    end
    % Rescale variables to impose equality constraints
    Xk_end = (Xk_nsc_end)./scaling.QsQdots';
    FTtildek_end = (FTtildek_nsc_end)./scaling.FTtilde';
    % Add equality constraints (next interval starts with end values of 
    % states from previous interval)
    g   = {g{:}, Xk_end-Xk{k+2,1}, FTtildek_end-FTtildek, ...
        ak_end-ak, a_ak_end-a_ak};
    lbg = [lbg; zeros(2*nq.all + NMuscle + NMuscle + nq.arms,1)];
    ubg = [ubg; zeros(2*nq.all + NMuscle + NMuscle + nq.arms,1)];  
end

% Additional path constraints
% Periodicity of the states

% Qs and Qdots
QsInvA = [jointi.pelvis.tilt:2*jointi.pelvis.tilt,...
    2*jointi.pelvis.tx,2*jointi.pelvis.ty-1:2*jointi.pelvis.ty,...
    2*jointi.hip_flex.l-1:2*jointi.trunk.ext,...
    2*jointi.sh_flex.l-1:2*jointi.elb.r]';
QsInvB = [jointi.pelvis.tilt:2*jointi.pelvis.tilt,...
    2*jointi.pelvis.tx,2*jointi.pelvis.ty-1:2*jointi.pelvis.ty,...
    2*jointi.hip_flex.r-1:2*jointi.hip_rot.r,...
    2*jointi.hip_flex.l-1:2*jointi.hip_rot.l,...
    2*jointi.knee.r-1:2*jointi.knee.r,...
    2*jointi.knee.l-1:2*jointi.knee.l,...
    2*jointi.ankle.r-1:2*jointi.ankle.r,...
    2*jointi.ankle.l-1:2*jointi.ankle.l,...
    2*jointi.subt.r-1:2*jointi.subt.r,...
    2*jointi.subt.l-1:2*jointi.subt.l,...
    2*jointi.trunk.ext-1,2*jointi.trunk.ext,...
    2*jointi.sh_flex.r-1:2*jointi.sh_rot.r,...
    2*jointi.sh_flex.l-1:2*jointi.sh_rot.l,...
    2*jointi.elb.r-1:2*jointi.elb.r,...
    2*jointi.elb.l-1:2*jointi.elb.l]';             
g   = {g{:}, Xk_end(QsInvA)-X0(QsInvB,1)};
lbg = [lbg; zeros(length(QsInvB),1)];
ubg = [ubg; zeros(length(QsInvB),1)];     
g   = {g{:}, Xk_end(orderQsOpp)+X0(orderQsOpp,1)};
lbg = [lbg; zeros(length(orderQsOpp),1)];
ubg = [ubg; zeros(length(orderQsOpp),1)]; 

% Muscle activations
g   = {g{:}, ak_end-a0(orderMusInv,1)};
lbg = [lbg; zeros(NMuscle,1)];
ubg = [ubg; zeros(NMuscle,1)];

% Muscle-tendon forces
g   = {g{:}, FTtildek_end-FTtilde0(orderMusInv,1)};
lbg = [lbg; zeros(NMuscle,1)];
ubg = [ubg; zeros(NMuscle,1)];

% Arm activations
g   = {g{:}, a_ak_end-a_a0(orderArmInv,1)};
lbg = [lbg; zeros(nq.arms,1)];
ubg = [ubg; zeros(nq.arms,1)];    

% Average speed
vel_aver_tot = dist_trav_tot/tf; 
g   = {g{:}, vel_aver_tot - v_tgt};
lbg = [lbg; 0];
ubg = [ubg; 0];

% Assert bounds / IG

% Lower bounds smaller than upper bounds    
assert_bw = isempty(find(lbw <= ubw == 0,1)); % Design variables
assert_bg = isempty(find(lbg <= ubg == 0,1)); % Constraints

% Design varibles between -1 and 1
assert_bwl = isempty(find(lbw < -1 == 1,1));
assert_bwu = isempty(find(1 < ubw == 1,1));   

% Initial guess within bounds
%%%%% ISSUE HERE WITH INITIAL GUESS NOT BEING WITHIN LOWER BOUNDS AT 1 POINT
%%%%% MANUALLY CHANGED IT SO IT WAS --- UNSURE WHAT IT IS THOUGH...
assert_w0_ubw = isempty(find(w0 <= ubw == 0,1));
assert_w0_lbw = isempty(find(lbw <= w0 == 0,1));

%Create NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.mu_strategy      = 'adaptive';
options.ipopt.max_iter = 10000;
options.ipopt.tol = 1*10^(-tol_ipopt);
solver = nlpsol('solver', 'ipopt', prob, options);

%Create and save diary
saveName = 'predictiveSim';
saveFileName = ['predictiveSim_',date];
pathResults = [mainPath,'\..\results'];
if ~(exist([pathResults,'\',saveName],'dir')==7)
    mkdir(pathResults,saveName);
end 
diary([pathResults,'\',saveName,'\',saveFileName]); 

%Solve problem
clc
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
    'lbg', lbg, 'ubg', ubg);

%Turn diary off
diary off

%Extract results
w_opt = full(sol.x);
g_opt = full(sol.g); 

%Create setup
setup.tolerance.ipopt = tol_ipopt;
setup.bounds = bounds;
setup.scaling = scaling;
setup.guess = guess;
setup.lbw = lbw;
setup.ubw = ubw;

%Save results and setup
save([pathResults,'\',saveName,'\',saveFileName,'_w.mat'],'w_opt');
save([pathResults,'\',saveName,'\',saveFileName,'_g.mat'],'g_opt');
save([pathResults,'\',saveName,'\',saveFileName,'_setup.mat'],'setup');

%% Extract results

%All optimized design variables are saved in a single column vector      

%Number of design variables    
NControls = NMuscle+NMuscle+nq.all+nq.arms;
NStates = NMuscle+NMuscle+2*nq.all+nq.arms;
NParameters = 1;
%In the loop
Nwl = NControls+d*NStates+NStates;
%In total
Nw = NParameters+NStates+N*Nwl;
%Before the variable corresponding to the first collocation point
Nwm = NParameters+NStates+NControls;

%Here we extract the results and re-organize them for analysis  

%Static parameters
tf_opt  = w_opt(1:NParameters);

%Mesh points

%Muscle activations and muscle-tendon forces
a_opt = zeros(N+1,NMuscle);
FTtilde_opt = zeros(N+1,NMuscle);
for i = 1:NMuscle
    a_opt(:,i) = w_opt(NParameters+i:Nwl:Nw);
    FTtilde_opt(:,i) = w_opt(NParameters+NMuscle+i:Nwl:Nw);
end

%Qs and Qdots
q_opt = zeros(N+1,nq.all);
qdot_opt = zeros(N+1,nq.all);
count = 0;
for i = 1:2:2*nq.all
    count = count +1;
    q_opt(:,count) = w_opt(NParameters+NMuscle+NMuscle+i:Nwl:Nw);
    qdot_opt(:,count) = w_opt(NParameters+NMuscle+NMuscle+i+1:Nwl:Nw);
end

%Arm activations
a_a_opt = zeros(N+1,nq.arms);
for i = 1:nq.arms
    a_a_opt(:,i) =w_opt(NParameters+NMuscle+NMuscle+2*nq.all+i:Nwl:Nw);
end   

%Time derivative of muscle activations and muscle-tendon forces
vA_opt = zeros(N,NMuscle);
dFTtilde_opt = zeros(N,NMuscle);
for i = 1:NMuscle
    vA_opt(:,i) = w_opt(NParameters+NStates+i:Nwl:Nw);
    dFTtilde_opt(:,i) = w_opt(NParameters+NStates+NMuscle+i:Nwl:Nw);
end

%Time derivative of joint velocities
qdotdot_opt = zeros(N,nq.all);
for i = 1:nq.all
    qdotdot_opt(:,i) = ...
        w_opt(NParameters+NStates+NMuscle+NMuscle+i:Nwl:Nw);
end

%Arm excitations
e_a_opt = zeros(N,nq.arms);
for i = 1:nq.arms
    e_a_opt(:,i) = ...
       w_opt(NParameters+NStates+NMuscle+NMuscle+nq.all+i:Nwl:Nw); 
end 

%Collocation points

%Muscle activations
a_opt_ext=zeros(N*(d+1)+1,NMuscle);
a_opt_ext(1:(d+1):end,:)= a_opt;
for nmusi=1:NMuscle
    a_opt_ext(2:(d+1):end,nmusi) = w_opt(Nwm+nmusi:Nwl:Nw);
    a_opt_ext(3:(d+1):end,nmusi) = w_opt(Nwm+NMuscle+nmusi:Nwl:Nw);
    a_opt_ext(4:(d+1):end,nmusi) = ...
        w_opt(Nwm+NMuscle+NMuscle+nmusi:Nwl:Nw);
end

%Muscle activations at collocation points only
a_opt_ext_col = zeros(N*d,NMuscle); 
for nmusi=1:NMuscle
    a_opt_ext_col(1:d:end,nmusi) = w_opt(Nwm+nmusi:Nwl:Nw);
    a_opt_ext_col(2:d:end,nmusi) = w_opt(Nwm+NMuscle+nmusi:Nwl:Nw);
    a_opt_ext_col(3:d:end,nmusi) = ...
        w_opt(Nwm+NMuscle+NMuscle+nmusi:Nwl:Nw);   
end  

%Muscle-tendon forces
FTtilde_opt_ext=zeros(N*(d+1)+1,NMuscle);
FTtilde_opt_ext(1:(d+1):end,:)= FTtilde_opt;
for nmusi=1:NMuscle
    FTtilde_opt_ext(2:(d+1):end,nmusi) = ...
        w_opt(Nwm+d*NMuscle+nmusi:Nwl:Nw);
    FTtilde_opt_ext(3:(d+1):end,nmusi) = ...
        w_opt(Nwm+d*NMuscle+NMuscle+nmusi:Nwl:Nw);
    FTtilde_opt_ext(4:(d+1):end,nmusi) = ...
        w_opt(Nwm+d*NMuscle+NMuscle+NMuscle+nmusi:Nwl:Nw);
end

%Qs and Qdots
q_opt_ext=zeros(N*(d+1)+1,nq.all);
q_opt_ext(1:(d+1):end,:)= q_opt;
q_dot_opt_ext=zeros(N*(d+1)+1,nq.all);
q_dot_opt_ext(1:(d+1):end,:)= qdot_opt;
nqi_col = 1:2:2*nq.all;
for nqi=1:nq.all
    nqi_q = nqi_col(nqi);
    q_opt_ext(2:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
        d*NMuscle+nqi_q:Nwl:Nw);   
    q_opt_ext(3:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
        d*NMuscle+2*nq.all+nqi_q:Nwl:Nw);  
    q_opt_ext(4:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
        d*NMuscle+2*nq.all+2*nq.all+nqi_q:Nwl:Nw);  
    q_dot_opt_ext(2:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
        d*NMuscle+nqi_q+1:Nwl:Nw);   
    q_dot_opt_ext(3:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
        d*NMuscle+2*nq.all+nqi_q+1:Nwl:Nw);  
    q_dot_opt_ext(4:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
        d*NMuscle+2*nq.all+2*nq.all+nqi_q+1:Nwl:Nw);
end

%Arm activations
a_a_opt_ext=zeros(N*(d+1)+1,nq.arms);
a_a_opt_ext(1:(d+1):end,:)= a_a_opt;
for nmusi=1:nq.arms
    a_a_opt_ext(2:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
        d*NMuscle+d*2*nq.all+nmusi:Nwl:Nw);
    a_a_opt_ext(3:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
        d*NMuscle+d*2*nq.all+nq.arms+nmusi:Nwl:Nw);
    a_a_opt_ext(4:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
        d*NMuscle+d*2*nq.all+nq.arms+nq.arms+nmusi:Nwl:Nw);
end    

%% Unscale results

%States at mesh points

%Qs (1:N-1)
q_opt_unsc.rad = q_opt(1:end-1,:).*repmat(...
    scaling.Qs,size(q_opt(1:end-1,:),1),1); 

%Convert in degrees
q_opt_unsc.deg = q_opt_unsc.rad;
q_opt_unsc.deg(:,[1:3,7:end]) = q_opt_unsc.deg(:,[1:3,7:end]).*180/pi;

%Qs (1:N)
q_opt_unsc_all.rad = q_opt(1:end,:).*repmat(...
    scaling.Qs,size(q_opt(1:end,:),1),1); 

%Convert in degrees
q_opt_unsc_all.deg = q_opt_unsc_all.rad;
q_opt_unsc_all.deg(:,[1:3,7:end]) = ...
    q_opt_unsc_all.deg(:,[1:3,7:end]).*180/pi;       

%Qdots (1:N-1)
qdot_opt_unsc.rad = qdot_opt(1:end-1,:).*repmat(...
    scaling.Qdots,size(qdot_opt(1:end-1,:),1),1);

%Convert in degrees
qdot_opt_unsc.deg = qdot_opt_unsc.rad;
qdot_opt_unsc.deg(:,[1:3,7:end]) = ...
    qdot_opt_unsc.deg(:,[1:3,7:end]).*180/pi;

%Qdots (1:N)
qdot_opt_unsc_all.rad = qdot_opt(1:end,:).*repmat(...
    scaling.Qdots,size(qdot_opt(1:end,:),1),1); 

%Muscle activations
a_opt_unsc = a_opt(1:end-1,:).*repmat(...
    scaling.a,size(a_opt(1:end-1,:),1),size(a_opt,2));

%Muscle-tendon forces
FTtilde_opt_unsc = FTtilde_opt(1:end-1,:).*repmat(...
    scaling.FTtilde,size(FTtilde_opt(1:end-1,:),1),1);

%Arm activations
a_a_opt_unsc = a_a_opt(1:end-1,:).*repmat(...
    scaling.a_a,size(a_a_opt(1:end-1,:),1),size(a_a_opt,2));
%Controls at mesh points

%Time derivative of Qdots
qdotdot_opt_unsc.rad = ...
    qdotdot_opt.*repmat(scaling.Qdotdots,size(qdotdot_opt,1),1);

%Convert in degrees
qdotdot_opt_unsc.deg = qdotdot_opt_unsc.rad;
qdotdot_opt_unsc.deg(:,[1:3,7:end]) = ...
    qdotdot_opt_unsc.deg(:,[1:3,7:end]).*180/pi;

%Time derivative of muscle activations (states)
vA_opt_unsc = vA_opt.*repmat(scaling.vA,size(vA_opt,1),size(vA_opt,2));
tact = 0.015;
tdeact = 0.06;

%Get muscle excitations from time derivative of muscle activations
e_opt_unsc = computeExcitationRaasch(a_opt_unsc,vA_opt_unsc,...
    ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);

%Time derivative of muscle-tendon forces
dFTtilde_opt_unsc = dFTtilde_opt.*repmat(...
    scaling.dFTtilde,size(dFTtilde_opt,1),size(dFTtilde_opt,2));

%Arm excitations
e_a_opt_unsc = e_a_opt.*repmat(scaling.e_a,size(e_a_opt,1),...
    size(e_a_opt,2));

%% Time grid    

%Mesh points
tgrid = linspace(0,tf_opt,N+1);
dtime = zeros(1,d+1);
for i=1:4
    dtime(i)=tau_root(i)*(tf_opt/N);
end

%Mesh points and collocation points
tgrid_ext = zeros(1,(d+1)*N+1);
for i=1:N
    tgrid_ext(((i-1)*4+1):1:i*4)=tgrid(i)+dtime;
end
tgrid_ext(end) = tf_opt;     
 
%% Joint torques and ground reaction forces at optimal solution

Xk_Qs_Qdots_opt             = zeros(N,2*nq.all);
Xk_Qs_Qdots_opt(:,1:2:end)  = q_opt_unsc.rad;
Xk_Qs_Qdots_opt(:,2:2:end)  = qdot_opt_unsc.rad;
Xk_Qdotdots_opt             = qdotdot_opt_unsc.rad;
out_res_opt = zeros(N,nq.all+NGRF+NcalcOrall);
for i = 1:N
    [res] = F1([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)']);
    out_res_opt(i,:) = full(res);    
end
GRF_opt_unsc = out_res_opt(:,GRFi.all);

%assertArmTmax should be 0
assertArmTmax = max(max(abs(out_res_opt(:,armsi)-(a_a_opt_unsc*...
    scaling.ArmTau - 0.1*Xk_Qs_Qdots_opt(:,armsi*2)))));    

%% Stride length and width 

%%%%% NOTE: these calculations might need altering for right leg first
%%%%% trial rather than left leg as per Falisse data

%For the stride length we also need the values at the end of the
%interval so N+1 where states but not controls are defined
Xk_Qs_Qdots_opt_all = zeros(N+1,2*size(q_opt_unsc_all.rad,2));
Xk_Qs_Qdots_opt_all(:,1:2:end)  = q_opt_unsc_all.rad;
Xk_Qs_Qdots_opt_all(:,2:2:end)  = qdot_opt_unsc_all.rad;

%We just want to extract the positions of the calcaneus origins so we
%do not really care about Qdotdot that we set to 0
Xk_Qdotdots_opt_all = zeros(N+1,size(q_opt_unsc_all.rad,2));
out_res_opt_all = zeros(N+1,nq.all+NGRF+NcalcOrall);
for i = 1:N+1
    [res] = F1([Xk_Qs_Qdots_opt_all(i,:)';Xk_Qdotdots_opt_all(i,:)']);
    out_res_opt_all(i,:) = full(res);    
end

%The stride length is the distance covered by the calcaneus origin
%Right leg
dist_r = sqrt(f_Jnn3(out_res_opt_all(end,calcOrall.r)-...
    out_res_opt_all(1,calcOrall.r)));

%Left leg
dist_l = sqrt(f_Jnn3(out_res_opt_all(end,calcOrall.l)-...
    out_res_opt_all(1,calcOrall.l)));

%The total stride length is the sum of the right and left stride 
%lengths after a half gait cycle, since we assume symmetry
StrideLength_opt = full(dist_r + dist_l);

%The stride width is the medial distance btween the calcaneus origins
StepWidth_opt = full(abs(out_res_opt_all(:,calcOrall.r(3)) - ...
    out_res_opt_all(:,calcOrall.l(3))));
stride_width_mean = mean(StepWidth_opt);
stride_width_std = std(StepWidth_opt);
    
%% Passive joint torques at optimal solution

Tau_pass_opt_all = zeros(N,15);
for i = 1:N    
    Tau_pass_opt.hip.flex.l    = f_PassiveMoments(k_pass.hip.flex,...
       theta.pass.hip.flex,Xk_Qs_Qdots_opt(i,jointi.hip_flex.l*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.hip_flex.l*2));
    Tau_pass_opt.hip.flex.r    = f_PassiveMoments(k_pass.hip.flex,...
       theta.pass.hip.flex,Xk_Qs_Qdots_opt(i,jointi.hip_flex.r*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.hip_flex.r*2));
    Tau_pass_opt.hip.add.l     = f_PassiveMoments(k_pass.hip.add,...
       theta.pass.hip.add,Xk_Qs_Qdots_opt(i,jointi.hip_add.l*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.hip_add.l*2));
    Tau_pass_opt.hip.add.r     = f_PassiveMoments(k_pass.hip.add,...
       theta.pass.hip.add,Xk_Qs_Qdots_opt(i,jointi.hip_add.r*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.hip_add.r*2));
    Tau_pass_opt.hip.rot.l     = f_PassiveMoments(k_pass.hip.rot,...
       theta.pass.hip.rot,Xk_Qs_Qdots_opt(i,jointi.hip_rot.l*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.hip_rot.l*2));
    Tau_pass_opt.hip.rot.r     = f_PassiveMoments(k_pass.hip.rot,...
       theta.pass.hip.rot,Xk_Qs_Qdots_opt(i,jointi.hip_rot.r*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.hip_rot.r*2));
    Tau_pass_opt.knee.l        = f_PassiveMoments(k_pass.knee,...
       theta.pass.knee,Xk_Qs_Qdots_opt(i,jointi.knee.l*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.knee.l*2));
    Tau_pass_opt.knee.r        = f_PassiveMoments(k_pass.knee,...
       theta.pass.knee,Xk_Qs_Qdots_opt(i,jointi.knee.r*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.knee.r*2));
    Tau_pass_opt.ankle.l       = f_PassiveMoments(k_pass.ankle,...
       theta.pass.ankle,Xk_Qs_Qdots_opt(i,jointi.ankle.l*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.ankle.l*2));
    Tau_pass_opt.ankle.r       = f_PassiveMoments(k_pass.ankle,...
       theta.pass.ankle,Xk_Qs_Qdots_opt(i,jointi.ankle.r*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.ankle.r*2));
    Tau_pass_opt.subt.l       = f_PassiveMoments(k_pass.subt,...
       theta.pass.subt,Xk_Qs_Qdots_opt(i,jointi.subt.l*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.subt.l*2));
    Tau_pass_opt.subt.r       = f_PassiveMoments(k_pass.subt,...
       theta.pass.subt,Xk_Qs_Qdots_opt(i,jointi.subt.r*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.subt.r*2));
    Tau_pass_opt.trunk.ext     = f_PassiveMoments(k_pass.trunk.ext,...
       theta.pass.trunk.ext,Xk_Qs_Qdots_opt(i,jointi.trunk.ext*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.trunk.ext*2));
    Tau_pass_opt.trunk.ben     = f_PassiveMoments(k_pass.trunk.ben,...
       theta.pass.trunk.ben,Xk_Qs_Qdots_opt(i,jointi.trunk.ben*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.trunk.ben*2));
    Tau_pass_opt.trunk.rot     = f_PassiveMoments(k_pass.trunk.rot,...
       theta.pass.trunk.rot,Xk_Qs_Qdots_opt(i,jointi.trunk.rot*2-1),...
       Xk_Qs_Qdots_opt(i,jointi.trunk.rot*2));        
    Tau_pass_opt_all(i,:) = full([Tau_pass_opt.hip.flex.l,...
       Tau_pass_opt.hip.add.l,Tau_pass_opt.hip.rot.l,...             
       Tau_pass_opt.hip.flex.r,Tau_pass_opt.hip.add.r,...
       Tau_pass_opt.hip.rot.r,Tau_pass_opt.knee.l,...
       Tau_pass_opt.knee.r,Tau_pass_opt.ankle.l,...
       Tau_pass_opt.ankle.r,Tau_pass_opt.subt.l,...
       Tau_pass_opt.subt.r,Tau_pass_opt.trunk.ext,...
       Tau_pass_opt.trunk.ben,Tau_pass_opt.trunk.rot]);
end  

%% Assert average speed    

%Distance travelled
dist_trav_opt = q_opt_ext(end,jointi.pelvis.tx)*...
    scaling.Qs(jointi.pelvis.tx) - q_opt_ext(1,jointi.pelvis.tx)*...
    scaling.Qs(jointi.pelvis.tx);

%Time elapsed
time_elaps_opt = tf_opt; % time elapsed

%Velocity
vel_aver_opt = dist_trav_opt/time_elaps_opt; 

%assert_v_tg should be 0
assert_v_tg = abs(vel_aver_opt-v_tgt);

%% Decompose optimal cost

%Starting values
J_opt           = 0;
E_cost          = 0;
A_cost          = 0;
Arm_cost        = 0;
Qdotdot_cost    = 0;
Pass_cost       = 0;
GRF_cost        = 0;
vA_cost         = 0;
dFTtilde_cost   = 0;
QdotdotArm_cost = 0;
count           = 1;
h_opt           = tf_opt/N;

for k=1:N      
    % Get muscle-tendon lengths, velocities, moment arms
    % Left leg
    qin_l_opt_all = [Xk_Qs_Qdots_opt(k,jointi.hip_flex.l*2-1),...
        Xk_Qs_Qdots_opt(k,jointi.hip_add.l*2-1), ...
        Xk_Qs_Qdots_opt(k,jointi.hip_rot.l*2-1), ...
        Xk_Qs_Qdots_opt(k,jointi.knee.l*2-1), ...
        Xk_Qs_Qdots_opt(k,jointi.ankle.l*2-1),...
        Xk_Qs_Qdots_opt(k,jointi.subt.l*2-1),...
        Xk_Qs_Qdots_opt(k,jointi.trunk.ext*2-1),...
        Xk_Qs_Qdots_opt(k,jointi.trunk.ben*2-1),...
        Xk_Qs_Qdots_opt(k,jointi.trunk.rot*2-1)];  
    qdotin_l_opt_all = [Xk_Qs_Qdots_opt(k,jointi.hip_flex.l*2),...
        Xk_Qs_Qdots_opt(k,jointi.hip_add.l*2),...
        Xk_Qs_Qdots_opt(k,jointi.hip_rot.l*2),...
        Xk_Qs_Qdots_opt(k,jointi.knee.l*2),...
        Xk_Qs_Qdots_opt(k,jointi.ankle.l*2),...
        Xk_Qs_Qdots_opt(k,jointi.subt.l*2),...
        Xk_Qs_Qdots_opt(k,jointi.trunk.ext*2),...
        Xk_Qs_Qdots_opt(k,jointi.trunk.ben*2),...
        Xk_Qs_Qdots_opt(k,jointi.trunk.rot*2)];  
    [lMTk_l_opt_all,vMTk_l_opt_all,~] = ...
        f_lMT_vMT_dM(qin_l_opt_all,qdotin_l_opt_all);    
    % Right leg
    qin_r_opt_all = [Xk_Qs_Qdots_opt(k,jointi.hip_flex.r*2-1),...
        Xk_Qs_Qdots_opt(k,jointi.hip_add.r*2-1),...
        Xk_Qs_Qdots_opt(k,jointi.hip_rot.r*2-1),...
        Xk_Qs_Qdots_opt(k,jointi.knee.r*2-1),...
        Xk_Qs_Qdots_opt(k,jointi.ankle.r*2-1),...
        Xk_Qs_Qdots_opt(k,jointi.subt.r*2-1),...
        Xk_Qs_Qdots_opt(k,jointi.trunk.ext*2-1),...
        Xk_Qs_Qdots_opt(k,jointi.trunk.ben*2-1),...
        Xk_Qs_Qdots_opt(k,jointi.trunk.rot*2-1)];  
    qdotin_r_opt_all = [Xk_Qs_Qdots_opt(k,jointi.hip_flex.r*2),...
        Xk_Qs_Qdots_opt(k,jointi.hip_add.r*2),...
        Xk_Qs_Qdots_opt(k,jointi.hip_rot.r*2),...
        Xk_Qs_Qdots_opt(k,jointi.knee.r*2),...
        Xk_Qs_Qdots_opt(k,jointi.ankle.r*2),...
        Xk_Qs_Qdots_opt(k,jointi.subt.r*2),...
        Xk_Qs_Qdots_opt(k,jointi.trunk.ext*2),...
        Xk_Qs_Qdots_opt(k,jointi.trunk.ben*2),...
        Xk_Qs_Qdots_opt(k,jointi.trunk.rot*2)];      
    [lMTk_r_opt_all,vMTk_r_opt_all,~] = ...
        f_lMT_vMT_dM(qin_r_opt_all,qdotin_r_opt_all);
    % Both legs
    lMTk_lr_opt_all = ...
        [lMTk_l_opt_all([1:43,47:49],1);lMTk_r_opt_all(1:46,1)];
    vMTk_lr_opt_all = ...
        [vMTk_l_opt_all([1:43,47:49],1);vMTk_r_opt_all(1:46,1)];        
    % Metabolic energy rate
    [~,~,Fce_opt_all,Fpass_opt_all,Fiso_opt_all,vMmax_opt_all,...
        massM_opt_all] = f_forceEquilibrium_FtildeState_all_tendon(...
            a_opt_unsc(k,:)',FTtilde_opt_unsc(k,:)',...
            dFTtilde_opt_unsc(k,:)',full(lMTk_lr_opt_all),...
            full(vMTk_lr_opt_all),tensions,aTendon,shift);                  
    [~,lMtilde_opt_all] = f_FiberLength_TendonForce_tendon(...
        FTtilde_opt_unsc(k,:)',full(lMTk_lr_opt_all),aTendon,shift);                
    [vM_opt_all,~] = ...
        f_FiberVelocity_TendonForce_tendon(FTtilde_opt_unsc(k,:)',...
        dFTtilde_opt_unsc(k,:)',full(lMTk_lr_opt_all),...
        full(vMTk_lr_opt_all),aTendon,shift);   

    if mE == 0 % Bhargava et al. (2004)
        [e_tot_all,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(...
            a_opt_unsc(k,:)',a_opt_unsc(k,:)',full(lMtilde_opt_all),...
            full(vM_opt_all),full(Fce_opt_all)',full(Fpass_opt_all)',...
            full(massM_opt_all)',pctsts,full(Fiso_opt_all)',...
            MTparameters_m(1,:)',body_mass,10); 
    elseif mE == 1 % Umberger et al. (2003)
        % vMtilde defined for this model as vM/lMopt
        vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
        [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2003all(...
            a_opt_unsc(k,:)',a_opt_unsc(k,:)',full(lMtilde_opt_all),...
            vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
            full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
            full(Fiso_opt_all)',body_mass,10);                
    elseif mE == 2 % Umberger (2010)
        % vMtilde defined for this model as vM/lMopt
        vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
        [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2010all(...
            a_opt_unsc(k,:)',a_opt_unsc(k,:)',full(lMtilde_opt_all),...
            vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
            full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
            full(Fiso_opt_all)',body_mass,10);            
    elseif mE == 3 % Uchida et al. (2016)
        % vMtilde defined for this model as vM/lMopt
        vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
        [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2016all(...
            a_opt_unsc(k,:)',a_opt_unsc(k,:)',full(lMtilde_opt_all),...
            vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
            full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
            full(Fiso_opt_all)',body_mass,10);      
    elseif mE == 4 % Umberger (2010) treating muscle lengthening 
        % heat rate as Umberger et al. (2003)
        % vMtilde defined for this model as vM/lMopt
        vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
        [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2010all_hl(...
            a_opt_unsc(k,:)',a_opt_unsc(k,:)',full(lMtilde_opt_all),...
            vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
            full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
            full(Fiso_opt_all)',body_mass,10);  
    elseif mE == 5 % Umberger (2010) treating negative mechanical 
        % work as Umberger et al. (2003)
        % vMtilde defined for this model as vM/lMopt
        vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
        [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2010all_neg(...
            a_opt_unsc(k,:)',a_opt_unsc(k,:)',full(lMtilde_opt_all),...
            vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
            full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
            full(Fiso_opt_all)',body_mass,10);  
    end
    e_tot_opt_all = full(e_tot_all)';
    for j=1:d                     
        J_opt = J_opt + 1/(dist_trav_opt)*(...
            W.E*B(j+1)      *...
                (f_J92exp(e_tot_opt_all,exp_E))/body_mass*h_opt + ...                   
            W.A*B(j+1)      *(f_J92(a_opt_ext_col(count,:)))*h_opt +...      
            W.ArmE*B(j+1)   *(f_J8(e_a_opt(k,:)))*h_opt +...                    
            W.Ak*B(j+1)     *...
                (f_J21(qdotdot_opt(k,residuals_noarmsi)))*h_opt + ...                     
            W.passMom*B(j+1)*(f_J15(Tau_pass_opt_all(k,:)))*h_opt + ...                    
            W.u*B(j+1)      *(f_J92(vA_opt(k,:)))*h_opt + ...        
            W.u*B(j+1)      *(f_J92(dFTtilde_opt(k,:)))*h_opt + ...                           
            W.u*B(j+1)      *(f_J8(qdotdot_opt(k,armsi)))*h_opt);                                   
            E_cost = E_cost + W.E*B(j+1)*...
                (f_J92exp(e_tot_opt_all,exp_E))/body_mass*h_opt;
            A_cost = A_cost + W.A*B(j+1)*...
                (f_J92(a_opt_ext_col(count,:)))*h_opt;
            Arm_cost = Arm_cost + W.ArmE*B(j+1)*...
                (f_J8(e_a_opt(k,:)))*h_opt;
            Qdotdot_cost = Qdotdot_cost + W.Ak*B(j+1)*...
                (f_J21(qdotdot_opt(k,residuals_noarmsi)))*h_opt;
            Pass_cost = Pass_cost + W.passMom*B(j+1)*...
                (f_J15(Tau_pass_opt_all(k,:)))*h_opt;
            vA_cost = vA_cost + W.u*B(j+1)*...
                (f_J92(vA_opt(k,:)))*h_opt;
            dFTtilde_cost = dFTtilde_cost + W.u*B(j+1)*...
                (f_J92(dFTtilde_opt(k,:)))*h_opt;
            QdotdotArm_cost = QdotdotArm_cost + W.u*B(j+1)*...
                (f_J8(qdotdot_opt(k,armsi)))*h_opt;                
            count = count + 1;  
    end
end

%Extract values
J_optf = full(J_opt);     
E_costf = full(E_cost);
A_costf = full(A_cost);
Arm_costf = full(Arm_cost);
Qdotdot_costf = full(Qdotdot_cost);
Pass_costf = full(Pass_cost);
vA_costf = full(vA_cost);
dFTtilde_costf = full(dFTtilde_cost);
QdotdotArm_costf = full(QdotdotArm_cost);

% assertCost should be 0 
assertCost = J_optf - 1/(dist_trav_opt)*(E_costf+A_costf+Arm_costf+...
    Qdotdot_costf+Pass_costf+vA_costf+dFTtilde_costf+...
    QdotdotArm_costf);

%% Reconstruct full gait cycle

%We reconstruct the full gait cycle from the simulated half gait cycle
%Identify heel strike
threshold = 20; % there is foot-ground contact above the threshold

% % % % We adjust some thresholds manually
% % % if ww == 17
% % %     threshold = 31;
% % % elseif ww == 133
% % %     threshold = 21;    
% % % elseif ww == 159
% % %     threshold = 42;
% % % elseif ww == 160
% % %     threshold = 44;
% % % elseif ww == 161
% % %     threshold = 57;
% % % elseif ww == 162
% % %     threshold = 69;
% % % end

%Cler any existing heel strike analysis (only relebvant if repeating sims)
if exist('HS1','var')
    clear HS1
end

%Check if heel strike is on the right side   
phase_tran_tgridi = find(GRF_opt_unsc(:,2)<threshold,1,'last');
if ~isempty(phase_tran_tgridi)        
    if phase_tran_tgridi == N
        temp_idx = find(GRF_opt_unsc(:,2)>threshold,1,'first');
        if ~isempty(temp_idx)
            if temp_idx-1 ~= 0 && ...
                    find(GRF_opt_unsc(temp_idx-1,2)<threshold)
                phase_tran_tgridi_t = temp_idx;             
                IC1i = phase_tran_tgridi_t;
                HS1 = 'r';
            end 
        else            
            IC1i = phase_tran_tgridi + 1; 
            HS1 = 'r';
        end
    else            
        IC1i = phase_tran_tgridi + 1; 
        HS1 = 'r';
    end        
end
if ~exist('HS1','var')
    % Check if heel strike is on the left side  
    phase_tran_tgridi = find(GRF_opt_unsc(:,5)<threshold,1,'last');       
    if phase_tran_tgridi == N
        temp_idx = find(GRF_opt_unsc(:,5)>threshold,1,'first');
        if ~isempty(temp_idx)  
            if temp_idx-1 ~= 0 && ...
                    find(GRF_opt_unsc(temp_idx-1,5)<threshold)
                phase_tran_tgridi_t = temp_idx;             
                IC1i = phase_tran_tgridi_t;
                HS1 = 'l';
            else
                IC1i = phase_tran_tgridi + 1; 
                HS1 = 'l';
            end 
        else
            IC1i = phase_tran_tgridi + 1; 
            HS1 = 'l';
        end
    else            
        IC1i = phase_tran_tgridi + 1; 
        HS1 = 'l';
    end        
end
if isempty(phase_tran_tgridi)
    disp('No heel strike detected, consider increasing the threshold');
end

%Joint kinematics
%Helper variables to reconstruct full gait cycle assuming symmetry
QsSymA = [jointi.pelvis.tilt,jointi.pelvis.ty,...
    jointi.hip_flex.l:jointi.trunk.ext,...
    jointi.sh_flex.l:jointi.elb.r];
QsSymB = [jointi.pelvis.tilt,jointi.pelvis.ty,...    
    jointi.hip_flex.r:jointi.hip_rot.r,...
    jointi.hip_flex.l:jointi.hip_rot.l,...
    jointi.knee.r,jointi.knee.l,...
    jointi.ankle.r,jointi.ankle.l,...
    jointi.subt.r,jointi.subt.l,...
    jointi.trunk.ext,...
    jointi.sh_flex.r:jointi.sh_rot.r,...
    jointi.sh_flex.l:jointi.sh_rot.l,...
    jointi.elb.r,jointi.elb.l]; 
QsOpp = [jointi.pelvis.list:jointi.pelvis.rot,jointi.pelvis.tz,...
    jointi.trunk.ben:jointi.trunk.rot];
QsSymA_ptx = [jointi.pelvis.tilt,jointi.pelvis.tx,...
    jointi.pelvis.ty,...
    jointi.hip_flex.l:jointi.trunk.ext,...
    jointi.sh_flex.l:jointi.elb.r];
QsSymB_ptx = [jointi.pelvis.tilt,jointi.pelvis.tx,...
    jointi.pelvis.ty,...    
    jointi.hip_flex.r:jointi.hip_rot.r,...
    jointi.hip_flex.l:jointi.hip_rot.l,...
    jointi.knee.r,jointi.knee.l,...
    jointi.ankle.r,jointi.ankle.l,...
    jointi.subt.r,jointi.subt.l,...
    jointi.trunk.ext,...
    jointi.sh_flex.r:jointi.sh_rot.r,...
    jointi.sh_flex.l:jointi.sh_rot.l,...
    jointi.elb.r,jointi.elb.l]; 

%Qs
Qs_opt = zeros(N*2,size(q_opt_unsc.deg,2));
Qs_opt(1:N-IC1i+1,:) = q_opt_unsc.deg(IC1i:end,:);   
Qs_opt(N-IC1i+2:N-IC1i+1+N,QsSymA) = q_opt_unsc.deg(1:end,QsSymB);
Qs_opt(N-IC1i+2:N-IC1i+1+N,QsOpp) = -q_opt_unsc.deg(1:end,QsOpp);    
Qs_opt(N-IC1i+2:N-IC1i+1+N,jointi.pelvis.tx) = ...
    q_opt_unsc.deg(1:end,jointi.pelvis.tx) + ...
    q_opt_unsc_all.deg(end,jointi.pelvis.tx);        
Qs_opt(N-IC1i+2+N:2*N,:) = q_opt_unsc.deg(1:IC1i-1,:);    
Qs_opt(N-IC1i+2+N:2*N,jointi.pelvis.tx) = ...
    q_opt_unsc.deg(1:IC1i-1,jointi.pelvis.tx) + ...
    2*q_opt_unsc_all.deg(end,jointi.pelvis.tx);

%If the first heel strike was on the left foot then we invert so that
%we always start with the right foot, for analysis purpose
if strcmp(HS1,'l')         
    Qs_opt(:,QsSymA_ptx)  = Qs_opt(:,QsSymB_ptx);
    Qs_opt(:,QsOpp)       = -Qs_opt(:,QsOpp);    
end
temp_Qs_opt_pelvis_tx = Qs_opt(1,jointi.pelvis.tx);
Qs_opt(:,jointi.pelvis.tx) = Qs_opt(:,jointi.pelvis.tx)-...
    temp_Qs_opt_pelvis_tx;

%For visualization in OpenSim GUI
q_opt_GUI_GC = zeros(2*N,1+nq.all+2);
q_opt_GUI_GC(1:N-IC1i+1,1) = tgrid(:,IC1i:end-1)';
q_opt_GUI_GC(N-IC1i+2:N-IC1i+1+N,1)  = tgrid(:,1:end-1)' + tgrid(end);
q_opt_GUI_GC(N-IC1i+2+N:2*N,1) = tgrid(:,1:IC1i-1)' + 2*tgrid(end);    
q_opt_GUI_GC(:,2:end-2) = Qs_opt;
q_opt_GUI_GC(:,end-1:end) = 1.51*180/pi*ones(2*N,2); % pro_sup (locked)
q_opt_GUI_GC(:,1) = q_opt_GUI_GC(:,1)-q_opt_GUI_GC(1,1);

%Create .mot file for OpenSim GUI

%Set joint angle labels
JointAngle.labels = {'time','pelvis_tilt','pelvis_list',...
    'pelvis_rotation','pelvis_tx','pelvis_ty','pelvis_tz',...
    'hip_flexion_l','hip_adduction_l','hip_rotation_l',...
    'hip_flexion_r','hip_adduction_r','hip_rotation_r',...
    'knee_angle_l','knee_angle_r','ankle_angle_l',...
    'ankle_angle_r','subtalar_angle_l','subtalar_angle_r',...
    'lumbar_extension','lumbar_bending','lumbar_rotation',...
    'arm_flex_l','arm_add_l','arm_rot_l',...
    'arm_flex_r','arm_add_r','arm_rot_r',...
    'elbow_flex_l','elbow_flex_r',...
    'pro_sup_l','pro_sup_r'};

%Two gait cycles
q_opt_GUI_GC_2 = [q_opt_GUI_GC;q_opt_GUI_GC];
q_opt_GUI_GC_2(2*N+1:4*N,1) = q_opt_GUI_GC_2(2*N+1:4*N,1) + ...
    q_opt_GUI_GC_2(end,1) + ...
    q_opt_GUI_GC_2(end,1)-q_opt_GUI_GC_2(end-1,1);
q_opt_GUI_GC_2(2*N+1:4*N,jointi.pelvis.tx+1) = ...
    q_opt_GUI_GC_2(2*N+1:4*N,jointi.pelvis.tx+1) + ...
    2*q_opt_unsc_all.deg(end,jointi.pelvis.tx);

%Set joint angle data
JointAngle.data = q_opt_GUI_GC_2;

%Set filename
filenameJointAngles = [pathResults,'\',saveName,'\',saveFileName,'_kinematics.mot'];

%Write motion file using helper function
write_motionFile(JointAngle, filenameJointAngles)

%Qdots
Qdots_opt = zeros(N*2,size(q_opt,2));
Qdots_opt(1:N-IC1i+1,:) = qdot_opt_unsc.deg(IC1i:end,:);
Qdots_opt(N-IC1i+2:N-IC1i+1+N,QsSymA_ptx) = ...
    qdot_opt_unsc.deg(1:end,QsSymB_ptx);
Qdots_opt(N-IC1i+2:N-IC1i+1+N,QsOpp) = ...
    -qdot_opt_unsc.deg(1:end,QsOpp);
Qdots_opt(N-IC1i+2+N:2*N,:) = qdot_opt_unsc.deg(1:IC1i-1,:);

%If the first heel strike was on the left foot then we invert so that
%we always start with the right foot, for analysis purpose
if strcmp(HS1,'l')
    Qdots_opt(:,QsSymA_ptx) = Qdots_opt(:,QsSymB_ptx);
    Qdots_opt(:,QsOpp) = -Qdots_opt(:,QsOpp);    
end

%Qdotdots
qdotdot_opt_GC = zeros(N*2,size(q_opt,2));
qdotdot_opt_GC(1:N-IC1i+1,:) = qdotdot_opt_unsc.deg(IC1i:end,:);
qdotdot_opt_GC(N-IC1i+2:N-IC1i+1+N,QsSymA_ptx) = ...
    qdotdot_opt_unsc.deg(1:end,QsSymB_ptx);
qdotdot_opt_GC(N-IC1i+2:N-IC1i+1+N,QsOpp) = ...
    -qdotdot_opt_unsc.deg(1:end,QsOpp);
qdotdot_opt_GC(N-IC1i+2+N:2*N,:) = qdotdot_opt_unsc.deg(1:IC1i-1,:);

%If the first heel strike was on the left foot then we invert so that
%we always start with the right foot, for analysis purpose
if strcmp(HS1,'l')
    qdotdot_opt_GC(:,QsSymA_ptx) = qdotdot_opt_GC(:,QsSymB_ptx);
    qdotdot_opt_GC(:,QsOpp) = -qdotdot_opt_GC(:,QsOpp);    
end

%Ground reaction forces
GRFs_opt = zeros(N*2,NGRF);
GRFs_opt(1:N-IC1i+1,:) = GRF_opt_unsc(IC1i:end,1:6);
GRFs_opt(N-IC1i+2:N-IC1i+1+N,:) = GRF_opt_unsc(1:end,[4:6,1:3]);
GRFs_opt(N-IC1i+2:N-IC1i+1+N,[3,6]) = ...
    -GRFs_opt(N-IC1i+2:N-IC1i+1+N,[3,6]);
GRFs_opt(N-IC1i+2+N:2*N,:) = GRF_opt_unsc(1:IC1i-1,1:6);
GRFs_opt = GRFs_opt./(body_weight/100);

%If the first heel strike was on the left foot then we invert so that
%we always start with the right foot, for analysis purpose
if strcmp(HS1,'l')
    GRFs_opt(:,[4:6,1:3]) = GRFs_opt(:,:);
    GRFs_opt(:,[3,6]) = -GRFs_opt(:,[3,6]);
end

%Joint torques
Ts_opt = zeros(N*2,size(q_opt,2));
Ts_opt(1:N-IC1i+1,1:nq.all) = ...
    out_res_opt(IC1i:end,1:nq.all)./body_mass;
Ts_opt(N-IC1i+2:N-IC1i+1+N,QsSymA_ptx) = ...
    out_res_opt(1:end,QsSymB_ptx)./body_mass;
Ts_opt(N-IC1i+2:N-IC1i+1+N,QsOpp) = ...
    -out_res_opt(1:end,QsOpp)./body_mass;
Ts_opt(N-IC1i+2+N:2*N,1:nq.all) = ...
    out_res_opt(1:IC1i-1,1:nq.all)./body_mass;

%If the first heel strike was on the left foot then we invert so that
%we always start with the right foot, for analysis purpose
if strcmp(HS1,'l')
    Ts_opt(:,QsSymA_ptx) = Ts_opt(:,QsSymB_ptx);
    Ts_opt(:,QsOpp) = -Ts_opt(:,QsOpp);    
end

%Muscle-Tendon Forces
orderMusInv = [NMuscle/2+1:NMuscle,1:NMuscle/2];
FTtilde_opt_GC = zeros(N*2,NMuscle);
FTtilde_opt_GC(1:N-IC1i+1,:) = FTtilde_opt_unsc(IC1i:end,:);
FTtilde_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = ...
    FTtilde_opt_unsc(1:end,orderMusInv);
FTtilde_opt_GC(N-IC1i+2+N:2*N,:) = FTtilde_opt_unsc(1:IC1i-1,:);

%If the first heel strike was on the left foot then we invert so that
%we always start with the right foot, for analysis purpose
if strcmp(HS1,'l')
    FTtilde_opt_GC(:,:) = FTtilde_opt_GC(:,orderMusInv);
end

%Muscle activations
Acts_opt = zeros(N*2,NMuscle);
Acts_opt(1:N-IC1i+1,:) = a_opt_unsc(IC1i:end,:);
Acts_opt(N-IC1i+2:N-IC1i+1+N,:) = a_opt_unsc(1:end,orderMusInv);
Acts_opt(N-IC1i+2+N:2*N,:) = a_opt_unsc(1:IC1i-1,:);

%If the first heel strike was on the left foot then we invert so that
%we always start with the right foot, for analysis purpose
if strcmp(HS1,'l')
    Acts_opt(:,:) = Acts_opt(:,orderMusInv);
end

%Time derivative of muscle-tendon force
dFTtilde_opt_GC = zeros(N*2,NMuscle);
dFTtilde_opt_GC(1:N-IC1i+1,:) = dFTtilde_opt_unsc(IC1i:end,:);
dFTtilde_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = ...
    dFTtilde_opt_unsc(1:end,orderMusInv);
dFTtilde_opt_GC(N-IC1i+2+N:2*N,:) = dFTtilde_opt_unsc(1:IC1i-1,:);

%If the first heel strike was on the left foot then we invert so that
%we always start with the right foot, for analysis purpose
if strcmp(HS1,'l')
    dFTtilde_opt_GC(:,:) = dFTtilde_opt_GC(:,orderMusInv);
end

%Muscle excitations
vA_opt_GC = zeros(N*2,NMuscle);
vA_opt_GC(1:N-IC1i+1,:) = vA_opt_unsc(IC1i:end,:);
vA_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = vA_opt_unsc(1:end,orderMusInv);
vA_opt_GC(N-IC1i+2+N:2*N,:) = vA_opt_unsc(1:IC1i-1,:);

%If the first heel strike was on the left foot then we invert so that
%we always start with the right foot, for analysis purpose
if strcmp(HS1,'l')
    vA_opt_GC(:,:) = vA_opt_GC(:,orderMusInv);
end
tact = 0.015;
tdeact = 0.06;
e_opt_GC = computeExcitationRaasch(Acts_opt,vA_opt_GC,...
    ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);   

%Arm activations
orderArmInv = [jointi.sh_flex.r:jointi.sh_rot.r,...
        jointi.sh_flex.l:jointi.sh_rot.l,...
        jointi.elb.r:jointi.elb.r,...
        jointi.elb.l:jointi.elb.l]-jointi.sh_flex.l+1; 
a_a_opt_GC = zeros(N*2,nq.arms);
a_a_opt_GC(1:N-IC1i+1,:) = a_a_opt_unsc(IC1i:end,:);
a_a_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = a_a_opt_unsc(1:end,orderArmInv);
a_a_opt_GC(N-IC1i+2+N:2*N,:) = a_a_opt_unsc(1:IC1i-1,:);

%If the first heel strike was on the left foot then we invert so that
%we always start with the right foot, for analysis purpose
if strcmp(HS1,'l')
    a_a_opt_GC(:,:) = a_a_opt_GC(:,orderArmInv);
end

% Arm excitations
e_a_opt_GC = zeros(N*2,nq.arms);
e_a_opt_GC(1:N-IC1i+1,:) = e_a_opt_unsc(IC1i:end,:);
e_a_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = e_a_opt_unsc(1:end,orderArmInv);
e_a_opt_GC(N-IC1i+2+N:2*N,:) = e_a_opt_unsc(1:IC1i-1,:);
% If the first heel strike was on the left foot then we invert so that
% we always start with the right foot, for analysis purpose
if strcmp(HS1,'l')
    e_a_opt_GC(:,:) = e_a_opt_GC(:,orderArmInv);
end  

%Passive joint torques
Tau_pass_opt_inv = [4:6,1:3,8,7,10,9,12,11,13:15];
Tau_pass_opt_GC = zeros(N*2,15);
Tau_pass_opt_GC(1:N-IC1i+1,:) = Tau_pass_opt_all(IC1i:end,:); 
Tau_pass_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = ...
    Tau_pass_opt_all(1:end,Tau_pass_opt_inv);
Tau_pass_opt_GC(N-IC1i+2+N:2*N,:) = Tau_pass_opt_all(1:IC1i-1,:); 

%If the first heel strike was on the left foot then we invert so that
%we always start with the right foot, for analysis purpose
if strcmp(HS1,'l')
    Tau_pass_opt_GC(:,Tau_pass_opt_inv) = Tau_pass_opt_GC(:,:);
end

%% Metabolic cost of transport

Qs_opt_rad = Qs_opt;
Qs_opt_rad(:,[1:3,7:end]) = Qs_opt_rad(:,[1:3,7:end]).*pi/180;
qdot_opt_GC_rad = Qdots_opt;
qdot_opt_GC_rad(:,[1:3,7:end])= qdot_opt_GC_rad(:,[1:3,7:end]).*pi/180;   

% Pre-allocations
e_mo_opt = zeros(2*N,1);     
e_mo_optb = zeros(2*N,1);    
vMtilde_opt_all = zeros(2*N, NMuscle);
lMtilde_opt_all = zeros(2*N, NMuscle);

for nn = 1:2*N
    % Get muscle-tendon lengths, velocities, moment arms
    % Left leg
    qin_l_opt = [Qs_opt_rad(nn,jointi.hip_flex.l),...
        Qs_opt_rad(nn,jointi.hip_add.l), ...
        Qs_opt_rad(nn,jointi.hip_rot.l), ...
        Qs_opt_rad(nn,jointi.knee.l), ...
        Qs_opt_rad(nn,jointi.ankle.l),...
        Qs_opt_rad(nn,jointi.subt.l),...
        Qs_opt_rad(nn,jointi.trunk.ext),...
        Qs_opt_rad(nn,jointi.trunk.ben),...
        Qs_opt_rad(nn,jointi.trunk.rot)];  
    qdotin_l_opt = [qdot_opt_GC_rad(nn,jointi.hip_flex.l),...
        qdot_opt_GC_rad(nn,jointi.hip_add.l),...
        qdot_opt_GC_rad(nn,jointi.hip_rot.l),...
        qdot_opt_GC_rad(nn,jointi.knee.l),...
        qdot_opt_GC_rad(nn,jointi.ankle.l),...
        qdot_opt_GC_rad(nn,jointi.subt.l),...
        qdot_opt_GC_rad(nn,jointi.trunk.ext),...
        qdot_opt_GC_rad(nn,jointi.trunk.ben),...
        qdot_opt_GC_rad(nn,jointi.trunk.rot)];  
    [lMTk_l_opt,vMTk_l_opt,~] = ...
        f_lMT_vMT_dM(qin_l_opt,qdotin_l_opt);    
    % Right leg
    qin_r_opt = [Qs_opt_rad(nn,jointi.hip_flex.r),...
        Qs_opt_rad(nn,jointi.hip_add.r),...
        Qs_opt_rad(nn,jointi.hip_rot.r),...
        Qs_opt_rad(nn,jointi.knee.r),...
        Qs_opt_rad(nn,jointi.ankle.r),...
        Qs_opt_rad(nn,jointi.subt.r),...
        Qs_opt_rad(nn,jointi.trunk.ext),...
        Qs_opt_rad(nn,jointi.trunk.ben),...
        Qs_opt_rad(nn,jointi.trunk.rot)];  
    qdotin_r_opt = [qdot_opt_GC_rad(nn,jointi.hip_flex.r),...
        qdot_opt_GC_rad(nn,jointi.hip_add.r),...
        qdot_opt_GC_rad(nn,jointi.hip_rot.r),...
        qdot_opt_GC_rad(nn,jointi.knee.r),...
        qdot_opt_GC_rad(nn,jointi.ankle.r),...
        qdot_opt_GC_rad(nn,jointi.subt.r),...
        qdot_opt_GC_rad(nn,jointi.trunk.ext),...
        qdot_opt_GC_rad(nn,jointi.trunk.ben),...
        qdot_opt_GC_rad(nn,jointi.trunk.rot)];      
    [lMTk_r_opt,vMTk_r_opt,~] = ...
        f_lMT_vMT_dM(qin_r_opt,qdotin_r_opt);
    % Both legs
    lMTk_lr_opt     = [lMTk_l_opt([1:43,47:49],1);lMTk_r_opt(1:46,1)];
    vMTk_lr_opt     = [vMTk_l_opt([1:43,47:49],1);vMTk_r_opt(1:46,1)];    
    [Hilldiff_optt,FT_optt,Fce_optt,Fpass_optt,Fiso_optt,...
        vMmax_optt,massM_optt] = ...
            f_forceEquilibrium_FtildeState_all_tendon(...
            Acts_opt(nn,:)',FTtilde_opt_GC(nn,:)',...
            dFTtilde_opt_GC(nn,:)',full(lMTk_lr_opt),...
            full(vMTk_lr_opt),tensions,aTendon,shift);             
    [~,lMtilde_opt] = f_FiberLength_TendonForce_tendon(...
        FTtilde_opt_GC(nn,:)',full(lMTk_lr_opt),aTendon,shift);         
    lMtilde_opt_all(nn,:) = full(lMtilde_opt)';   
    [vM_opt,vMtilde_opt] = ...
        f_FiberVelocity_TendonForce_tendon(FTtilde_opt_GC(nn,:)',...
        dFTtilde_opt_GC(nn,:)',full(lMTk_lr_opt),full(vMTk_lr_opt),...
        aTendon,shift);        
    vMtilde_opt_all(nn,:) = full(vMtilde_opt)';    
    if mE == 0 % Bhargava et al. (2004)
    [~,~,~,~,~,e_mot] = ...
        fgetMetabolicEnergySmooth2004all(Acts_opt(nn,:)',...
            Acts_opt(nn,:)',full(lMtilde_opt),full(vM_opt),...
            full(Fce_optt),full(Fpass_optt),full(massM_optt),pctsts,...
            full(Fiso_optt),MTparameters_m(1,:)',body_mass,10);           
    elseif mE == 1 % Umberger et al. (2003)
        % vMtilde defined for this model as vM/lMopt
        vMtildeUmbk_opt = full(vM_opt)./(MTparameters_m(2,:)');
        [~,~,~,~,e_mot] = fgetMetabolicEnergySmooth2003all(...
            Acts_opt(nn,:)',Acts_opt(nn,:)',full(lMtilde_opt),...
            vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
            full(massM_optt)',pctsts,full(vMmax_optt)',...
            full(Fiso_optt)',body_mass,10); 
    elseif mE == 2 % Umberger (2010)
        % vMtilde defined for this model as vM/lMopt
        vMtildeUmbk_opt = full(vM_opt)./(MTparameters_m(2,:)');
        [~,~,~,~,e_mot] = fgetMetabolicEnergySmooth2010all(...
            Acts_opt(nn,:)',Acts_opt(nn,:)',full(lMtilde_opt),...
            vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
            full(massM_optt)',pctsts,full(vMmax_optt)',...
            full(Fiso_optt)',body_mass,10); 
    elseif mE == 3 % Uchida et al. (2016)
        % vMtilde defined for this model as vM/lMopt
        vMtildeUmbk_opt = full(vM_opt)./(MTparameters_m(2,:)');
        [~,~,~,~,e_mot] = fgetMetabolicEnergySmooth2016all(...
            Acts_opt(nn,:)',Acts_opt(nn,:)',full(lMtilde_opt),...
            vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
            full(massM_optt)',pctsts,full(vMmax_optt)',...
            full(Fiso_optt)',body_mass,10); 
    elseif mE == 4 % Umberger (2010) treating muscle lengthening 
        % heat rate as Umberger et al. (2003)
        % vMtilde defined for this model as vM/lMopt
        vMtildeUmbk_opt = full(vM_opt)./(MTparameters_m(2,:)');
        [~,~,~,~,e_mot] = fgetMetabolicEnergySmooth2010all_hl(...
            Acts_opt(nn,:)',Acts_opt(nn,:)',full(lMtilde_opt),...
            vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
            full(massM_optt)',pctsts,full(vMmax_optt)',...
            full(Fiso_optt)',body_mass,10); 
    elseif mE == 5 % Umberger (2010) treating negative mechanical 
        % work as Umberger et al. (2003)
        % vMtilde defined for this model as vM/lMopt
        vMtildeUmbk_opt = full(vM_opt)./(MTparameters_m(2,:)');
        [~,~,~,~,e_mot] = fgetMetabolicEnergySmooth2010all_neg(...
            Acts_opt(nn,:)',Acts_opt(nn,:)',full(lMtilde_opt),...
            vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
            full(massM_optt)',pctsts,full(vMmax_optt)',...
            full(Fiso_optt)',body_mass,10);             
    end   
    e_mo_opt(nn,:) = full(e_mot)';   
    % Bhargava et al. (2004)
    [~,~,~,~,~,e_motb] = ...
        fgetMetabolicEnergySmooth2004all(Acts_opt(nn,:)',...
        Acts_opt(nn,:)',full(lMtilde_opt),full(vM_opt),...
        full(Fce_optt),full(Fpass_optt),full(massM_optt),pctsts,...
        full(Fiso_optt),MTparameters_m(1,:)',body_mass,10);
    e_mo_optb(nn,:) = full(e_motb)'; 
end

%Get COT
dist_trav_opt_GC = Qs_opt_rad(end,jointi.pelvis.tx) - ...
    Qs_opt_rad(1,jointi.pelvis.tx); % distance traveled    
time_GC = q_opt_GUI_GC(:,1);
e_mo_opt_tr = trapz(time_GC,e_mo_opt);
e_mo_opt_trb = trapz(time_GC,e_mo_optb);

%Cost of transport: J/kg/m

%Energy model used for optimization
COT_opt_InOpt = e_mo_opt_tr/body_mass/dist_trav_opt_GC; 

%Energy model from Bhargava et al. (2004)
COT_opt = e_mo_opt_trb/body_mass/dist_trav_opt_GC; 

%% Optimal cost and CPU time

pathDiary = [pathResults,'\',saveName,'\',saveFileName];
[CPU_IPOPT,CPU_NLP,~,Cost,~,~,~,~,OptSol] = readDiary(pathDiary);

%% Save results       

%Structure results

%%%%%% TODO: results structures will be slightly different to Falisse
%%%%%% structures with the different options etc.

