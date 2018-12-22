% Simulate_Stochastic_Swim
% by David Ehrlich, Jan 12, 2016

% Associated with "Control of Movement Initiation Underlies the Development of Balance" 
% Ehrlich DE and Schoppik D, Curr Biol. 2017 Feb 6;27(3):334-344.

% To fit simulated distribution of IEIs to empirical values

%% Input and Output
% 'coeffs' is the solutions for the coefficients of the relative bout likelihood function (alpha, zeta, beta, gamma).
% 'psi' is the coefficient of the Bayesian prior (psi), which scales the relative bout likelihood
% 'IEI_dist' is the probability distribution of IEIs for the simulation.


function IEI_dist = Simulate_Stochastic_Swim (psi, coeffs)

%% Some switches
Prior_Short = 0; % If set to 1, psi is fit with Relative Bout Likelihood held constant.
% Be sure Prior_Short is set identically in any simulations run with 'Psi' calculated here.

Pitch_Correction_Switch = 1; % If set to 1, bout pitch rotation is correlated with pre-bout pitch. if 0,
% bout pitch rotation is randomly drawn from Gaussian distribution
% approximating observed values.

AngVel_Correction_Switch = 1; % If set to 1, bout angular acceleration is correlated with pre-bout angular velocity. if 0,
% bout angular acceleration is randomly drawn from Gaussian distribution
% approximating observed values.

%% Some parameters
Model_Duration = 15000; %duration of simulation in sec
num_Iterations = 20; % how many times is the simulation run?
num_Bins = 25; %number of bins for output (probability distribution of IEIs)
samplerate = 40; %in Hz
SwimBout_Duration = 0.15; % duration of swim bouts in sec (how far is simulation advanced upon bout initiation?)
t_offset = 6; % absolute refractory period for bout initiation, from initiation of previous bout, in samples. Currently set to match 'SwimBout_Duration'

% Relative Bout Likelihood function coefficients
alpha = coeffs(1);
zeta = coeffs(2);
beta = coeffs(3);
gamma = coeffs(4);

%% Calculate constant angular acceleration
%load variable for angular acceleration observed during IEIs
% % % % load('PropBoutIEIangAcc.mat')
% % % % load('PropBoutIEItime.mat')
wdname = pwd;
load(['Group_' wdname(end-4:end) '.mat'])

% isolate IEIs during light phase
[day,~]=lightdarksplit(PropBoutIEItime);
dayIEIangAcc = PropBoutIEIangAcc(day);
% calculate median angular acceleration to apply during simulation
passive_AngAcc = median(dayIEIangAcc(~isnan(dayIEIangAcc)));

%% Define bout pitch correction from correlation of pre-bout pitch and net bout rotation
% % % % % load('PropBoutInitPitch.mat');
% % % % % load('PropBoutAlignedTime.mat');
% % % % % load('PropBoutNetPitchChg.mat');
% isolate frames during light phase
[day,~]=lightdarksplit(PropBoutAlignedTime);
% statistics of net bout rotation
Net_Bout_Rotation_Mean = mean(PropBoutNetPitchChg(day));
Net_Bout_Rotation_SD = std(PropBoutNetPitchChg(day));
% find slope and intercept of correlation of pre-bout pitch and net bout rotation
Net_Bout_Rotation_Fit = polyfit(PropBoutInitPitch(day), PropBoutNetPitchChg(day), 1);
% find correlation coefficient of correlation of pre-bout pitch and net bout rotation
Net_Bout_Rotation_Corr_Coeff = corrcoef(PropBoutInitPitch(day), PropBoutNetPitchChg(day));

%% Define bout angular velocity correction from correlation of pre-bout angular velocity and net angular acceleration
% % % % load('PropBoutAlignedAngVel.mat')
% % % % load('PropBoutAlignedTime.mat')
% isolate frames during light phase
[day,~]=lightdarksplit(PropBoutAlignedTime);
% Bouts were aligned such that fastest translation speed falls on element
% 31.  Pre-bout angular velocity is centered at 175 msec preceding, and
% post-bout angular velocity begins 250 msec following peak translation.
Prebout_AV = mean(PropBoutAlignedAngVel(day,22:26),2);
Postbout_AV = mean(PropBoutAlignedAngVel(day,41:46),2);
% statistics of net bout angular acceleration
Net_Bout_AngAcc_Mean = mean(Postbout_AV-Prebout_AV);
Net_Bout_AngAcc_SD = std(Postbout_AV-Prebout_AV);
% Variables defining correlation of pre-bout angular velocity and net angular acceleration,
% distinct for nose-up and nose-down pre-bout angular velocity
up_inds = find(Prebout_AV > 0);
dn_inds = find(Prebout_AV < 0);
up_AngAcc_Fit = polyfit(Prebout_AV(up_inds), Postbout_AV(up_inds)-Prebout_AV(up_inds),1);
up_AngAcc_Rsquared = corrcoef(Prebout_AV(up_inds), Postbout_AV(up_inds)-Prebout_AV(up_inds));
dn_AngAcc_Fit = polyfit(Prebout_AV(dn_inds), Postbout_AV(dn_inds)-Prebout_AV(dn_inds),1);
dn_AngAcc_Rsquared = corrcoef(Prebout_AV(dn_inds), Postbout_AV(dn_inds)-Prebout_AV(dn_inds));

%% Simulate swimming
for model_Iteration = 1:num_Iterations
    % Initialize variables for appending simulated bout properties
    Net_Bout_AngAcc = []; Net_Bout_Rotation = [];
    Pre_Bout_Pitch = []; Pre_Bout_AngVel = [];
    % Initialize pitch and angular velocity variables (each simulation
    % begins at horiziontal and with no angular velocity).
    Pitch = zeros(1);
    AngVel = zeros(1);

    % The simulation advances time-steps using a while loop.
    t=1; % 't' is the time index.

    % The simulation begins as though a bout was just terminated, 
    % required for determining time-variant bout initiation. 'Bout_Index' is appended
    % with the initiation time of each bout and is used to calculated IEIs.
    Bout_Index = -1;  
    
    % advance until the simulation reaches the set duration
    while t < (Model_Duration * samplerate)                                
        t = t+1;       
        % Each time-step, angular acceleration is used to calculate new angular velocity 
        AngVel(t) = AngVel(t-1) + cosd(Pitch(t-1))*passive_AngAcc/samplerate;
        % Angular velocity is used to calculate new pitch
        Pitch(t) = Pitch(t-1) + AngVel(t)/samplerate;
        
        % limit pitches to +/- 180 deg
        if Pitch(t) < -180
            Pitch(t) = Pitch(t) + 360;
        elseif Pitch(t) > 180
            Pitch(t) = Pitch(t) -360;
        end 
        
        %% Bout initiation mechanism
        if Prior_Short==1
            % Calculate a time-variant but posture-invariant P_bout from
            % the Bayesian prior only
            P_bout(t) = psi*sqrt((t-Bout_Index(end) - t_offset)/samplerate);
        else
            % Calculate a time- and posture-variant P_bout, from the
            % product of the Bayesian prior and the Relative Bout Likelihood
            RelativeBoutLikelihood(t) = (zeta + alpha*abs(Pitch(t)) + beta*abs(AngVel(t)) - gamma*(AngVel(t)));
            P_bout(t) = psi * sqrt((t-Bout_Index(end) - t_offset)/samplerate) * RelativeBoutLikelihood(t);
        end
        
        %% Initiate bout if a random number is smaller than P_bout
        if rand(1) < P_bout(t)            
            % Calculate net pitch change across bout
            if Pitch_Correction_Switch == 0
                % If bouts are not corrective for pitch, Net Bout Rotation
                % is randomly drawn from Gaussian approxmiate of observed distribution.
                Net_Bout_Rotation(end+1) = normrnd( Net_Bout_Rotation_Mean, Net_Bout_Rotation_SD);
            else %Net bout rotation is correlated with pre-bout pitch
                Net_Bout_Rotation(end+1) = normrnd( Net_Bout_Rotation_Fit(1)*Pitch(t) + Net_Bout_Rotation_Fit(2), Net_Bout_Rotation_SD * (1-Net_Bout_Rotation_Corr_Coeff(1,2)^2));        
            end
            % Calculat net angular velocity change across bout
            if AngVel_Correction_Switch == 0
                % If bouts are not corrective for angular velocity, Net
                % Bout Angular Acceleration is randomly drawn from Gaussian approxmiate of observed distribution.
                Net_Bout_AngAcc(end+1) = normrnd( Net_Bout_AngAcc_Mean, Net_Bout_AngAcc_SD);
            else    %Bout angular acceleration is correlated with pre-bout angular velocity
                if AngVel(t) > 0
                    Net_Bout_AngAcc(end+1) = normrnd( up_AngAcc_Fit(1)*AngVel(t) + up_AngAcc_Fit(2), Net_Bout_AngAcc_SD * (1-up_AngAcc_Rsquared(1,2)^2));
                else
                    Net_Bout_AngAcc(end+1) = normrnd( dn_AngAcc_Fit(1)*AngVel(t) + dn_AngAcc_Fit(2), Net_Bout_AngAcc_SD * (1-dn_AngAcc_Rsquared(1,2)^2));
                end                
            end
            
            % Append bout initiation time
            Bout_Index(end+1) = t;
            % Append pre-bout posture for diagnostics
            Pre_Bout_Pitch(end+1) = Pitch(t);
            Pre_Bout_AngVel(end+1) = AngVel(t);

            % Advance by bout duration
            t = t + SwimBout_Duration * samplerate;
            % Apply net bout rotation and net bout angular acceleration
            % across bout
            Pitch(t) = Pitch(t-(SwimBout_Duration * samplerate)) + Net_Bout_Rotation(end);
            AngVel(t) = AngVel(t-(SwimBout_Duration * samplerate)) + Net_Bout_AngAcc(end); 
            
            %interpolate pitch during bout (skipped time-steps)
            Pitch(t-(SwimBout_Duration * samplerate)+1:t-1) = linspace(Pitch(t-(SwimBout_Duration * samplerate)),Pitch(t), (SwimBout_Duration * samplerate)-1);                       
        end
    end

    % Calculate IEIs from bout initiation times
    modelIEIs = diff(Bout_Index)/samplerate;

    % Crop pitch and angular velocity if bout advanced past Model_Duration (for plotting)
    if length(Pitch) > Model_Duration * samplerate
        Pitch(Model_Duration * samplerate + 1: end) = [];
    end
    if length(AngVel) > Model_Duration * samplerate
        AngVel(Model_Duration * samplerate + 1: end) = [];
    end
    
    % Calculate IEI histogram
    edges = linspace(0,20,num_Bins+1);
    IEI_dist = histc(modelIEIs,edges);
    % Make probability distribution
    IEI_dist = IEI_dist(1:end-1)./sum(IEI_dist);
    
    % Store IEI probability distribution across simulations
    All_Counts_IEI(model_Iteration, :) = IEI_dist;
end

% Average IEI probability distribution across simulations
if num_Iterations > 1
    IEI_dist = mean(All_Counts_IEI);
end

end