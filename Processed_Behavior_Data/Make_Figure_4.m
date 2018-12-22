% Generate_Figure_4.m
% by David Ehrlich, April 14, 2016

% To generate Figure 4 from "Control of Movement Initiation Underlies the Development of Balance" 
% Ehrlich DE and Schoppik D, Curr Biol. 2017 Feb 6;27(3):334-344.

% Be sure to first run preprocessing scripts.
% Run this in the provided 'Raw_Behavior_Data' folder containing one subfolder 
% corresponding to each group.

clear all
%initialize figure
hhh = figure;
set(hhh, 'Position',[50 50 1200 600])

CI{1} = 0; %initialize confidence interval array

%list of ages for plotting
Ages = [4 7 14 21]; %dpf

%iterate through each time-point folder
d = dir;
isub = [d(:).isdir]; %# returns logical vector
ageFolders = {d(isub).name}';

for ageFolder=1:4
    cd(ageFolders{ageFolder+2})  
    
    %% Figure 3A
% % %     mat_name = dir('*.mat');
% % %     load(mat_name.name);
    wdname = pwd;
    load(['Group_' wdname(end-4:end) '.mat'])

    load('PitchEdges');
    load('AVedges');
    load('BoutProbability');
    load('CoeffsBinFit.mat')
    load('RbinFit.mat')
    load('CovBbinFit.mat')

    %% Plot empirical Relative Bout Likelihood across the pitch and angular velocity bins
    subplot(3,12,[3*(ageFolder-1)+1 3*(ageFolder-1)+2 3*(ageFolder-1)+3])
    hhhh = surf(PitchEdges,AVedges,zeros(length(PitchEdges)), BoutProbability');
    set(hhhh, 'edgecolor', 'none')
    xlabel('Pitch (deg)'); ylabel('Angular Velocity (deg/sec)');
    colorbar; colormap(parula(256)); view(0,90);
    axis([-15 30 -10 10]);
    title([ageFolders{ageFolder+2} ' Obs. RBL vs Posture']);
    caxis([0 4]);
    axis square

    %% Plot Relative Bout Likelihood function evaluated for a range of pitches and angular velocities
    subplot(3,12,[3*(ageFolder-1)+13 3*(ageFolder-1)+14 3*(ageFolder-1)+15])
    clear Pbout;
    % generate linearly spaced pairs of pitch and angular velocity to
    % evaluate Relative Bout Likelihood function for plotting
    Pitches = linspace(PitchEdges(1), PitchEdges(end), (length(PitchEdges)-1) * 20);
    AngVels = linspace(AVedges(1), AVedges(end), (length(PitchEdges)-1) * 250);
    for i=1:length(Pitches)
        for j=1:length(AngVels)
            Pbout(i,j) = coeffs(2) + coeffs(1)*abs(Pitches(i)) + coeffs(3)*abs(AngVels(j)) - coeffs(4)*(AngVels(j));
        end
    end
    h = surf(Pitches - (Pitches(2)-Pitches(1)), AngVels - (AngVels(2)-AngVels(1)), zeros(length(Pitches),length(AngVels))', Pbout');
    set(h, 'edgecolor', 'none')
    xlabel('Pitch (deg)'); ylabel('Angular Vel. (deg/sec)'); 
    title('Continuously eval. RBL');
    colorbar; colormap(parula(256)); view(0,90);
    axis([-15 30 -10 10]);
    caxis([0 4]);
    axis square
    
    %% Load parameter solutions for plotting
% % % %     load('CoeffsBinFit.mat')
% % % %     load('RbinFit.mat')
% % % %     load('CovBbinFit.mat')
    
    % Returns 99% confidence intervals to correct for multiple (6) comparisons across age
    CI{ageFolder} = nlparci(coeffs(1:3),R,'covar', CovB, 'alpha', 0.0083);
    % save alpha solution and error for plotting confidence intervals
    alpha(ageFolder) = coeffs(1);
    alpha_error(ageFolder) = coeffs(1) - CI{ageFolder}(1,1);
    % save zeta solution and error for plotting confidence intervals
    zeta(ageFolder) = coeffs(2);
    zeta_error(ageFolder) = coeffs(2) - CI{ageFolder}(2,1);
    % save beta solution and error for plotting confidence intervals 
    beta(ageFolder) = coeffs(3);
    beta_error(ageFolder) = coeffs(3) - CI{ageFolder}(3,1);
        
    %% Calculate R^2 of fit of estimated to observed Relative Bout Likelihood
    FitPbouts = Fit_bout_likelihood(coeffs, vertcat(PitchEdges, AVedges));
    tmp = corrcoef(BoutProbability(:), FitPbouts);
    Relative_Bout_Likelihood_Fit(ageFolder) = tmp(1,2)^2;
    
    %exit subfolder
    cd .. 
end
    
%% Figure 4B: plot alpha vs. age
subplot(3,12,[25 26 27 28]); hold on
errorbar(Ages, alpha, alpha_error)
ylabel('alpha (deg^-^1)')
xlabel('age (dpf)')
title('Pitch sensitivity (best & 99CI)')
axis([0 22 -0.02 0.12])
axis square

%% Figure 4C: plot beta vs. age
subplot(3,12,[29 30 31 32]); hold on
errorbar(Ages, beta, beta_error)
ylabel('beta (sec x deg^-^1)')
xlabel('age (dpf)')
title('Angular vel sensitivity')
axis([0 22 0 0.25])
axis square

%% Figure 4D: plot zeta vs. age
subplot(3,12,[33 34 35 36]); hold on
errorbar(Ages, zeta, zeta_error)
ylabel('zeta (dimensionless)')
xlabel('age (dpf)')
title('Baseline')
axis([0 22 0 1])
axis square

