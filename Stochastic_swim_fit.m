% Stochastic_swim_fit.m
% by David Ehrlich, Jan 18, 2016

% Associated with "Control of Movement Initiation Underlies the Development of Balance" 
% Ehrlich DE and Schoppik D, Curr Biol. 2017 Feb 6;27(3):334-344.

% Run this in an age-specific subfolder of the 'Raw_Behavior_Data' directory 
% to fit empirical relative bout likelihood and then solve 
% for the coefficient of the Bayesion prior that best matches the empirical 
% distribution of inter-bout intervals. 
% Please note: the solver in this script may take a while to run.

% First run 'Compute_bout_likelihood.m' to generate necessary 'BoutProbability', 'PitchEdges', and 'AVedges' variables.

clear all

%% Some switches
LoadCoeffs = 0; % if 1, Relative bout likelihood function coeffs are loaded instead of solved using nlinfit
plotsOn = 1; % set to 1 to plot output

%% Some parameters
samplerate = 40; %Hz
numBins = 25; % How many bins for comparing IEI histograms when evaluating simulations?

%% Fit relative bout likelihood
load('PitchEdges.mat')
load('AVedges.mat')
load('BoutProbability.mat')
% this variable conveys quantiles to FitRelativeLikelihood
InputSets = vertcat(PitchEdges, AVedges);

% Initial values for parameters of RelativeBoutLikelihood function
coeffInit = [.01  0.5  0.01];
stepSize = [.002  0.1   0.002];

% fitting options
opts = statset('nlinfit');
% opts.Display = 'iter'; %uncomment to display iterative fit improvement
opts.DerivStep = stepSize;

if LoadCoeffs == 1
    % load solutions of previous fit to relative bout probability
    load('CoeffsBinFit.mat')
    load('RbinFit.mat')
    load('CovBbinFit.mat')
else
    disp('Fitting relative bout likelihood')
    % run @nlinfit, using BoutProbability as target data (Y) 
    [coeffs,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(InputSets, BoutProbability(:), @Fit_bout_likelihood, coeffInit, opts);
    % Gamma is fixed at 0.035
    coeffs(4) = 0.035;
    % save coefficient solutions
    save('CoeffsBinFit','coeffs')
    % save additional variables for confidence interval calculation
    save('RbinFit','R')
    save('CovBbinFit','CovB')
end

%% Identify best psi in Bayesian prior based on above solutions to relative bout likelihood
% Calculate target IEI histogram for fitting
edges = linspace(0,20,numBins+1);

% % % % mat_name = dir('*.mat');
% % % % load(mat_name.name);
wdname = pwd;
load(['Group_' wdname(end-4:end) '.mat'])
% % % % % load('PropBoutIEI.mat')
% % % % % load('PropBoutIEItime.mat')

% find frames during light phase
[day,~]=lightdarksplit(PropBoutIEItime);
targetIEIdist = histc(PropBoutIEI(day),edges);
targetIEIdist = targetIEIdist(1:end-1)./sum(targetIEIdist); %make probability distribution

coeffInit = 0.04; %initial psi guess    
stepSize =  0.01;   

% fitting options
opts = statset('nlinfit');
opts.Display = 'iter'; %'final'; %uncomment to display iterative fit improvement
opts.DerivStep = stepSize;

disp('Matching IEI distribution')
%run @nlinfit, using targetIEIdist as target data (Y) 
[coeffsIEI,Riei,Jiei,CovBiei,MSEiei,ErrorModelInfoIEI] = nlinfit(coeffs, targetIEIdist, @Simulate_Stochastic_Swim, coeffInit, opts);
save('PsiFit','coeffsIEI')
save('PsiR','Riei')
save('PsiCovB','CovBiei')

if plotsOn == 1
    fh = figure;
    set(fh, 'Position', [100, 100, 1049, 195]);  

    % Plot observed Relative Bout Likelihood 
    subplot(1,4,1); hold on;
    % plot surface with color defined by observed relative likelihood
    surf(PitchEdges,AVedges,zeros(length(PitchEdges)), BoutProbability')
    xlabel('Pitch (deg)'); ylabel('Angular Velocity (deg/sec)');
    colorbar; colormap(parula); view(0,90);
    axis([-20 40 -10 10]);
    title('Observed Relative Likelihood by Posture');
    caxis([0.3 5]);
    caxisSP1 = caxis;

    % Plot Relative Bout Likelihood for same bins as plot 1 but calculated using continuous function evaluated in each bin     
    subplot(1,4,2); hold on;
    % evaluate Relative Likelihood function using parameter solutions across each posture bin 
    BinPboutEval = Fit_bout_likelihood(coeffs, InputSets);
    for i=1: length(PitchEdges) - 1
        for j=1: length(PitchEdges) - 1
            % convert linear function output back to matrix for plotting
            BinPboutEvalMatrix(i,j) = BinPboutEval((j-1)*(length(PitchEdges)-1) + i);
        end
    end
    % plot surface with color defined by evaluated relative likelihood
    surf(PitchEdges, AVedges, zeros(length(PitchEdges)), BinPboutEvalMatrix');
    xlabel('Pitch (deg)'); ylabel('Angular Velocity (deg/sec)'); 
    title('Evaluated Relative Likelihood by Posture');
    colorbar; colormap(parula); view(0,90);
    axis([-20 40 -10 10]);
    caxis(caxisSP1);

    % plot continuous Relative Bout Likelihood function
    subplot(1,4,3); hold on;
    clear Pbout;
    % generate linearly spaced pairs of pitch and angular velocity to evaluate Relative Bout Likelihood function
    Pitches = linspace(PitchEdges(1), PitchEdges(end), (length(PitchEdges)-1) * 20);
    AngVels = linspace(AVedges(1), AVedges(end), (length(PitchEdges)-1) * 250);
    for i=1:length(Pitches)
        for j=1:length(AngVels)
            Pbout(i,j) = coeffs(2) + coeffs(1)*abs(Pitches(i)) + coeffs(3)*abs(AngVels(j)) - 0.035*(AngVels(j));
        end
    end
    h = surf(Pitches - (Pitches(2)-Pitches(1)), AngVels - (AngVels(2)-AngVels(1)), zeros(length(Pitches),length(AngVels))', Pbout');
    set(h, 'edgecolor', 'none')
    xlabel('Pitch (deg)'); ylabel('Angular Velocity (deg/sec)'); 
    title('Continuously evaluated Relative Likelihood by Posture');
    colorbar; colormap(parula); view(0,90);
    axis([-20 40 -10 10]);
    caxis(caxisSP1);

    % Plot evaluated vs. observed Relative Bout Likelihood for each bin
    subplot(1,4,4); hold on;
    plot(BoutProbability(:), BinPboutEval, 'b.', 'MarkerSize', 5);
    xlabel('Observed Relative Likelihood'); ylabel('Evaluated Relative Likelihood');
    sub4axis = axis;
    plot([min(sub4axis) max(sub4axis)],[min(sub4axis) max(sub4axis)],'k-')   
    mTextBox = uicontrol('style','text');
    set(mTextBox,'String',['alpha= ' num2str(coeffs(1)) ', zeta= ' num2str(coeffs(2)) ', beta= ' num2str(coeffs(3)) ', gamma= ' num2str(coeffs(4))])
    set(mTextBox,'Units','characters')
    set(mTextBox,'Position',[138 12 25 4])
end  
 
% Calculate R^2 of observed vs. evaluated relative bout likelihood 
tmp = corrcoef(BoutProbability(:), BinPboutEval);
Relative_Bout_Likelihood_Fit = tmp(1,2)^2




