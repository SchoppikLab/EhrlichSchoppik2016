% Fit_bout_likelihood.m
% by David Ehrlich, April 13, 2016

% Associated with "Control of Movement Initiation Underlies the Development of Balance" 
% Ehrlich DE and Schoppik D, Curr Biol. 2017 Feb 6;27(3):334-344.

% For fitting function to empirical Relative Bout Likelihood vs. binned by pitch and angular velocity

% Inputs and outputs
% coeffs contains values of coefficients for the relative bout likelihood function
% Quantile_Edges contains edges for pitch and angular velocity bins.
% binned_Output is a vector containing Relative Bout Likelihood evaluated for each of the bins.

function binned_Output = Fit_bout_likelihood(coeffs, Quantile_Edges)
alpha = coeffs(1);
zeta = coeffs(2);
beta = coeffs(3);
%Gamma is fixed
gamma = 0.035;

% bin edges taken from input
PitchEdges = Quantile_Edges(1,:);
AngVelEdges = Quantile_Edges(2,:);

% number of bins across pitch and angular velocity
numBins = length(Quantile_Edges)-1;

%% Evaluate Relative Bout Likelihood for every pair of observed pitch and angular velocity during an IEI
wdname = pwd;
load(['Group_' wdname(end-4:end) '.mat'])
% % % % load(mat_name.name);
% % % % load('PropBoutIEIpitch.mat')
% % % % load('PropBoutIEIangVel.mat')
% % % % load('PropBoutIEItime.mat')

% isolate frames during light phase
[day,~]=lightdarksplit(PropBoutIEItime);
Pitches = PropBoutIEIpitch(day);
AngVels = PropBoutIEIangVel(day);
% remove NaNs
Pitches(isnan(AngVels)) = [];
AngVels(isnan(AngVels)) = [];

% Evaluate bout likelihood at each observed pair of pitch and angular
% velocity during IEIs.
for i=1:length(Pitches)
    %zeta (baseline coefficient) is restricted to be positive
    P_bout(i) = abs(zeta) + alpha*(abs(Pitches(i))) + beta*(abs(AngVels(i))) - gamma*((AngVels(i)));
end

% average bout likelihood for each posture bin
for PitchInd = 1:length(PitchEdges)-1
    %Find indices of observed IEIs occupying this pitch bin
    PitchHit = intersect(find(Pitches>=PitchEdges(PitchInd)),find(Pitches<PitchEdges(PitchInd+1)));
    for AVind = 1:length(AngVelEdges)-1
        %Find indices of observed IEIs occupying this angular velocity bin
        AVhit = intersect(find(AngVels>=AngVelEdges(AVind)),find(AngVels<AngVelEdges(AVind+1)));
        % average Pbout values for every pair of pitch and angular velocity
        % occupying this bin
        bin_P_bout(PitchInd, AVind) = mean(P_bout(intersect(PitchHit,AVhit)));
    end
end

%Make output vector from binned bout likelihoods
binned_Output = bin_P_bout(:);

end

