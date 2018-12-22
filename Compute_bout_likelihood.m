% Compute_bout_likelihood.m
% by David Ehrlich, Jan 21, 2016

% Associated with "Control of Movement Initiation Underlies the Development of Balance" 
% Ehrlich DE and Schoppik D, Curr Biol. 2017 Feb 6;27(3):334-344.

% Creates a lookup table for relative bout likelihood as a function of pitch and angular velocity

clear all

%% Some parameters
numBins = 8; % number of quantiles for splitting pitch and angular velocity probability distributions
samplerate = 40; %Hz

%% Some switches
% if plotsOn equals 1, diagnostic plots will be produced
plotsOn = 1;
% if saveOn equals 1, output variables will be saved
saveOn = 1;

Dev_Folders = {'04dpf','07dpf','14dpf','21dpf'};
for dev_group = 1:4
    cd(Dev_Folders{dev_group})  

    clear ProbPreBout ProbAllObs BoutProbability
    
    %% Find quantiles of pitch and angular velocity during IEIs 
    % load inter-event interval variables
    % % % % mat_name = dir('*.mat');
    % % % % load(mat_name.name);
    wdname = pwd;
    load(['Group_' wdname(end-4:end) '.mat'])

    % find frames during light phase
    [day,~]=lightdarksplit(PropBoutIEItime);
    % [~,day]=lightdarksplit(PropBoutIEItime);

    % initialize variables to append pitch and angular velocity during IEI
    AllPitches = []; AllAngVels = []; 

    % iterate through each cell of the arrays for pitch and angular velocity during IEI, ignoring the first (empty) 
    for IEIcell = 1: (length(PropBoutIEIalignedAngVel) - 1)
        % if IEI occurred during light phase
        if ~isempty(intersect(day,IEIcell))
            % append IEI pitch and angular velocity to variables
            AllPitches( end+1: end+length(PropBoutIEIalignedPitch{IEIcell+1}) - 1) = PropBoutIEIalignedPitch{IEIcell+1}(1:end-1);       
    %         AllPitches( end+1: end+length(PropBoutIEIalignedAngVel{IEIcell+1}) - 1) = PropBoutIEIalignedPitch{IEIcell+1}(1:end-1);       
            AllAngVels( end+1: end+length(PropBoutIEIalignedAngVel{IEIcell+1}) - 1) = PropBoutIEIalignedAngVel{IEIcell+1}(1:end-1);
        end
    end

    % calculate quantiles and define bin edges
    tmp=100/numBins*(1:numBins);
    PitchEdges = prctile(AllPitches, [0 tmp]);
    AVedges = prctile(AllAngVels, [0 tmp]);

    %% Find Observation probability for pre-bout values
    % load bout-aligned posture variables
    % % % % load('PropBoutAlignedPitch.mat')
    % % % % load('PropBoutAlignedAngVel.mat')
    % % % % load('PropBoutAlignedTime.mat')
    % find frames during light phase
    [day,~]=lightdarksplit(PropBoutAlignedTime);

    % find mean pre-bout pitch and angular velocity 200 msec before fastest translation (column 31)
    prePitch = mean(PropBoutAlignedPitch(day,22:24),2);
    preAngVel = mean(PropBoutAlignedAngVel(day,22:24),2);

    % iterate through each quantile of pitch and angular velocity to find which
    % elements belong to which bin to generate a 2-d probability distribution
    for PitchInd = 1:length(PitchEdges)-1
        PitchHit = intersect(find(prePitch>PitchEdges(PitchInd)),find(prePitch<PitchEdges(PitchInd+1)));
        for AVind = 1:length(AVedges)-1
            AVhit = intersect(find(preAngVel>AVedges(AVind)),find(preAngVel<AVedges(AVind+1)));
            %probability distribution (sums to 1) of pre-bout posture occupying a given bin
            ProbPreBout(PitchInd, AVind) = length(intersect(PitchHit,AVhit)) / length(prePitch);
        end
    end

    %% Find Observation probability for all observed pitches and angular velocities (excluding refractory period)
    PostRefractPitches = []; PostRefractAngVels = []; 
    % iterate through each IEI and append pitches and angular velocities after refractory pd
    for IEIcell=1:length(PropBoutIEIalignedAngVel) - 1
        if ~isempty(intersect(day,IEIcell))
            % PropBoutIEIalignedAngVel and Pitch already exclude 50 msec after the
            % frame in which speed dips below 5 mm/sec.  After a subsequent 100 msec refractory period,
            % each pitch and angular velocity is appended.
            try PostRefractAngVels( end+1: end+length(PropBoutIEIalignedAngVel{IEIcell+1}) - 4) = PropBoutIEIalignedAngVel{IEIcell+1}(5:end); end
            try PostRefractPitches( end+1: end+length(PropBoutIEIalignedAngVel{IEIcell+1}) - 4) = PropBoutIEIalignedPitch{IEIcell+1}(5:end); end        
        end
    end
    % generate probability distributions of pitch and angular velocity across bins
    for PitchInd = 1:length(PitchEdges)-1
        PitchHit = intersect(find(PostRefractPitches>=PitchEdges(PitchInd)),find(PostRefractPitches<PitchEdges(PitchInd+1)));
        for AVind = 1:length(AVedges)-1
            AVhit = intersect(find(PostRefractAngVels>=AVedges(AVind)),find(PostRefractAngVels<AVedges(AVind+1)));
            ProbAllObs(PitchInd, AVind) = length(intersect(PitchHit,AVhit)) / length(PostRefractPitches);
        end
    end

    %% Calculate Bayesian posterior
    for PitchInd = 1:length(PitchEdges)-1
        for AVind = 1:length(AVedges)-1                               
            BoutProbability(PitchInd, AVind) = ProbPreBout(PitchInd, AVind) / ProbAllObs(PitchInd, AVind);
        end
    end

    %% Plot output
    if plotsOn == 1   
        figure;
        %Plot probability of observing a given pitch and angular velocity
        %at bout initiation
        subplot(2,1,1);
        surf(PitchEdges,AVedges,zeros(length(PitchEdges)), BoutProbability')
        xlabel('Pitch (deg)'); ylabel('Angular Velocity (deg/sec)'); zlabel('P(bout)'); title('P(pitch,AV | bout)');
        colorbar; colormap('parula'); view(0,90);
        axis([-20 30 -10 10])
        caxis([0 4])
        
        %Plot joint probability distribution of pitch and angular velocity
        subplot(2,1,2);
        surf(PitchEdges,AVedges,zeros(length(PitchEdges)), ProbAllObs')
        xlabel('Pitch (deg)'); ylabel('Angular Velocity (deg/sec)'); zlabel('P(obs)'); title('P(pitch,AV)');
        colorbar; colormap('parula'); caxis([0 0.02]); view(0,90);
        axis([-20 30 -10 10])
        caxis([0 4])

    end

    if saveOn == 1
        save('PitchEdges','PitchEdges');
        save('AVedges','AVedges');
        save('BoutProbability','BoutProbability')
    end
    
    cd ..
end