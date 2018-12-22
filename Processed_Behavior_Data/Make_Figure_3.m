% Make_Figure_3.m
% by David Ehrlich, April 14, 2016

% To generate Figure 3 from "Control of Movement Initiation Underlies the Development of Balance" 
% Ehrlich DE and Schoppik D, Curr Biol. 2017 Feb 6;27(3):334-344.

% Be sure to first run preprocessing scripts.
% Run this in the provided 'Raw_Behavior_Data' directory containing one folder corresponding to
% each group.


hhh = figure; %initialize figure

%iterate through both conditions
d = dir;
isub = [d(:).isdir]; %# returns logical vector
groupFolders = {d(isub).name}'; 

for condition=[5 6]
    cd(groupFolders{condition+2})
    
    clear proportion_nose_up pause_ang_vel IEIs_by_clutch
    all_IEIs = [];

    % Iterate through each subfolder to get these variables for each group
    d = dir;
    isub = [d(:).isdir]; %# returns logical vector
    clutchFolders = {d(isub).name}';
    for clutch = 1:length(clutchFolders)-2
        cd(clutchFolders{clutch+2})    

        mat_name = dir('*.mat');
        load(mat_name.name);

% % % %         load('PropBoutPeakAngVel.mat')
% % % %         load('PropBoutTime.mat')
        %analyze only bouts in light
        [day,~] = lightdarksplit(PropBoutTime);
        %find proportion of bouts with nose-up peak angular velocity 
        proportion_nose_up(clutch) = sum(PropBoutPeakAngVel(day)>0)./length(day);

        %average angular velocities between bouts
% % % % %         load('PropBoutIEIangVel.mat')
% % % % %         load('PropBoutIEItime.mat')
        [day,~]=lightdarksplit(PropBoutIEItime);
        pause_ang_vel(clutch) = nanmean(PropBoutIEIangVel(day)); 
        
        %load inter-event intervals
% % % %         load('PropBoutIEI.mat')
% % % %         load('PropBoutIEItime.mat')
        %analyze only events in light
        [day,~] = lightdarksplit(PropBoutIEItime);
        %concatenate
        all_IEIs(end+1:end+length(day)) = PropBoutIEI(day);
        IEIs_by_clutch(clutch) = mean(PropBoutIEI(day));
        
% % % %         load('BodyAngles.mat')
% % % %         load('GrabbedTimes.mat')
        [day,~]=lightdarksplit(GrabbedTimes);
        pitch_by_clutch(clutch) = mean(BodyAngles(day));
        
        cd ..
    end
    
    if condition==5 %controls
        color='k';
    else %oil
        color='r';    
    end
    
    %Plot percent of bouts nose-up vs angular velocity during pauses for
    %each group
    subplot(1,2,1); hold on
    plot(pause_ang_vel, proportion_nose_up, ['.' color])
    axis square; axis([-50 10 0.3 0.9])
    ylabel('% bouts nose-up')
    xlabel('Inter-bout angular velocity (deg/sec)')
    
    %Calculate IEI histogram
    edges = 0:0.1:20;
    centers = edges(1:end-1)+diff(edges(1:2))/2;
    counts = histc(all_IEIs,edges);
    
    %Plot figure 3E
    subplot(1,2,2);
    semilogy(centers, counts(1:end-1)./sum(counts), color)
    hold on;
    axis square;
    axis([0 3 0.0003  0.3])
    xlabel('IEI (sec)')
    ylabel('Probability')

    cd ..
end

legend('air swimbladder','oil swimbladder')

