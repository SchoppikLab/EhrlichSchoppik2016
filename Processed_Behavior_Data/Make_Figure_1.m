% Make_Figure_1.m
% by David Ehrlich, July 23, 2016

% To generate Figures 1D and 1E from "Control of Movement Initiation Underlies the Development of Balance" 
% Ehrlich DE and Schoppik D, Curr Biol. 2017 Feb 6;27(3):334-344.

% Be sure to first run preprocessing scripts
% Run this in the provided 'Raw_Behavior_Data' folder containing one folder corresponding to
% each age.

clear all

%initialize figure
hhh = figure;
set(hhh, 'Position',[50 50 1200 600])

%list of ages for plotting
Ages = [4 7 14 21]; %days post-fertilization

%iterate through each time-point folder
d = dir;
isub = [d(:).isdir]; %# returns logical vector
ageFolders = {d(isub).name}';

for age=1:4
    cd(ageFolders{age+2})  
    
    d = dir;
    isub = [d(:).isdir]; % returns logical vector
    clutchFolders = {d(isub).name}';
    
    %iterate through each clutch folder
    for clutch=1:4
        cd(clutchFolders{clutch+2}) 
        
        mat_name = dir('*.mat');
        load(mat_name.name);
    
% % % %         load('PropBoutPeakAngVel.mat')
% % % %         load('PropBoutTime.mat')
        %analyze only bouts during light phase
        [day,~]=lightdarksplit(PropBoutTime);
        %find proportion of bouts with nose-up peak angular velocity 
        proportion_nose_up(age,clutch) = sum(PropBoutPeakAngVel(day)>0)./length(day);
        
        %load instantaneous speeds, angular velocities, and times of day
% % % %         load('HeadingMatchedSpeeds.mat')
% % % %         load('HeadingMatchedTimes.mat')
% % % %         load('HeadingMatchedAngVels.mat')
        %analyze only bouts during light phase
        [day,~]=lightdarksplit(HeadingMatchedTimes);
        %translation speeds and angular velocities during light phase
        dayspd=HeadingMatchedSpeeds(day);
        dayvel=HeadingMatchedAngVels(day);
        %compute mean angular velocity when speed is below 1 mm/sec
        pause_ang_vel(age,clutch) = mean(dayvel(dayspd<1));                      
        
        cd ..
    end
        
    cd ..
end

%Plot results for each clutch
for clutch=1:4
    %Plot % of bouts nose-up
    subplot(1,4,1); hold on;
    plot(Ages, proportion_nose_up(:,clutch)*100, 'g-')
    %Plot average angular velocity during pauses
    subplot(1,4,2); hold on;
    plot(Ages, pause_ang_vel(:,clutch), '-', 'Color',[1 .5 0]);
    %Plot paired % of nose-up bouts and angular velocity during pauses by
    %clutch
    subplot(1,4,3); hold on
    plot(0,mean(proportion_nose_up(:,clutch))*100,'.','MarkerSize',15)
    
    subplot(1,4,4); hold on
    plot(0,mean(pause_ang_vel(:,clutch)),'.','MarkerSize',15)
end
%Plot average across clutches
subplot(1,4,1); hold on;
plot(Ages, mean(proportion_nose_up')*100, 'g-','LineWidth',5)
title('% bouts nose-up vs age by clutch')
ylabel('% bouts nose-up')
xlabel('Age (dpf)')
axis([0 24 7 85])

subplot(1,4,2); hold on;
plot(Ages, mean(pause_ang_vel'), '-', 'Color',[1 .5 0], 'LineWidth',5)
title('Angular velocity during pauses vs age by clutch')
ylabel('Angular velocity (deg/sec)')
xlabel('Age (dpf)')
axis([0 24 -5.2 4.2])

subplot(1,4,3); hold on;
plot(0,mean(mean(proportion_nose_up))*100,'.','MarkerSize',15)
title('% Bouts nose up by clutch')
ylabel('% bouts nose-up')
legend('Clutch 1','C2','C3','C4','Mean')
axis([-1 1 7 85])

subplot(1,4,4); hold on;
plot(0,mean(mean(pause_ang_vel)),'.','MarkerSize',15)
title('Angular velocity during pauses by clutch')
ylabel('Angular velocity (deg/sec)')
legend('Clutch 1','C2','C3','C4','Mean')
axis([-1 1 -5.2 4.2])





    


