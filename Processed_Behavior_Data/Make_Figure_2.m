% Make_Figure_2.m
% by David Ehrlich, July 25, 2016

% To generate Figure 2 from "Control of Movement Initiation Underlies the Development of Balance" 
% Ehrlich DE and Schoppik D, Curr Biol. 2017 Feb 6;27(3):334-344.

% Be sure to first run preprocessing scripts
% Run this in the provided 'Raw_Behavior_Data' folder containing one folder corresponding to
% each age.

clear nose_down_gain nose_up_gain

%initialize figure
hhh = figure;
set(hhh, 'Position',[0 0 1200 700])

%list of ages for plotting
Ages = [4 7 14 21]; %dpf

%iterate through each time-point folder
d = dir;
isub = [d(:).isdir];
ageFolders = {d(isub).name}'; %name of folders in directory

for age=1:4
    cd(ageFolders{age+2}) %switch to age folder 
    %get names of folders for each clutch at this age
    d = dir;
    isub = [d(:).isdir]; %# returns logical vector
    clutchFolders = {d(isub).name}';
    %iterate through each clutch subfolder
    for clutch=1:4
        cd(clutchFolders{clutch+2}) %switch to clutch subfolder
        %load matrices for the time, translation speed, and angular velocity of each
        %recorded swim bout, aligned to peak translation speed. Each bout
        %occupies a row and each column corresponds to a frame, aligned with peak speed
        %at sample 31.

        mat_name = dir('*.mat');
        load(mat_name.name);
% % % % % 
% % % % %         load('PropBoutAlignedAngVel.mat')
% % % % %         load('PropBoutAlignedSpeed.mat')
% % % % %         load('PropBoutAlignedTime.mat')
        %Use 'lightdarksplit' to identify bouts that occurred during the
        %light phase of the light-dark cycle.
        [day,~]=lightdarksplit(PropBoutAlignedTime);
        %make variable for angular velocity of bouts that occurred during the day
        day_ang_vels = PropBoutAlignedAngVel(day,:);
        %calculate pre-bout angular velocity for each bout
        pre_bout_ang_vels = mean(day_ang_vels(:,22:26),2);
        %calculate post-bout angular velocity for each bout
        post_bout_ang_vels = mean(day_ang_vels(:,41:46),2);
        %calculate net change in angular velocity across each bout
        net_bout_angular_acceleration = post_bout_ang_vels - pre_bout_ang_vels;
        %calculate slope and intercept for best-fit line of net change in angular velocity vs pre-bout
        %angular velocity, separately for nose-down and nose-up pre-bout
        %angular velocities.
        nose_down_coeffs = polyfit(pre_bout_ang_vels(pre_bout_ang_vels <= 0), net_bout_angular_acceleration(pre_bout_ang_vels <= 0), 1);
        nose_up_coeffs = polyfit(pre_bout_ang_vels(pre_bout_ang_vels > 0), net_bout_angular_acceleration(pre_bout_ang_vels > 0), 1);
        %calculate gain of angular velocity correction as negative slope of best-fit lines,
        %separate for nose-up and nose-down pre-bout angular velocity 
        nose_down_gain(age,clutch) = -nose_down_coeffs(1);
        nose_up_gain(age,clutch) = -nose_up_coeffs(1);
        
        %% Figure 2E, individual clutch data
        %after the last age is processed, plot gains of angular velocity correction as a function of age
        %and clutch for nose-up and nose-down pre-bout angular velocities
        if age==4
            subplot(5,6,[23 29]); hold on;
            plot(Ages,nose_down_gain(:,clutch),'Color',[.5 .5 .5])
            subplot(5,6,[24 30]); hold on;
            plot(Ages,nose_up_gain(:,clutch),'Color',[.5 .5 .5])
        end
        %exit clutch subfolder
        cd ..
    end
       
    %load matrices, collected for all clutches at a given age, 
    %for the time, translation speed, and angular velocity of each
    %recorded swim bout. Bouts are aligned to peak translation speed. Each bout
    %occupies a row and each column corresponds to a frame, aligned with peak speed
    %at sample 31.

% % % %     mat_name = dir('*.mat');
% % % %     load(mat_name.name);
    wdname = pwd;
    load(['Group_' wdname(end-4:end) '.mat'])
       
% % % % 
% % % %     load('PropBoutAlignedAngVel.mat')
% % % %     load('PropBoutAlignedSpeed.mat')
% % % %     load('PropBoutAlignedTime.mat')
    %Use 'lightdarksplit' to identify bouts that occurred during the
    %light phase of the light-dark cycle.
    [day,~]=lightdarksplit(PropBoutAlignedTime);
    day_ang_vels = PropBoutAlignedAngVel(day,:);
    
    %calculate pre-bout angular velocity for each bout
    pre_bout_ang_vels = mean(day_ang_vels(:,22:26),2);
    %calculate post-bout angular velocity for each bout
    post_bout_ang_vels = mean(day_ang_vels(:,41:46),2);
    %calculate net change in angular velocity across each bout
    net_bout_angular_acceleration = post_bout_ang_vels - pre_bout_ang_vels;

    %for grouped data, calculate slope and intercept for best-fit line of net change in angular velocity vs pre-bout
    %angular velocity, separately for nose-down and nose-up pre-bout
    %angular velocities.
    nose_down_coeffs = polyfit(pre_bout_ang_vels(pre_bout_ang_vels <= 0), net_bout_angular_acceleration(pre_bout_ang_vels <= 0), 1);
    nose_up_coeffs = polyfit(pre_bout_ang_vels(pre_bout_ang_vels > 0), net_bout_angular_acceleration(pre_bout_ang_vels > 0), 1);
    %for grouped data, calculate gain of angular velocity correction as negative slope of best-fit lines,
    %separate for nose-up and nose-down pre-bout angular velocity 
    nose_down_gain(age,5) = -nose_down_coeffs(1);
    nose_up_gain(age,5) = -nose_up_coeffs(1);
    
    %make representative plots of data at 7 and 21 dpf
    if age==2
        %% Plot Figure 2C
        %for all individual bouts at 7dpf with nose-down pre-bout angular velocity,
        %plot net change in angular velocity vs. pre-bout angular velocity.
        subplot(5,6,[19 25]); hold on;
        plot(pre_bout_ang_vels(pre_bout_ang_vels <= 0), net_bout_angular_acceleration(pre_bout_ang_vels <= 0), '.', 'Color',[.5 .5 .5])
        %superimpose average of evenly-spaced bins calculated with 'makeEvenHistogram' function
        [binned_pre, binned_net]= makeEvenHistogram(pre_bout_ang_vels, net_bout_angular_acceleration,100);
        plot(binned_pre, binned_net,'k-')
        axis([-20 0 -20 20]);    
        ylabel('Net change in angular velocity (deg/sec)')
        xlabel('Pre-bout ang. vel. (deg/sec)')
        title('7 dpf, nose-down')
        %plot best-fit line to nose-down data
        plot([-20 0], polyval(nose_down_coeffs,[-20 0]), '-b', 'LineWidth',3);

        %for all individual bouts at 7dpf with nose-up pre-bout angular velocity,
        %plot net change in angular velocity vs. pre-bout angular velocity.
        subplot(5,6,[20 26]); hold on;
        plot(pre_bout_ang_vels(pre_bout_ang_vels > 0), net_bout_angular_acceleration(pre_bout_ang_vels > 0), '.', 'Color',[.5 .5 .5])
        %superimpose average of evenly-spaced bins calculated with 'makeEvenHistogram' function
        [binned_pre, binned_net]= makeEvenHistogram(pre_bout_ang_vels, net_bout_angular_acceleration,100);
        plot(binned_pre, binned_net,'k-');
        axis([0 20 -20 20]);    
        ylabel('Net change in angular velocity (deg/sec)')
        xlabel('Pre-bout ang. vel. (deg/sec)')
        title('7 dpf, nose-up')
        %plot best-fit line to nose-up data
        plot([0 20], polyval(nose_up_coeffs,[0 20]), '-b', 'LineWidth',3);

    elseif age==4
        %% Plot Figure 2D
        %for all individual bouts at 21dpf with nose-down pre-bout angular velocity,
        %plot net change in angular velocity vs. pre-bout angular velocity.
        subplot(5,6,[21 27]); hold on;
        plot(pre_bout_ang_vels(pre_bout_ang_vels <= 0), net_bout_angular_acceleration(pre_bout_ang_vels <= 0), '.', 'Color',[.5 .5 .5])
        %superimpose average of evenly-spaced bins calculated with 'makeEvenHistogram' function
        [binned_pre, binned_net]= makeEvenHistogram(pre_bout_ang_vels, net_bout_angular_acceleration,100);
        plot(binned_pre, binned_net,'k-');
        axis([-20 0 -20 20]);    
        ylabel('Net change in angular velocity (deg/sec)')
        xlabel('Pre-bout ang. vel. (deg/sec)')
        title('21 dpf, nose-down')
        %plot best-fit line to nose-down data
        plot([-20 0], polyval(nose_down_coeffs,[-20 0]), '-m', 'LineWidth',3);

        %for all individual bouts at 21dpf with nose-up pre-bout angular velocity,
        %plot net change in angular velocity vs. pre-bout angular velocity.
        subplot(5,6,[22 28]); hold on;
        plot(pre_bout_ang_vels(pre_bout_ang_vels > 0), net_bout_angular_acceleration(pre_bout_ang_vels > 0), '.', 'Color',[.5 .5 .5])
        %superimpose average of evenly-spaced bins calculated with 'makeEvenHistogram' function
        [binned_pre, binned_net]= makeEvenHistogram(pre_bout_ang_vels, net_bout_angular_acceleration,100);
        plot(binned_pre, binned_net,'k-');
        axis([0 20 -20 20]);    
        ylabel('Net change in angular velocity (deg/sec)')
        xlabel('Pre-bout ang. vel. (deg/sec)')
        title('21 dpf, nose-up')
        %plot best-fit line to nose-up data
        plot([0 20], polyval(nose_up_coeffs,[0 20]), '-m', 'LineWidth',3);

        %% Plot grouped data in Figure 2E
        %corrective gains by age and clutch for nose-up and
        %nose-down pre-bout angular velocities
        subplot(5,6,[23 29]); hold on;
        %Plot gain of angular velocity correction for bouts grouped across
        %clutches with nose-down pre-bout angular velocity
        plot(Ages,nose_down_gain(:,5),'k-','LineWidth',2)
        ylabel('Gain'); xlabel('Age (dpf)');
        title('Nose-down pre-bout ang. vel.')
        axis([0 22 0 1.5]);
        %Plot gain of angular velocity correction for bouts grouped across
        %clutches with nose-up pre-bout angular velocity
        subplot(5,6,[24 30]); hold on;
        plot(Ages,nose_up_gain(:,5),'k-','LineWidth',2)
        ylabel('Gain'); xlabel('Age (dpf)');
        title('Nose-up pre-bout ang. vel.')
        axis([0 22 0 1.5]);
    end
    
    
    %% Figure 2A & 2B, sort bout angular velocities into quintiles
    %sort pre-bout angular velocities by index
    [~,sortIndices] = sort(pre_bout_ang_vels);
    numBands = 5; %number of quantiles (5 for quintiles)
    tileIndices(1) = 0; %vector of indices at edges of quantiles. first quantile starts at 0.
    %iterate for each quantile
    for i=1:numBands
        %Calculate and plot bout angular velocities for this quantile
        subplot(5,6,[5 6 11 12 17 18]); hold on
        %calculate index corresponding to edge of quantile
        tileIndices(i+1) = floor(i/numBands*(length(pre_bout_ang_vels)));
        %average bout angular velocity for this quantile
        quintile_means(i,:) = mean(day_ang_vels(sortIndices(tileIndices(i)+1:tileIndices(i+1)),:));
        %plot pre-bout and post-bout angular velocity for this quantile in
        %Figure 2B, in age-specific column
        plot([age age+0.4],[mean(quintile_means(i,22:26)) mean(quintile_means(i,41:46))],'k.-','MarkerSize',15)
        
        %Plot angular velocity for this quantile in Fig 2A
        if age==2 %if 7dpf
            subplot(5,6,[7 8 13 14]); hold on;
            plot(linspace(-0.45,0.45,37),quintile_means(i,(13:49)),'b','LineWidth',2);
            axis([-0.45 0.45 -20 60]);    
            ylabel('Angular Velocity (deg/sec)')
            xlabel('Time (sec)')
        elseif age==4 %if 21dpf
            subplot(5,6,[9 10 15 16]); hold on;
            plot(linspace(-0.45,0.45,37),quintile_means(i,(13:49)),'m','LineWidth',2);
            axis([-0.45 0.45 -20 60]);    
            ylabel('Angular Velocity (deg/sec)')
            xlabel('Time (sec)')
        end
    
    end
    %Add labels, legends, and titles to Fig 2B
    subplot(5,6,[5 6 11 12 17 18]); hold on
    text(age, -14, ageFolders{age+2})
    ylabel('Angular Velocity (deg/sec)')
    xlabel('Age [4, 7, 14, 21dpf]')
    title('Pre- and post-bout angular velocity by quintile and age')
    axis([0.5 5 -15 10])
    set(gca,'XColor','w')
    
    %% Plot bout translation speed (mean +/- standard deviation) in Fig 2A
    if age==2 %if 7dpf
        subplot(5,6,[1 2]); hold on; title('7 dpf')
        plot(linspace(-0.45,0.45,37),mean(PropBoutAlignedSpeed(day,(13:49))),'k','LineWidth',2)
        plot(linspace(-0.45,0.45,37),mean(PropBoutAlignedSpeed(day,(13:49)))+std(PropBoutAlignedSpeed(day,(13:49))),'b')
        plot(linspace(-0.45,0.45,37),mean(PropBoutAlignedSpeed(day,(13:49)))-std(PropBoutAlignedSpeed(day,(13:49))),'b')
        ylabel('Speed (mm/sec)');
        axis([-0.45 0.45 -5 20]);
        
    elseif age==4 %if 21dpf
        subplot(5,6,[3 4]); hold on; title('21 dpf')
        plot(linspace(-0.45,0.45,37),mean(PropBoutAlignedSpeed(day,(13:49))),'k','LineWidth',2)
        plot(linspace(-0.45,0.45,37),mean(PropBoutAlignedSpeed(day,(13:49)))+std(PropBoutAlignedSpeed(day,(13:49))),'m')
        plot(linspace(-0.45,0.45,37),mean(PropBoutAlignedSpeed(day,(13:49)))-std(PropBoutAlignedSpeed(day,(13:49))),'m')
        ylabel('Speed (mm/sec)');
        axis([-0.45 0.45 -5 20]); 
    end
    
    %exit age folder
    cd ..
end







    


