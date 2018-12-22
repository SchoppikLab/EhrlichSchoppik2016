% Preprocess_LabView_output.m
% by David Ehrlich, April 21, 2015

% Associated with "Control of Movement Initiation Underlies the Development of Balance" 
% Ehrlich DE and Schoppik D, Curr Biol. 2017 Feb 6;27(3):334-344.

% Analyzes raw LabView swimming data to extract epochs, 
% continuous series of frames that had one fish-sized particle.

function output = Preprocess_LabView_output(raw,filenum)
filenameList = dir('*.dlm');
filename = filenameList(filenum).name;
% loop through and get rid of all epochs that are too short
minDur = 100; %in frames, at 40 Hz
% all epochs that have more that one fish
maxFish = 1;
% epochs where there was more than one fish but appear as one fish will
% have an improbably large instantaneous displacement
maxInstDisp = 35;

% or an improbably large angular velocity
maxAngVel = 100;
% smoothing for this detection
smWindow = 10;

% epochs with inexplicably large gaps between frame timestamps
maxDeltaT = .05;
% avoid edge effects
epochBuf = 3;

% scaling of video frames
scale = 60; %pixels/mm

% first, find the epochs
epochStop = [1;find(diff(raw(:,9)))];

% epoch index
ind = 1;

% check each epoch and write data
for i = 1:length(epochStop)-1
    % define the indices corresponding to each epoch, avoiding edge effects
    epochDex = epochStop(i):epochStop(i+1);
    epochDex = epochDex(epochBuf:end-epochBuf);
    
    % determine whether epoch meets our criteria
    if length(epochDex) >= minDur  % long enough?
        if max(raw(epochDex,2)) < maxFish    % only has one fish          
            x = raw(epochDex,4)-raw(epochDex(1),4); % horizontal position of fish centroid            
            y = - (raw(epochDex,5)-raw(epochDex(1),5)); %vertical position, with sign flipped so positive deltaY corresponds to upward motion            
            
            %get head position for determining orientation
            headx = (raw(epochDex,6)-raw(epochDex(1),6));
            heady = (raw(epochDex,7)-raw(epochDex(1),7));
            
            displacement = diff(sqrt(x.^2 + y.^2));            
            %if displacements are reasonable
            if isempty(find(abs(displacement) > maxInstDisp,1))
                %time
                t = diff(raw(epochDex,1));
                %if timesteps are reasonable
                if isempty(find(t > maxDeltaT, 1))
                    %if angular velocities are reasonables
                    angvel = smooth(diff(raw(epochDex,3))./diff(raw(epochDex,1)),smWindow);
                    if isempty(find(abs(angvel) > maxAngVel,1))                                              
                        % extract the things we care about
                        output(ind).absOrientation = raw(epochDex,3);
                        output(ind).centeredOrientation = raw(epochDex,3)-raw(epochDex(1),3);
                        % return centered and scaled (mm) coordinates for position
                        output(ind).x = x./scale;
                        output(ind).y = y./scale;
                        % and translation velocity
                        output(ind).xvel = diff(output(ind).x) ./ diff(raw(epochDex,1));
                        output(ind).yvel = diff(output(ind).y) ./ diff(raw(epochDex,1));
                        % angular and translational velocity 
                        output(ind).angularVelocity = diff(raw(epochDex,3)) ./ diff(raw(epochDex,1));
                        output(ind).velocity = (displacement./scale) ./ diff(raw(epochDex,1));
                        % timestamps
                        output(ind).t = raw(epochDex(2:end),1)-raw(epochDex(1),1);
                      
                        %absolute time in hours, mins, secs
                        hoursElapsed = floor(raw(epochDex,1)/3600);
                        minsElapsed = floor(raw(epochDex,1)/60)-60*hoursElapsed;
                        secsElapsed = raw(epochDex,1)-60*minsElapsed-3600*hoursElapsed;                     
                        output(ind).absHours = hoursElapsed+str2num(filename(8:9));
                        output(ind).absMins = minsElapsed+str2num(filename(11:12));
                        output(ind).absSecs = secsElapsed+str2num(filename(14:15));
                        if output(ind).absSecs >= 60
                            output(ind).absMins = output(ind).absMins + 1;
                            output(ind).absSecs = output(ind).absSecs - 60;
                        end
                        if output(ind).absMins >= 60
                            output(ind).absHours = output(ind).absHours + 1;
                            output(ind).absMins = output(ind).absMins - 60;
                        end
                        output(ind).timeInHrs = output(ind).absHours+(output(ind).absMins)./60+(output(ind).absSecs)./3600;

                        %output head position
                        output(ind).headx = headx./scale;
                        output(ind).heady = heady./scale;
                                               
                        ind = ind + 1;                        
                    end                    
                end
            end
        end
    end
end

end
