function r = extract_cage_contact_points_4classes(ds,pstring)

% pstring = 'p', filter positions with p-values, 'nop' = do not filter
% using p.

% extracting the instants where sample ends, choice ends and reward starts
% (when mouse just touched the screen or the reward port). There is a
% subfield for each dataset.
mydata = ds;

% list of good trials
Ntrials = max(mydata.trialNumber);

% "touching" pad 1
% cue 1->5
sample1 = (mydata.trialType==1).*(mydata.taskPhase==2);
% target 5->1
target5 = (mydata.trialType==5).*(mydata.taskPhase==6);

% touching pad 5
% cue 5->1
sample5 = (mydata.trialType==5).*(mydata.taskPhase==2);
% target 1->5
target1 = (mydata.trialType==1).*(mydata.taskPhase==6);

% touching reward port
padR = (mydata.taskPhase==8);

% reward for each trial (we only keep correct trials)
R = zeros(1,Ntrials);
for tr=1:Ntrials
    ttime = find(mydata.trialNumber==tr);
    R(tr) = sum(mydata.reward(ttime));
    % clean up sample and target
    if R(tr)==-1
        % pad 1
        sample1(ttime) = 0;
        target5(ttime) = 0;
        % pad 5
        sample5(ttime) = 0;
        target1(ttime) = 0;
        % reward port
        padR(ttime) = 0;
    end
end

% Use diff to find instants when sample, target and reward pads are
% touched (= when taskPhase decreases by 1): time1, time5 and timeR.

% bring both together: all pad 1 touches during cue and choice periods of
% correct trials
pad1 = sample1 + target5;
pad1 = [0; diff(pad1)];
% <0 = end of sample/choice = when screen is touched
pad1 = +(pad1<0);
% mask out the points in the ill-sampled trials
pad1 = pad1.*(1-ds.illSampled);
% mask the instants where the mouse could not be tracked.
if strcmp(pstring,'p')
    pad1 = pad1.*mydata.headPosition.p;
end
% extract the corresponding frames
time1 = find(pad1==1);

% bring both together: all pad 5 touches during cue and choice periods of
% correct trials
pad5 = sample5 + target1;
pad5 = [0; diff(pad5)];
% <0 = end of sample/choice = when screen is touched
pad5 = +(pad5<0);
% mask out the points in the ill-sampled trials
pad5 = pad5.*(1-ds.illSampled);
% mask the instants where the mouse could not be tracked.
if strcmp(pstring,'p')
    pad5 = pad5.*mydata.headPosition.p;
end
% extract the corresponding frames
time5 = find(pad5==1);

% and times when the mouse touches the reward port
padR = [0; diff(padR)];
padR = +(padR>0);
% mask out the points in the ill-sampled trials
padR = padR.*(1-ds.illSampled);
% mask the instants where the mouse could not be tracked.
if strcmp(pstring,'p')
    padR = padR.*mydata.headPosition.p;
end

% The mouse moving in a straight line between touchpads 1 and 5 and the
% reward port, the positions of the mouse when receiving reward depends on
% the trial type.
% points closest to pad 1 (so when cue=5)
padR1 = padR.*(mydata.trialType==5);
% points closest to pad 5 (so when cue=1)
padR5 = padR.*(mydata.trialType==1);

% extract the corresponding frames
timeR1 = find(padR1==1);
timeR5 = find(padR5==1);

% extract the corresponding positions
traj = [mydata.headPosition.x mydata.headPosition.y];
pos1 = traj(time1,:);
pos5 = traj(time5,:);
posR1 = traj(timeR1,:);
posR5 = traj(timeR5,:);

% display
% clf
% hold on
% plot(traj(:,1),traj(:,2),'-','Color',[1 1 1]*0.75)
% plot(pos1(:,1),pos1(:,2),'or')
% plot(pos5(:,1),pos5(:,2),'ob')
% plot(posR1(:,1),posR1(:,2),'og')
% plot(posR5(:,1),posR5(:,2),'oc')

positions = [];

positions.pos1 = pos1;
positions.time1 = time1;
positions.pos5 = pos5;
positions.time5 = time5;
positions.posR1 = posR1;
positions.timeR1 = timeR1;
positions.posR5 = posR5;
positions.timeR5 = timeR5;

r = positions;

end