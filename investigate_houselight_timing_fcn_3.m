function r = investigate_houselight_timing_fcn(fol,dataset)

% create rectangular kernel with width 150 (approx 5s of dt = 1/30). Uses
% to detect houslight pulses (which are 5s long).
kernel = zeros(1,200+150);
kernel(1,100+(1:150)) = 1;

fig = figure('Position',[327 313 2030 1022]);

datadir = fol;

disp(['Session ' datadir])

% read the houselight onset/offset straight from the schedules
if exist([datadir '/schedules.csv'],'file')==2
    
    fid = fopen([datadir '/schedules.csv']);
    tline = fgetl(fid);
    onset_schedules = [];
    offset_schedules = [];
    ttl = [];
    
    % csv file line number
    while ischar(tline)
        
        tline = fgetl(fid);
        
        if not(isequal(tline,-1))
            
            % analyse the content of each line. "," is the separator
            separator_pos = find(tline==',');
            % read columns 1 (time) and 4 (house light state)
            col1 = tline(1:separator_pos(1)-1);
            col3 = tline(separator_pos(2)+1:separator_pos(3)-1);
            col4 = tline(separator_pos(3)+1:separator_pos(4)-1);
            
            if strcmp(col4,'TO Houselight On')
                onset_schedules = [onset_schedules; str2num(col1)];
            end
            if strcmp(col4,'Reset Houselight off')
                offset_schedules = [offset_schedules; str2num(col1)];
            end
            
            % also look for ttl while we're at it
            if strcmp(col3,'Output On Event') && strcmp(col4,'TTL #1')
                ttl = str2num(col1);
            end
        end
        
    end
    
    fclose(fid);
    
end

% only keep the part of the schedules on/off-set that are part of
% dataset: reject all onset/offset that take place after last
% experiment instant in the dataset. This will be used just in the plot
% below.
lasttime = dataset.experimentTimes(end);
pos = find(onset_schedules>lasttime);
onset_schedules(pos) = [];
pos = find(offset_schedules>lasttime);
offset_schedules(pos) = [];

% remove onset/offsets that are in ill sampled trials (label them with NaN
% and remove them all at once).
for trial=1:max(dataset.trialNumber)
    tframes = find(dataset.trialNumber==trial);
    ttimes = dataset.experimentTimes([tframes(1) tframes(end)]);
    if unique(dataset.illSampled(tframes))==1
        for p=1:size(onset_schedules,1)
            if onset_schedules(p)>=ttimes(1) && onset_schedules(p)<ttimes(2)
                onset_schedules(p) = NaN;
            end
        end
        for p=1:size(offset_schedules,1)
            if offset_schedules(p)>=ttimes(1) && offset_schedules(p)<ttimes(2)
                offset_schedules(p) = NaN;
            end
        end
    end
end
onset_schedules(find(isnan(onset_schedules))) = [];
offset_schedules(find(isnan(offset_schedules))) = [];


% We will compute the lag between the bh videos and the cage using only the
% good trial (illSampled==0) = on/offset_dataset - on/offset_videos.

% A) read the houselight onset/offset from dataset.houseLight variable (derivative
% is 1 or -1).
onset_dataset = dataset.experimentTimes(find([0 diff(dataset.houseLight')]==1));
offset_dataset = dataset.experimentTimes(find([0 diff(dataset.houseLight')]==-1));

% remove the onsets/offsets that fall in ill sampled trials (label them with NaN
% and remove them all at once).
for trial=1:max(dataset.trialNumber)
    tframes = find(dataset.trialNumber==trial);
    ttimes = dataset.experimentTimes([tframes(1) tframes(end)]);
    if unique(dataset.illSampled(tframes))==1
        for p=1:size(onset_dataset,1)
            if onset_dataset(p)>=ttimes(1) && onset_dataset(p)<ttimes(2)
                onset_dataset(p) = NaN;
            end
        end
        for p=1:size(offset_dataset,1)
            if offset_dataset(p)>=ttimes(1) && offset_dataset(p)<ttimes(2)
                offset_dataset(p) = NaN;
            end
        end
    end
end
onset_dataset(find(isnan(onset_dataset))) = [];
offset_dataset(find(isnan(offset_dataset))) = [];

% B) read the houselight onset/offset from the mean video using the
% synchronization frames -  already done earlier. Simply read the adequate
% .mat file.
hl_videos = [];
load([fol '/houselight_bh_videos.mat'],'hl_videos','meansig');

houselight = hl_videos;

% houselight was created from the video signal, so we need to synchronize
% it with the miniscope, exactly as we did for the trajectory. We do this
% with the synchronization matrix:
if exist([datadir '/msTouchSync_new.mat'],'file')
    load([datadir '/msTouchSync_new.mat']);
else
    load([datadir '/msTouchSync.mat']);
end
n = length(dataset.experimentTimes);
houselight = houselight(synchronization.miniscopeMaster.slaveFrames(1:n));

% Remove houselight pulses that fall in ill sampled trials
houselight = houselight.*(1-dataset.illSampled');

% extract the frame rate (which we will plot at the bottom).
fr_rate = [0 diff(synchronization.miniscopeMaster.slaveTimes)];
fr_rate = fr_rate(1:n);

% extract the onset and offset values from the videos (which will be
% compared with those from the schedule file).
onset_videos = dataset.experimentTimes(find([0 diff(houselight)]==1));
offset_videos = dataset.experimentTimes(find([0 diff(houselight)]==-1));

% this is from this that we
% will evaluate the lag between the video and the schedule
vd_dtst_delays = [onset_dataset-onset_videos offset_dataset-offset_videos];

% return the delays between videos and dataset. If there is no well-sampled
% failed trials, then simply return [0 0]
if isempty(vd_dtst_delays)
    r = [0 0];
else
    r = vd_dtst_delays;
end

% Make sure that nb of houselight switches is equal to the nb of failed
% trials
failed = [];
Ntrials = max(dataset.trialNumber);
for tr=1:Ntrials
    trframes = find(dataset.trialNumber==tr);
    R = max(dataset.reward(trframes));
    if (R==0) && (unique(dataset.illSampled(trframes))==0)
        failed = [failed tr];
    end
end
Nfailed = length(failed);

% print out a summary if number of onset/offset is not equal to the number
% of failed trails.
if length(onset_videos)~=Nfailed || length(offset_videos)~=Nfailed
    disp('Problem:')
    disp(['Nb of failed trials: ' num2str(Nfailed)])
    disp(['Nb of house light onset: ' num2str(length(onset_videos))])
    disp(['Nb of house light offset: ' num2str(length(offset_videos))])
    pause(0.1)
    problem = 1;
    
    for i1=1:10
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    end
    pause
else
    problem = 0;
end

% ================================= display ===============================

% compare schedules vs synchronized vs synchronized video
clf(fig)

% number of columns in the summary
ncols = 5;

% -------------------------------------------------------------------------
% display the video luminosity and the houselight on/off events from the
% schedule.
subplot(5,ncols,unique([1:ncols-1 ncols+(1:ncols-1)]))
hold on
% from the schedules
for i1=1:Nfailed
    if (i1<=length(onset_schedules)) && (i1<=length(offset_schedules))
        line(onset_schedules(i1)*[1 1],[0 1.25],'Color',[0 1 0])
        line(offset_schedules(i1)*[1 1],[0 1.25],'Color',[1 0 0])
    end
end

% from the synchronized data
plot(dataset.experimentTimes,dataset.houseLight,'b')

% and from the videos
N = length(dataset.trialNumber);
smeansig = meansig(synchronization.miniscopeMaster.slaveFrames(1:N));
temp = smeansig - min(smeansig);
temp = temp/max(temp);
plot(dataset.experimentTimes,temp,'m')
xlim([0 dataset.experimentTimes(end)])

% and the ill-sampled data
plot(dataset.experimentTimes,dataset.illSampled,'c-')

% the houselight variable extracted from the videos
plot(dataset.experimentTimes,houselight,'--k')


% -------------------------------------------------------------------------
% Display the frame rate and illsampled trials
subplot(5,ncols,2*ncols+(1:ncols-1))
% the frame rate
plot(dataset.experimentTimes,fr_rate/1000,'-','Color',[0.5 0.5 0])
ylabel('frame rate')
xlim([0 dataset.experimentTimes(end)])
ylim([0 1.1])

% -------------------------------------------------------------------------
% display the lag in a table
subplot(5,ncols,unique([3*ncols+(1:ncols-1) 4*ncols+(1:ncols-1)]))
cla
set(gca,'XTick',[])
set(gca,'YTick',[])
box on
% display onsets
ylim([-Nfailed 3])
xlim([-0.05 1])
dx = 0.1;
dy = 1;

xlabel(datadir,'Interpreter','None')
if length(onset_dataset)==length(onset_videos) && length(offset_dataset)==length(offset_videos)
    
    x = 0;
    y = 2;
    text(x,y,'Onsets:','HorizontalAlignment','left')
    y = y - dy;
    text(x+dx,y,'schedules','HorizontalAlignment','left')
    for i1=1:Nfailed
        y = y - dy;
        text(x+dx,y,num2str(onset_schedules(i1)),'HorizontalAlignment','left')
    end
    
    x = x+dx;
    y = 2;
    y = y - dy;
    text(x+dx,y,'dataset','HorizontalAlignment','left')
    for i1=1:Nfailed
        y = y - dy;
        text(x+dx,y,num2str(onset_dataset(i1)),'HorizontalAlignment','left')
    end
    
    x = x+dx;
    y = 2;
    y = y - dy;
    text(x+dx,y,'videos','HorizontalAlignment','left')
    for i1=1:Nfailed
        y = y - dy;
        text(x+dx,y,num2str(onset_videos(i1)),'HorizontalAlignment','left')
    end
    
    x = x+dx;
    y = 2;
    y = y - dy;
    text(x+dx,y,'dataset-videos','HorizontalAlignment','left')
    for i1=1:Nfailed
        y = y - dy;
        text(x+dx,y,num2str(onset_dataset(i1)-onset_videos(i1)),'HorizontalAlignment','left')
    end
    
    % display offsets
    x = x + 2*dx;
    y = 2;
    text(x,y,'Offsets:','HorizontalAlignment','left')
    y = y - dy;
    text(x+dx,y,'schedules','HorizontalAlignment','left')
    for i1=1:Nfailed
        y = y - dy;
        text(x+dx,y,num2str(offset_schedules(i1)),'HorizontalAlignment','left')
    end
    
    x = x+dx;
    y = 2;
    y = y - dy;
    text(x+dx,y,'dataset','HorizontalAlignment','left')
    for i1=1:Nfailed
        y = y - dy;
        text(x+dx,y,num2str(offset_dataset(i1)),'HorizontalAlignment','left')
    end
    
    x = x+dx;
    y = 2;
    y = y - dy;
    text(x+dx,y,'videos','HorizontalAlignment','left')
    for i1=1:Nfailed
        y = y - dy;
        text(x+dx,y,num2str(offset_videos(i1)),'HorizontalAlignment','left')
    end
    
    x = x+dx;
    y = 2;
    y = y - dy;
    text(x+dx,y,'dataset-videos','HorizontalAlignment','left')
    for i1=1:Nfailed
        y = y - dy;
        text(x+dx,y,num2str(offset_dataset(i1)-offset_videos(i1)),'HorizontalAlignment','left')
    end
    
    xlabel(datadir,'Interpreter','None')
    
    % ---------------------------------------------------------------------
    % Display stats of onsets et offsets for the session
    subplot(5,ncols,ncols*(1:5))
    hold on
    onset = vd_dtst_delays(:,1);
    for i1=1:length(onset)
        line(1/3*[1 1],[min(onset) max(onset)],'Color',[0 0 0])
        plot(1/3,onset,'ro')
        plot(1/3,mean(onset),'ob')
    end
    offset = vd_dtst_delays(:,2);
    for i1=1:length(offset)
        line(2/3*[1 1],[min(offset) max(offset)],'Color',[0 0 0])
        plot(2/3,offset,'ro')
        plot(2/3,mean(offset),'ob')
    end
    ylabel('mean shift of onset/offset of houselight')
    xlim([0 1])
    if not(isempty(vd_dtst_delays))
        ylim([min(vd_dtst_delays(:)) max(vd_dtst_delays(:))]);
    end
    
else
    disp('Could not match the number of failed trials with the houselight pulses. Saving the figure anyways.');
end

drawnow

% Save the result
img = getframe(gcf);
imwrite(img.cdata,[datadir '/qualitycheck/houselight_alignment','.png']);

if problem==1
    imwrite(img.cdata,['bad_houselight_alignment','.png']);
end

%close all

end