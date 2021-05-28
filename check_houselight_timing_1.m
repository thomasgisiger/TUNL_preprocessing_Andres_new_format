clear;

% data directory
datadirs{1} = 'D:/Andres/0403E/TUNL/STAGE 1/S3/2 SEC/AUGUST_29_19_H9_M33_S18';

% average sampling interval (s)
dt = 1/30;

% boundaries to pick the house light pulses from the bh video mean
good_max = [80 140];


% create rectangular kernel with width 150 (approx 5s of dt = 1/30). Uses
% to detect houslight pulses (which are 5s long).
kernel = zeros(1,200+150);
kernel(1,100+(1:150)) = 1;

for d=1:length(datadirs)
    
    datadir = datadirs{d};
    
    mkdir([datadir '/qualitycheck']);
    
    % =========================================================================
    % read the mean signal from behavioral videos/mean files (used to visualize
    % the timing of the houselight events).
    
    % make video list
    viddir = [datadir '/BehavCam_0'];
    file_pattern = [viddir '/*.avi'];
    files = dir(file_pattern);
    Nvids = length(files);
    nums = [];
    Vnums = [];
    
    if Nvids>0
        nums = [];
        for i1=1:Nvids
            temp = files(i1).name;
            temp(end-3:end) = [];
            if not(isempty(temp)) && isempty(find(temp==' '))
                nums = [nums; i1 str2num(temp)];
            end
        end
        
        % sort the file numbers
        [u,v] = sort(nums(:,2),'ascend');
        Vnums = nums(v,:);
    end
    
    % make mean signal file list
    meandir = [datadir '/qualitycheck'];
    file_pattern = [meandir '/bh_mean*.mat'];
    files = dir(file_pattern);
    Nmeans = length(files);
    nums = [];
    Mnus = [];
    
    if Nmeans>0
        nums = [];
        for i1=1:Nmeans
            temp = files(i1).name;
            temp(1:length('bh_mean')) = [];
            temp(end-3:end) = [];
            if not(isempty(temp)) && isempty(find(temp==' '))
                nums = [nums; i1 str2num(temp)];
            end
        end
        
        % sort the file numbers
        [u,v] = sort(nums(:,2),'ascend');
        Mnums = nums(v,:);
    end
    
    % merge the signals from the .mat and (if necessary).avi files
    Nfiles = max(size(Vnums,1),size(Mnums,1));
    
    mean_bh = [];
    
    for f=1:Nfiles
        
        avifile = '';
        if Nvids>0
            numv = Vnums(f,2);
            avifile = [viddir '/' num2str(numv) '.avi'];
        end
        
        meanfile = '';
        if Nmeans>0
            numm = Mnums(f,2);
            meanfile = [meandir '/bh_mean' num2str(numm) '.mat'];
        end
        
        % load the mat file if it exists, the avi file if not, or stop
        % otherwise
        if exist(meanfile,'file')==2
            disp(['Loading mean file ' meanfile '...']);
            load(meanfile);
            mean_bh = [mean_bh bh_mean_signal];
        else
            if exist(avifile,'file')==2
                msvidObj = VideoReader(avifile);
                disp(['Loading video ' avifile '...']);
                video=msvidObj.read();
                bh_mean_signal = squeeze(mean(mean(mean(video,1),2),3));
                mean_bh = [mean_bh; bh_mean_signal];
            else
                disp('Cannot find mean signal file. Stopping here.')
            end
        end
    end
    
    
    % =========================================================================
    % extract the houselight events from the schedules file
    
    schedules = [datadir '/schedules.csv'];
    
    % do some schedule.csv file cleanup
    
    % remove all the '"'
    fid = fopen(schedules,'rt');
    X = fread(fid);
    fclose(fid);
    X = char(X.');
    % replace string S1 with string S2
    Y = strrep(X, '"', '');
    
    % remove the file's header which finishes at 'Arg5_Value'
    pos_start = strfind(Y,'Arg5_Value');
    if not(isempty(pos_start))
        pos_start = pos_start + length('Arg5_Value');
        Y(1:pos_start) = [];
    end
    
    fid2 = fopen(schedules,'wt');
    fwrite(fid2,Y);
    fclose (fid2);
   
    
    % declare arrays to speed up reading files
    % 1 = time, 2 = on(1)/off(-1)/otherwise(0),
    houselight = zeros(100000,2);
    
    index = 0;
    
    % read the trajectory and the detection probability
    fid = fopen(schedules);
    tline = fgetl(fid);
    
    % important information starts after line that starts with Evnt_Time
    good_stuff = 0;
    
    while ischar(tline)
        
        tline = fgetl(fid);
        
        % make sure that we are in the relevant part of the file
        if length(tline)>length('Evnt_Time')
            % if line is the header, skip it and start reading after that
            if contains(tline,'Evnt_Time')
                good_stuff = 1;
                tline = fgetl(fid);
            else
                % if line starts with a 0, we are in the data that needs to
                % be read
                if str2num(tline(1))==0
                    good_stuff = 1;
                end
            end
        end
        
        if not(isequal(tline,-1)) && good_stuff
            
            % analyse the content of each line. "," is the separator
            separator_pos = find(tline==',');
            
            % extract 4th field
            field4 = tline(separator_pos(3)+1:separator_pos(4)-1);
            
            % extract third field
            field3 = tline(separator_pos(2)+1:separator_pos(3)-1);

            % extract first field
            field1 = tline(1:separator_pos(1)-1);

            % check if this is a houselight event
            index = index + 1;
            
            % store the time
            houselight(index,1) = str2num(field1);
            
            % store flags if this is a houselight event
            if strcmp(field4,'HouseLight #1')
                % check if this is a on/off
                if strcmp(field3,'Output On Event')
                    houselight(index,2) = 1;
                else
                    if strcmp(field3,'Output Off Event')
                        houselight(index,2) = -1;
                    else
                        disp('Strange Houselight field in schedules')
                        pause
                    end
                end
            end
            
        end
    end
    
    fclose(fid);
    
    % remove the excess
    houselight(index+1:end,:) = [];
    
    % =====================================================================
    
    % read the houselight onset/offset from the mean video using the
    % synchronization frames.
    meansig = mean_bh;
    nsig = meansig - min(meansig);
    nsig = nsig/max(nsig);
    
    % 2) compute its ddt
    dsig = [0 diff(nsig)];
    
    % 3) convolve the signal with a square kernel: maxes should be
    % where kernel is sitting over houselight pulse
    tempc = conv(nsig,kernel,'same');
    
    % 4) look for local maxima of convolved signal above some threshold value
    [x,y]=findpeaks(tempc,'MinPeakProminence',40);
    
    % filter the maxima using boundaries (see start of script)
    pos = [y' x'];
    good_pos = (pos(:,2)>=good_max(1)).*(pos(:,2)<=good_max(2));
    pos = pos(find(good_pos),:);
    
    % plot the bh signal and recreated houselight
    figure('Position',[484 793 1724 550]);
    subplot(3,1,2)
    hold on
    plot(tempc)
    plot(y,x,'or')
    axis tight
    xlabel('time (frames)')
    ylabel('-d/dt video signal');
    
    
    % 5) extract from it the corresponding approximate onset and offset
    % values
    aonsets = [];
    aoffsets = [];
    % NB: house light pulse is 5s, so 150 dt, half of which is 75.
    for j1=1:size(pos,1)
        aonsets = [aonsets pos(j1,1)-75];
        aoffsets = [aoffsets pos(j1,1)+75];
    end
    
    % 6) tweak the onset and offset positions so as to maximize the value at
    % onset, and the dddt at offset
    for j1=1:length(aonsets)
        % adjust onsets
        posa = aonsets(j1)+(-10:10);
        scan = dsig(posa);
        [m,p] = max(scan);
        aonsets(j1) = posa(p);
        
        % adjust offset
        posb = aoffsets(j1)+(-10:10);
        scan = dsig(posb);
        [m,p] = min(scan);
        aoffsets(j1) = posb(p);
    end
    
    % recreate binary house light signal and label it houselight.
    Nt = length(meansig);
    hl_videos = zeros(1,Nt);
    for t=2:Nt
        if hl_videos(t-1)==0 && any(t==aonsets)
            hl_videos(t) = 1;
        else
            if hl_videos(t-1)==1 && any(t==aoffsets)
                hl_videos(t) = 0;
            else
                hl_videos(t) = hl_videos(t-1);
            end
        end
    end
    
    % save it: we will need it when we merge all data in a single dataset
    save([datadir '/houselight_bh_videos.mat'],'hl_videos','meansig');

    % =========================================================================
    
    % make summary figure
    subplot(3,1,1)
    hold on
    % recreated houselight from bh videos
    plot(hl_videos*max(mean_bh))
    title(datadir,'Interpreter','None')
    % mean bh video signal
    plot(mean_bh)
    axis tight
    xlabel('Time (frames)')
    ylabel('Mean bh video signal')
    
    % plot the part of schedules that was recorded by the bh camera
    subplot(3,1,3)
    Nf = length(mean_bh);
    upper = Nf*dt;
    houselight(find(houselight(:,1)>upper),:) = [];
    plot(houselight(:,1),houselight(:,2))
    axis tight
    xlabel('Time (s)')
    ylabel('houselight from schedules')
    
    % save image
    img = getframe(gcf);
    imwrite(img.cdata, [datadir '/qualitycheck/houselight_verification.png']);
pause
    close all
    
end















