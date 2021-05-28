clear;

% data directroy
datadirs{1} = 'D:/Andres/0403E/TUNL/STAGE 1/S3/2 SEC/AUGUST_29_19_H9_M33_S18';

% mean sampling frequency inverse
dt = 1/30;

% color palette
palette = [1 0 0;0 1 0;0 0 1;1 1/2 0;1 0 1;0 1 1;0 0 0;1/2 1/2 1/2];

for d=1:length(datadirs)
    
    datadir = datadirs{d};
    
    mkdir([datadir '/qualitycheck']);
    
    disp(' ');
    disp('------------------------------------------------------------------')
    disp(['Reading DLC output from ' datadir]);
    
    % read the csv files and plot the trajectory
    trajdir = [datadir '/BehavCam_0'];
    
    file_pattern = [trajdir '/*.csv'];
    files = dir(file_pattern);
    Nfiles = length(files);
    
    % remove timestamp from list
    to_remove = -100;
    for f=1:Nfiles
        if strcmp(files(f).name,'timeStamps.csv')
            to_remove = f;
        end
    end
    if to_remove>0
        files(to_remove) = [];
    end
    Nfiles = length(files);
    
    % look for common name pattern, then remove it to extract the file number
    max_length = -10;
    for f=1:Nfiles
        max_length = max(max_length,length(files(f).name));
    end
    Nfiles = length(files);
    
    disp(['Found ' num2str(Nfiles) ' DLC .csv output files....']);
    
    % look for common pattern from beginning of file name
    common_begin = zeros(1,max_length);
    for c=1:max_length
        column = [];
        for f=1:Nfiles
            if c<=length(files(f).name)
                column = [column; files(f).name(c)];
            end
        end
        verdict = (column(1)==column);
        if sum(verdict)==Nfiles
            common_begin(c) = 1;
        end
    end
    lower = find(common_begin==0,1,'first');
    
    % look for common pattern from ending of file name
    common_end = zeros(1,max_length);
    for c=max_length:-1:1
        column = [];
        for f=1:Nfiles
            % right justify the name
            str = [blanks(max_length-length(files(f).name)) files(f).name];
            if c<=length(str)
                column = [column; str(c)];
            end
        end
        verdict = (column(1)==column);
        if sum(verdict)==Nfiles
            common_end(c) = 1;
        end
    end
    upper = find(common_end==0,1,'last');
    
    % extract the file number for each rank f and make list of the csv files read
    file_nums = zeros(Nfiles,2);
    for f=1:Nfiles
        str = files(f).name;
        % remove end
        n = max_length - upper;
        str(end-n+1:end) = [];
        % remove beginning
        str(1:lower-1) = [];
        file_nums(f,:) = [f str2num(str)];
    end
    
    % done! :)
    
    
    % read the csv files in the order of the recording
    [u,v] = sort(file_nums(:,2),'ascend');
    snums = file_nums(v,:);
    
    % declare arrays to speed up reading files
    traj = zeros(100000,18);
    prob = zeros(100000,9);
    index = 0;
    
    for f=1:Nfiles
        
        my_file = [files(snums(f,1)).folder '/' files(snums(f,1)).name];
        
        disp(['Processing ' my_file])
        
        % read the trajectory and the detection probability
        fid = fopen(my_file);
        tline = fgetl(fid);
        trial = 0;
        
        % csv file line number
        line_number = 0;
        
        while ischar(tline)
            
            tline = fgetl(fid);
            
            % count the lines
            line_number = line_number + 1;
            
            if not(isequal(tline,-1))
                
                % analyse the content of each line. "," is the separator
                separator_pos = find(tline==',');
                
                % number of points that were tracked
                nb_tracked = length(separator_pos)/3;
                
                if line_number>3
                    
                    index = index + 1;
                    
                    x = zeros(1,nb_tracked );
                    y = zeros(1,nb_tracked );
                    p = zeros(1,nb_tracked );
                    
                    xtemp = [];
                    ptemp = [];
                    
                    for i1=1:nb_tracked
                        
                        if i1<nb_tracked
                            x(i1) = str2num(tline(separator_pos((i1-1)*3+1)+1:separator_pos((i1-1)*3+2)-1));
                            y(i1) = str2num(tline(separator_pos((i1-1)*3+2)+1:separator_pos((i1-1)*3+3)-1));
                            p(i1) = str2num(tline(separator_pos((i1-1)*3+3)+1:separator_pos((i1-1)*3+4)-1));
                        else
                            x(i1) = str2num(tline(separator_pos((i1-1)*3+1)+1:separator_pos((i1-1)*3+2)-1));
                            y(i1) = str2num(tline(separator_pos((i1-1)*3+2)+1:separator_pos((i1-1)*3+3)-1));
                            p(i1) = str2num(tline(separator_pos((i1-1)*3+3)+1:end));
                        end
                        
                        xtemp = [xtemp x(i1) y(i1)];
                        ptemp = [ptemp p(i1)];
                    end
                    
                    traj(index,1:2*nb_tracked) = xtemp;
                    prob(index,1:nb_tracked) = ptemp;
                    
                end
                
            end
            
        end
        
        fclose(fid);
        
    end
    
    % remove the lines and columns full of 0s
    empty_cols = find(sum(abs(traj),1)==0);
    empty_lines = find(sum(abs(traj),2)==0);
    traj = traj(1:(empty_lines(1)-1),1:(empty_cols(1)-1));
    
    empty_cols = find(sum(abs(prob),1)==0);
    empty_lines = find(sum(abs(prob),2)==0);
    prob = prob(1:(empty_lines(1)-1),1:(empty_cols(1)-1));
    
    
    % =================== visual summary of the tracking ======================
    
    figure('Position',[716 221 1558 1076])
    
    Npoints = size(traj,2)/2;
    Nt = size(traj,1);
    
    % plot trajectory in 2d
    subplot(4,4,[1 2 5 6])
    hold on
    for point=Npoints:-1:1
        plot(traj(:,(point-1)*2+1),traj(:,(point-1)*2+2),'-','Color',palette(point,:))
    end
    xlabel('mouse position')
    title(datadir,'Interpreter','None')
    
    % plot x(t)
    subplot(4,4,[3 4])
    temp = traj(:,1:2:end);
    range = max(temp(:)) - min(temp(:));
    hold on
    for point=Npoints:-1:1
        plot(traj(:,(point-1)*2+1),'.','Color',palette(point,:))
    end
    axis([0 Nt min(temp(:))-range/10 max(temp(:))+range/10]);
    ylabel('x(t)')
    
    % plot y(t)
    subplot(4,4,[7 8])
    temp = traj(:,2:2:end);
    range = max(temp(:)) - min(temp(:));
    hold on
    for point=Npoints:-1:1
        plot(traj(:,(point-1)*2+2),'.','Color',palette(point,:))
    end
    axis([0 Nt min(temp(:))-range/10 max(temp(:))+range/10]);
    ylabel('y(t)')
    xlabel('time (frames)')
    
    % velocity distributions
    variations = diff(traj);
    var_per_point = [];
    
    subplot(4,4,[9 10 13 14])
    hold on
    for point=1:Npoints
        temp = [variations(:,(point-1)*2+1) variations(:,(point-1)*2+2)];
        var_per_point(point,:) = sqrt(sum(temp.*temp,2))/dt;
        
        [u,v] = hist(var_per_point(point,:),100);
        
        plot(v,u,'-','Color',palette(point,:))
        % plot y on a log scale
        set(gca,'YScale','log')
        ylabel(['Point ' num2str(point)])
    end
    
    xlabel('velocity (px/s)')
    axis tight
    
    labels = {};
    for point=1:Npoints
        labels{point} = ['point ' num2str(point)];
    end
    legend(labels);
    
    % mouse position distribution
    tempx = zeros(Npoints*Nt,1);
    tempy = zeros(Npoints*Nt,1);
    for point=1:Npoints
        tempx((point-1)*Nt+(1:Nt),1) = traj(:,(point-1)*2+1);
        tempy((point-1)*Nt+(1:Nt),1) = traj(:,(point-1)*2+2);
    end
    subplot(4,4,[11 12 15 16])
    hist3([tempx tempy],[100 100]);
    xlabel('x')
    ylabel('y')
    
    drawnow
    
    %pause
    
    % save image
    img = getframe(gcf);
    imwrite(img.cdata, [datadir '/qualitycheck/DLC_traj_summary.png']);
    
    close all
    
    % save the trajectory
    save([datadir '/trajectory.mat'],'traj','prob');
    
end















































