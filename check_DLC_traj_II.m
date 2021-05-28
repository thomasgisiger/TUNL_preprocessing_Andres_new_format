clear;

% data directroy
datadirs{1} = 'D:/Andres/0403E/TUNL/STAGE 1/S3/2 SEC/AUGUST_29_19_H9_M33_S18';


% color palette
palette = [1 0 0;0 1 0;0 0 1;1 1/2 0;1 0 1;0 1 1;0 0 0;1/2 1/2 1/2];

% read number of elements from current matlab's version colormap
temp = colormap;
nb_colors = size(temp,1);


for d=1:length(datadirs)
    
    datadir = datadirs{d};
    
    mkdir([datadir '/qualitycheck']);
    
    disp(' ');
    disp('------------------------------------------------------------------')
    disp(['Reading DLC output from ' datadir]);
    
    % Part 2: verify the relative timing between the trajectory and videos.
    viddir = [datadir '/BehavCam_0'];
    file_pattern = [viddir '/*.avi'];
    files = dir(file_pattern);
    Nvids = length(files);
    
    if Nvids==0
        disp('There are no videos to read in this folder. Moving on to the next session.');
    else
        
        nums = [];
        for i1=1:Nvids
            temp = files(i1).name;
            temp(end-3:end) = [];
            temp(1:length('behavCam')) = [];
            if not(isempty(temp)) && isempty(find(temp==' '))
                nums = [nums; i1 str2num(temp)];
            end
        end
        
        % sort the file numbers
        [u,v] = sort(nums(:,2),'ascend');
        Vnums = nums(v,:);
        
        % load the trajectories
        load([datadir '/trajectory.mat']);
        Nf = size(traj,1);
        Npoints = size(traj,2)/2;
        
        % make a visual summary
        lines = floor(sqrt(Nvids));
        columns = ceil(Nvids/lines);
        
        figure('Position',[177 129 2276 1222]);
        
        for v = 1:Nvids
            
            subplot(lines,columns,v)
            
            num = Vnums(v,2);
            avi_name = [datadir '/BehavCam_0/' num2str(num) '.avi'];
            msvidObj = VideoReader(avi_name);
            disp(['Loading video ' avi_name '...']);
            video=msvidObj.read();
            
            % while we are at it, compute the mean signal and save it in
            % .mat file.
            if exist([datadir '/qualitycheck/bh_mean' num2str(num) '.mat'])==0
                bh_mean_signal = squeeze(mean(mean(mean(video,1),2),3))';
                save([datadir '/qualitycheck/bh_mean' num2str(num) '.mat'],'bh_mean_signal');
            end
            
            f = 1+floor(rand*size(video,4));
            frame = nb_colors*mean(video(:,:,:,f),3)/255;
            image(frame)
            colormap(gray)
            hold on
            for point=1:Npoints
                fnum = (v-1)*1000+f;
                plot(traj(fnum,(point-1)*2+1),traj(fnum,(point-1)*2+2),'o',...
                    'MarkerEdgeColor',palette(point,:),'MarkerFaceColor',palette(point,:));
            end
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            xlabel(['Video ' num2str(v)]);
            
            if v==1
                title(datadir,'Interpreter','None');
            end
            
            drawnow
            
        end
        
        % save in a file
        img = getframe(gcf);
        imwrite(img.cdata, [datadir '/qualitycheck/DLC_traj_delay.png']);
        
    end
    
    close all
    
end








