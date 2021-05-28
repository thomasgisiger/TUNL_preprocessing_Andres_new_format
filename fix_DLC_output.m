clear;

session_folders = {};

d = 1;

session_folders{end+1} = 'AUGUST_29_19_H9_M33_S18';

session_path = ['D:\Andres\0403E\TUNL\STAGE 1\S3\2 SEC' '/' session_folders{d}];

disp([num2str(d) ') Fixing trajectory from session ' session_path '...']);

% % load the output from dlc
% load([session_path '/trajectory.mat'],'traj','prob');
% 
% % store it in files traj_raw.mat and traj_fixed.mat
% % raw
% traj_raw = [];
% traj_raw.x = traj(:,1:2:end);
% traj_raw.y = traj(:,2:2:end);
% traj_raw.p = prob;
% traj_raw.n = size(traj,1);
% save([session_path '/traj_raw.mat'],'traj_raw');
% 
% % fixed
% traj_fixed = [];
% traj_fixed.x = traj(:,1:2:end);
% traj_fixed.y = traj(:,2:2:end);
% traj_fixed.p = prob;
% traj_fixed.n = size(traj,1);
% save([session_path '/traj_fixed.mat'],'traj_fixed');

% if we are up to raw and fixed already, skip the first steps:
load([session_path '/traj_raw.mat'],'traj_raw');
load([session_path '/traj_fixed.mat'],'traj_fixed');

% =========================================================================
%                              CODE OUTLINE
% Goal: take fixed (interpolated and smoothed) trajectories (points 1 to 8) 
% from DLC and extract from it a trajectory for the miniscope/head of the
% mouse.

% We follow this list of steps:
% 1': extract delay between raw and fixed trajectories. Used to compensate
% for primitive low pass filtering.

% 2: extract miniscope/head position using point 8 (when p(8) = 1) or by 
% inferring its position using a straight line or circle fit onto points on
% the body of the mouse.

% 3: replace by NaN the points where velocity takes on huge values (this is
% an artefact of the process in 2).

% 3': replace points whose positions are clearly outliers by NaN.

% 4: merge trajectories in 3 and 3'

% 4': linear interpolation over all NaNs.

% 5: low-pass filtering of the trajectory.

% 6: mask out parts of the trajectory where the mouse is untractable for
% more than 5 time steps, or covers more than 10cm.

% 7: add in the delay computed in 2 and save the trajectory and bad data
% mask.

% =========================================================================



% used in inference of the position of the miniscope (yields distance
% between given tracking point and miniscope position).
% proportions of the dots on the mouse
props = [5 4 3 2.75 2] - 0.5;
% gives the distance between point # x and the miniscope position
fprops = @(x) ( (x==2)*props(1) + (x==3)*props(2) + (x==4)*props(3) + (x==57)*props(4) + (x==6)*props(5) );

% used to extract the list of points for making miniscope position
% inference (produces label that is stored in traj.track).
signature = @(x) ( 10.^(0:length(x)-1)*x );

% maximum mouse velocity (m/s)
global vth;
vth = 0.5;

% time resolution
global dt;
dt = 1/30;

% conversion factor pixel to mm (roughly 170mm = 303 pixels)
global fac;
fac = 0.6;

% sample freq is about 30 fps.
Fs = 1/dt;

% maximum of 5 outlier frames (= half wave of period equal to 10 frames then)
global limit_freq;
limit_freq = 10/Fs;

% IR camera frame dimensions
global IR_dim;
IR_dim = [303 422];

% threshold for rejecting series of frames where p is not 1
% consecutive frames where the mouse was untrackable
global consec_frames;
consec_frames = 5;

% covered distance by the mouse while it was untrackable (in mm)
global covered_dist;
covered_dist = 100;

% ------------------------------ variables --------------------------------


% 1) load the raw and fixed trajectories
[traj_raw,traj_fixed,Nraw,Nfixed] = load_trajs(session_path);
Nr = length(traj_raw.x(:,1));

% 1') check the relative timing between traj_raw and traj_fixed. A delay is
% possible because of temporal convolutoin with gaussian (low-pass filter)
% by Moh. opt_shift used below to have a trajectory that follows the videos 
% as closely as possible.
opt_shift = fix_delay_lp(traj_raw,traj_fixed);

% 2) extract the miniscope position from traj_fixed, either from the position
% of point 8 (miniscope), or by fitting a straight line or a circle on
% points 1 to 4 (spine), 5, 7 (ears) and 6 (head). If the position cannot
% be extracted, position is set to (NaN,NaN) and p and track to 0.
traj = extract_miniscope_position(traj_fixed,fprops,signature);

% 3) The inference performed in 2 often produces artefacts (jumps in
% position) some of which we can detect using the velocity: if it exceeds a
% certain threshold. We replace these points by NaNs.
nantraj1 = detect_velocity_artefacts(traj);

% 3') replace with NaNs points with x and y that are outliers, i.e. outside
% compact range of x values with non-zero frequencies.
nantraj2 = detect_position_artefacts(traj);

% 4) merge the trajectories corrected for velocity (nantraj1) and position
% (nantraj2) (essentially union of NaNs)
nantraj = merge_artefacts(nantraj1,nantraj2);

% 4') Interpolate over the NaNs created above
interptraj = interpolate_artefacts(nantraj);

% 5) low pass filter the interpolated trajectory
filttraj = lowpassfilt(interptraj);
filttraj.betas = traj.betas;
filttraj.params = traj.params;
filttraj.track = traj.track;
filttraj.line_inferred = traj.line_inferred;
filttraj.circle_inferred = traj.circle_inferred;
filttraj.p = traj.p;

% 6) mask out parts where the mouse is untrackable for too long or over too long
% distances
final_traj = mask_bad_parts(filttraj);

% 7) implement the delay from 2) and save the final processed trajectory
save_traj(final_traj,session_path,opt_shift);

% 8) save summary file for the trajectory's quality before and after miniscope
% position inference
save_summary(traj_fixed,final_traj,session_path);


% % 6) ====================== visual inspection =============================
% my_fig = figure('Position',[271 176 1820 1145]);
% 
% % check 1% of frames where point 8 (miniscope) has p=1, and all frames
% % where p(8) < 1.
% miniscope = traj_fixed.p(:,8);
% 
% % good frames: p(miniscope) = 1
% good_frames = find(miniscope==1);
% Ngood = length(good_frames);
% disp(['There are ' num2str(Ngood) ' (' num2str(round(100*Ngood/Nfixed)) '%) frames with p(miniscope) = 1.']);
% 
% % bad frames: p(miniscope) < 1
% bad_frames = find(miniscope<1);
% Nbad = length(bad_frames);
% disp(['There are ' num2str(Nbad) ' (' num2str(round(100*Nbad/Nfixed)) '%) frames with p(miniscope) < 1.']);
% 
% mkdir('good_frames')
% mkdir('bad_frames') 
% 
% checked = 1; %0.01;
% sample_good = good_frames(find(rand(1,length(good_frames))<checked));
% if Nbad/Nfixed>checked
%     % adjust so that bad frames verified represent 1% of whole dataset
%     prop = checked/(Nbad/Nfixed);
%     sample_bad = bad_frames(find(rand(1,length(bad_frames))<prop));
% else
%     sample_bad = bad_frames;
% end
% 
% % labeled list of frames (2nd column = quality of the frame: 1 = good, 2 = bad)
% labeled = zeros(Nfixed,2);
% labeled(1:Ngood,:) = [good_frames ones(Ngood,1)];
% labeled(Ngood+(1:Nbad),:) = [bad_frames 2*ones(Nbad,1)];
% [u,v] = sort(labeled(:,1),'ascend');
% labeled = labeled(v,:);
% 
% % number of movies to inspect
% movie_numbers = unique(ceil(union(sample_good,sample_bad)/1000));
% 
% % loop on movies
% for m=11:length(movie_numbers)
%     
%     % read movie file
%     my_file = [session_path '/behavCam' num2str(movie_numbers(m)) '.avi'];
%     disp(['Reading ' my_file]);
%     obj = VideoReader(my_file);
%     one_movie = double(obj.read());
%     disp('done')
%     
%     % inspect good and then bad frames
%     for quality=3
%         
%         % 1: display good frames, 2: display bad frames, 3: display all
%         % frames
%         switch quality
%             case 1
%                 % check good frames
%                 label = 'Good frames';
%                 disp(label)
%                 my_frames = sample_good;
%                 my_dir = 'good_frames';
%             case 2
%                 % check bad frames
%                 label = 'Bad frames';
%                 disp(label)
%                 my_frames = sample_bad;
%                 my_dir = 'bad_frames';
%             case 3
%                 % look at all frames
%                 label = 'all frames';
%                 disp(label)
%                 my_frames = labeled(:,1);
%                 my_dir = '';
%         end
%         
%         % frames for this video
%         myframes = my_frames(find((my_frames>1000*(m-1)).*(my_frames<=1000*m)));
%         frames = myframes - 1000*(m-1);
%         
%         for f=1:length(frames)
%             
%             clf
%             image(squeeze(one_movie(:,:,:,frames(f)))/255)
%             hold on
%             
%             if quality<3
%                 display_inferences(traj_fixed,filttraj,myframes(f),session_path,quality);
%             else
%                 display_inferences(traj_fixed,filttraj,myframes(f),session_path,labeled(myframes(f),2));
%             end
%             
%             drawnow
%             
% %             img = getframe(gcf);
% %             imwrite(img.cdata,[my_dir '/frame' num2str(temp(f)) '.png']);
%             
%             pause
%         end
%         
%     end
%     
% end


% ============================ functions ==================================

% read raw and fixed trajectories
function [traj_raw,traj_fixed,Nraw,Nfixed] = load_trajs(pth)

disp('Reading the trajectories...')

raw_path = [pth '/traj_raw.mat'];
fixed_path = [pth '/traj_fixed.mat'];

load(raw_path,'traj_raw');
load(fixed_path,'traj_fixed');

Nraw = traj_raw.n;
Nfixed = traj_fixed.n;

disp('Done')
disp(' ');

end

% extract delay between raw and fixed trajectories
function optsh = fix_delay_lp(trajr,trajf)

padding = 20;

% pad both time series to evaluate optimal timing
Nr = length(trajr.x);
tsr = [ones(1,padding)*trajr.x(1,8) trajr.x(:,8)' ones(1,padding)*trajr.x(end,8)];
tsf = [ones(1,padding)*trajf.x(1,8) trajf.x(:,8)' ones(1,padding)*trajf.x(end,8)];

% compute correlation with various shifts
rs = zeros(1,1+length(padding));
shifts = -padding/2:padding/2;
n =0;
for i1=shifts
    tsrc = trajr.x(:,8)';
    tsfc = tsf(padding+i1+(1:Nr));
    n = n + 1;
    rs(n) = corr(tsrc',tsfc');
end

% pick the relative delay that maximizes the correlation
[mmax,mpos] = max(rs);
optsh = shifts(mpos);

end

% extract miniscope position from traj_fixed
function traj = extract_miniscope_position(traj_fixed,fpr,sgn)

disp('Extracting/extrapolating the mouse position...');

Nfixed = length(traj_fixed.x);

% miniscope position
traj.x = zeros(Nfixed,1);
traj.y = zeros(Nfixed,1);
% tracking probability
traj.p = zeros(Nfixed,1);
% point(s) used for tracking
traj.track = zeros(Nfixed,1);

% straight line fitting
% coefficients
traj.betas = zeros(Nfixed,2);
% miniscope position inferred from the line
traj.line_inferred = zeros(Nfixed,2);
% residual of the fit
traj.line_res = zeros(Nfixed,1);

% circle fitting
% circle parameters
traj.params = zeros(Nfixed,4);
% miniscope position inferred from the circle
traj.circle_inferred = zeros(Nfixed,2);
% residual of the fit
traj.circle_res = zeros(Nfixed,1);

% handles to fprops and signature
fprops = fpr;
signature = sgn;

% extract the trajectory
for f=1:Nfixed
    
    % contains relevant data for the fitting, testing and inference.
    data = [];
    data.num = [];
    data.x = [];
    data.y = [];
    r = Inf;
    
    % points to be fit
    % mouse spine
    for i1=1:4
        if traj_fixed.p(f,i1)==1
            data.num = [data.num; i1];
            data.x = [data.x; traj_fixed.x(f,i1)];
            data.y = [data.y; traj_fixed.y(f,i1)];
        end
    end
    % ears (5 and 7, 5 for convenience)
    if traj_fixed.p(f,5)==1 && traj_fixed.p(f,7)==1
        data.num = [data.num; 5];
        data.x = [data.x; mean(traj_fixed.x(f,[5 7]))];
        data.y = [data.y; mean(traj_fixed.y(f,[5 7]))];
    end
    % head
    if traj_fixed.p(f,6)==1
        data.num = [data.num; 6];
        data.x = [data.x; traj_fixed.x(f,6)];
        data.y = [data.y; traj_fixed.y(f,6)];
    end
    Np = size(data.num,1);
    
    % fitting straight line
    if Np>=2
        % ls fit
        M = [ones(Np,1) data.x];
        % solution
        betas = pinv(M)*data.y;
        traj.betas(f,:) = betas;
        
        % closest points on the line to the fitted data
        myline = [0:450; betas(1)+betas(2)*(0:450)];
        for i1=1:Np
            my_dist = myline' - [data.x(i1) data.y(i1)];
            my_dist = sqrt(sum(my_dist.^2,2));
            [mm,mypos] = min(my_dist);
            data.line_closest_x(i1) = myline(1,mypos);
            data.line_closest_y(i1) = myline(2,mypos);
        end
        
        % deviation of the line to the original points
        temp = [data.x - data.line_closest_x' data.y - data.line_closest_y'];
        traj.line_res(f) = mean(sqrt(sum(temp.^2,2)));
        
        % mean distance between fitted points
        d = [data.x data.y];
        d = dist(d');
        d = triu(d,1)-triu(d,2);
        d = d(find(d~=0));
        d = mean(d);
        
        % director vector extracted from a straight line only defined up to a
        % factor, so we choose the convention that points in the same
        % direction as the vector pointing from 1 to 4:
        bodyvect = [data.x(end)-data.x(1) data.y(end)-data.y(1)];
        a = betas(2);
        vect = [1 a]/sqrt(1+a^2);
        if bodyvect(1)<0
            vect = - vect;
        end
        
        % inference of the miniscope position
        traj.line_inferred(f,:) = [data.line_closest_x(end) data.line_closest_y(end)] + d*fprops(data.num(end))*vect;
        
    end
    
    % fitting circle
    if Np>=3
        % circle fit
        [x0,y0,r,orientation,myangles] = circle_fit_or(data.x,data.y);
        % infer the head/miniscope position
        traj.params(f,:) = [x0 y0 r orientation];
        
        % closest points on the circle to the fitted data
        circle = [x0+r*cos(0:pi/50:2*pi); y0+r*sin(0:pi/50:2*pi)];
        for i1=1:Np
            my_dist = circle' - [data.x(i1) data.y(i1)];
            my_dist = sqrt(sum(my_dist.^2,2));
            [mm,mypos] = min(my_dist);
            data.circle_closest_x(i1) = circle(1,mypos);
            data.circle_closest_y(i1) = circle(2,mypos);
        end
        
        % deviation of the line to the original points
        temp = [data.x - data.circle_closest_x' data.y - data.circle_closest_y'];
        traj.circle_res(f) = mean(sqrt(sum(temp.^2,2)));
        
        % mean angle between points fitted on
        d = mean(diff(myangles));
        
        % infer the position of the miniscope
        circle_inferred = [x0 + cos(d*fprops(data.num(end)))*(data.x(end)-x0) - sin(d*fprops(data.num(end)))*(data.y(end)-y0) ...
            y0 + cos(d*fprops(data.num(end)))*(data.y(end)-y0) + sin(d*fprops(data.num(end)))*(data.x(end)-x0)];
        traj.circle_inferred(f,:) = circle_inferred;
    end
    
    % mouse position extraction: we want to use the miniscope as position
    % of the head (to follow the same convention as for the visible light
    % tracking).
    if traj_fixed.p(f,8)==1
        % use the miniscope
        traj.x(f) = traj_fixed.x(f,8);
        traj.y(f) = traj_fixed.y(f,8);
        traj.p(f) = 1;
        traj.track(f) = 8;
    else
        % select the trajectory according to the number of points that are reliable, thefitted straight
        % line, and the circle
        if length(data.num)<2
            % at most one point is reliable: nothing much to do here.
            traj.x(f) = NaN;
            traj.y(f) = NaN;
            traj.p(f) = 0;
            traj.track(f) = 0;
        else
            if length(data.num)<=3
                % two points are reliable: use the straight line
                % fit
                traj.x(f) = traj.line_inferred(f,1);
                traj.y(f) = traj.line_inferred(f,2);
                traj.p(f) = 1;
                traj.track(f) = signature(data.num);
            else
                if length(data.num)>=4
                    % at least 4 points are reliable and we can
                    % infer the positoin of the miniscope with a
                    % straight line or a circle. We use the one
                    % with the smallest residual.
                    if (traj.line_res(f)<traj.circle_res(f)) || (r<=40)
                        % straight line gives better inference
                        traj.x(f) = traj.line_inferred(f,1);
                        traj.y(f) = traj.line_inferred(f,2);
                        traj.p(f) = 1;
                        traj.track(f) = signature(data.num);
                    else
                        % straight line gives better inference
                        traj.x(f) = traj.circle_inferred(f,1);
                        traj.y(f) = traj.circle_inferred(f,2);
                        traj.p(f) = 1;
                        traj.track(f) = signature(data.num);
                    end
                end
            end
        end
    end
    
end

disp('Done')
disp(' ');

end

% detect points with unrealistic velocity
function nantraj = detect_velocity_artefacts(traj)
global vth dt fac;

disp('Setting velocity artefacts to NaN...')

% velocity in m/s
% converts pixels to m
traj_m.x = fac*traj.x/1000;
traj_m.y = fac*traj.y/1000;

vx = [0; diff(traj_m.x)]/dt;
vy = [0; diff(traj_m.y)]/dt;
v = sqrt(vx.*vx + vy.*vy);

% detect excessive velocities
Nf = length(v);
nantraj.x = traj.x;
nantraj.y = traj.y;
for t=1:Nf
    if v(t)>vth 
        nantraj.x(t) = NaN;
        nantraj.y(t) = NaN;
    end
end

disp('Done')
disp(' ');

end

% detect points with outlying positions
function outtraj = detect_position_artefacts(nantraj)

global IR_dim;

disp('Setting position artefacts to NaN...')

temp = nantraj;

% We find outliers by computing the distributions of x and y, and then
% looking for the values that give the tightest value lumps. For this, we
% find the center of the (x,y) lump and find the range along x and y from
% each side.

center = [nanmean(nantraj.x) nanmean(nantraj.y)];

% extract compact range of x values with non-zero frequencies
[ux,vx] = hist(nantraj.x,IR_dim(2));
[uy,vy] = hist(nantraj.y,IR_dim(1));

% find vx and vy closest to center.
% x
my_dist = vx - center(1);
[mm,mp] = min(abs(my_dist));
middlex = mp;
% y
my_dist = vy - center(2);
[mm,mp] = min(abs(my_dist));
middley = mp;

% extract compact range of x values with non-zero frequencies
rangex = [];
% left
index = find(ux(1:middlex)==0,1,'last');
if isempty(index)
    rangex(1) = vx(1);
else
    rangex(1) = vx(index);
end
% right
index = find(ux(end:-1:middlex)==0,1,'last');
if isempty(index)
    rangex(2) = vx(end);
else
    rangex(2) = vx(IR_dim(2)-index);
end

% extract compact range of y values with non-zero frequencies
rangey = [];
% left
index = find(uy(1:middley)==0,1,'last');
if isempty(index)
    rangey(1) = vy(1);
else
    rangey(1) = vy(index);
end
% right
index = find(uy(end:-1:middley)==0,1,'last');
if isempty(index)
    rangey(2) = vy(end);
else
    rangey(2) = vy(IR_dim(1)-index);
end

% flag outliers with NaNs
for f=1:length(nantraj.x)
    if nantraj.x(f)<rangex(1) || nantraj.x(f)>rangex(2)
        temp.x(f) = NaN;
    end
    if nantraj.y(f)<rangey(1) || nantraj.y(f)>rangey(2)
        temp.y(f) = NaN;
    end
end

outtraj = temp;

disp('Done')
disp(' ')

end

% merge trajectories with NaNs
function ntraj = merge_artefacts(ntraj1,ntraj2)

N = length(ntraj1.x);
ntraj.x = zeros(1,N);
ntraj.y = zeros(1,N);

for f=1:N
    % any point with NaN in traj1 or traj2, x or y yield NaN 
    if isnan(ntraj1.x(f)) || isnan(ntraj1.y(f)) || isnan(ntraj2.x(f)) || isnan(ntraj2.y(f))
        ntraj.x(f) = NaN;
        ntraj.y(f) = NaN;
    else
        % otherwise keep the position
        ntraj.x(f) = ntraj1.x(f);
        ntraj.y(f) = ntraj1.y(f);
    end
end

end

% interpolate over the artefacts
function interptraj = interpolate_artefacts(nantraj)

Nf = size(nantraj.x,1);

disp('Interpolating over NaNs...')

% linear interpolation
% interpolated over NaNs that are in the trajectories
% .x
nan_traj = isnan(nantraj.x);
traj_t = 1:numel(nantraj.x);
interptraj.x = nantraj.x;
interptraj.x(nan_traj) = round(interp1(traj_t(~nan_traj), interptraj.x(~nan_traj), traj_t(nan_traj)));

% .y
nan_traj = isnan(nantraj.y);
traj_t = 1:numel(nantraj.y);
interptraj.y = nantraj.y;
interptraj.y(nan_traj) = round(interp1(traj_t(~nan_traj), interptraj.y(~nan_traj), traj_t(nan_traj)));

% fix NaNs at both ends of all_traj
% .x
nanlistx = find(isnan(interptraj.x));
if not(isempty(nanlistx))
    % check if led invisible at the start
    if nanlistx(1)==1
        temp = find(not(isnan(interptraj.x)),1,'first');
        interptraj.x(1:temp-1) = interptraj.x(temp);
    end
    
    % check if led invisible at the start
    if nanlistx(end)==Nf
        temp = find(not(isnan(interptraj.x)),1,'last');
        interptraj.x(temp+1:end) = interptraj.x(temp);
    end
end

% .y
nanlisty = find(isnan(interptraj.y));
if not(isempty(nanlisty))
    % check if led invisible at the start
    if nanlisty(1)==1
        temp = find(not(isnan(interptraj.y)),1,'first');
        interptraj.y(1:temp-1) = interptraj.y(temp);
    end
    
    % check if led invisible at the start
    if nanlisty(end)==Nf
        temp = find(not(isnan(interptraj.y)),1,'last');
        interptraj.y(temp+1:end) = interptraj.y(temp);
    end
end

disp('Done')
disp(' ');

end

% low-pass filter trajectory
function filttraj = lowpassfilt(interptraj)

disp('Low-pass filtering the trajectories...')

% next, do a low pass filtering using filtfilt and butterworth filter
global limit_freq;
d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',limit_freq,'DesignMethod','butter');
filttraj.x = filtfilt(d1,interptraj.x);
filttraj.y = filtfilt(d1,interptraj.y);

disp('Done')
disp(' ');

end

% function that displays the tracked points, the fitted straight line and 
% circle, and the inferred miniscope positions.
function display_inferences(traj_fixed,traj,f,session_path,quality)

    colors = [1 1 113; 3 17 231; 14 164 249; 76 248 184; 185 248 83; 245 183 18; 239 49 11; 113 3 7];

    % display what was extracted
   
    % draw all tracked points with p = 1
    for p=1:8
        if traj_fixed.p(f,p)==1
            plot(traj_fixed.x(f,p),traj_fixed.y(f,p),'o', ...
                'MarkerFaceColor',colors(p,:)/255, ...
                'MarkerEdgeColor',[1 1 1])
            text(traj_fixed.x(f,p)+10,traj_fixed.y(f,p),num2str(p),'Color',[1 0 0])
        end
    end
    
    % indicate with green crosshair the extracted mouse position
    if traj.track(f)>0
        line(traj.x(f)*[1 1],[1 300],'Color',[0 1 0])
        line([1 450],traj.y(f)*[1 1],'Color',[0 1 0])
    end
    
    % display fitted straight line and circle if miniscope position was
    % infered
    if quality==2
        
        % draw straight line passing through the points tracking the mouse's spine
        b = traj.betas(f,1);
        a = traj.betas(f,2);
        line([1 422],[a*1+b a*422+b],'Color',[0 1 1])
        
        % draw inferred position of head point
        % from straight line
        plot(traj.line_inferred(f,1),traj.line_inferred(f,2),'dr');
        % from circle
        plot(traj.circle_inferred(f,1),traj.circle_inferred(f,2),'*m');
        
        % draw the circle spanning the spine points
        x0 = traj.params(f,1);
        y0 = traj.params(f,2);
        r = traj.params(f,3);
        circle = [x0+r*cos(0:pi/50:2*pi); y0+r*sin(0:pi/50:2*pi)];
        plot(circle(1,:),circle(2,:),'Color',[1 0 0])
        
    end
    
    xlabel(['frame ' num2str(f) ', session ' session_path],'Interpreter','None');
    
    xlim([-100 550])
    ylim([-100 400])
    
    drawnow

end

% function to mask out parts of the trajectory where the mouse can't be
% tracked for more than a certain number of frames, or over more than a
% certain distance.
function ftraj = mask_bad_parts(mytraj)

global consec_frames covered_dist fac

missing = 1-mytraj.p;
% divide in blocks of 0's bordered by 1's
[labeled_blocks,numRegions] = bwlabel(missing,4);

% compute the duration and distance spanned by each block in labeled_blocks
% compute duration of the block where the led disappears
durations = zeros(1,numRegions);
for i1=1:numRegions
    durations(i1) = length(find(labeled_blocks==i1));
end

% compute the distance between the beginning and end of the block
distances = zeros(1,numRegions);
Nf = length(mytraj.x);

for s=1:numRegions
    temp = find(labeled_blocks==s);
    
    % led position immediately before and after the block
    initial = temp(1)-1;
    if initial==0
        initial = 1;
    end
    final = temp(end)+1;
    if final>Nf
        final = Nf;
    end
    p1 = round([mytraj.x(initial) mytraj.y(initial)]);
    p2 = round([mytraj.x(final) mytraj.y(final)]);
    my_dist = p1-p2;
    my_dist = sqrt(sum(my_dist.^2));
    my_dist = round(my_dist);
    
    distances(s) = my_dist;
end

% finally make a mask for the data we keep: durations<5 and
% distances<10cm/fact
bad_data = zeros(1,Nf);
for s=1:numRegions
    my_block = find(labeled_blocks==s);
    if durations(s)>consec_frames || distances(s)>(covered_dist/fac)
        bad_data(my_block) = 1;
    end
end

% save the preprocessed trajectory and a variable labeling frames where
% the mouse's position could not be determined.
ftraj.x = mytraj.x;
ftraj.y = mytraj.y;
ftraj.p = 1 - bad_data;

end

% function to save the preprocessed trajectory into a format adequate for
% synchronization
function save_traj(traj,session_path,opt_sh)

disp('Saving the trajectory...')

% create final trajectory
position = [];
position.bodypart1.x = traj.x;
position.bodypart1.y = traj.y;
position.bodypart1.bad_data = 1 - traj.p;

% shift the trajectory by value computed in 2.
N = length(position.bodypart1.x);
n = 2*abs(opt_sh);
tsx = [ones(1,n)*position.bodypart1.x(1) position.bodypart1.x(:)' ones(1,n)*position.bodypart1.x(end)];
tsy = [ones(1,n)*position.bodypart1.y(1) position.bodypart1.y(:)' ones(1,n)*position.bodypart1.y(end)];
tsbd = [ones(1,n)*position.bodypart1.bad_data(1) position.bodypart1.bad_data(:)' ones(1,n)*position.bodypart1.bad_data(end)];

position.bodypart1.x = tsx(n+opt_sh+(1:N));
position.bodypart1.y = tsy(n+opt_sh+(1:N));
position.bodypart1.bad_data = tsbd(n+opt_sh+(1:N));

% Finally, save the trajectory
save([session_path '/mouse_traj_corrected.mat'],'position');

disp('Done');

end

% function to make trajectory quality summary
function save_summary(traj_fixed,fintraj,session_path)

disp('Saving trajectory summary file...')

fig = figure('Position',[912 663 980 672]);

% traj after DLC
good = find(traj_fixed.p(:,8)==1);
Nf = length(traj_fixed.p(:,8));
bad = setdiff(1:Nf,good);

subplot(2,2,1)
plot(traj_fixed.x(:,8),traj_fixed.y(:,8),'-b')
hold on
plot(traj_fixed.x(bad,8),traj_fixed.y(bad,8),'.r')
xlabel('miniscope x')
ylabel('miniscope y')
axis tight

subplot(2,2,2)
hold on
bar(-1,length(good)/Nf,'FaceColor',[0 0 1])
bar(1,length(bad)/Nf,'FaceColor',[1 0 0])
set(gca,'XTick',[-1 1],'XTickLabel',{'=1','<1'})
xlabel('trajectory from DLC')
ylabel('Proportion of frames with p(miniscope)')

% traj after fixing by inference
good = find(fintraj.p==1);
bad = setdiff(1:Nf,good);

subplot(2,2,3)
plot(fintraj.x(:),fintraj.y(:),'-b')
hold on
plot(fintraj.x(bad),fintraj.y(bad),'.r')
xlabel('miniscope x')
ylabel('miniscope y')
axis tight

subplot(2,2,4)
hold on
bar(-1,length(good)/Nf,'FaceColor',[0 0 1])
bar(1,length(bad)/Nf,'FaceColor',[1 0 0])
set(gca,'XTick',[-1 1],'XTickLabel',{'=1','<1'})
xlabel('trajectory after inference')
ylabel('Proportion of frames with p(miniscope)')

drawnow

img = getframe(gcf);
imwrite(img.cdata,[session_path '/qualitycheck/summary_miniscope_position.png']);

close(fig);

disp('Done')
disp(' ');

end


































