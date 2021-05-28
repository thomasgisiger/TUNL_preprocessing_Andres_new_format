function r = correct_dataset(ds,lags,fol)

% cage dimensions
cage_dims = [240 332];
cage_prop = cage_dims(1)/cage_dims(2);

% shift to be performed
dt = 1/30;
delay = round(mean(lags(:)/dt));

disp(['Shifting by ' num2str(delay) ' frames.']);

% shift trajectory and dcs
Nf = size(ds.trialNumber,1);
% traj
traj = [ds.headPosition.x ds.headPosition.y ds.headPosition.p];
if delay>0
    temp = [traj(1,1:3).*ones(delay,3); traj(1:(Nf-delay),:)];
elseif delay<0
    temp = [traj((abs(delay)+1):Nf,:); traj(Nf,1:3).*ones(abs(delay),3)];
else
    temp = traj;
end
ds.headPosition.x = temp(:,1);
ds.headPosition.y = temp(:,2);
ds.headPosition.p = temp(:,3);

% dcs
dcs = ds.dcs;
Ncells = size(dcs,2);
if delay>0
    temp = [dcs(1,1:Ncells).*ones(delay,Ncells); dcs(1:(Nf-delay),1:Ncells)];
elseif delay<0
    temp = [dcs((abs(delay)+1):Nf,1:Ncells); dcs(Nf,1:Ncells).*ones(abs(delay),Ncells)];
else
    temp = dcs;
end

% points that will be used for the registration
temp = extract_cage_contact_points_4classes(ds,'p');
ds.registration.pos1 = temp.pos1;
ds.registration.pos5 = temp.pos5;
ds.registration.posR1 = temp.posR1;
ds.registration.posR5 = temp.posR5;

% do the same thing for the frames kept from the original videos
kf = ds.keptFrames';
if delay>0
    temp = [kf(1).*ones(delay,1); kf(1:(Nf-delay))];
elseif delay<0
    temp = [kf((abs(delay)+1):Nf); kf(Nf).*ones(abs(delay),1)];
else
    temp = kf;
end
ds.keptFrames = temp;


% -----------------------------------------------------------------

% Here, we used to scale the trajectory so that x in [-1,1] and y in [0,1],
% and also so that pad 1 would be top left, pad 5, top right, and R at the
% bottom. However, instead I will leave it as is with just a diagram just
% to show where the mean positions for 1, 5 and R are and let the user do
% this manually afterwards as we want to be able to match the trajectory
% carefully with the videos.

% (no correction)
dataset_corrected = ds;
dataset_corrected.headPosition.x = ds.headPosition.x;
dataset_corrected.headPosition.y = ds.headPosition.y;
% and the contact points
dataset_corrected.registration.pos1 = ds.registration.pos1; 
dataset_corrected.registration.pos5 = ds.registration.pos5; 
dataset_corrected.registration.posR1 = ds.registration.posR1; 
dataset_corrected.registration.posR5 = ds.registration.posR5; 

% -----------------------------------------------------------------

% check that the pad 1 points are on the left of the points for pad 5, and
% above points for reward
mean_pos1 = mean(dataset_corrected.registration.pos1,1);
mean_pos5 = mean(dataset_corrected.registration.pos5,1);
mean_posR1 = mean(dataset_corrected.registration.posR1,1);
mean_posR5 = mean(dataset_corrected.registration.posR5,1);

figure
plot(dataset_corrected.headPosition.x,dataset_corrected.headPosition.y)
hold on
text(mean_pos1(1),mean_pos1(2),'1','BackgroundColor',[1 1 1]);
text(mean_pos5(1),mean_pos5(2),'5','BackgroundColor',[1 1 1]);
text(mean_posR1(1),mean_posR1(2),'R1','BackgroundColor',[1 1 1]);
text(mean_posR5(1),mean_posR5(2),'R5','BackgroundColor',[1 1 1]);
xlabel('x')
xlabel('y')

drawnow

img = getframe(gcf);
imwrite(img.cdata,[fol '\qualitycheck\corrected_trajectory_markers.png']);

r = dataset_corrected;

end