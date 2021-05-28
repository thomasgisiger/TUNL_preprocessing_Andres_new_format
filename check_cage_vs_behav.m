function check_cage_behav(ds)

W = 1000;
H = 800;

f = figure('Position',[500 500 W H]);

% list of trials
Ntrials = max(ds.trialNumber);
nt = length(ds.trialNumber);

% shift between cage and behavioral time series.
shift = 0;

% index quantifiying the quality of the synchronization of the cage and
% trajectory
dt = 1/30;
my_shifts = (-2/dt):(2/dt);
Ns = length(my_shifts);
my_index = zeros(1,Ns);

% extracting the instants where sample ends, choice ends and reward starts
% (when mouse just touched the screen or the reward port). 
mydata = ds;

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

% bring both together: all pad 1 touches during cue and choice periods of
% correct trials
pad1 = sample1 + target5;
pad1 = [0; diff(pad1)];

% mask out the points in the ill-sampled trials
pad1 = pad1.*(1-ds.illSampled);

% <0 = end of sample/choice = when screen is touched
pad1 = +(pad1<0);
% time1 = find(pad1==1);

% bring both together: all pad 5 touches during cue and choice periods of
% correct trials
pad5 = sample5 + target1;
pad5 = [0; diff(pad5)];
pad5 = +(pad5<0);

% mask out the points in the ill-sampled trials
pad5 = pad5.*(1-ds.illSampled);

% and times when the mouse touches the reward port
padR = [0; diff(padR)];

% mask out the points in the ill-sampled trials
padR = padR.*(1-ds.illSampled);

padR = +(padR>0);
timeR = find(padR==1);

time1 = [];
time5 = [];
timeR = [];

compute_index(ds);
show_me(0);


% ------------------------ build the gui ----------------------------------

Ygui = 20;
Xgui = 290;

% buttons to change shift
inc = uicontrol('Style','pushbutton','String','>','Position',[400+Xgui Ygui 50 25],'Callback',@F_increase);
decr = uicontrol('Style','pushbutton','String','<','Position',[10+Xgui Ygui 50 25],'Callback',@F_decrease);

% value of shift
txt = uicontrol('Style','text','String','Shift value (1/30 s)','Position',[180+Xgui Ygui+20 100 25]);
str = uicontrol('Style','edit','String',num2str(shift),'Enable','Off','Position',[200+Xgui Ygui 50 25]);



% ----------------------- function definitions ----------------------------
% vary shift values with the two buttons
    function F_increase(source,event)
        shift = shift + 1;
        set(str,'String',num2str(shift));
        show_me(shift);
    end

    function F_decrease(source,event)
        shift = shift - 1;
        set(str,'String',num2str(shift));
        show_me(shift);
    end

    function compute_index(ds)

        for i1=1:Ns
            my_shift = my_shifts(i1);
            
            % pad 1
            time1 = find(pad1==1);
            time1 = time1 +  my_shift;
            nogood = union(find(time1<1),find(time1>nt));
            time1(nogood) = [];
            % mouse positions at times when cage detected pad 1 contact
            pos1 = [ds.headPosition.x(time1) ds.headPosition.y(time1)];

            % pad 5
            time5 = find(pad5==1);
            time5 = time5 +  my_shift;
            nogood = union(find(time5<1),find(time5>nt));
            time5(nogood) = [];
            % mouse positions at times when cage detected pad 5 contact
            pos5 = [ds.headPosition.x(time5) ds.headPosition.y(time5)];
            
            % pad R
            timeR = find(padR==1);
            timeR = timeR +  my_shift;
            nogood = union(find(timeR<1),find(timeR>nt));
            timeR(nogood) = [];
            % mouse positions at times when cage detected pad R contact
            posR = [ds.headPosition.x(timeR) ds.headPosition.y(timeR)];
            
            % cluster the points at positions pos1 and pos5 and posR
            my_set = [pos1; pos5; posR];
            [id,c,sumd] = kmeans(my_set,3);
            
            % index to characterize the distribution of clusters and their relative
            % distance
            nb1 = length(find(id==1));
            nb2 = length(find(id==2));
            nb3 = length(find(id==3));
            d12 = c(1,:) - c(2,:);
            d12 = sqrt(d12*d12');
            d23 = c(2,:) - c(3,:);
            d23 = sqrt(d23*d23');
            d31 = c(3,:) - c(1,:);
            d31 = sqrt(d31*d31');
            
            % index is maximum if centroids are far apart and each contain lots of
            % points.
            my_index(i1) = nb1*nb2*d12 + nb2*nb3*d23 + nb3*nb1*d31;
        end

    end

    function show_me(shift)
        
        % clustering index
        subplot(10,4,1:12)
        plot(my_shifts*dt,my_index,'-ok')
        line(shift*[1 1]*dt,[min(my_index) 1.1*max(my_index)],'Color',[0 0 1]);
        line([my_shifts(1) my_shifts(end)]*dt,my_index(find(my_shifts==shift))*[1 1],'Color',[0 0 1]);
        ylim([min(my_index) 1.1*max(my_index)]);
        xlabel('shift (1/30s)')
        ylabel('index')
        
        
        % cage view
        h = subplot(10,4,union(18:4:34,19:4:35));
        cla(h);
        hold on
        % mouse traj in light gray
        plot(ds.headPosition.x,ds.headPosition.y,'Color',0.75*[1 1 1])
        
        % shifted position extraction
        % pad 1
        time1 = find(pad1==1);
        time1 = time1 + shift;
        nogood = union(find(time1<1),find(time1>nt));
        time1(nogood) = [];
        % mouse positions at times when cage detected pad 1 contact
        pos1 = [ds.headPosition.x(time1) ds.headPosition.y(time1)];
        
        % pad 5
        time5 = find(pad5==1);
        time5 = time5 + shift;
        nogood = union(find(time5<1),find(time5>nt));
        time5(nogood) = [];
        % mouse positions at times when cage detected pad 5 contact
        pos5 = [ds.headPosition.x(time5) ds.headPosition.y(time5)];
        
        % Reward pad
        timeR = find(padR==1);
        timeR = timeR + shift;
        nogood = union(find(timeR<1),find(timeR>nt));
        timeR(nogood) = [];
        % mouse positions at times when cage detected pad 5 contact
        posR = [ds.headPosition.x(timeR) ds.headPosition.y(timeR)];
        
        % mouse positions at time1 and 5
        plot(pos1(:,1),pos1(:,2),'or');
        plot(pos5(:,1),pos5(:,2),'ob');
        plot(posR(:,1),posR(:,2),'og');
        xlabel('x')
        ylabel('y')
        hold off
    end

end