function r = convergence(fol)

% new features:
% produces trial with NaN taskPhase variable for trials with unsolvable
% problems (e.g. frame rate drops to huge levels).

% ------------------------ folder/file names ------------------------------
% data folder
datafolder = fol;

% ------------------- variables -------------------------------------------
% trajectory (will contain the raw, synched and scaled trajectories)
% trajectory.bodypart1.raw.x/y = raw imported trajectory
% trajectory.sync.x/y = after being synchronized with miniscope
% trajectory.scaled.x/y = after being normalized in [-1 1] for x, and [0 1]
% 
% for y.
% trajectory.p = probability that the mouse position is correct.
% (contains both normalized and original trajectory).
trajectory = [];


% miniscope data: miniscope = deconvolved calcium signal
% .original = from ms. and .cleaned = after removal last part after
% synchronization
miniscope = [];
synchro_matrix = [];

% cage data
% ---------
% cage.events (session events and times in cage time):
% .trial_times/correction/start_delay_times/delay_end_times/samples/sample_times/display_sample_times/ ...
% targets/target_times/touch_down_times/touch_down_position/touch_up_times/touch_up_position/break_BIR_times/ ...
% poke_tray_times/tray_light_on_times/tray_light_pulse_times/tray_light_pulse_durations/tray_light_off_times/ ...
% sound1_times/sound1_durations/correction_trial/display_image/display_image_times/hide_image/hide_image_times/ ...
% reward/original_reward/reward_times/ITI_start_times/house_light_on_times/house_light_off_times/feeder_times/ ...
% feeder_durations/incentives_times

% cage.timeseries (event times series in synchronized time):
% trial_number_ts/trial_type_ts/sample_times_dt/sample_ts/target_times_dt/target_ts/hide_image_times_dt/break_BIR_ts/...
% poke_tray_ts/reward_ts/incentives_ts/touch_down_ts/touch_up_ts/sound_ts/light_ts/feeder_ts/houselight_ts/epochs_ts
cage = [];
cage.events = [];
cage.timeseries = [];

% deconvolved calcium signal
Ncells = [];

% time increment (in seconds)
global dt;  % global is ok here as this is essentially a constant.
dt = 1/30;

Ntrials = [];
Nnew = [];
Ttotal = [];

% function's final output
dataset = [];

% total duration of the recording (in cage time)

% threshold radial distance used to classify if the mouse is still eating
% reward (in normalized units x in [-1,1] and y in [0,1]);
% or not
global R_dist_thresh   % global is ok here as this is essentially a constant.
R_dist_thresh = 1/5;


% ----------------- main body of the script -------------------------------

% read the files for the cage, synchronization, trajectory and miniscope
% in datafolder, and store the mouse position in trajectory, the
% deconvolved calcium signal in miniscope, and the synchronization matrix in
% synchro_matrix. It also calls function read_variables which goes through
% the content of schedules.csv and stores it in cage.events and Ttotal.
[trajectory,miniscope,synchro_matrix,cage,Ttotal] = readfiles(datafolder);

% synchronize the trajectory and cage events with the miniscope time
% series, and save everything in a single structure.
% It also calls
[dataset,Nnew,Ntrials] = synchro(synchro_matrix,trajectory,cage,miniscope,Ttotal,fol);

% simple gui that allows to display to what extent mouse is near pads 1 and
% 5 when the cage detects contacts there.
checkcagebehav(datafolder,dataset)

% return the dataset for further use
r = dataset;


%  ======================= functions =============================

% function to read event timing and durations in the cage structure
% some information about which information is stored where:
% Time is read in the first columne of the file (A).

% The sample/target identities are a bit trickier since in correction
% trials, where the sample/target of the preceding trial are used, the
% sample/target are not explicitly displayed. The sample and target
% identities in non-correcion trials can be read by looking in column D for
% Sample_position (=sample) and Correct_position (=target) in column I
% (=arg1_values). Correction trials can be spotted by looking for
% "Start_correction_trial" in column D, and we simply repeat the sample and
% target from the previous trial.

% The "Whisker - Display Image" (column C) and "Bussey Mouse Operant Mode 5 x 1 x low" in
% column D mark when the sample/target/distractor images are displayed (column
% E = "image") and hidden (column E = "background") once for the sample, and once for the target
% and another for the distractor. Note: "Bussey Mouse Operant Mode 5 x 1 x low" is also
% a flag to spot when the mouse touches the screen (see below).

% The presentation onset time of the sample is simpler to detect as it is labeled
% by "Display_sample" in column D for correction and non-correction trials.

% The delay boundaries can be easily spotted for all trials using flags
% "Start_delay" and "End_delay" in column D (again for all trial types).

% The times when the mouse breaks the back infra-red beam can be read using
% flags "Input Event" or "Input Transition On Event" or "Input Transition Off Event"
% in column C and "BIRBeam #1" in column D. There appears to be different nomenclatures
% depending on the version of the cage softward, which separates the event into
% on and off. My guess is that the upgrade was meant to extract from the
% cage the BIR break not as a binary event, but rather also making its
% duration available. We have to cover both cases below.
%
% The times when the mouse pokes into the reward hole are marked by "Input Event" or
% "Input Transition On Event" or "Input Transition Off Event" in column C
% and "Tray #1" in column D (again there appears to be different nomenclatures
% depending on the version of the cage softward, which separates the event into
% on and off).

% Contact of the mouse onto the response wall are read through variables
% "Touch Down Event" or "Touch Up Event" in column C, and
% "Bussey Mouse Operant Mode 5 x 1 x low" in column D, and nothing in column E.

% The tray light can have three states: on (= on until some mouse action
% turns it off), pulse (= on during at most during a specific duration) and
% off. This is read with flags "Output On Event", "Pulse Output Event" and
% "Output Off Event" in column C and "TrayLight #1" in column D.

% Sounds in the experiment are only pulses of specific length labeled by
% "Pulse Output Event" in column C, "Sound_On #1" in column D and the pulse
% duration in column I.

% There are several flags that can be used to read if a trial is correct or
% not. We chose looking for flag "Reward Collected Start ITI" in column D
% for correct trial, and "Incorrect Trails to Correction Trial" or
% incorrect ones. We collect the reward (1/-1) and the time of collection.

% The presence of an incentive during the delay can be spotted with flag
% "Sample Touch rewarded" or "Sample Touch not rewarded" in column D. We
% collect if the incentive is present or not, and the time when the cage
% decides whether to dispense it. (this is not perfect but we cannot use
% nose poke times as if the animal does not poke, then we get a really
% wrong value - stay with this.).

% The start of the ITI period is marked by "Reward Collected Start ITI" in
% column D for rewarded trials, and "Incorrect Trails to Correction Trial"
% for non-rewarded trials. The former is just 1ms after the tray signals
% that is was contacted.

% Houselight activity is marked by flags "Output On Event" or "Output Off Event"
% in column C and "HouseLight #1" in column D.

% The feeder is activated in pulses of .4 s for the incentives, and .8 for
% trial reward. This is marked by "Pulse Output Event" in column C and
% "Feeder #1" in column D. The duration is in column I.

    function [cg,ttot] = read_variables(folder)
        
        fid = fopen([folder '/schedules.csv']);
        tline = fgetl(fid);
        trial = 0;
        
        % reinitialize the cage variable
        cg.events = [];
        
        % flags the beginning of the task
        task_started = 0;
        
        % csv file line number
        line_number = 0;
        
        while ischar(tline)
            
            %disp(tline)
            tline = fgetl(fid);
            
            % count the lines
            line_number = line_number + 1;
            
            if not(isequal(tline,-1))
                
                % analyse the content of each line. "," is the separator
                separator_pos = find(tline==',');
                
                % first column contains the time stamp
                time = str2num(tline(1:separator_pos(1)));
                
                % fourth column contains a lot of the event stamps.
                
                % read the item name
                item_name = tline(separator_pos(3)+1:separator_pos(4)-1);
                
                % read the event name
                evnt_name = tline(separator_pos(2)+1:separator_pos(3)-1);
                
                % read alias name
                alias_name = tline(separator_pos(4)+1:separator_pos(5)-1);
                
                % read arg1 value (contains the sample, target, nose poked position,
                % etc.)
                arg1_value = tline(separator_pos(8)+1:separator_pos(9)-1);
                
                % count the trials: locate when the first trial begins: "Tray Entry starts first trial"
                % or locate each new trial: "Next trial"
                if strcmp(item_name,'Tray Entry starts first trial') || strcmp(item_name,'Next trial')
                    
                    % detect beginnin of first trial
                    if task_started==0
                        task_started = 1;
                    end
                    
                    trial = trial + 1;
                    if isfield(cg.events,'trial_times')==0
                        cg.events.trial_times = [];
                    end
                    cg.events.trial_times(end+1) = time;
                    disp(' ');
                    disp(['Trial ' num2str(trial) ' starts at time ' num2str(time)]);
                end
                
                % locate trial sample, target and correction type:
                % -------------------------------------------------
                % 1) non-correction trials (the first trial can't be a
                % correction trial)
                if  task_started && strcmp(item_name,'Tray Entry starts first trial') || strcmp(item_name,'Next trial non corrected Trial')
                    if isfield(cg.events,'correction')==0
                        cg.events.correction = [];
                    end
                    cg.events.correction(end+1) = 0;
                    disp('Trial is not type correction.');
                end
                
                % read the sample/distractor time and position (if not correction
                % trial)
                if task_started && strcmp(item_name,'Sample_Position')
                    sample_label = str2num(arg1_value);
                    if isfield(cg.events,'samples')==0
                        cg.events.samples = [];
                    end
                    cg.events.samples(end+1) = sample_label;
                    disp(['Sample ' num2str(sample_label) ' presented at time ' num2str(time)]);
                end
                
                % 2) correction trials
                % if the current trial is a correction trial, just repeat the same
                % sample/distractor and target (first trial can't be a correction trial
                % so cage.events.correction, cage.events.samples and cage.events.targets are not empty).
                if task_started && strcmp(item_name,'Start Correction Trial')
                    cg.events.correction(end+1) = 1;
                    disp('Trial is type correction');
                    % we use the same sample/target as in the previous trial
                    cg.events.samples(end+1) = cg.events.samples(end);
                    cg.events.targets(end+1) = cg.events.targets(end);
                end
                
                % read the sample presentation time during correction/non
                % correction trials
                if task_started && strcmp(item_name,'Display Sample') && strcmp(evnt_name,'Condition Event')
                    if isfield(cg.events,'sample_times')==0
                        cg.events.sample_times = [];
                    end
                    cg.events.sample_times(end+1) = time;
                    sample_time = time;
                    disp(['Sample displayed at time ' num2str(time)]);
                end
                
                % read the target position
                if task_started && strcmp(item_name,'Correct_Position')
                    target_label = str2num(arg1_value);
                    if isfield(cg.events,'targets')==0
                        cg.events.targets = [];
                    end
                    cg.events.targets(end+1) = target_label;
                    disp(['Target ' num2str(target_label)]);
                end
                
                % locate the start and end of the delay
                if task_started && strcmp(item_name,'Start Delay')
                    if isfield(cg.events,'start_delay_times')==0
                        cg.events.start_delay_times = [];
                    end
                    cg.events.start_delay_times(end+1) = time;
                    disp(['Start delay at time ' num2str(time)]);
                end
                if task_started && strcmp(item_name,'Delay End')
                    if isfield(cg.events,'delay_end_times')==0
                        cg.events.delay_end_times = [];
                    end
                    cg.events.delay_end_times(end+1) = time;
                    disp(['End delay at time ' num2str(time)]);
                end
                
                % read the time where the target and distractor are presented: look for
                % whisker display image - Bussey Mouse Operant Mode 5 x 1 x low - Image
                % after the sample has been presented
                if task_started && strcmp(evnt_name,'Whisker - Display Image') && ...
                        strcmp(item_name,'Bussey Mouse Operant Mode 5 x 1 x low') && ...
                        strcmp(alias_name,'Image') && ...
                        time>sample_time
                    if isfield(cg.events,'target_times')==0
                        cg.events.target_times = [];
                    end
                    cg.events.target_times(end+1) = time;
                    disp(['Target/distractor presented at time ' num2str(time)]);
                end
                
                % look for when images are hidden
                if task_started && strcmp(evnt_name,'Whisker - Display Image') && ...
                        strcmp(item_name,'Bussey Mouse Operant Mode 5 x 1 x low') && ...
                        strcmp(alias_name,'Background')
                    if isfield(cg.events,'hide_image')==0
                        cg.events.hide_image = [];
                    end
                    cg.events.hide_image(end+1) = str2num(arg1_value);
                    if isfield(cg.events,'hide_image_times')==0
                        cg.events.hide_image_times = [];
                    end
                    cg.events.hide_image_times(end+1) = time;
                end
                
                % mouse breaks the back IR beam: two cases = 1) event is
                % punctual, and 2) event has a finite duration. We use two
                % variables: cage.events.break_BIR_start_times, and if
                % needed cage.events.break_BIR_end_times
                if task_started
                    % start of event
                    if (strcmp(evnt_name,'Input Event') || ...
                            strcmp(evnt_name,'Input Transition On Event')) && ...
                            strcmp(item_name,'BIRBeam #1')
                        % record start of BIR break
                        if isfield(cg.events,'break_BIR_start_times')==0
                            cg.events.break_BIR_start_times = [];
                        end
                        cg.events.break_BIR_start_times(end+1) = time;
                        disp(['Back IR beam broken at time ' num2str(time)]);
                    end
                    
                    % end of event (if it was recorded)
                    if strcmp(evnt_name,'Input Transition Off Event') && ...
                            strcmp(item_name,'BIRBeam #1')
                        % record end of BIR break
                        if isfield(cg.events,'break_BIR_end_times')==0
                            cg.events.break_BIR_end_times = [];
                        end
                        cg.events.break_BIR_end_times(end+1) = time;
                        disp(['Back IR beam released at time ' num2str(time)]);
                    end
                end
                
                % mouse pokes in the hole (same thing as before: event is
                % considered punctual in older cage softward versions, and
                % with a finite duration in the newer versions.
                if task_started
                    % start of event
                    if (strcmp(evnt_name,'Input Event') || ...
                            strcmp(evnt_name,'Input Transition On Event')) && ...
                            strcmp(item_name,'Tray #1')
                        % record start of poke event
                        if isfield(cg.events,'poke_tray_start_times')==0
                            cg.events.poke_tray_start_times = [];
                        end
                        cg.events.poke_tray_start_times(end+1) = time;
                        disp(['Mouse poked in the tray at time ' num2str(time)]);
                    end
                    
                    % end of event (if it was recorded)
                    if strcmp(evnt_name,'Input Transition Off Event') && ...
                            strcmp(item_name,'Tray #1')
                        % record end of poke event
                        if isfield(cg.events,'poke_tray_end_times')==0
                            cg.events.poke_tray_end_times = [];
                        end
                        cg.events.poke_tray_end_times(end+1) = time;
                        disp(['Mouse left tray at time ' num2str(time)]);
                    end
                end
                
                % mouse presses on the screen
                if task_started && strcmp(evnt_name,'Touch Down Event') && ...
                        strcmp(item_name,'Bussey Mouse Operant Mode 5 x 1 x low')
                    if isfield(cg.events,'touch_down_times')==0
                        cg.events.touch_down_times = [];
                    end
                    cg.events.touch_down_times(end+1) = time;
                    if isfield(cg.events,'touch_down_position')==0
                        cg.events.touch_down_position = [];
                    end
                    cg.events.touch_down_position(end+1) = str2num(arg1_value);
                    if isfield(cg.events,'touch_down_lines')==0
                        cg.events.touch_down_lines = [];
                    end
                    cg.events.touch_down_lines(end+1) = line_number;
                    disp(['Nose pressed against the screens at position ' arg1_value ' at time ' num2str(time)]);
                end
                
                % mouse leaves the screen
                if task_started && strcmp(evnt_name,'Touch Up Event') && ...
                        strcmp(item_name,'Bussey Mouse Operant Mode 5 x 1 x low')
                    if isfield(cg.events,'touch_up_times')==0
                        cg.events.touch_up_times = [];
                    end
                    cg.events.touch_up_times(end+1) = time;
                    if isfield(cg.events,'touch_up_position')==0
                        cg.events.touch_up_position = [];
                    end
                    cg.events.touch_up_position(end+1) = str2num(arg1_value);
                    if isfield(cg.events,'touch_up_lines')==0
                        cg.events.touch_up_lines = [];
                    end
                    cg.events.touch_up_lines(end+1) = line_number;
                    disp(['Nose leaves the screen at position ' arg1_value ' at time ' num2str(time)]);
                end
                
                %  tray light turns on
                if task_started && strcmp(evnt_name,'Output On Event') && ...
                        strcmp(item_name,'TrayLight #1')
                    if isfield(cg.events,'tray_light_on_times')==0
                        cg.events.tray_light_on_times = [];
                    end
                    if isfield(cg.events,'tray_light_on_lines')==0
                        cg.events.tray_light_on_lines = [];
                    end
                    cg.events.tray_light_on_times(end+1) = time;
                    cg.events.tray_light_on_lines(end+1) = line_number;
                    disp(['Tray light turns on at time ' num2str(time)]);
                end
                
                %  tray light pulses - this has to be treated differently
                %  as in some cases, there are no light pulses in the
                %  experiment I am not sure why. So we create the light
                %  pulse fields in all cases, and in some cases, they will
                %  just be empty.
                if isfield(cg.events,'tray_light_pulse_times')==0
                    cg.events.tray_light_pulse_times = [];
                end
                if isfield(cg.events,'tray_light_pulse_lines')==0
                    cg.events.tray_light_pulse_lines = [];
                end
                if isfield(cg.events,'tray_light_pulse_durations')==0
                    cg.events.tray_light_pulse_durations = [];
                end
                if task_started && strcmp(evnt_name,'Pulse Output Event') && ...
                        strcmp(item_name,'TrayLight #1')
                    cg.events.tray_light_pulse_times(end+1) = time;
                    cg.events.tray_light_pulse_durations(end+1) = str2num(arg1_value);
                    cg.events.tray_light_pulse_lines(end+1) = line_number;
                end
                
                %  tray light turns off
                if task_started && strcmp(evnt_name,'Output Off Event') && ...
                        strcmp(item_name,'TrayLight #1')
                    if isfield(cg.events,'tray_light_off_times')==0
                        cg.events.tray_light_off_times = [];
                    end
                    cg.events.tray_light_off_times(end+1) = time;
                    if isfield(cg.events,'tray_light_off_lines')==0
                        cg.events.tray_light_off_lines = [];
                    end
                    cg.events.tray_light_off_lines(end+1) = line_number;
                    disp(['Tray light turns off at time ' num2str(time)]);
                end
                
                % sound pulses
                if task_started && strcmp(evnt_name,'Pulse Output Event') && ...
                        strcmp(item_name,'Sound_On #1')
                    if isfield(cg.events,'sound1_times')==0
                        cg.events.sound1_times = [];
                    end
                    cg.events.sound1_times(end+1) = time;
                    if isfield(cg.events,'sound1_lines')==0
                        cg.events.sound1_lines = [];
                    end
                    cg.events.sound1_lines(end+1) = line_number;
                    if isfield(cg.events,'sound1_durations')==0
                        cg.events.sound1_durations = [];
                    end
                    cg.events.sound1_durations(end+1) = str2num(arg1_value);
                    disp(['Sound 1 pulses at time ' num2str(time) ' with duration ' arg1_value]);
                end
                
                % reward collection
                if task_started && strcmp(item_name,'Reward Collected Start ITI')
                    if isfield(cg.events,'reward')==0
                        cg.events.reward = [];
                    end
                    cg.events.reward(end+1) = 1;
                    if isfield(cg.events,'reward_times')==0
                        cg.events.reward_times = [];
                    end
                    cg.events.reward_times(end+1) = time;
                    disp(['Reward collected at time ' num2str(time)]);
                end
                
                % reward witheld
                if task_started && strcmp(item_name,'Incorrect Trails to Correction Trial')
                    if isfield(cg.events,'reward')==0
                        cg.events.reward = [];
                    end
                    cg.events.reward(end+1) = -1;
                    if isfield(cg.events,'reward_times')==0
                        cg.events.reward_times = [];
                    end
                    cg.events.reward_times(end+1) = time;
                    disp(['Reward witheld at time ' num2str(time)]);
                end
                
                % incentive collection
                % if the sample is rewarded
                if task_started && strcmp(item_name,'Sample Touch rewarded')
                    if isfield(cg.events,'incentives')==0
                        cg.events.incentives = [];
                    end
                    cg.events.incentives(end+1) = 1;
                    if isfield(cg.events,'incentives_times')==0
                        cg.events.incentives_times = [];
                    end
                    cg.events.incentives_times(end+1) = time;
                    disp(['Incentive received at time ' num2str(time)]);
                end
                
                % there is no incentive
                if task_started && strcmp(item_name,'Sample Touch not rewarded')
                    if isfield(cg.events,'incentives')==0
                        cg.events.incentives = [];
                    end
                    cg.events.incentives(end+1) = -1;
                    if isfield(cg.events,'incentives_times')==0
                        cg.events.incentives_times = [];
                    end
                    cg.events.incentives_times(end+1) = time;
                    disp(['Incentive withheld at time ' num2str(time)]);
                end
                
                % measure ITI boundaries
                if task_started && strcmp(item_name,'Reward Collected Start ITI') || ...
                        strcmp(item_name,'Incorrect Trails to Correction Trial')
                    if isfield(cg.events,'ITI_start_times')==0
                        cg.events.ITI_start_times = [];
                    end
                    cg.events.ITI_start_times(end+1) = time;
                    disp(['ITI started at time ' num2str(time)]);
                end
                
                % house light on/off
                if task_started && strcmp(item_name,'HouseLight #1') && strcmp(evnt_name,'Output On Event')
                    if isfield(cg.events,'house_light_on_times')==0
                        cg.events.house_light_on_times = [];
                    end
                    cg.events.house_light_on_times(end+1) = time;
                    disp(['House light on at time ' num2str(time)]);
                end
                if task_started && strcmp(item_name,'HouseLight #1') && strcmp(evnt_name,'Output Off Event')
                    if isfield(cg.events,'house_light_off_times')==0
                        cg.events.house_light_off_times = [];
                    end
                    cg.events.house_light_off_times(end+1) = time;
                    disp(['House light off at time ' num2str(time)]);
                end
                
                % feeder
                if task_started && strcmp(item_name,'Feeder #1') && strcmp(evnt_name,'Pulse Output Event')
                    if isfield(cg.events,'feeder_times')==0
                        cg.events.feeder_times = [];
                    end
                    cg.events.feeder_times(end+1) = time;
                    if isfield(cg.events,'feeder_durations')==0
                        cg.events.feeder_durations = [];
                    end
                    cg.events.feeder_durations(end+1) = str2num(arg1_value);
                    disp(['feeder on at time ' num2str(time) ' with duration ' arg1_value]);
                end
                
            end
        end
        
        fclose(fid);
        
        % total duration of the recordings (will give the right value on the
        % last execution of the loop): Ttotal = t.
        ttot = time;
        
        %         disp('trial times')
        %         cg.events.trial_times
        %         disp('samples')
        %         cg.events.samples
        %         disp('sample times')
        %         cg.events.sample_times
        %         disp('targets')
        %         cg.events.targets
        %         disp('corrections')
        %         cg.events.correction
        %         disp('delay start times')
        %         cg.events.start_delay_times
        %         disp('delay end times')
        %         cg.events.delay_end_times
        %         disp('target times')
        %         cg.events.target_times
        %         disp('hide_image_times')
        %         cg.events.hide_image_times
        %         disp('break bir start times')
        %         cg.events.break_BIR_start_times
        %         disp('break bir end times')
        %         cg.events.break_BIR_end_times
        %         disp('poke tray times')
        %         cg.events.poke_tray_start_times
        %         disp('screen press down')
        %         cg.events.touch_down_position
        %         disp('screen press up')
        %         cg.events.touch_up_position
        %         disp('light on')
        %         cg.events.tray_light_on_times
        %         disp('light pulse')
        %         cg.events.tray_light_pulse_times
        %         disp('light off')
        %         cg.events.tray_light_off_times
        %         disp('sound times')
        %         cg.events.sound1_times
        %         disp('sound durations')
        %         cg.events.sound1_durations
        %         disp('reward')
        %         cg.events.reward
        %         disp('reward times')
        %         cg.events.reward_times
        %         disp('incentives')
        %         cg.events.incentives
        %         disp('incentive times')
        %         cg.events.incentives_times
        %         disp('start of ITI')
        %         cg.events.ITI_start_times
        %         disp('house light on times')
        %         cg.events.house_light_on_times
        %         disp('house light off times')
        %         cg.events.house_light_off_times
        %         disp('feeder times')
        %         cg.events.feeder_times
        %         disp('feeder durations')
        %         cg.events.feeder_durations

    end

% read the files for the cage, synchronization, trajectory and miniscope
    function [tr,mn,syn,cg,ttot] = readfiles(folder)
        
        % cage: read events timing and duration
        disp(' ');
        disp('Reading the content of the cage output file...')
        [cg,ttot] = read_variables(folder);
        disp('Done!')
        
        % synchronization: read the matrix
        disp(' ');
        disp('Reading the synchronization matrix with miniscope as master and behavior as slave...');
        synchronization = [];
        load([folder '/msTouchSync_new.mat']);
        syn  = synchronization;
        
        % check if the first frame is 0 or 1. In the former case, shift the
        % whole column forward by 1
        if syn.miniscopeMaster.masterFrames(1)==0
            disp('Miniscope frames started at 0. Shifted them down one unit');
            syn.miniscopeMaster.masterFrames = syn.miniscopeMaster.masterFrames + 1;
        end
        if syn.miniscopeMaster.slaveFrames(1)==0
            disp('Behavior frames started at 0. Shifted them down one unit');
            syn.miniscopeMaster.slaveFrames = syn.miniscopeMaster.slaveFrames + 1;
        end

        disp('Done!');
        
        % trajectory: read the vectors (in lines and columns)
        disp('Reading the mouse trajectory...');
        position = [];
        load([folder '/mouse_traj_corrected.mat'],'position');
        % We need the original values to lay the mouse position onto the
        % video, and also a normalized version that we can use for the rest
        % of the analysis:
        tr = [];
        tr.x = position.bodypart1.x;
        tr.y = position.bodypart1.y;
        tr.norm.x = -1+2*(tr.x-min(tr.x))/(max(tr.x)-min(tr.x));
        tr.norm.y = (tr.y-min(tr.y))/(max(tr.y)-min(tr.y));
        disp('Done!');
        
        % trajectory: read the mask for the bad data (when the mouse
        % position was impossible to determine)
        disp('Reading the mouse trajectory mask...');
        tr.bad_data = position.bodypart1.bad_data;
        disp('Done!');
        
        % read the deconvolved calcium signal
        disp('Reading the deconvolved calcium signal...');
        ms = [];
        
        load([folder '/Miniscope_2/msDeconvolved.mat']);

        if isfield(ms,'deconvolution')
            mn.original = ms.deconvolution.deconvolvedSig;
        else
            mn.original = ms.deconvolvedSig;
        end
        
        % add the raw and filtered traces if they are present in the ms
        % file.
        mn.raw = ms.rawTraces;
        mn.filt = ms.filtTraces;
        
        clear('ms','environment','calcium','processed','properties');
        disp('Done!');
        
    end


% synchronize the datastreams ---------------------------------------------

% function to convert event time list into a time series synchronized with
% the miniscope. Times in my_events can be larger than the largest
% miniscope time value, making the same frame repeat itself at the end.
% This is removed last.
    function my_ts = event2timeseries(my_event,sync_matrix)
        my_ts = zeros(size(my_event));
        for i1=1:length(my_ts)
            temp = my_event(i1);
            my_dist = abs(sync_matrix.miniscopeMaster.clean_cageTimes-temp);
            [mmin,pos] = min(my_dist);
            my_ts(i1) = sync_matrix.miniscopeMaster.clean_masterFrames(pos);
        end
        
        % check that the number of entries equals that of results
        if length(my_ts)~=length(my_event)
            disp('Original and synchronized time series do not have the same length... ');
            pause(0.1) 
        end
    end

% function to display repeating elements in a series of numbers
    function display_repeats(x)
        uniqueitems = unique(x);
        for i1=1:length(uniqueitems)
            nboccur = length(find(x==uniqueitems(i1)));
            if nboccur>1
                disp(['Time ' num2str(uniqueitems(i1)) ' appeared ' num2str(nboccur) ' times.']);
            end
        end
    end

% function to translate the list of times created above into vectors coding for
% each event, condition, etc. in frames.
    function [cg,nnew] = create_variables(cg,syn,tr,ntrials,fol)
        
        % number of time points of our data
        nnew = syn.miniscopeMaster.clean_N;
        
        % trial number (using trial_times) --------------------------------
        disp('Making trial number time series');
        
        % register trial beginnings and ends to miniscope frames.
        cg.timeseries.trial_times_dt = event2timeseries(cg.events.trial_times,syn);
        
        % Check if there are duplicate entries
        if isequal(cg.timeseries.trial_times_dt,unique(cg.timeseries.trial_times_dt))==0
            disp('Problem in registration...')
            display_repeats(cg.events.trial_times);
            pause(0.1)
        end
        % Make a time series with current trial value at each frame
        cg.timeseries.trial_number_ts = zeros(nnew,1);
        for i1=1:length(cg.timeseries.trial_times_dt)-1
            temp1 = cg.timeseries.trial_times_dt(i1);
            temp2 = cg.timeseries.trial_times_dt(i1+1);
            cg.timeseries.trial_number_ts(temp1:temp2-1) = i1;
        end
        % and the last point
        cg.timeseries.trial_number_ts(nnew) = ntrials;

        % trial type (using samples) ----------------------------------------------
        disp('Making trial type time series')
        cg.timeseries.trial_type_ts = zeros(nnew,1);
        for i1=1:ntrials
            ttime = find(cg.timeseries.trial_number_ts==i1);
            cg.timeseries.trial_type_ts(ttime) = cg.events.samples(i1);
        end
        
        % -------------------------------------------------------------------------
        % make sample and target/distractor ts (using sample_times, target_times and
        % hide_image_times):
        disp('Making sample and target time series');
        cg.timeseries.sample_times_dt = event2timeseries(cg.events.sample_times,syn);
        % Check if there are duplicate entries
        if isequal(cg.timeseries.sample_times_dt,unique(cg.timeseries.sample_times_dt))==0
            disp('Problem in registration...')
            display_repeats(cg.events.sample_times);
            pause(0.1)
        end
        cg.timeseries.target_times_dt = event2timeseries(cg.events.target_times,syn);
        % Check if there are duplicate entries
        if isequal(cg.timeseries.target_times_dt,unique(cg.timeseries.target_times_dt))==0
            disp('Problem in registration...');
            display_repeats(cg.events.target_times);
            pause(0.1)
        end
        cg.timeseries.hide_image_times_dt = event2timeseries(cg.events.hide_image_times,syn);
        % Check if there are duplicate entries
        if isequal(cg.timeseries.hide_image_times_dt,unique(cg.timeseries.hide_image_times_dt))==0
            disp('Problem in registration...')
            display_repeats(cg.events.hide_image_times);
            pause(0.1)
        end
        
        % Now, check if the frames when sample or target appear overlap
        % with frames when images disappear.
        if not(isempty(intersect(cg.timeseries.sample_times_dt,cg.timeseries.hide_image_times_dt)))
            disp('Overlap with sample on/off');
            pause(0.1)
        end
        if not(isempty(intersect(cg.timeseries.target_times_dt,cg.timeseries.hide_image_times_dt)))
            disp('Overlap with target on/off');
            pause(0.1)
        end
        
        % now, use these to construct time series that is 1/0 when the
        % sample(distractor) is presented/hidden.
        cg.timeseries.sample_ts = zeros(nnew,1);
        cg.timeseries.target_ts = zeros(nnew,1);
        
        % 0: after choice or before sample, 1: in sample, 2: delay before target, 3: in target
        in_stimulus = 0;
        for t=1:nnew
            % in the sample
            if in_stimulus==0 && any(t==cg.timeseries.sample_times_dt)
                in_stimulus = 1;
            end
            % in the delay
            if in_stimulus==1 && any(t==cg.timeseries.hide_image_times_dt)
                in_stimulus = 2;
            end
            if in_stimulus==1
                cg.timeseries.sample_ts(t) = 1;
            end
            
            % in the target/distractor
            if in_stimulus==2 && any(t==cg.timeseries.target_times_dt)
                in_stimulus = 3;
            end
            % after the choice period
            if in_stimulus==3 && any(t==cg.timeseries.hide_image_times_dt)
                in_stimulus = 0;
            end
            if in_stimulus==3
                cg.timeseries.target_ts(t) = 1;
            end
        end
        
%         hold on
%         plot(cg.timeseries.sample_ts)
%         plot(cg.timeseries.target_ts)
%         pause
        
        % BIR break -----------------------------------------------------
        % if cage.events.break_BIR_end_times is empty,then events are
        % punctual, otherwise, they have a finite duration.
        % register times to cage times. Again, here we are not too precise
        % as we use the light as task period separator, not the bir break.
        disp('Making BIR break time series');
        disp('bir break start')
        cg.timeseries.break_BIR_start_times_dt = event2timeseries(cg.events.break_BIR_start_times,syn);
        % Check if there are duplicate entries
        if isequal(cg.timeseries.break_BIR_start_times_dt,unique(cg.timeseries.break_BIR_start_times_dt))==0
            disp('Problem in frame registration...');
            display_repeats(cg.timeseries.break_BIR_start_times_dt);
            disp('Removing the duplicate values and continuing...');
            cg.timeseries.break_BIR_start_times_dt = unique(cg.timeseries.break_BIR_start_times_dt);
        end
        % event with finite duration
        if isfield(cg.events,'break_BIR_end_times')
            disp('bir break end');
            cg.timeseries.break_BIR_end_times_dt = event2timeseries(cg.events.break_BIR_end_times,syn);
            
            % Check if there are duplicate entries
            if isequal(cg.timeseries.break_BIR_end_times_dt,unique(cg.timeseries.break_BIR_end_times_dt))==0
                disp('Problem in frame registration...');
                display_repeats(cg.timeseries.break_BIR_end_times_dt);
                disp('Removing the duplicate values and continuing...');
                pause(0.1)
                cg.timeseries.break_BIR_end_times_dt = unique(cg.timeseries.break_BIR_end_times_dt);
            end
            
            % check for simultaneous on/off
            if not(isempty(intersect(cg.timeseries.break_BIR_end_times_dt,cg.timeseries.break_BIR_start_times_dt)))
                disp('Overlap with bir break on/off');
                pause(0.1)
            end
        end
        
        % Encode these in a time series.
        cg.timeseries.break_BIR_ts = zeros(nnew,1);
        if isfield(cg.events,'break_BIR_end_times')==0
            % punctual
            cg.timeseries.break_BIR_ts(cg.timeseries.break_BIR_start_times_dt) = 1;
        else
            % finite duration
            in_BIRbreak = 0;
            for t=1:nnew
                if in_BIRbreak==0 && any(t==cg.timeseries.break_BIR_start_times_dt)
                    in_BIRbreak = 1;
                end
                if in_BIRbreak==1 && any(t==cg.timeseries.break_BIR_end_times_dt)
                    in_BIRbreak = 0;
                end
                if in_BIRbreak==1
                    cg.timeseries.break_BIR_ts(t) = 1;
                end
            end
        end
        
        
        % nose poke -----------------------------------------------------
        % if cage.events.poke_tray_end_times is empty,then events are
        % punctual, otherwise, they have a finite duration.
        disp('Making nose poke time series');
        disp('nose poke start')
        cg.timeseries.poke_tray_start_times_dt = event2timeseries(cg.events.poke_tray_start_times,syn);
        % Check if there are duplicate entries
        if isequal(cg.timeseries.poke_tray_start_times_dt,unique(cg.timeseries.poke_tray_start_times_dt))==0
            disp('Problem in registration...');
            display_repeats(cg.events.poke_tray_start_times);
            pause(0.1)
        end
        % event with finite duration
        if isfield(cg.events,'poke_tray_end_times')
            disp('nose poke end');
            cg.timeseries.poke_tray_end_times_dt = event2timeseries(cg.events.poke_tray_end_times,syn);
            
            % Check if there are duplicate entries
            if isequal(cg.timeseries.poke_tray_end_times_dt,unique(cg.timeseries.poke_tray_end_times_dt))==0
                disp('Problem in registration...');
                display_repeats(cg.events.poke_tray_end_times);
                pause(0.1)
            end
            
            % check for simultaneous on/off
            if not(isempty(intersect(cg.timeseries.poke_tray_end_times_dt,cg.timeseries.poke_tray_start_times_dt)))
                disp('>>>>>>>>>>>>>> Overlap with nose poke on/off');
                disp('I keep going...');
            end
        end
        
        % Encode these in a time series.
        cg.timeseries.poke_tray_ts = zeros(nnew,1);
        if isfield(cg.events,'poke_tray_end_times')==0
            % punctual
            cg.timeseries.poke_tray_ts(cg.timeseries.poke_tray_start_times_dt) = 1;
        else
            % finite duration
            in_nose_poke = 0;
            for t=1:nnew
                if in_nose_poke==0 && any(t==cg.timeseries.poke_tray_start_times_dt)
                    in_nose_poke = 1;
                end
                if in_nose_poke==1 && any(t==cg.timeseries.poke_tray_end_times_dt)
                    in_nose_poke = 0;
                end
                if in_nose_poke==1
                    cg.timeseries.poke_tray_ts(t) = 1;
                end
            end
        end
        
        % reward -----------------------------------------------------
        disp('Making reward time series');
        % register times to cage time.
        cg.timeseries.reward_times_dt = event2timeseries(cg.events.reward_times,syn);
        % Check if there are duplicate entries
        if isequal(cg.timeseries.reward_times_dt,unique(cg.timeseries.reward_times_dt))==0
            disp('Problem in registration...')
            display_repeats(cg.events.reward_times);
            pause(0.1)
        end
        % Convert into a time series.
        cg.timeseries.reward_ts = zeros(nnew,1);
        cg.timeseries.reward_ts(cg.timeseries.reward_times_dt) = cg.events.reward(1:ntrials);
        % reward duration will be extracted using the mouse trajectory
        % later on.
        
        % incentives -----------------------------------------------------
        disp('Making incentives time series');
        % register times to cage time.
        cg.timeseries.incentives_times_dt = event2timeseries(cg.events.incentives_times,syn);
        % Check if there are duplicate entries
        if isequal(cg.timeseries.incentives_times_dt,unique(cg.timeseries.incentives_times_dt))==0
            disp('Problem in registration...');
            display_repeats(cg.events.incentives_times);
            pause(0.1)
        end
        % Convert to a time series.
        cg.timeseries.incentives_ts = zeros(nnew,1);
        cg.timeseries.incentives_ts(cg.timeseries.incentives_times_dt) = cg.events.incentives(1:ntrials);
        % (their duration will be set using the mouse trajectory if need be).
        
        % screen contact -----------------------------------------------------
        disp('Screen contact');
        % We will make separate time series for touch down and touch up, so
        % there is nothing too complicated (like for sound and tray light) to do here...
        % Also, no need to be too precise (like we were for the tray light)
        % as this is not specify task phases. So we can just remove events
        % that fall in identical frames.
        disp('touch');
        % Check if there are duplicate entries
        if isequal(cg.events.touch_down_times,unique(cg.events.touch_down_times))==0
            disp('Problem in native time registration...');
            display_repeats(cg.events.touch_down_times);
            disp('Removing the duplicate values and continuing...');
            cg.events.touch_down_times = unique(cg.events.touch_down_times);
        end
        % Register time to cage time
        cg.timeseries.touch_down_times_dt = event2timeseries(cg.events.touch_down_times,syn);
        % Check if there are duplicate entries
        if isequal(cg.timeseries.touch_down_times_dt,unique(cg.timeseries.touch_down_times_dt))==0
            disp('Problem in registered time registration...');
            display_repeats(cg.timeseries.touch_down_times_dt);
            disp('Removing the duplicate values and continuing...');
            cg.timeseries.touch_down_times_dt = unique(cg.timeseries.touch_down_times_dt);
        end
        % Convert to time series
        cg.timeseries.touch_down_ts = zeros(nnew,1);
        cg.timeseries.touch_down_ts(cg.timeseries.touch_down_times_dt) = cg.events.touch_down_position(1:length(cg.timeseries.touch_down_times_dt));
        
        disp('release')
        % Check if there are duplicate entries
        if isequal(cg.events.touch_up_times,unique(cg.events.touch_up_times))==0
            disp('Problem in native time registration...');
            display_repeats(cg.events.touch_up_times);
            disp('Removing the duplicate values and continuing...');
            cg.events.touch_up_times = unique(cg.events.touch_up_times);
        end
        % Register time to cage time
        cg.timeseries.touch_up_times_dt = event2timeseries(cg.events.touch_up_times,syn);
        % Check if there are duplicate entries
        if isequal(cg.timeseries.touch_up_times_dt,unique(cg.timeseries.touch_up_times_dt))==0
            disp('Problem in registered time registration...');
            display_repeats(cg.timeseries.touch_up_times_dt);
            disp('Removing the duplicate values and continuing...');
            cg.timeseries.touch_up_times_dt = unique(cg.timeseries.touch_up_times_dt);
        end
        % Convert to time series
        cg.timeseries.touch_up_ts = zeros(nnew,1);
        cg.timeseries.touch_up_ts(cg.timeseries.touch_up_times_dt) = cg.events.touch_up_position(1:length(cg.timeseries.touch_up_times_dt));
        
        % sound -----------------------------------------------------
        disp('sound')
        
        % bringing together all the times when the sound is on.
        % sound_on_off is stuctured as follows:
        % line number, 1/-1, time and frame after discretization.
        sound_on_off = [];
        
        % sound on
        n = length(cg.events.sound1_times);
        sound_on_off = [sound_on_off;cg.events.sound1_lines' ones(n,1) cg.events.sound1_times' cg.events.sound1_durations'];
        
        % now discretize the times in frames
        temp = event2timeseries(sound_on_off(:,3),syn);
        sound_on_off(:,5) = temp;
        
        % remove events that take place at the same frame (and potentially
        % at same time)
        for ii1=2:size(sound_on_off,1)
            if sound_on_off(ii1,5)==sound_on_off(ii1-1,5)
                sound_on_off(ii1,:) = NaN;
            end
        end
        no_good = isnan(sound_on_off(:,5));
        sound_on_off(no_good,:) = [];
        
        % Next, we need to introduce lines corresponding to when the pulse
        % turns itself off.
        temp = [];
        for ii1=1:size(sound_on_off,1)
            % add a line with 0's
            temp = [temp; sound_on_off(ii1,:);0 0 0 0 0];
        end
        sound_on_off = temp;
        
        % add sound pulse off signal (-1) where there are 0's.
        for ii1=1:size(sound_on_off,1)
            if isequal(sound_on_off(ii1,:),zeros(1,5))
                % fill in the line: line number remains 0
                % 2n column: -1 (turn off line)
                sound_on_off(ii1,2) = -1;
                % 3rd column: time when pulse was turned on (time on
                % previous line) plus specific duration
                sound_on_off(ii1,3) = sound_on_off(ii1-1,3) + sound_on_off(ii1-1,4);
                % 4th column: duration - no need to write anything
                % 5th column: convert to frame
                sound_on_off(ii1,5) = event2timeseries(sound_on_off(ii1,3),syn);
            end
        end
        
        % Due to various reasons, several events can fall onto the same
        % frame. One of them is that the beginning and ending of 0.077ms pulses
        % can be registered to a single frame (so we have a pair (1,-1) on
        % the same frame) if the acquisition rate of the miniscope drops below 1/0.77
        % = about 15 fps. And to make things even simpler, sometimes, two such pulses
        % can land next to each other (e.g. frames n n n+1 n+1) and we need to correct
        % all four...
        for index=2:size(sound_on_off,1)
            
            % case 1: pair to fix are the last two items. Subtract 1 to
            % frame of the first
            if index==size(sound_on_off,1)
                if sound_on_off(index-1,5)==sound_on_off(index,5)
                    disp('Fixed last pulse')
                    sound_on_off(index-1,5) = sound_on_off(index-1,5) - 1;
                end
            end
            
            % case 1': pair is the first two elements. add 1 to frame of
            % second element
            if index==2
                if sound_on_off(1,5)==sound_on_off(2,5)
                    disp('Fixed fist pulse')
                    sound_on_off(2,5) = sound_on_off(1,5) + 1;
                end
            end
            
            if index>3
                
                % case 2: pair between top and next to last item: add 1 to
                % frame of the second
                if sound_on_off(index-1,5)==sound_on_off(index,5) && ...
                        sound_on_off(index+1,5)~=(sound_on_off(index,5)+1) && ...
                        sound_on_off(index-2,5)~=(sound_on_off(index,5)-1)
                    disp('Fixed single pulse in the middle of the array');
                    sound_on_off(index,5) = sound_on_off(index,5) + 1;
                end
                
                % case 3: two consecutive pairs in the middle of the array.
                % Add 0 1 2 3 to their frame number
                if (sound_on_off(index-3,5)==sound_on_off(index-2,5)) && (sound_on_off(index-1,5)==sound_on_off(index,5))
                    disp('Fixed two consecutive pulses')
                    sound_on_off(index-2,5) = sound_on_off(index-2,5) + 1;
                    sound_on_off(index-1,5) = sound_on_off(index-2,5) + 1;
                    sound_on_off(index,5) = sound_on_off(index-2,5) + 2;
                end
                
                % case 4: two different frames for +1 are matched to an
                % identical frame for -1. Just shift down the last one
                if (sound_on_off(index-3,5)~=sound_on_off(index-1,5)) && (sound_on_off(index-2,5)==sound_on_off(index,5))
                    disp('Fixed two interlaced pulses')
                    sound_on_off(index,5) = sound_on_off(index,5) + 1;
                end
                
            end
            
        end
        
        %         % visual check
        %         for ii1=1:size(sound_on_off,1)
        %             disp(num2str(sound_on_off(ii1,:)));
        %         end
        
        % now, make the time series for the sound (1 = on, 0 = off)
        cg.timeseries.sound_ts = zeros(nnew,1);
        for fr=1:nnew
            
            % look if the current frame is part of light_on_off (i.e. if
            % light changes state during that frame.
            temp = find(fr==sound_on_off(:,5));
            
            if length(temp)==0
                % nope: sound stays as is
                cg.timeseries.sound_ts(fr) = cg.timeseries.sound_ts(fr-1);
            else
                if length(temp)==1
                    % sound changes
                    if sound_on_off(temp,2)<0
                        % sound is turned off
                        cg.timeseries.sound_ts(fr) = 0;
                    end
                    if sound_on_off(temp,2)>0
                        % sound is turned on
                        cg.timeseries.sound_ts(fr) = 1;
                    end
                else
                    % something else is happening
                    disp('>>>>>>>>>>>> Problem when computing the sound activity');
                    disp('I keep going...');
                end
            end
        end
        
        % tray light/pulse -----------------------------------------------------
        % This is important to get right as it defines the various parts of
        % the delay amont other things. Note that the light pulse at the
        % beginning of the delay last at most 2s, so the mouse does not
        % have to break the IR beam to turn it off.
        disp('Tray light')
        
        % bringing together all the times when the light state changes.
        % light_on_off is stuctured as follows:
        % line number, 1/-1/2/-2, time and frame after discretization.
        light_on_off = [];
        
        % light on
        n = length(cg.events.tray_light_on_times);
        light_on_off = [light_on_off;cg.events.tray_light_on_lines' ones(n,1) cg.events.tray_light_on_times'];
        
        % light off
        n = length(cg.events.tray_light_off_times);
        light_on_off = [light_on_off;cg.events.tray_light_off_lines' -1*ones(n,1) cg.events.tray_light_off_times'];
        
        % pulse on
        n = length(cg.events.tray_light_pulse_times);
        light_on_off = [light_on_off;cg.events.tray_light_pulse_lines' 2*ones(n,1) cg.events.tray_light_pulse_times' ];
        
        % sort it in increasing line number
        [u,v] = sort(light_on_off(:,1),'ascend');
        slight_on_off = light_on_off(v,:);
        
        % now discretize the times in frames
        temp = event2timeseries(slight_on_off(:,3),syn);
        slight_on_off(:,4) = temp;
        
        % Next, we need to introduce lines corresponding to when the pulse
        % turns itself off.
        % Look for instances where +2 is not followed by -1 less than 2
        % seconds after (= when the pulse turns itself out by itself). If
        % so, add a line with the off signal (-2)
        temp = [];
        for ii1=1:size(slight_on_off,1)-1
            if slight_on_off(ii1,2)==2 && slight_on_off(ii1+1,2)==1
                % add a line with 0's
                %disp('Added line with 0s');
                temp = [temp; slight_on_off(ii1,:);0 0 0 0];
            else
                % copy the current line
                temp = [temp; slight_on_off(ii1,:)];
            end
        end
        temp = [temp; slight_on_off(end,:)];
        slight_on_off = temp;
        
        % add light pulse off signal (-2) where there are 0's.
        for ii1=1:size(slight_on_off,1)
            if isequal(slight_on_off(ii1,:),zeros(1,4))
                % fill in the line: line number remains 0
                % 2n column: -2 (turn off line)
                slight_on_off(ii1,2) = -2;
                % 3rd column: time when pulse was turned on (time on previous line) plus 2 seconds
                slight_on_off(ii1,3) = slight_on_off(ii1-1,3) + 2;
                % 4th column: convert to frame
                slight_on_off(ii1,4) = event2timeseries(slight_on_off(ii1,3),syn);
            end
        end
        
        % Due to various reasons, several events can fall onto the same
        % frame. However, as we used the original cage order for events,
        % we can just move the second event of the pair to the following frame
        % as we still respect the real order of events.
        for index=2:size(slight_on_off,1)
            if slight_on_off(index-1,4) == slight_on_off(index,4)
                % move the second one down one frame if it is not already
                % used
                new_index = slight_on_off(index,4) + 1;
                if not(ismember(new_index,slight_on_off(:,4)))
                    disp(['Disambiguated frame ' num2str(slight_on_off(index,:))]);
                    slight_on_off(index,4) = new_index;
                else
                    disp('Disambiguation problem...')
                    pause(0.1)
                end
            end
        end
        
        % check if there remain any problems
        if length(unique(slight_on_off(:,4)))~=length(slight_on_off(:,4))
            disp('Duplicate frame on/off/pulse on/pulse off entries in cage schedule...');
            display_repeats(slight_on_off(:,4));
            pause(0.1)
        end
        
        % If we reach this point, all should be clean now. :)
        if length(find(abs(slight_on_off(:,2))==2))==0
            disp('>>>>>>>>> There are no light pulses in the session.');
        end
        
        % now, make the time series for the light (1 = on, 0 = off)
        cg.timeseries.light_ts = zeros(nnew,1);
        cg.timeseries.light_ts(1) = 0;
        for fr=1:nnew
            
            % look if the current frame is part of light_on_off (i.e. if
            % light changes state during that frame.
            temp = find(fr==slight_on_off(:,4));
            
            if length(temp)==0
                % nope: light stays as is
                if fr>1
                    cg.timeseries.light_ts(fr) = cg.timeseries.light_ts(fr-1);
                end
            else
                if length(temp)==1
                    % light changes
                    if slight_on_off(temp,2)<0
                        % light is turned off
                        cg.timeseries.light_ts(fr) = 0;
                    end
                    if slight_on_off(temp,2)>0
                        % light is turned on
                        cg.timeseries.light_ts(fr) = 1;
                    end
                else
                    % something else is happening
                    disp('Problem when computing the light activity');
                    temp
                    slight_on_off(temp)
                    
                    pause(0.1)
                end
            end
        end
        
        %         % visual check
        %         for ii1=1:size(slight_on_off,1)
        %             disp([num2str(slight_on_off(ii1,1)) ' ' num2str(slight_on_off(ii1,2)) ' ' num2str(slight_on_off(ii1,3)) ' ' num2str(slight_on_off(ii1,4))]);
        %         end
        %
        %         figure;
        %         plot(cg.timeseries.light_ts);
        %         hold on;
        %         plot(cg.timeseries.trial_number_ts/30);
        %         plot(cg.timeseries.sample_ts/2,'*');
        %         plot(cg.timeseries.target_ts/2,'s');
        %         plot(cg.timeseries.reward_ts/2,'^');
        %         plot(cg.timeseries.sound_ts/2,'o');
        %         pause(0.1)
        
        
        % feeder ------------------------------------------------------------------
        disp('Feeder');
        % Register time to cage time
        cg.timeseries.feeder_times_dt = event2timeseries(cg.events.feeder_times,syn);
        % Check if there are duplicate entries
        if isequal(cg.timeseries.feeder_times_dt,unique(cg.timeseries.feeder_times_dt))==0
            disp('Problem in registration...');
            display_repeats(cg.events.feeder_times);
            pause(0.1)
        end
        % and build the times series
        cg.timeseries.feeder_ts = zeros(nnew,1);
        for i1=1:length(cg.timeseries.feeder_times_dt)
            width = round(cg.events.feeder_durations(i1)/dt);
            cg.timeseries.feeder_ts(cg.timeseries.feeder_times_dt(i1)+(0:width)) = 1;
        end
        
        % house light -------------------------------------------------------------
        disp('house light');
        disp('on')
        % Register time to cage time
        cg.timeseries.house_light_on_times_dt = event2timeseries(cg.events.house_light_on_times,syn);
        % Check if there are duplicate entries
        if isequal(cg.timeseries.house_light_on_times_dt,unique(cg.timeseries.house_light_on_times_dt))==0
            disp('Problem in registration...');
            display_repeats(cg.events.house_light_on_times);
            pause(0.1)
        end
        disp('off')
        % Register time to cage time
        cg.timeseries.house_light_off_times_dt = event2timeseries(cg.events.house_light_off_times,syn);
        % Check if there are duplicate entries
        if isequal(cg.timeseries.house_light_off_times_dt,unique(cg.timeseries.house_light_off_times_dt))==0
            disp('Problem in registration...');
            display_repeats(cg.events.house_light_off_times);
            pause(0.1)
        end
        % and build the time series
        cg.timeseries.houselight_ts = zeros(nnew,1);
        light_on = 0;
        for t=1:nnew
            if light_on==0 && any(t==cg.timeseries.house_light_on_times_dt)
                light_on = 1;
            end
            if light_on==1 && any(t==cg.timeseries.house_light_off_times_dt)
                light_on = 0;
            end
            cg.timeseries.houselight_ts(t) = light_on;
        end
        
        % task epochs -------------------------------------------------------------
        disp('Making epoch time series');
        cg.timeseries.epochs_ts = zeros(nnew,1);
        
        % = 2 during sample
        cg.timeseries.epochs_ts(logical(cg.timeseries.sample_ts)) = 2;
        
        % = 6 during choice
        cg.timeseries.epochs_ts(logical(cg.timeseries.target_ts)) = 6;

        % Crtierion 1: check that all trials have a least a cue and choice period. If
        % not (potentially because of large fluctuation in the frame rate),
        % just set cg.timeseries.epochs_ts to NaN for that trial(s).
        bad_trials = [];
        for trial=1:ntrials
            ttime = find(cg.timeseries.trial_number_ts==trial);
            nbsample = length(find(cg.timeseries.epochs_ts(ttime)==2));
            nbtarget = length(find(cg.timeseries.epochs_ts(ttime)==6));
            if nbsample==0 || nbtarget==0
               bad_trials = [bad_trials trial]; 
               disp(['Trial ' num2str(trial) ' was removed as we could not locate the cue or choice period.']);
            end
        end
        
        % Criterion 2: Another criterion is whether the mean frame rate in each trial is
        % high enough. It should be on average 30fps (or 1/30 per frame).
        for trial=1:ntrials
            ttime = find(cg.timeseries.trial_number_ts==trial);
            
            if isempty(ttime)
                % empty trial (can happen with enough undersampling)
                bad_trials = [bad_trials trial];
                disp(['Trial ' num2str(trial) ' was removed because the trial was empty.']);
            else
                mean_duration = (syn.miniscopeMaster.cageTimes(ttime(end)) - ...
                    syn.miniscopeMaster.cageTimes(ttime(1)))/length(ttime);
                if mean_duration>=0.1
                    % low frame rate
                    bad_trials = [bad_trials trial];
                    disp(['Trial ' num2str(trial) ' was removed because the frame rate was too low.']);
                end
            end
        end
        
        % remove repetitions
        bad_trials = unique(bad_trials);
        
        % assign 1,3,5 and 7
        % 1 is just any time interval before the cue in the trial
        for trial=1:ntrials
            ttime = find(cg.timeseries.trial_number_ts==trial);
            tepochs = cg.timeseries.epochs_ts(ttime);
            % look for indices to the left of the sample in the current
            % trial
            if not(ismember(trial,bad_trials))
                find2s = find(tepochs==2);
                if find2s(1)>1
                    cg.timeseries.epochs_ts(ttime(1:find2s(1)-1)) = 1;
                end
            end
        end
        
        % partition the light activity into blocks
        for trial=1:ntrials
            
            ttime = find(cg.timeseries.trial_number_ts==trial);
            trial_light = cg.timeseries.light_ts(ttime);
            
            if not(ismember(trial,bad_trials))
                
                % divide the light activation in consecutive blocks of 1's
                % separated by 0's.
                [labeled_blocks,numRegions] = bwlabel(trial_light,4);
                
                % reward
                my_reward = cg.events.reward(trial);
                
                % add 10 to all blocks
                for t=1:length(ttime)
                    if labeled_blocks(t)>0
                        labeled_blocks(t) = labeled_blocks(t) + 10;
                    end
                end
                
                % label each block according to the epoch it corresponds to: 1, 3, 5
                % and 7.
                % Note that the number of light tray block depends on the
                % trial number and the reward:
                % trial 1, R = 3 blocks
                % trial 1, no R = 2 blocks
                % trial>1, R = 4 blocks
                % trial>1, no R = 3 blocks.
                % Note, also that the first light pulse in the delay lasts at
                % most 2s: even if the mouse does not break the BIR, the light
                % will go off and d2 will start. So this case has also been
                % taken into account. :)
                % To take into account the various cases, I am using a tree
                % approach with R value and trial value: see below
                for n=1:numRegions
                    my_epoch = 0;
                    % tree starts here
                    if my_reward==1
                        if trial==1
                            % 1,R
                            switch n
                                case 1
                                    % block 1 = 3
                                    my_epoch = 3;
                                case 2
                                    % block 2 = 5
                                    my_epoch = 5;
                                case 3
                                    % block 3 = 7
                                    my_epoch = 7;
                            end
                        else
                            % >1,R
                            switch n
                                case 1
                                    % block 1 = 1
                                    my_epoch = 1;
                                case 2
                                    % block 2 = 3
                                    my_epoch = 3;
                                case 3
                                    % block 3 = 5
                                    my_epoch = 5;
                                case 4
                                    % block 4 = 7
                                    my_epoch = 7;
                            end
                        end
                        
                    else
                        if trial==1
                            % 1,no R
                            switch n
                                case 1
                                    % block 1 = 3
                                    my_epoch = 3;
                                case 2
                                    % block 2 = 5
                                    my_epoch = 5;
                            end
                        else
                            % >1, no R
                            switch n
                                case 1
                                    % block 1 = 1
                                    my_epoch = 1;
                                case 2
                                    % block 2 = 3
                                    my_epoch = 3;
                                case 3
                                    % block 3 = 5
                                    my_epoch = 5;
                            end
                        end
                        
                    end
                    % then put the correct trial period in the blocks that are
                    % present (a bit cumbersome but clear and hopefully works
                    % in all cases).
                    labeled_blocks(labeled_blocks==(10+n)) = my_epoch;
                end
                % put region numbers in epochs_ts.
                temp = find(labeled_blocks>0);
                temp1 = cg.timeseries.epochs_ts(ttime);
                temp1(temp) = labeled_blocks(temp);
                cg.timeseries.epochs_ts(ttime) = temp1;
                
            end
        end
        
        % Next up we add period 4 (between periods 3 and 5, duh!)
        for trial=1:ntrials
            % locate the end of 3 (=d1) and start of 5 (=d3)
            ttime = find(cg.timeseries.trial_number_ts==trial);
            tepochs = cg.timeseries.epochs_ts(ttime);
            if not(ismember(trial,bad_trials))
                % find where 3s and 5s are
                pos3 = find(tepochs==3);
                pos5 = find(tepochs==5);
                tepochs(pos3(end)+1:pos5(1)-1) = 4;
                % make the change
                cg.timeseries.epochs_ts(ttime) = tepochs;
            end
        end
        
        % reward (period 6, only if R=+1). We use reward_ts to detect when
        % the reward was contacted.
        near_reward = 0;
        reward_pt = [];
        my_dist = 1000000;
        
        for t=1:nnew
            
            trial = cg.timeseries.trial_number_ts(t);
            if not(ismember(trial,bad_trials))
                
                % corresponding frame in behavioral time
                bh_frame = syn.miniscopeMaster.slaveFrames(t);
                
                % wait until reward has been dispensed
                if near_reward==0 && cg.timeseries.reward_ts(t)==1
                    % keep in memory where the mouse was at that instant
                    reward_pt = [tr.norm.x(bh_frame) tr.norm.y(bh_frame)];
                    near_reward = 1;
                end
                
                % compute the distance of the mouse to the reward point
                if near_reward==1 && not(isempty(reward_pt))
                    my_dist = reward_pt - [tr.norm.x(bh_frame) tr.norm.y(bh_frame)];
                    my_dist = sqrt(my_dist*my_dist');
                end
                
                % end the reward period when the mouse moves away from the reward point
                if near_reward==1 && my_dist>R_dist_thresh
                    near_reward = 0;
                    reward_pt = [];
                    my_dist = 1000000;
                end
                
                % if mouse stays near the reward point, then epoch = 8
                if near_reward==1
                    cg.timeseries.epochs_ts(t) = 8;
                end
                
            end
        end
        
        % the rest of the trial is filled with the ITI period: should simply replace the 0s that are left.
        cg.timeseries.epochs_ts(find(cg.timeseries.epochs_ts==0)) = 9;
        
%         % visual illustration
%         subplot(3,1,1)
%         plot(cg.timeseries.epochs_ts,'-*')
%         xlabel('epochs (1 = wait, 2 = cue, etc.')
%         xlim([1 length(cg.timeseries.epochs_ts)])
%         temp = [0 diff(cg.timeseries.trial_number_ts)'];
%         trbnd = find(temp>0);
%         for l1=1:length(trbnd)
%            line(trbnd(l1)*[1 1],[0 10],'Color','r') 
%         end
%         
%         subplot(3,1,2)
%         plot(syn.miniscopeMaster.cageTimes,'-*')
%         xlabel('Frame instants (s)')
%         xlim([1 length(cg.timeseries.epochs_ts)])
%         for l1=1:length(trbnd)
%            line(trbnd(l1)*[1 1],max(syn.miniscopeMaster.cageTimes)*[0 1],'Color','r') 
%         end
%        
%         subplot(3,1,3)
%         plot(1./diff(syn.miniscopeMaster.cageTimes),'-*')
%         xlabel('frame rate (Hz)')
%         xlim([1 length(cg.timeseries.epochs_ts)])
%         for l1=1:length(trbnd)
%            line(trbnd(l1)*[1 1],[0 max(1./(diff(syn.miniscopeMaster.cageTimes)))],'Color','r') 
%         end
%         
%         pause

        % More triage is needed: 
        % Criterion 3: remove trials where there are not all
        % required phases, and if they are not in the proper order
        for trial=1:ntrials
            
            ttime = find(cg.timeseries.trial_number_ts==trial);
            temp = cg.timeseries.epochs_ts(ttime);
            if not(isempty(temp))
                
                % remove repetition but without changing the order
                stemp = temp(1);
                for i1=2:length(temp)
                    if temp(i1)~=stemp(end)
                        stemp(end+1) = temp(i1);
                    end
                end
                
                % is the trial rewarded or not
                R = max(cg.timeseries.reward_ts(ttime));
                if R==1
                    % rewarded trial
                    if trial==1
                        % no phase 1
                        if not(isequal(stemp,2:9))
                            bad_trials = [bad_trials trial];
                        end
                    else
                        % phase 1 is present
                        if not(isequal(stemp,1:9))
                            bad_trials = [bad_trials trial];
                        end
                    end
                else
                    % unrewarded trial
                    if trial==1
                        % no phase 1
                        if not(isequal(stemp,[2 3 4 5 6 9]))
                            bad_trials = [bad_trials trial];
                        end
                    else
                        % phase 1 is present
                        if not(isequal(stemp,[1 2 3 4 5 6 9]))
                            bad_trials = [bad_trials trial];
                        end
                    end
                end
                
            end
            
        end
        bad_trials = unique(bad_trials);
        
        % Criterion 4: Actually, we need one more test: computing the duration of the
        % delay (d1+d2) and making sure that it has the required number of
        % dt (do not look at the experiment time for this)
        disp('Delay durations from the recordings:')
        for trial=1:ntrials
            
            ttime = find(cg.timeseries.trial_number_ts==trial);
            
            % length of d1
            nbd1 = length(find(cg.timeseries.epochs_ts(ttime)==3));
            
            % length of d2
            nbd2 = length(find(cg.timeseries.epochs_ts(ttime)==4));
            
            % delay duration
            delayDuration = (nbd1 + nbd2)*dt;
            
            % read the delay duration from the folder name
            if contains(fol,'2 SEC','IgnoreCase',true) || contains(fol,'2s','IgnoreCase',true)
                delay = 2;
            elseif contains(fol,'4 SEC','IgnoreCase',true) || contains(fol,'4s','IgnoreCase',true)
                delay = 4;
            elseif contains(fol,'6 SEC','IgnoreCase',true) || contains(fol,'6s','IgnoreCase',true)
                delay = 6;
            elseif contains(fol,'8 SEC','IgnoreCase',true) || contains(fol,'8s','IgnoreCase',true)
                delay = 8;
            else
               disp('Unable to read the delay duration...')
               pause
            end

            % mark the trial as illSampled if the number of frames is
            % inferior to 90% of what it should be
            if delayDuration<0.9*delay
                bad_trials = [bad_trials trial];
                disp(['Trial ' num2str(trial) ': ' num2str(delayDuration) 's <<<<< ill sampled delay']);
            else
                disp(['Trial ' num2str(trial) ': ' num2str(delayDuration) 's']);
            end
        end
        
        % Criterion 5: more tests: make sure that experimentTimes is monotonous
        % increasing. If it is not, remove the corresponding trial
        for trial=1:ntrials
            ttime = find(cg.timeseries.trial_number_ts==trial);
            times = syn.miniscopeMaster.clean_cageTimes(ttime);
            increments = [0 diff(times)];
            if min(increments)<0
                bad_trials = [bad_trials trial];
                disp(['Trial ' num2str(trial) ' removed as experiment times are not monotonously increasing...']);
            end
        end
        
        % Criterion 6: the houselight, if it is turned on, must be on for 5
        % seconds. Trials where this is not the case should be removed from
        % analysis.
        for trial=1:ntrials
            tframes = find(cg.timeseries.trial_number_ts==trial);
            light = cg.timeseries.houselight_ts(tframes);
            if sum(light)*dt>1.1*5
                bad_trials = [bad_trials trial];
                disp(['Trial ' num2str(trial) ' removed as house Light on for ' num2str(sum(light)*dt) 's, i.e. more than 5s...']);
            end
        end

        
        disp(' ');
        
        % final list
        bad_trials = unique(bad_trials);
        
        % label ill-sampled trials with a new time series, set bad trial(s)
        % to 1
        cg.timeseries.illsampled_ts = zeros(size(cg.timeseries.trial_number_ts));
        for trial=1:ntrials
            ttime = find(cg.timeseries.trial_number_ts==trial);
            if ismember(trial,bad_trials)
                cg.timeseries.illsampled_ts(ttime) = 1;
                disp(['Note: trial ' num2str(trial) ' was not usable: I am setting the illSampled variable to 1 during this trial.']);
            end
        end

%         % visual illustration after cleaning up undersampled trials.
%         subplot(3,1,1)
%         plot(cg.timeseries.epochs_ts,'-*')
%         xlabel('epochs (1 = wait, 2 = cue, etc.')
%         xlim([1 length(cg.timeseries.epochs_ts)])
%         temp = [0 diff(cg.timeseries.trial_number_ts)'];
%         trbnd = find(temp>0);
%         for l1=1:length(trbnd)
%            line(trbnd(l1)*[1 1],[0 10],'Color','r') 
%         end
%         
%         subplot(3,1,2)
%         plot(syn.miniscopeMaster.cageTimes,'-*')
%         xlabel('Frame instants (s)')
%         xlim([1 length(cg.timeseries.epochs_ts)])
%         for l1=1:length(trbnd)
%            line(trbnd(l1)*[1 1],max(syn.miniscopeMaster.cageTimes)*[0 1],'Color','r') 
%         end
%        
%         subplot(3,1,3)
%         plot(1./diff(syn.miniscopeMaster.cageTimes),'-*')
%         xlabel('frame rate (Hz)')
%         xlim([1 length(cg.timeseries.epochs_ts)])
%         for l1=1:length(trbnd)
%            line(trbnd(l1)*[1 1],[0 max(1./(diff(syn.miniscopeMaster.cageTimes)))],'Color','r') 
%         end
%         
%         pause
        
        % The undersampled trials are labeled by NaN in taskPhase, so we do
        % not need a separate field in dataset to output them. However, a
        % visual summary would be nice.
        sumfig = figure('Position',[555 814 1377 553]);
        
        % percentage of undersampled trials
        subplot(1,2,1)
        temp = [ntrials-length(bad_trials) length(bad_trials)];
        pie(temp)
        legend('ok trials','undersampled trials')
        xlabel('Trials')
        
        % percentage of undersampled frames
        subplot(1,2,2)
        temp1 = length(find(cg.timeseries.illsampled_ts==1));
        temp = [length(cg.timeseries.epochs_ts)-temp1 temp1];
        pie(temp)
        legend('ok frames','undersampled frames')
        xlabel('Frames')
        drawnow
        
        img = getframe(gcf);
        imwrite(img.cdata,[fol '/qualitycheck/summary_undersampling.png']);
        close(sumfig);
        
    end

% synchronize the trajectory and cage events with the miniscope time
% series, and save everything in a single structure.
    function [ds,nnew,ntrials] = synchro(syn,tr,cg,mn,ttot,fol)
        
        % synchronization = build a set of data for instants for
        % trajectory, calcium signal and cage events. This implies that the
        % result will only be defined over the shortest time interval where
        % these three information exist. We add "clean" fields to store
        % these values:
        syn.miniscopeMaster.clean_N = syn.miniscopeMaster.N;
        syn.miniscopeMaster.clean_masterFrames = syn.miniscopeMaster.masterFrames;
        syn.miniscopeMaster.clean_masterTimes = syn.miniscopeMaster.masterTimes;
        syn.miniscopeMaster.clean_slaveFrames = syn.miniscopeMaster.slaveFrames;
        syn.miniscopeMaster.clean_slaveTimes = syn.miniscopeMaster.slaveTimes;
        syn.miniscopeMaster.clean_cageTimes = syn.miniscopeMaster.cageTimes;
        
        % We now truncate synchronization.miniscopeMaster in steps
        
        % 1) miniscope frames:
        % By definition, original matrix synchronization.miniscopeMaster
        % has the same number of frames as the miniscope recording, so no
        % need to trim it here.
        
        % 2) behavior video frames:
        % Also by definition, the largest frame in the behavioral data series
        % cannot be larger than the total number of frames in videos behavCam[].avi.
        % So,nothing to trim here either.
        
        % 3) trajectory
        % Some of the video can correspond to something that is not usable (e.g.
        % at the end, when the experiment is over). So the trajectory can
        % have fewer frames than the behavCam videos. We trim these here.
        Ntr = length(tr.x);
        temp = find(syn.miniscopeMaster.clean_slaveFrames>Ntr);
        disp(['Removing ' num2str(length(temp)) ' frames that exceed the trajectory length.']);
        syn.miniscopeMaster.clean_masterFrames(temp) = [];
        syn.miniscopeMaster.clean_masterTimes(temp) = [];
        syn.miniscopeMaster.clean_slaveFrames(temp) = [];
        syn.miniscopeMaster.clean_slaveTimes(temp) = [];
        syn.miniscopeMaster.clean_cageTimes(temp) = [];
        syn.miniscopeMaster.clean_N = length(syn.miniscopeMaster.clean_masterFrames);
        
        % 4) session duration
        % In case the cage's record of events is shorter than the videos,
        % we need to trim the synchronization matrix some more.
        temp = find(syn.miniscopeMaster.clean_cageTimes>ttot);
        disp(['Removing ' num2str(length(temp)) ' frames that exceed the cage events recording length.']);
        syn.miniscopeMaster.clean_masterFrames(temp) = [];
        syn.miniscopeMaster.clean_masterTimes(temp) = [];
        syn.miniscopeMaster.clean_slaveFrames(temp) = [];
        syn.miniscopeMaster.clean_slaveTimes(temp) = [];
        syn.miniscopeMaster.clean_cageTimes(temp) = [];
        syn.miniscopeMaster.clean_N = length(syn.miniscopeMaster.clean_masterFrames);
        
        % 5) trial boundaries
        % Most likely, the end of the synchronization matrix falls
        % somewhere inside a trial, so the rest of that trial is probably
        % unusable. So, we just want to remove lines that fall after the
        % last complete trial.
        last_time = syn.miniscopeMaster.clean_cageTimes(end);
        temp = find(last_time>cg.events.trial_times);
        ntrials = temp(end)-1;
        time_boundary = cg.events.trial_times(ntrials+1);
        temp = find(syn.miniscopeMaster.clean_cageTimes>time_boundary);
        disp(['Removing ' num2str(length(temp)) ' frames that exceeded the last complete trial.']);
        syn.miniscopeMaster.clean_masterFrames(temp) = [];
        syn.miniscopeMaster.clean_masterTimes(temp) = [];
        syn.miniscopeMaster.clean_slaveFrames(temp) = [];
        syn.miniscopeMaster.clean_slaveTimes(temp) = [];
        syn.miniscopeMaster.clean_cageTimes(temp) = [];
        syn.miniscopeMaster.clean_N = length(syn.miniscopeMaster.clean_masterFrames);
        
        % Output the label of frames that were kept (to help Mohammad with
        % his analysis).
        kept_frames = 1:syn.miniscopeMaster.clean_N;
        
        disp('---------------------------------------------------------------------------------------------------------');
        disp('Final synchronization matrix size:');
        disp(['There are ' num2str(syn.miniscopeMaster.clean_N) ' miniscope frames.']);
        disp(['The last behavior frame is number ' num2str(syn.miniscopeMaster.clean_slaveFrames(end)) '.']);
        disp(['The last cage time is ' num2str(syn.miniscopeMaster.clean_cageTimes(end)) '.']);
        disp(['The synchronization matrix spans ' num2str(ntrials) ' trials.']);
        disp('---------------------------------------------------------------------------------------------------------');
        
        
        % 6) To simplify the synchronization below, we can also trim from
        % all the event time lists instants that fall outside of
        % the ntrials.
        
        disp('Trimming event variables to match the synchronization matrix:');
        Tmax = syn.miniscopeMaster.cageTimes(syn.miniscopeMaster.clean_N+1);
        
        % trial times
        disp('trial times')
        temp = find(cg.events.trial_times>=Tmax);
        cg.events.trial_times(temp) = [];
        
        % trial correction
        disp('trial correction')
        cg.events.correction(ntrials+1:end) = [];
        
        % sample
        disp('sample')
        cg.events.samples(ntrials+1:end) = [];
        cg.events.sample_times(ntrials+1:end) = [];
        
        % target
        disp('target')
        cg.events.targets(ntrials+1:end) = [];
        
        % target times are counted double as the cage reports target and
        % distractor presentation
        % remove double entries
        cg.events.target_times = unique(cg.events.target_times);
        % and only keep entries for the trials we consider
        cg.events.target_times(ntrials+1:end) = [];
        
        % start and end of delay
        disp('start and end of delay')
        cg.events.start_delay_times(ntrials+1:end) = [];
        cg.events.delay_end_times(ntrials+1:end) = [];
        
        % hiding images in order (sample, target, distractor) and repeat for
        % each trial.
        disp('hiding images')
        cg.events.hide_image(3*ntrials+1:end) = [];
        cg.events.hide_image_times(3*ntrials+1:end) = [];
        % remove repeated entries for target and distractor
        temp = [];
        for row=length(cg.events.hide_image)-1:-1:1
            if cg.events.hide_image_times(row)==cg.events.hide_image_times(row+1)
                temp = [temp row];
            end
        end
        cg.events.hide_image(temp) = [];
        cg.events.hide_image_times(temp) = [];
        
        % BIR beam break
        disp('BIR beam break')
        % remove events for trials beyond Ntrials
        temp = find(cg.events.break_BIR_start_times>Tmax);
        cg.events.break_BIR_start_times(temp) = [];
        
        % do the same thing for end times
        if isfield(cg.events,'break_BIR_end_times')
            % remove events for trials beyond Ntrials
            temp = find(cg.events.break_BIR_end_times>Tmax);
            cg.events.break_BIR_end_times(temp) = [];
        end
        
        % nose poke in reward hole
        disp('nose poke in reward hole')
        % remove events for trials beyond Ntrials
        temp = find(cg.events.poke_tray_start_times>Tmax);
        cg.events.poke_tray_start_times(temp) = [];
        % do the same thing for end times
        if isfield(cg.events,'poke_tray_end_times')
            % remove events for trials beyond Ntrials
            temp = find(cg.events.poke_tray_end_times>Tmax);
            cg.events.poke_tray_end_times(temp) = [];
        end
        
        % mouse contact with the touch screen
        disp('mouse contact with the touch screen')
        % pressing
        % remove events for trials beyond Ntrials
        temp = find(cg.events.touch_down_times>Tmax);
        cg.events.touch_down_times(temp) = [];
        cg.events.touch_down_position(temp) = [];
        cg.events.touch_down_lines(temp) = [];
        % releasing
        % remove events for trials beyond Ntrials
        temp = find(cg.events.touch_up_times>Tmax);
        cg.events.touch_up_times(temp) = [];
        cg.events.touch_up_position(temp) = [];
        cg.events.touch_up_lines(temp) = [];
        
        % tray light on/pulse/off
        disp('tray light on/pulse/off')
        % remove events for trials beyond Ntrials
        % on
        temp = find(cg.events.tray_light_on_times>Tmax);
        cg.events.tray_light_on_times(temp) = [];
        cg.events.tray_light_on_lines(temp) = [];
        % pulse
        temp = find(cg.events.tray_light_pulse_times>Tmax);
        cg.events.tray_light_pulse_times(temp) = [];
        cg.events.tray_light_pulse_lines(temp) = [];
        cg.events.tray_light_pulse_durations(temp) = [];
        % off
        temp = find(cg.events.tray_light_off_times>Tmax);
        cg.events.tray_light_off_times(temp) = [];
        cg.events.tray_light_off_lines(temp) = [];
        
        % sound pulses
        disp('sound pulses')
        % remove events for trials beyond Ntrials
        temp = find(cg.events.sound1_times>Tmax);
        cg.events.sound1_times(temp) = [];
        cg.events.sound1_lines(temp) = [];
        cg.events.sound1_durations(temp) = [];
        
        % reward
        disp('reward')
        % remove events for trials beyond Ntrials
        temp = find(cg.events.reward_times>Tmax);
        cg.events.reward_times(temp) = [];
        cg.events.reward(temp) = [];
        
        % incentives
        disp('incentives')
        % remove events for trials beyond Ntrials
        temp = find(cg.events.incentives_times>Tmax);
        cg.events.incentives_times(temp) = [];
        cg.events.incentives(temp) = [];
        
        % ITI boundaries
        disp('ITI boundaries')
        % remove events for trials beyond Ntrials
        temp = find(cg.events.ITI_start_times>Tmax);
        cg.events.ITI_start_times(temp) = [];
        
        % house light
        disp('house light')
        % remove events for trials beyond Ntrials
        % on
        temp = find(cg.events.house_light_on_times>Tmax);
        cg.events.house_light_on_times(temp) = [];
        % off
        temp = find(cg.events.house_light_off_times>Tmax);
        cg.events.house_light_off_times(temp) = [];
        
        % feeder
        disp('feeder')
        % remove events for trials beyond Ntrials
        temp = find(cg.events.feeder_times>Tmax);
        cg.events.feeder_times(temp) = [];
        cg.events.feeder_durations(temp) = [];
        
        disp('Done');
        
        % ===================== perform the actual synchronization ========
        
        disp('Synching the trajectory...');
        
        % trajectory:
        % -----------
        % Synchronizing the mouse trajectory (in the behav frame) to that in the
        % miniscope frame. We do this by going sequentially through all time points
        % in the miniscope time frame, look up the corresponding behavior time
        % frame point, and use the trajectory at that instant.
        Nmini = syn.miniscopeMaster.clean_N;
        
        % Now, we can synch the original trajectory with miniscope time
        tr.sync.x = zeros(Nmini,1);
        tr.sync.y = zeros(Nmini,1);
        % loop on miniscope time
        for t=1:Nmini
            behav_frame = syn.miniscopeMaster.clean_slaveFrames(t);
            tr.sync.x(t) = tr.x(behav_frame);
            tr.sync.y(t) = tr.y(behav_frame);
        end
        disp('Done!');
        
        disp('Synching bad data mask ...');
        % and the bad data mask (should be done in the same way)
        bad_data_sync = zeros(Nmini,1);
        % loop on miniscope time
        bad_data = tr.bad_data;
        for t=1:Nmini
            behav_frame = syn.miniscopeMaster.clean_slaveFrames(t);
            bad_data_sync(t) = bad_data(behav_frame);
        end
        disp('Done!');
        
        % cage:
        % -----
        disp('Creating the corresponding variable vectors...')
        [cg,nnew] = create_variables(cg,syn,tr,ntrials,fol);
        disp('Done!');
        
        % save all the results in a single structure
        ds.trialNumber = cg.timeseries.trial_number_ts;
        ds.trialType = cg.timeseries.trial_type_ts;
        ds.taskPhase = cg.timeseries.epochs_ts;
        
        % also it would be nice to have the time in the experiment as well for each
        % frame
        ds.experimentTimes = syn.miniscopeMaster.clean_cageTimes';
        
        ds.sample = cg.timeseries.sample_ts;
        ds.target = cg.timeseries.target_ts;
        ds.BIRBreak = cg.timeseries.break_BIR_ts;
        ds.nosePoke = cg.timeseries.poke_tray_ts;
        ds.reward = cg.timeseries.reward_ts;
        ds.incentives = cg.timeseries.incentives_ts;
        ds.touchUp = cg.timeseries.touch_up_ts;
        ds.touchDown = cg.timeseries.touch_down_ts;
        ds.sound = cg.timeseries.sound_ts;
        ds.light = cg.timeseries.light_ts;
        ds.feeder = cg.timeseries.feeder_ts;
        ds.houseLight = cg.timeseries.houselight_ts;
        ds.illSampled = cg.timeseries.illsampled_ts;
        
        % mouse trajectory
        ds.headPosition.x = tr.sync.x;
        ds.headPosition.y = tr.sync.y;
        % we use bad_data_sync as a probability of tracking the mouse.
        ds.headPosition.p = 1-bad_data_sync;

        % Miniscope frames that were kept in the final data set
        ds.keptFrames = kept_frames;
        
        % labels
        ds.labels{1} = 'trialNumber';
        ds.labels{end+1} = 'trialType';
        ds.labels{end+1} = 'taskPhase';
        ds.labels{end+1} = 'experimentTimes';
        ds.labels{end+1} = 'sample';
        ds.labels{end+1} = 'target';
        ds.labels{end+1} = 'BIRBreak';
        ds.labels{end+1} = 'nosePoke';
        ds.labels{end+1} = 'reward';
        ds.labels{end+1} = 'incentives';
        ds.labels{end+1} = 'touchUp';
        ds.labels{end+1} = 'touchDown';
        ds.labels{end+1} = 'sound';
        ds.labels{end+1} = 'light';
        ds.labels{end+1} = 'feeder';
        ds.labels{end+1} = 'houseLight';
        ds.labels{end+1} = 'keptFrames';
        ds.labels{end+1} = 'illSampled';
        
        % Finally the calcium signal
        % only keep what we need
        % deconvolved signal
        mn.dcs_cleaned = mn.original(1:nnew,:);
        ds.dcs = mn.dcs_cleaned;
        % signal before deconvolution
        mn.rs_cleaned = mn.raw(1:nnew,:);
        ds.rawTraces = mn.rs_cleaned;
        % filtered signal
        mn.filt_cleaned = mn.filt(1:nnew,:);
        ds.filtTraces = mn.filt_cleaned;
        
        % Potential problem: last element of dataset.trialNumber ,might be
        % 0 instead of the actual value. Nothing serious: just remove the last point
        % in this case.
        if ds.trialNumber(end)==0 && ...
                ds.trialType(end)==0 && ...
                ds.taskPhase(end)==0
            for i1=1:15
                eval(['ds.' ds.labels{i1} '(end) = [];']);
            end
        end
        
        disp('Structure Dataset is ready to be verified and saved.');
        
    end


% check the synchronized data ---------------------------------------------

% simple gui that allows to display to what extent mouse is near pads 1 and
% 5 when the cage detects contacts there.
    function checkcagebehav(fol,ds)
        check_cage_vs_behav(ds);
        drawnow
        img = getframe(gcf);
        imwrite(img.cdata,[fol '/qualitycheck/summary_cage_traj_synch.png']);
        close all
    end


end
