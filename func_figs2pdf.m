function func_figs2pdf(mydir)

import mlreportgen.report.*
import mlreportgen.dom.*

pos = find(mydir=='/');
name = mydir(pos(end)+1:end);

report = Report([mydir '/session_summary_' name],'pdf');

files = dir([mydir '/*.png']);
Nfiles = length(files);

report.Layout.Landscape = 1;

% preferred order to save figures
order = {};

order{end+1} = 'summary_bad_data.png';
order{end+1} = 'DLC_traj_summary.png';
order{end+1} = 'DLC_traj_delay.png';
order{end+1} = 'houselight_verification.png';
order{end+1} = 'summary_masking.png';
order{end+1} = 'summary_miniscope_position.png';
order{end+1} = 'summary_timestamp.png';
order{end+1} = 'summary_videos.png';
order{end+1} = 'summary_synchronization.png';
order{end+1} = 'summary_undersampling.png';
order{end+1} = 'houselight_alignment.png';
order{end+1} = 'summary_cage_traj_synch_no_lag_shift.png';
order{end+1} = 'summary_cage_traj_synch.png';
order{end+1} = 'summary_cage_traj_synch_with_lag_shift.png';
order{end+1} = 'corrected_trajectory_markers.png';

Norder = length(order);

for o=1:Norder
    my_order = order{o};
    for f=1:Nfiles
        if strcmp(my_order,files(f).name)
            temp = Image([files(f).folder '/' files(f).name]);
            temp.Width = '9.5in';
            temp.Height = [];
            add(report,temp);
        end
    end
end

close(report);

end