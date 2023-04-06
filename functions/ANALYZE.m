%% MAIN Function: Analyze
function [DLC,EXPERIMENT_PARAMS] = ANALYZE(DLC,EXPERIMENT_PARAMS, params, funcs)
    [DLC,EXPERIMENT_PARAMS] = extract_scale_shift(DLC,EXPERIMENT_PARAMS, params, funcs);
    [DLC,EXPERIMENT_PARAMS] = detect_licks_and_bouts(DLC,EXPERIMENT_PARAMS, params, funcs);
    [DLC,EXPERIMENT_PARAMS] = calculate_lick_kinematics(DLC,EXPERIMENT_PARAMS, params, funcs);
    [DLC,EXPERIMENT_PARAMS] = geometrization(DLC,EXPERIMENT_PARAMS, params, funcs);
    [DLC,EXPERIMENT_PARAMS] = sort_licks(DLC,EXPERIMENT_PARAMS, params, funcs);
    [DLC,EXPERIMENT_PARAMS] = quantify_food(DLC,EXPERIMENT_PARAMS, params, funcs);
    [DLC,EXPERIMENT_PARAMS] = detect_harvest(DLC,EXPERIMENT_PARAMS, params, funcs);
    %PGH_gif_maker(DLC,200:203,'png');
    %PGH_plot_sample_trace(DLC);
end

%% Function: Extract, scale, and shift DLC data
function [DLC, EXPERIMENT_PARAMS] = extract_scale_shift(DLC,EXPERIMENT_PARAMS, params, funcs)
fprintf('Extracting, scaling, and shifting DLC data ...')
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs;

% specify scale for pixel to mm conversion
webcam_scale = 7/70;
blackfly_scale = 7/47;
scale = [webcam_scale blackfly_scale];
EXPERIMENT_PARAMS.scale = scale;
FPS = EXPERIMENT_PARAMS.FPS;

% load, interpolate (1K), scale, and rotate DLC data from table into variables
data_FPS = table2array(DLC.data);
time_vid = ((1/FPS) : (1/FPS) : size(data_FPS,1)/FPS)';

% interpolation from FPS to 1K
DLC.FILE.interp = 1;
if DLC.FILE.interp == 1
    time_1K = (time_vid(1) : 0.001 : time_vid(end))';
    for counter_fields = 1 : size(data_FPS,2)
        data_1K(:,counter_fields) = interp1(time_vid,data_FPS(:,counter_fields),time_1K);
    end
else
    time_1K = time_vid; % for making gifs
    data_1K = data_FPS;
end

tip_tongue_x = scale(2)*data_1K(:,1);
tip_tongue_y = scale(2)*data_1K(:,2);
r_tongue_x = scale(2)*data_1K(:,3);
r_tongue_y = scale(2)*data_1K(:,4);
l_tongue_x = scale(2)*data_1K(:,5);
l_tongue_y = scale(2)*data_1K(:,6);
mid_tongue_x = scale(2)*data_1K(:,7);
mid_tongue_y = scale(2)*data_1K(:,8);
r_nose_x = scale(2)*data_1K(:,9);
r_nose_y = scale(2)*data_1K(:,10);
l_nose_x = scale(2)*data_1K(:,11);
l_nose_y = scale(2)*data_1K(:,12);
r_food_x = scale(2)*data_1K(:,13);
r_food_y = scale(2)*data_1K(:,14);
l_food_x = scale(2)*data_1K(:,15);
l_food_y = scale(2)*data_1K(:,16);
r_tube_r_x = scale(2)*data_1K(:,17);
r_tube_r_y = scale(2)*data_1K(:,18);
r_tube_l_x = scale(2)*data_1K(:,19);
r_tube_l_y = scale(2)*data_1K(:,20);
l_tube_r_x = scale(2)*data_1K(:,21);
l_tube_r_y = scale(2)*data_1K(:,22);
l_tube_l_x = scale(2)*data_1K(:,23);
l_tube_l_y = scale(2)*data_1K(:,24);

% find origin and shift DLC data
% clustering approach
% [~, cent_cluster_x] = kmeans(tip_tongue_x, 3);
% [~, cent_cluster_y] = kmeans(tip_tongue_y, 3);
% [x0, ind_min] = min(cent_cluster_x);
% if x0 < (mean(r_nose_x) + mean(l_nose_x))/2
%     [~, ind_max] = max(cent_cluster_x);
%     ind_mid = 1:3;
%     ind_mid(ind_max) = [];
%     ind_mid(ind_min) = [];
%     x0 = cent_cluster_x(ind_mid);
% end
%
% midpoint = (mean(r_nose_y) + mean(l_nose_y)) / 2;
% if abs(cent_cluster_y(1) - midpoint) < abs(cent_cluster_y(2) - midpoint) && abs(cent_cluster_y(1) - midpoint) < abs(cent_cluster_y(3) - midpoint)
%     y0 = cent_cluster_y(1);
% elseif abs(cent_cluster_y(2) - midpoint) < abs(cent_cluster_y(1) - midpoint) && abs(cent_cluster_y(2) - midpoint) < abs(cent_cluster_y(3) - midpoint)
%     y0 = cent_cluster_y(2);
% else
%     y0 = cent_cluster_y(3);
% end

% clustering approach
% nose marker approach
y0 = (mean(r_nose_y) + mean(l_nose_y)) / 2;
x0 = (mean(r_nose_x) + mean(l_nose_x)) / 2 + 5;

EXPERIMENT_PARAMS.duration_video = time_1K(end) - time_1K(1) ;
DLC.POINTS.tip_tongue_x = tip_tongue_x - x0;
DLC.POINTS.tip_tongue_y = tip_tongue_y - y0;
DLC.POINTS.r_tongue_x = r_tongue_x - x0;
DLC.POINTS.r_tongue_y = r_tongue_y - y0;
DLC.POINTS.l_tongue_x = l_tongue_x - x0;
DLC.POINTS.l_tongue_y = l_tongue_y - y0;
DLC.POINTS.mid_tongue_x = mid_tongue_x - x0;
DLC.POINTS.mid_tongue_y = mid_tongue_y - y0;
DLC.POINTS.r_nose_x = r_nose_x - x0;
DLC.POINTS.r_nose_y = r_nose_y - y0;
DLC.POINTS.l_nose_x = l_nose_x - x0;
DLC.POINTS.l_nose_y = l_nose_y - y0;
DLC.POINTS.l_food_x = l_food_x - x0;
DLC.POINTS.l_food_y = l_food_y - y0;
DLC.POINTS.l_tube_r_x = l_tube_r_x - x0;
DLC.POINTS.l_tube_r_y = l_tube_r_y - y0;
DLC.POINTS.l_tube_l_x = l_tube_l_x - x0;
DLC.POINTS.l_tube_l_y = l_tube_l_y - y0;
DLC.POINTS.r_food_x = r_food_x - x0;
DLC.POINTS.r_food_y = r_food_y - y0;
DLC.POINTS.r_tube_r_x = r_tube_r_x - x0;
DLC.POINTS.r_tube_r_y = r_tube_r_y - y0;
DLC.POINTS.r_tube_l_x = r_tube_l_x - x0;
DLC.POINTS.r_tube_l_y = r_tube_l_y - y0;
DLC.POINTS.x0 = x0;
DLC.POINTS.y0 = y0;
DLC.TIME.time_1K = time_1K;

fprintf(' --> Completed. \n')
end

%% Function: Detect licks and bouts
function [DLC, EXPERIMENT_PARAMS] = detect_licks_and_bouts(DLC,EXPERIMENT_PARAMS, params, funcs)
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs;
fprintf('Detecting licks and bouts ... ');
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs;
tip_tongue_x = DLC.POINTS.tip_tongue_x;
tip_tongue_y = DLC.POINTS.tip_tongue_y;

% set windows
if DLC.FILE.interp == 1
    smooth_window = 20;
    detection_window = 100;
    window_time = 500;
    bout_threshold = 1000;
else
    smooth_window = 2;
    detection_window = 10;
    window_time = 50;
    bout_threshold = 100;
end
prom_thresh = 0.85;

% compute and smooth displacement trace d_tip
d_tip = movmean(sqrt(tip_tongue_x.^2 + tip_tongue_y.^2),smooth_window)';

% compute negative signal and subtract within bout mins to take licks to 0
% while tongue in mouth
[peak_neg,ind_peak_neg,~,~]= findpeaks(-d_tip, 'MinPeakProminence', 1, 'MinPeakDistance',detection_window,'Annotate', 'extent');
is_del = isoutlier(peak_neg);
peak_neg(is_del) = [];
ind_peak_neg(is_del) = [];
d_tip_neg_interp = interp1(ind_peak_neg,peak_neg,1: length(d_tip));
d_tip = d_tip + d_tip_neg_interp;
d_tip(d_tip < 0) = 0;
d_tip(isnan(d_tip)) = 0;
if DLC.FILE.interp == 1
    v_tip = [0 diff(d_tip)]*1000;
else
    v_tip = [0 diff(d_tip)]*100;
end

% detect onset and offset of licks
[peak,ind_peak,~,prom]= findpeaks(d_tip, 'MinPeakProminence', 3, 'MinPeakDistance',detection_window,'Annotate', 'extent');
ind_lick_onset = nan(length(peak),1);
ind_lick_offset = nan(length(peak),1);
for counter_peak = 1 : length(peak)
    if ind_peak(counter_peak)-window_time < 1
        window_time = window_time + (ind_peak(counter_peak)-window_time) - 1;
    elseif ind_peak(counter_peak)+window_time > length(d_tip)
        window_time = window_time - (ind_peak(counter_peak) + window_time - length(d_tip));
    else
        if DLC.FILE.interp == 1
            window_time = 500;
        else
            window_time = 50;
        end
    end
    ind_lick_onset_ = find(d_tip(ind_peak(counter_peak)-window_time : ind_peak(counter_peak)) <= peak(counter_peak) - prom(counter_peak)*prom_thresh,1, 'last' ) + ind_peak(counter_peak)-window_time - 1;
    ind_lick_offset_ =  find(d_tip(ind_peak(counter_peak) : ind_peak(counter_peak)+window_time)  <= peak(counter_peak) - prom(counter_peak)*prom_thresh,1, 'first' ) + ind_peak(counter_peak);

    %recompute onset/offset using v_tip and values from d_tip
    ind_lick_onset_ = (find(v_tip(ind_lick_onset_-detection_window/2 : ind_lick_onset_+detection_window/2 ) <= 30,1, 'last' ) + ind_lick_onset_ - detection_window/2-1)-1;
    ind_lick_offset_ = (find(v_tip(ind_lick_offset_-detection_window/2 : ind_lick_offset_+detection_window/2 ) <= -30,1, 'last' ) + ind_lick_offset_ - detection_window/2-1) + 1;

    if isempty(ind_lick_onset_) || isempty(ind_lick_offset_)
        ind_lick_onset_ = nan;
        ind_lick_offset_ = nan;
    elseif d_tip(ind_lick_onset_)<prom(counter_peak) && d_tip(ind_lick_onset_) > 3
        ind_lick_onset_ = find(v_tip(ind_lick_onset_-detection_window : ind_lick_onset_) >= 30,1, 'first' ) + ind_lick_onset_ - detection_window-1;
    elseif d_tip(ind_lick_onset_)>prom(counter_peak) || d_tip(ind_lick_offset_)>prom(counter_peak)
        ind_lick_onset_ = nan;
        ind_lick_offset_ = nan;
    end

    if isempty(ind_lick_onset_) || isempty(ind_lick_offset_)
        ind_lick_onset_ = nan;
        ind_lick_offset_ = nan;
    end
    ind_lick_onset(counter_peak,1) = ind_lick_onset_;
    ind_lick_offset(counter_peak,1) =  ind_lick_offset_;
end
is_nan = isnan(ind_lick_onset) | isnan(ind_lick_offset) | (ind_lick_offset - ind_lick_onset)<=0;
ind_lick_onset(is_nan) = [];
ind_lick_offset(is_nan) = [];

is_outlier = isoutlier(d_tip(ind_lick_onset)) | isoutlier(d_tip(ind_lick_offset));
ind_lick_onset(is_outlier) = [];
ind_lick_offset(is_outlier) = [];

time_1K = DLC.TIME.time_1K;
time_lick_onset = time_1K(ind_lick_onset);
time_lick_offset =  time_1K(ind_lick_offset);
time_lick_duration = time_lick_offset - time_lick_onset;

validity = time_lick_duration <= 500;
ind_lick_onset(~validity) = [];
ind_lick_offset(~validity) = [];
time_lick_onset(~validity) = [];
time_lick_offset(~validity) = [];
time_lick_duration(~validity) = [];

num_lick = length(ind_lick_onset);

% 1s threshold for lick bouts
ind_lick_onset_str_bout_ = [1; find(diff(ind_lick_onset) >bout_threshold) + 1];
ind_lick_onset_str_bout = ind_lick_onset(ind_lick_onset_str_bout_);
ind_lick_onset_end_bout_ = [find(diff(ind_lick_onset) >bout_threshold); length(ind_lick_onset)];
ind_lick_onset_end_bout = ind_lick_onset(ind_lick_onset_end_bout_);

num_lick_bout = (ind_lick_onset_end_bout_ - ind_lick_onset_str_bout_) + 1;

validity = num_lick_bout>=3;
num_lick_bout(~validity) = [];
ind_lick_onset_str_bout(~validity) = [];
ind_lick_onset_end_bout(~validity) = [];

time_lick_onset_str_bout = time_1K(ind_lick_onset_str_bout);
time_lick_onset_end_bout = time_1K(ind_lick_onset_end_bout);
time_bout_duration = time_lick_onset_end_bout - time_lick_onset_str_bout;
num_bout = length(ind_lick_onset_str_bout);

DLC.KINEMATIC.d_tip = d_tip;
DLC.IND.num_lick = num_lick;
DLC.IND.num_bout = num_bout;
DLC.IND.num_lick_bout = num_lick_bout;
DLC.IND.ind_lick_onset = ind_lick_onset;
DLC.IND.ind_lick_offset = ind_lick_offset;
DLC.IND.ind_lick_onset_str_bout = ind_lick_onset_str_bout;
DLC.IND.ind_lick_onset_end_bout = ind_lick_onset_end_bout;
DLC.TIME.time_lick_onset_str_bout = time_lick_onset_str_bout;
DLC.TIME.time_lick_onset_end_bout = time_lick_onset_end_bout;
DLC.TIME.time_lick_onset = time_lick_onset;
DLC.TIME.time_lick_offset = time_lick_offset;
DLC.TIME.time_lick_duration = time_lick_duration;
DLC.TIME.time_bout_duration = time_bout_duration;

if EXPERIMENT_PARAMS.flag_figure_debug == 1
    figure
    hold on;
    plot(time_1K,d_tip,'k');
    plot(time_1K(ind_lick_onset_str_bout), d_tip(ind_lick_onset_str_bout), 'or', 'MarkerSize',10 );
    plot(time_1K(ind_lick_onset_end_bout), d_tip(ind_lick_onset_end_bout), 'ob', 'MarkerSize',10 );
    plot(time_1K(ind_lick_onset), d_tip(ind_lick_onset), '.r', 'MarkerSize',10 );
    plot(time_1K(ind_lick_offset), d_tip(ind_lick_onset), '.b', 'MarkerSize',10 );

    yyaxis right;
    hold on
    plot(time_1K,v_tip,'c');
    plot(time_1K(ind_lick_onset_str_bout), v_tip(ind_lick_onset_str_bout), 'or', 'MarkerSize',10 );
    plot(time_1K(ind_lick_onset_end_bout), v_tip(ind_lick_onset_end_bout), 'ob', 'MarkerSize',10 );
    plot(time_1K(ind_lick_onset), v_tip(ind_lick_onset), '.r', 'MarkerSize',10 );
    plot(time_1K(ind_lick_offset), v_tip(ind_lick_offset), '.b', 'MarkerSize',10 );

    xlabel('Time (s)');
    ylabel('Distance (mm)');
    title([EXPERIMENT_PARAMS.file_name ': ' num2str(num_lick) ' licks | ' num2str(num_bout) ' bouts'], 'interpreter', 'none')
    ESN_Beautify_Plot(gcf, [20 10])
end
fprintf(' --> Completed. \n')
end

%% Function: Calculate lick kinematics: d, v, a, angle, ILI(inter lick interval), ILR(instantaneous lick rate = 1/ILI)
function [DLC, EXPERIMENT_PARAMS] = calculate_lick_kinematics(DLC,EXPERIMENT_PARAMS, params, funcs)
fprintf('Calculating lick kinematics ...');
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs;
params.lick.length_trace = 500;
length_trace = params.lick.length_trace;
d_tip = DLC.KINEMATIC.d_tip;
time_lick_duration = DLC.TIME.time_lick_duration;

num_lick = DLC.IND.num_lick;
num_bout = DLC.IND.num_bout;
time_1K = DLC.TIME.time_1K;
ind_lick_onset = DLC.IND.ind_lick_onset;
ind_lick_offset = DLC.IND.ind_lick_offset;
ind_lick_onset_str_bout = DLC.IND.ind_lick_onset_str_bout;
ind_lick_onset_end_bout = DLC.IND.ind_lick_onset_end_bout;
tip_tongue_x = DLC.POINTS.tip_tongue_x;
tip_tongue_y = DLC.POINTS.tip_tongue_y;
r_tongue_x = DLC.POINTS.r_tongue_x;
r_tongue_y = DLC.POINTS.r_tongue_y;
l_tongue_x = DLC.POINTS.l_tongue_x;
l_tongue_y = DLC.POINTS.l_tongue_y;
mid_tongue_x = DLC.POINTS.mid_tongue_x;
mid_tongue_y = DLC.POINTS.mid_tongue_y;
l_tube_r_x = DLC.POINTS.l_tube_r_x;
l_tube_r_y = DLC.POINTS.l_tube_r_y;
l_tube_l_x = DLC.POINTS.l_tube_l_x;
l_tube_l_y = DLC.POINTS.l_tube_l_y;
l_food_x = DLC.POINTS.l_food_x;
l_food_y = DLC.POINTS.l_food_y;
r_tube_r_x = DLC.POINTS.r_tube_r_x;
r_tube_r_y = DLC.POINTS.r_tube_r_y;
r_tube_l_x = DLC.POINTS.r_tube_l_x;
r_tube_l_y = DLC.POINTS.r_tube_l_y;
r_food_x = DLC.POINTS.r_food_x;
r_food_y = DLC.POINTS.r_food_y;
r_nose_x = DLC.POINTS.r_nose_x;
r_nose_y = DLC.POINTS.r_nose_y;
l_nose_x = DLC.POINTS.l_nose_x;
l_nose_y = DLC.POINTS.l_nose_y;
midtip_y = tip_tongue_y - mid_tongue_y;
midtip_x = tip_tongue_x - mid_tongue_x;
time_lick_onset = DLC.TIME.time_lick_onset;

% build kinematicks _lick
for counter_lick =  1 : 1 : num_lick
    inds_ = (ind_lick_onset(counter_lick) ) : ind_lick_offset(counter_lick);
    if length(inds_) > 500
        inds_ = inds_(1:500);
        pad = [];
    elseif length(inds_) < 500
        % inds_ = padarray(inds_,[0 500-length(inds_)], inds_(end), 'post');
        pad = nan(1,500-length(inds_));
    end

    % d_lick
    d_lick(counter_lick, :) = [d_tip(inds_) pad];
    [d_lick_max(counter_lick,1), ind_d_lick_max_local] = max(d_lick(counter_lick, :));
    ind_d_lick_max(counter_lick,1) = ind_lick_onset(counter_lick) + ind_d_lick_max_local - 1;

    % v_lick
    v_lick(counter_lick, :) = [0 (diff(d_lick(counter_lick,:))./ (time_1K(2) - time_1K(1)))];
    [v_lick_max(counter_lick,1), ind_v_lick_max_local] = max(v_lick(counter_lick,:));
    ind_v_lick_max(counter_lick,1) = ind_lick_onset(counter_lick) + ind_v_lick_max_local - 1;

    [v_lick_min(counter_lick,1), ind_v_lick_min_local] = min(v_lick(counter_lick,:));
    ind_v_lick_min(counter_lick,1) = ind_lick_onset(counter_lick) + ind_v_lick_min_local - 1;

    % angle_lick
    angle_lick(counter_lick, :) = [rad2deg(atan2(midtip_y(inds_),midtip_x(inds_)))' pad];

    % points
    ind_onset_(counter_lick) = ind_lick_onset(counter_lick);
    ind_vmax_(counter_lick) = ind_v_lick_max(counter_lick);
    ind_dmax_(counter_lick) = ind_d_lick_max(counter_lick);
    ind_vmin_(counter_lick) = ind_v_lick_min(counter_lick);
    ind_offset_(counter_lick) = ind_lick_offset(counter_lick);

    tip_tongue_x_lick(counter_lick, :) = [tip_tongue_x(inds_)' pad];
    tip_tongue_y_lick(counter_lick, :) = [tip_tongue_y(inds_)' pad];
    tip_tongue_x_onset(counter_lick) = tip_tongue_x(ind_onset_(counter_lick));
    tip_tongue_y_onset(counter_lick) = tip_tongue_y(ind_onset_(counter_lick));
    tip_tongue_x_vmax(counter_lick) = tip_tongue_x(ind_vmax_(counter_lick));
    tip_tongue_y_vmax(counter_lick) = tip_tongue_y(ind_vmax_(counter_lick));
    tip_tongue_x_dmax(counter_lick) = tip_tongue_x(ind_dmax_(counter_lick));
    tip_tongue_y_dmax(counter_lick) = tip_tongue_y(ind_dmax_(counter_lick));
    tip_tongue_x_vmin(counter_lick) = tip_tongue_x(ind_vmin_(counter_lick));
    tip_tongue_y_vmin(counter_lick) = tip_tongue_y(ind_vmin_(counter_lick));
    tip_tongue_x_offset(counter_lick) = tip_tongue_x(ind_offset_(counter_lick));
    tip_tongue_y_offset(counter_lick) = tip_tongue_y(ind_offset_(counter_lick));

    r_tongue_x_lick(counter_lick, :) = [r_tongue_x(inds_)' pad];
    r_tongue_y_lick(counter_lick, :) = [r_tongue_y(inds_)' pad];
    r_tongue_x_onset(counter_lick) = r_tongue_x(ind_onset_(counter_lick));
    r_tongue_y_onset(counter_lick) = r_tongue_y(ind_onset_(counter_lick));
    r_tongue_x_vmax(counter_lick) = r_tongue_x(ind_vmax_(counter_lick));
    r_tongue_y_vmax(counter_lick) = r_tongue_y(ind_vmax_(counter_lick));
    r_tongue_x_dmax(counter_lick) = r_tongue_x(ind_dmax_(counter_lick));
    r_tongue_y_dmax(counter_lick) = r_tongue_y(ind_dmax_(counter_lick));
    r_tongue_x_vmin(counter_lick) = r_tongue_x(ind_vmin_(counter_lick));
    r_tongue_y_vmin(counter_lick) = r_tongue_y(ind_vmin_(counter_lick));
    r_tongue_x_offset(counter_lick) = r_tongue_x(ind_offset_(counter_lick));
    r_tongue_y_offset(counter_lick) = r_tongue_y(ind_offset_(counter_lick));

    l_tongue_x_lick(counter_lick, :) = [l_tongue_x(inds_)' pad];
    l_tongue_y_lick(counter_lick, :) = [l_tongue_y(inds_)' pad];
    l_tongue_x_onset(counter_lick) = l_tongue_x(ind_onset_(counter_lick));
    l_tongue_y_onset(counter_lick) = l_tongue_y(ind_onset_(counter_lick));
    l_tongue_x_vmax(counter_lick) = l_tongue_x(ind_vmax_(counter_lick));
    l_tongue_y_vmax(counter_lick) = l_tongue_y(ind_vmax_(counter_lick));
    l_tongue_x_dmax(counter_lick) = l_tongue_x(ind_dmax_(counter_lick));
    l_tongue_y_dmax(counter_lick) = l_tongue_y(ind_dmax_(counter_lick));
    l_tongue_x_vmin(counter_lick) = l_tongue_x(ind_vmin_(counter_lick));
    l_tongue_y_vmin(counter_lick) = l_tongue_y(ind_vmin_(counter_lick));
    l_tongue_x_offset(counter_lick) = l_tongue_x(ind_offset_(counter_lick));
    l_tongue_y_offset(counter_lick) = l_tongue_y(ind_offset_(counter_lick));

    mid_tongue_x_lick(counter_lick, :) = [mid_tongue_x(inds_)' pad];
    mid_tongue_y_lick(counter_lick, :) = [mid_tongue_y(inds_)' pad];
    mid_tongue_x_onset(counter_lick) = mid_tongue_x(ind_onset_(counter_lick));
    mid_tongue_y_onset(counter_lick) = mid_tongue_y(ind_onset_(counter_lick));
    mid_tongue_x_vmax(counter_lick) = mid_tongue_x(ind_vmax_(counter_lick));
    mid_tongue_y_vmax(counter_lick) = mid_tongue_y(ind_vmax_(counter_lick));
    mid_tongue_x_dmax(counter_lick) = mid_tongue_x(ind_dmax_(counter_lick));
    mid_tongue_y_dmax(counter_lick) = mid_tongue_y(ind_dmax_(counter_lick));
    mid_tongue_x_vmin(counter_lick) = mid_tongue_x(ind_vmin_(counter_lick));
    mid_tongue_y_vmin(counter_lick) = mid_tongue_y(ind_vmin_(counter_lick));
    mid_tongue_x_offset(counter_lick) = mid_tongue_x(ind_offset_(counter_lick));
    mid_tongue_y_offset(counter_lick) = mid_tongue_y(ind_offset_(counter_lick));

    r_food_x_lick(counter_lick, :) = [r_food_x(inds_)' pad];
    r_food_y_lick(counter_lick, :) = [r_food_y(inds_)' pad];
    r_food_x_onset(counter_lick) = r_food_x(ind_onset_(counter_lick));
    r_food_y_onset(counter_lick) = r_food_y(ind_onset_(counter_lick));
    r_food_x_vmax(counter_lick) = r_food_x(ind_vmax_(counter_lick));
    r_food_y_vmax(counter_lick) = r_food_y(ind_vmax_(counter_lick));
    r_food_x_dmax(counter_lick) = r_food_x(ind_dmax_(counter_lick));
    r_food_y_dmax(counter_lick) = r_food_y(ind_dmax_(counter_lick));
    r_food_x_vmin(counter_lick) = r_food_x(ind_vmin_(counter_lick));
    r_food_y_vmin(counter_lick) = r_food_y(ind_vmin_(counter_lick));
    r_food_x_offset(counter_lick) = r_food_x(ind_offset_(counter_lick));
    r_food_y_offset(counter_lick) = r_food_y(ind_offset_(counter_lick));

    l_food_x_lick(counter_lick, :) = [l_food_x(inds_)' pad];
    l_food_y_lick(counter_lick, :) = [l_food_y(inds_)' pad];
    l_food_x_onset(counter_lick) = l_food_x(ind_onset_(counter_lick));
    l_food_y_onset(counter_lick) = l_food_y(ind_onset_(counter_lick));
    l_food_x_vmax(counter_lick) = l_food_x(ind_vmax_(counter_lick));
    l_food_y_vmax(counter_lick) = l_food_y(ind_vmax_(counter_lick));
    l_food_x_dmax(counter_lick) = l_food_x(ind_dmax_(counter_lick));
    l_food_y_dmax(counter_lick) = l_food_y(ind_dmax_(counter_lick));
    l_food_x_vmin(counter_lick) = l_food_x(ind_vmin_(counter_lick));
    l_food_y_vmin(counter_lick) = l_food_y(ind_vmin_(counter_lick));
    l_food_x_offset(counter_lick) = l_food_x(ind_offset_(counter_lick));
    l_food_y_offset(counter_lick) = l_food_y(ind_offset_(counter_lick));

    r_nose_x_lick(counter_lick, :) = [r_nose_x(inds_)' pad];
    r_nose_y_lick(counter_lick, :) = [r_nose_y(inds_)' pad];
    r_nose_x_onset(counter_lick) = r_nose_x(ind_onset_(counter_lick));
    r_nose_y_onset(counter_lick) = r_nose_y(ind_onset_(counter_lick));
    r_nose_x_vmax(counter_lick) = r_nose_x(ind_vmax_(counter_lick));
    r_nose_y_vmax(counter_lick) = r_nose_y(ind_vmax_(counter_lick));
    r_nose_x_dmax(counter_lick) = r_nose_x(ind_dmax_(counter_lick));
    r_nose_y_dmax(counter_lick) = r_nose_y(ind_dmax_(counter_lick));
    r_nose_x_vmin(counter_lick) = r_nose_x(ind_vmin_(counter_lick));
    r_nose_y_vmin(counter_lick) = r_nose_y(ind_vmin_(counter_lick));
    r_nose_x_offset(counter_lick) = r_nose_x(ind_offset_(counter_lick));
    r_nose_y_offset(counter_lick) = r_nose_y(ind_offset_(counter_lick));

    l_nose_x_lick(counter_lick, :) = [l_nose_x(inds_)' pad];
    l_nose_y_lick(counter_lick, :) = [l_nose_y(inds_)' pad];
    l_nose_x_onset(counter_lick) = l_nose_x(ind_onset_(counter_lick));
    l_nose_y_onset(counter_lick) = l_nose_y(ind_onset_(counter_lick));
    l_nose_x_vmax(counter_lick) = l_nose_x(ind_vmax_(counter_lick));
    l_nose_y_vmax(counter_lick) = l_nose_y(ind_vmax_(counter_lick));
    l_nose_x_dmax(counter_lick) = l_nose_x(ind_dmax_(counter_lick));
    l_nose_y_dmax(counter_lick) = l_nose_y(ind_dmax_(counter_lick));
    l_nose_x_vmin(counter_lick) = l_nose_x(ind_vmin_(counter_lick));
    l_nose_y_vmin(counter_lick) = l_nose_y(ind_vmin_(counter_lick));
    l_nose_x_offset(counter_lick) = l_nose_x(ind_offset_(counter_lick));
    l_nose_y_offset(counter_lick) = l_nose_y(ind_offset_(counter_lick));

    r_tube_r_x_lick(counter_lick, :) = [r_tube_r_x(inds_)' pad];
    r_tube_r_y_lick(counter_lick, :) = [r_tube_r_y(inds_)' pad];
    r_tube_r_x_onset(counter_lick) = r_tube_r_x(ind_onset_(counter_lick));
    r_tube_r_y_onset(counter_lick) = r_tube_r_y(ind_onset_(counter_lick));
    r_tube_r_x_vmax(counter_lick) = r_tube_r_x(ind_vmax_(counter_lick));
    r_tube_r_y_vmax(counter_lick) = r_tube_r_y(ind_vmax_(counter_lick));
    r_tube_r_x_dmax(counter_lick) = r_tube_r_x(ind_dmax_(counter_lick));
    r_tube_r_y_dmax(counter_lick) = r_tube_r_y(ind_dmax_(counter_lick));
    r_tube_r_x_vmin(counter_lick) = r_tube_r_x(ind_vmin_(counter_lick));
    r_tube_r_y_vmin(counter_lick) = r_tube_r_y(ind_vmin_(counter_lick));
    r_tube_r_x_offset(counter_lick) = r_tube_r_x(ind_offset_(counter_lick));
    r_tube_r_y_offset(counter_lick) = r_tube_r_y(ind_offset_(counter_lick));

    r_tube_l_x_lick(counter_lick, :) = [r_tube_l_x(inds_)' pad];
    r_tube_l_y_lick(counter_lick, :) = [r_tube_l_y(inds_)' pad];
    r_tube_l_x_onset(counter_lick) = r_tube_l_x(ind_onset_(counter_lick));
    r_tube_l_y_onset(counter_lick) = r_tube_l_y(ind_onset_(counter_lick));
    r_tube_l_x_vmax(counter_lick) = r_tube_l_x(ind_vmax_(counter_lick));
    r_tube_l_y_vmax(counter_lick) = r_tube_l_y(ind_vmax_(counter_lick));
    r_tube_l_x_dmax(counter_lick) = r_tube_l_x(ind_dmax_(counter_lick));
    r_tube_l_y_dmax(counter_lick) = r_tube_l_y(ind_dmax_(counter_lick));
    r_tube_l_x_vmin(counter_lick) = r_tube_l_x(ind_vmin_(counter_lick));
    r_tube_l_y_vmin(counter_lick) = r_tube_l_y(ind_vmin_(counter_lick));
    r_tube_l_x_offset(counter_lick) = r_tube_l_x(ind_offset_(counter_lick));
    r_tube_l_y_offset(counter_lick) = r_tube_l_y(ind_offset_(counter_lick));

    l_tube_r_x_lick(counter_lick, :) = [l_tube_r_x(inds_)' pad];
    l_tube_r_y_lick(counter_lick, :) = [l_tube_r_y(inds_)' pad];
    l_tube_r_x_onset(counter_lick) = l_tube_r_x(ind_onset_(counter_lick));
    l_tube_r_y_onset(counter_lick) = l_tube_r_y(ind_onset_(counter_lick));
    l_tube_r_x_vmax(counter_lick) = l_tube_r_x(ind_vmax_(counter_lick));
    l_tube_r_y_vmax(counter_lick) = l_tube_r_y(ind_vmax_(counter_lick));
    l_tube_r_x_dmax(counter_lick) = l_tube_r_x(ind_dmax_(counter_lick));
    l_tube_r_y_dmax(counter_lick) = l_tube_r_y(ind_dmax_(counter_lick));
    l_tube_r_x_vmin(counter_lick) = l_tube_r_x(ind_vmin_(counter_lick));
    l_tube_r_y_vmin(counter_lick) = l_tube_r_y(ind_vmin_(counter_lick));
    l_tube_r_x_offset(counter_lick) = l_tube_r_x(ind_offset_(counter_lick));
    l_tube_r_y_offset(counter_lick) = l_tube_r_y(ind_offset_(counter_lick));

    l_tube_l_x_lick(counter_lick, :) = [l_tube_l_x(inds_)' pad];
    l_tube_l_y_lick(counter_lick, :) = [l_tube_l_y(inds_)' pad];
    l_tube_l_x_onset(counter_lick) = l_tube_l_x(ind_onset_(counter_lick));
    l_tube_l_y_onset(counter_lick) = l_tube_l_y(ind_onset_(counter_lick));
    l_tube_l_x_vmax(counter_lick) = l_tube_l_x(ind_vmax_(counter_lick));
    l_tube_l_y_vmax(counter_lick) = l_tube_l_y(ind_vmax_(counter_lick));
    l_tube_l_x_dmax(counter_lick) = l_tube_l_x(ind_dmax_(counter_lick));
    l_tube_l_y_dmax(counter_lick) = l_tube_l_y(ind_dmax_(counter_lick));
    l_tube_l_x_vmin(counter_lick) = l_tube_l_x(ind_vmin_(counter_lick));
    l_tube_l_y_vmin(counter_lick) = l_tube_l_y(ind_vmin_(counter_lick));
    l_tube_l_x_offset(counter_lick) = l_tube_l_x(ind_offset_(counter_lick));
    l_tube_l_y_offset(counter_lick) = l_tube_l_y(ind_offset_(counter_lick));
end

% ILI bout
for counter_bout = 1 : 1 : num_bout
    ILI_bout(counter_bout,1) = mean(diff(time_lick_onset(find(ind_lick_onset >= ind_lick_onset_str_bout(counter_bout) & ind_lick_onset <= ind_lick_onset_end_bout(counter_bout)))));
end



v_tip = ([0 (diff(d_tip)./ (time_1K(2) - time_1K(1)))]);
angle_midtip = rad2deg(atan2(midtip_y, midtip_x))';
angle_lick_max = angle_midtip(ind_d_lick_max)';

ILR_bout = 1./ILI_bout;
ILI_bout((ILR_bout==Inf)) = nan;
ILR_bout((ILR_bout==Inf)) = nan;

ILI_lick = diff(time_lick_onset);
ILI_lick = [ILI_lick; nan];
ILR_lick = 1./ILI_lick;

time_d_lick_max_abs = time_1K(ind_d_lick_max);
time_d_lick_max_rel = (time_d_lick_max_abs - time_lick_onset);
time_v_lick_max_abs = time_1K(ind_v_lick_max);
time_v_lick_max_rel = (time_v_lick_max_abs - time_lick_onset);
time_v_lick_min_abs = time_1K(ind_v_lick_min);
time_v_lick_min_rel = (time_v_lick_min_abs - time_lick_onset);

DLC.KINEMATIC.d_tip = d_tip;
DLC.KINEMATIC.d_lick = d_lick;
DLC.KINEMATIC.d_lick_max = d_lick_max;
DLC.IND.ind_d_lick_max = ind_d_lick_max;
DLC.KINEMATIC.v_tip = v_tip;
DLC.KINEMATIC.angle_lick = angle_lick;
DLC.KINEMATIC.angle_midtip = angle_midtip;
DLC.KINEMATIC.angle_lick_max = angle_lick_max;
DLC.KINEMATIC.v_lick = v_lick;
DLC.KINEMATIC.v_lick_max = v_lick_max;
DLC.IND.ind_v_lick_max = ind_v_lick_max;
DLC.KINEMATIC.v_lick_min = v_lick_min;
DLC.IND.ind_v_lick_min = ind_v_lick_min;
DLC.TIME.time_d_lick_max_abs = time_d_lick_max_abs;
DLC.TIME.time_d_lick_max_rel = time_d_lick_max_rel;
DLC.TIME.time_v_lick_max_abs = time_v_lick_max_abs;
DLC.TIME.time_v_lick_max_rel = time_v_lick_max_rel;
DLC.TIME.time_v_lick_min_abs = time_v_lick_min_abs;
DLC.TIME.time_v_lick_min_rel = time_v_lick_min_rel;
DLC.KINEMATIC.ILI_bout = ILI_bout;
DLC.KINEMATIC.ILR_bout = ILR_bout;
DLC.KINEMATIC.ILI_lick = ILI_lick;
DLC.KINEMATIC.ILR_lick = ILR_lick;

DLC.KINEMATIC.tip_tongue_x_lick = tip_tongue_x_lick;
DLC.KINEMATIC.tip_tongue_y_lick = tip_tongue_y_lick ;
DLC.KINEMATIC.r_tongue_x_lick = r_tongue_x_lick ;
DLC.KINEMATIC.r_tongue_y_lick = r_tongue_y_lick ;
DLC.KINEMATIC.l_tongue_x_lick = l_tongue_x_lick ;
DLC.KINEMATIC.l_tongue_y_lick = l_tongue_y_lick ;
DLC.KINEMATIC.mid_tongue_x_lick = mid_tongue_x_lick ;
DLC.KINEMATIC.mid_tongue_y_lick = mid_tongue_y_lick ;
DLC.KINEMATIC.r_nose_x_lick = r_nose_x_lick ;
DLC.KINEMATIC.r_nose_y_lick = r_nose_y_lick ;
DLC.KINEMATIC.l_nose_x_lick = l_nose_x_lick ;
DLC.KINEMATIC.l_nose_y_lick = l_nose_y_lick ;
DLC.KINEMATIC.l_food_x_lick = l_food_x_lick ;
DLC.KINEMATIC.l_food_y_lick = l_food_y_lick ;
DLC.KINEMATIC.l_tube_r_x_lick = l_tube_r_x_lick ;
DLC.KINEMATIC.l_tube_r_y_lick = l_tube_r_y_lick ;
DLC.KINEMATIC.l_tube_l_x_lick = l_tube_l_x_lick ;
DLC.KINEMATIC.l_tube_l_y_lick = l_tube_l_y_lick ;
DLC.KINEMATIC.r_food_x_lick = r_food_x_lick ;
DLC.KINEMATIC.r_food_y_lick = r_food_y_lick ;
DLC.KINEMATIC.r_tube_r_x_lick = r_tube_r_x_lick ;
DLC.KINEMATIC.r_tube_r_y_lick = r_tube_r_y_lick ;
DLC.KINEMATIC.r_tube_l_x_lick = r_tube_l_x_lick ;
DLC.KINEMATIC.r_tube_l_y_lick = r_tube_l_y_lick ;

DLC.KINEMATIC.tip_tongue_x_onset = tip_tongue_x_onset;
DLC.KINEMATIC.tip_tongue_y_onset = tip_tongue_y_onset ;
DLC.KINEMATIC.r_tongue_x_onset = r_tongue_x_onset ;
DLC.KINEMATIC.r_tongue_y_onset = r_tongue_y_onset ;
DLC.KINEMATIC.l_tongue_x_onset = l_tongue_x_onset ;
DLC.KINEMATIC.l_tongue_y_onset = l_tongue_y_onset ;
DLC.KINEMATIC.mid_tongue_x_onset = mid_tongue_x_onset ;
DLC.KINEMATIC.mid_tongue_y_onset = mid_tongue_y_onset ;
DLC.KINEMATIC.r_nose_x_onset = r_nose_x_onset ;
DLC.KINEMATIC.r_nose_y_onset = r_nose_y_onset ;
DLC.KINEMATIC.l_nose_x_onset = l_nose_x_onset ;
DLC.KINEMATIC.l_nose_y_onset = l_nose_y_onset ;
DLC.KINEMATIC.l_food_x_onset = l_food_x_onset ;
DLC.KINEMATIC.l_food_y_onset = l_food_y_onset ;
DLC.KINEMATIC.l_tube_r_x_onset = l_tube_r_x_onset ;
DLC.KINEMATIC.l_tube_r_y_onset = l_tube_r_y_onset ;
DLC.KINEMATIC.l_tube_l_x_onset = l_tube_l_x_onset ;
DLC.KINEMATIC.l_tube_l_y_onset = l_tube_l_y_onset ;
DLC.KINEMATIC.r_food_x_onset = r_food_x_onset ;
DLC.KINEMATIC.r_food_y_onset = r_food_y_onset ;
DLC.KINEMATIC.r_tube_r_x_onset = r_tube_r_x_onset ;
DLC.KINEMATIC.r_tube_r_y_onset = r_tube_r_y_onset ;
DLC.KINEMATIC.r_tube_l_x_onset = r_tube_l_x_onset ;
DLC.KINEMATIC.r_tube_l_y_onset = r_tube_l_y_onset ;

DLC.KINEMATIC.tip_tongue_x_vmax = tip_tongue_x_vmax;
DLC.KINEMATIC.tip_tongue_y_vmax = tip_tongue_y_vmax ;
DLC.KINEMATIC.r_tongue_x_vmax = r_tongue_x_vmax ;
DLC.KINEMATIC.r_tongue_y_vmax = r_tongue_y_vmax ;
DLC.KINEMATIC.l_tongue_x_vmax = l_tongue_x_vmax ;
DLC.KINEMATIC.l_tongue_y_vmax = l_tongue_y_vmax ;
DLC.KINEMATIC.mid_tongue_x_vmax = mid_tongue_x_vmax ;
DLC.KINEMATIC.mid_tongue_y_vmax = mid_tongue_y_vmax ;
DLC.KINEMATIC.r_nose_x_vmax = r_nose_x_vmax ;
DLC.KINEMATIC.r_nose_y_vmax = r_nose_y_vmax ;
DLC.KINEMATIC.l_nose_x_vmax = l_nose_x_vmax ;
DLC.KINEMATIC.l_nose_y_vmax = l_nose_y_vmax ;
DLC.KINEMATIC.l_food_x_vmax = l_food_x_vmax ;
DLC.KINEMATIC.l_food_y_vmax = l_food_y_vmax ;
DLC.KINEMATIC.l_tube_r_x_vmax = l_tube_r_x_vmax ;
DLC.KINEMATIC.l_tube_r_y_vmax = l_tube_r_y_vmax ;
DLC.KINEMATIC.l_tube_l_x_vmax = l_tube_l_x_vmax ;
DLC.KINEMATIC.l_tube_l_y_vmax = l_tube_l_y_vmax ;
DLC.KINEMATIC.r_food_x_vmax = r_food_x_vmax ;
DLC.KINEMATIC.r_food_y_vmax = r_food_y_vmax ;
DLC.KINEMATIC.r_tube_r_x_vmax = r_tube_r_x_vmax ;
DLC.KINEMATIC.r_tube_r_y_vmax = r_tube_r_y_vmax ;
DLC.KINEMATIC.r_tube_l_x_vmax = r_tube_l_x_vmax ;
DLC.KINEMATIC.r_tube_l_y_vmax = r_tube_l_y_vmax ;

DLC.KINEMATIC.tip_tongue_x_dmax = tip_tongue_x_dmax;
DLC.KINEMATIC.tip_tongue_y_dmax = tip_tongue_y_dmax ;
DLC.KINEMATIC.r_tongue_x_dmax = r_tongue_x_dmax ;
DLC.KINEMATIC.r_tongue_y_dmax = r_tongue_y_dmax ;
DLC.KINEMATIC.l_tongue_x_dmax = l_tongue_x_dmax ;
DLC.KINEMATIC.l_tongue_y_dmax = l_tongue_y_dmax ;
DLC.KINEMATIC.mid_tongue_x_dmax = mid_tongue_x_dmax ;
DLC.KINEMATIC.mid_tongue_y_dmax = mid_tongue_y_dmax ;
DLC.KINEMATIC.r_nose_x_dmax = r_nose_x_dmax ;
DLC.KINEMATIC.r_nose_y_dmax = r_nose_y_dmax ;
DLC.KINEMATIC.l_nose_x_dmax = l_nose_x_dmax ;
DLC.KINEMATIC.l_nose_y_dmax = l_nose_y_dmax ;
DLC.KINEMATIC.l_food_x_dmax = l_food_x_dmax ;
DLC.KINEMATIC.l_food_y_dmax = l_food_y_dmax ;
DLC.KINEMATIC.l_tube_r_x_dmax = l_tube_r_x_dmax ;
DLC.KINEMATIC.l_tube_r_y_dmax = l_tube_r_y_dmax ;
DLC.KINEMATIC.l_tube_l_x_dmax = l_tube_l_x_dmax ;
DLC.KINEMATIC.l_tube_l_y_dmax = l_tube_l_y_dmax ;
DLC.KINEMATIC.r_food_x_dmax = r_food_x_dmax ;
DLC.KINEMATIC.r_food_y_dmax = r_food_y_dmax ;
DLC.KINEMATIC.r_tube_r_x_dmax = r_tube_r_x_dmax ;
DLC.KINEMATIC.r_tube_r_y_dmax = r_tube_r_y_dmax ;
DLC.KINEMATIC.r_tube_l_x_dmax = r_tube_l_x_dmax ;
DLC.KINEMATIC.r_tube_l_y_dmax = r_tube_l_y_dmax ;

DLC.KINEMATIC.tip_tongue_x_vmin = tip_tongue_x_vmin;
DLC.KINEMATIC.tip_tongue_y_vmin = tip_tongue_y_vmin ;
DLC.KINEMATIC.r_tongue_x_vmin = r_tongue_x_vmin ;
DLC.KINEMATIC.r_tongue_y_vmin = r_tongue_y_vmin ;
DLC.KINEMATIC.l_tongue_x_vmin = l_tongue_x_vmin ;
DLC.KINEMATIC.l_tongue_y_vmin = l_tongue_y_vmin ;
DLC.KINEMATIC.mid_tongue_x_vmin = mid_tongue_x_vmin ;
DLC.KINEMATIC.mid_tongue_y_vmin = mid_tongue_y_vmin ;
DLC.KINEMATIC.r_nose_x_vmin = r_nose_x_vmin ;
DLC.KINEMATIC.r_nose_y_vmin = r_nose_y_vmin ;
DLC.KINEMATIC.l_nose_x_vmin = l_nose_x_vmin ;
DLC.KINEMATIC.l_nose_y_vmin = l_nose_y_vmin ;
DLC.KINEMATIC.l_food_x_vmin = l_food_x_vmin ;
DLC.KINEMATIC.l_food_y_vmin = l_food_y_vmin ;
DLC.KINEMATIC.l_tube_r_x_vmin = l_tube_r_x_vmin ;
DLC.KINEMATIC.l_tube_r_y_vmin = l_tube_r_y_vmin ;
DLC.KINEMATIC.l_tube_l_x_vmin = l_tube_l_x_vmin ;
DLC.KINEMATIC.l_tube_l_y_vmin = l_tube_l_y_vmin ;
DLC.KINEMATIC.r_food_x_vmin = r_food_x_vmin ;
DLC.KINEMATIC.r_food_y_vmin = r_food_y_vmin ;
DLC.KINEMATIC.r_tube_r_x_vmin = r_tube_r_x_vmin ;
DLC.KINEMATIC.r_tube_r_y_vmin = r_tube_r_y_vmin ;
DLC.KINEMATIC.r_tube_l_x_vmin = r_tube_l_x_vmin ;
DLC.KINEMATIC.r_tube_l_y_vmin = r_tube_l_y_vmin ;

DLC.KINEMATIC.tip_tongue_x_offset = tip_tongue_x_offset;
DLC.KINEMATIC.tip_tongue_y_offset = tip_tongue_y_offset ;
DLC.KINEMATIC.r_tongue_x_offset = r_tongue_x_offset ;
DLC.KINEMATIC.r_tongue_y_offset = r_tongue_y_offset ;
DLC.KINEMATIC.l_tongue_x_offset = l_tongue_x_offset ;
DLC.KINEMATIC.l_tongue_y_offset = l_tongue_y_offset ;
DLC.KINEMATIC.mid_tongue_x_offset = mid_tongue_x_offset ;
DLC.KINEMATIC.mid_tongue_y_offset = mid_tongue_y_offset ;
DLC.KINEMATIC.r_nose_x_offset = r_nose_x_offset ;
DLC.KINEMATIC.r_nose_y_offset = r_nose_y_offset ;
DLC.KINEMATIC.l_nose_x_offset = l_nose_x_offset ;
DLC.KINEMATIC.l_nose_y_offset = l_nose_y_offset ;
DLC.KINEMATIC.l_food_x_offset = l_food_x_offset ;
DLC.KINEMATIC.l_food_y_offset = l_food_y_offset ;
DLC.KINEMATIC.l_tube_r_x_offset = l_tube_r_x_offset ;
DLC.KINEMATIC.l_tube_r_y_offset = l_tube_r_y_offset ;
DLC.KINEMATIC.l_tube_l_x_offset = l_tube_l_x_offset ;
DLC.KINEMATIC.l_tube_l_y_offset = l_tube_l_y_offset ;
DLC.KINEMATIC.r_food_x_offset = r_food_x_offset ;
DLC.KINEMATIC.r_food_y_offset = r_food_y_offset ;
DLC.KINEMATIC.r_tube_r_x_offset = r_tube_r_x_offset ;
DLC.KINEMATIC.r_tube_r_y_offset = r_tube_r_y_offset ;
DLC.KINEMATIC.r_tube_l_x_offset = r_tube_l_x_offset ;
DLC.KINEMATIC.r_tube_l_y_offset = r_tube_l_y_offset ;

if EXPERIMENT_PARAMS.flag_figure_debug == 1
    figure;
    hold on;
    subplot(1,6,1)
    boxplot(d_lick_max, {'D max'});
    ylabel('mm');
    ylim([0 inf])
    subplot(1,6,2)
    boxplot(v_lick_max, {'V max'});
    ylabel('mm/s');
    ylim([0 inf])
    subplot(1,6,3)
    boxplot(v_lick_min, {'V min'});
    ylabel('mm/s');
    ylim([-inf 0])
    subplot(1,6,4)
    boxplot(time_lick_duration*1000, {'Duration'});
    ylabel('ms');
    ylim([0 inf])
    subplot(1,6,5)
    boxplot(ILI_bout*1000, {'ILI'});
    ylabel('ms');
    ylim([0 inf])
    subplot(1,6,6)
    boxplot(ILR_bout, {'ILR'});
    ylabel('hz');
    ylim([0 inf])
    sgtitle([EXPERIMENT_PARAMS.file_name], 'interpreter', 'none')
    ESN_Beautify_Plot(gcf, [20 10])
end

fprintf(' --> Completed. \n')
end

%% Function: Geometrization
function [DLC, EXPERIMENT_PARAMS] = geometrization(DLC,EXPERIMENT_PARAMS, params, funcs)
fprintf('Geometrizing frames ...');

tip_tongue_x = DLC.POINTS.tip_tongue_x;
tip_tongue_y = DLC.POINTS.tip_tongue_y;
r_tongue_x = DLC.POINTS.r_tongue_x;
r_tongue_y = DLC.POINTS.r_tongue_y;
l_tongue_x = DLC.POINTS.l_tongue_x;
l_tongue_y = DLC.POINTS.l_tongue_y;
mid_tongue_x = DLC.POINTS.mid_tongue_x;
mid_tongue_y = DLC.POINTS.mid_tongue_y;
l_tube_r_x = DLC.POINTS.l_tube_r_x;
l_tube_r_y = DLC.POINTS.l_tube_r_y;
l_tube_l_x = DLC.POINTS.l_tube_l_x;
l_tube_l_y = DLC.POINTS.l_tube_l_y;
l_food_x = DLC.POINTS.l_food_x;
l_food_y = DLC.POINTS.l_food_y;
r_tube_r_x = DLC.POINTS.r_tube_r_x;
r_tube_r_y = DLC.POINTS.r_tube_r_y;
r_tube_l_x = DLC.POINTS.r_tube_l_x;
r_tube_l_y = DLC.POINTS.r_tube_l_y;
r_food_x = DLC.POINTS.r_food_x;
r_food_y = DLC.POINTS.r_food_y;
num_lick = DLC.IND.num_lick;
ind_d_lick_max = DLC.IND.ind_d_lick_max;

for counter_lick = 1:1:num_lick
    counter_frame = ind_d_lick_max(counter_lick);

    geo_tongue = polyshape([tip_tongue_x(counter_frame),l_tongue_x(counter_frame),...
        mid_tongue_x(counter_frame),r_tongue_x(counter_frame)],[tip_tongue_y(counter_frame),...
        l_tongue_y(counter_frame),mid_tongue_y(counter_frame),r_tongue_y(counter_frame)], 'Simplify', false);

    geo_r_tube_empty = polyshape([r_food_x(counter_frame),r_tube_r_x(counter_frame),...
        r_tube_r_x(counter_frame),r_tube_l_x(counter_frame), r_tube_l_x(counter_frame)],[r_food_y(counter_frame),...
        r_food_y(counter_frame), r_tube_r_y(counter_frame),r_tube_l_y(counter_frame), r_food_y(counter_frame)], 'Simplify', false);

    geo_r_tube_full = polyshape([r_tube_r_x(counter_frame), r_tube_r_x(counter_frame),...
        r_tube_l_x(counter_frame), r_tube_l_x(counter_frame), r_food_x(counter_frame)],...
        [r_food_y(counter_frame), r_tube_r_y(counter_frame) + abs(max(r_food_y)-min(r_food_y)), ...
        r_tube_l_y(counter_frame) + abs(max(r_food_y)-min(r_food_y)), r_food_y(counter_frame), r_food_y(counter_frame)], 'Simplify', false);

    geo_inter_tongue_r_tube_empty = intersect(geo_tongue, geo_r_tube_empty);

    geo_inter_tongue_r_tube_full = intersect(geo_tongue, geo_r_tube_full);

    geo_all = [geo_tongue geo_r_tube_empty geo_r_tube_full ...
        geo_inter_tongue_r_tube_empty geo_inter_tongue_r_tube_full];

    geo_l_tube_empty = polyshape([l_food_x(counter_frame),l_tube_r_x(counter_frame),...
        l_tube_r_x(counter_frame),l_tube_l_x(counter_frame), l_tube_l_x(counter_frame)],[l_food_y(counter_frame),...
        l_food_y(counter_frame), l_tube_r_y(counter_frame),l_tube_l_y(counter_frame), l_food_y(counter_frame)], 'Simplify', false);

    geo_l_tube_full = polyshape([l_tube_r_x(counter_frame), l_tube_r_x(counter_frame),...
        l_tube_l_x(counter_frame), l_tube_l_x(counter_frame), l_food_x(counter_frame)],...
        [l_food_y(counter_frame), l_tube_r_y(counter_frame) - abs(max(l_food_y)-min(l_food_y)), ...
        l_tube_l_y(counter_frame) - abs(max(l_food_y)-min(l_food_y)), l_food_y(counter_frame), l_food_y(counter_frame)], 'Simplify', false);

    geo_inter_tongue_l_tube_empty = intersect(geo_tongue, geo_l_tube_empty);

    geo_inter_tongue_l_tube_full = intersect(geo_tongue, geo_l_tube_full);

    geo_all = [geo_tongue geo_r_tube_empty geo_r_tube_full geo_l_tube_empty geo_l_tube_full...
        geo_inter_tongue_r_tube_empty geo_inter_tongue_r_tube_full geo_inter_tongue_l_tube_empty...
        geo_inter_tongue_l_tube_full ];


    [cent_tongue_r_tube_empty_x(counter_lick, 1), cent_tongue_r_tube_empty_y(counter_lick, 1) ] = centroid(geo_inter_tongue_r_tube_empty);
    [cent_tongue_r_tube_full_x(counter_lick, 1), cent_tongue_r_tube_full_y(counter_lick, 1)] = centroid(geo_inter_tongue_r_tube_full);
    bool_overlaps_all = overlaps(geo_all);
    bool_tongue_r_tube_empty(counter_lick, 1) = bool_overlaps_all(1,2);
    bool_tongue_r_tube_full(counter_lick, 1) = bool_overlaps_all(1,3);
    area_tongue(counter_lick, 1) = area(geo_tongue);
    area_r_tube_empty(counter_lick, 1) = area(geo_r_tube_empty);
    area_r_tube_full(counter_lick, 1) = area(geo_r_tube_full);
    area_inter_tongue_r_tube_empty(counter_lick, 1) = area(geo_inter_tongue_r_tube_empty);
    area_inter_tongue_r_tube_full(counter_lick, 1) = area(geo_inter_tongue_r_tube_full);
    [cent_tongue_l_tube_empty_x(counter_lick, 1),cent_tongue_l_tube_empty_y(counter_lick, 1)]  = centroid(geo_inter_tongue_l_tube_empty);
    [cent_tongue_l_tube_full_x(counter_lick, 1),cent_tongue_l_tube_full_y(counter_lick, 1)] = centroid(geo_inter_tongue_l_tube_full);
    bool_tongue_l_tube_empty(counter_lick, 1)= bool_overlaps_all(1,4);
    bool_tongue_l_tube_full(counter_lick, 1) = bool_overlaps_all(1,5);
    area_l_tube_empty(counter_lick, 1) = area(geo_l_tube_empty);
    area_l_tube_full(counter_lick, 1) = area(geo_l_tube_full);
    area_inter_tongue_l_tube_empty(counter_lick, 1) = area(geo_inter_tongue_l_tube_empty);
    area_inter_tongue_l_tube_full(counter_lick, 1) = area(geo_inter_tongue_l_tube_full);

end

DLC.GEO.area_tongue = area_tongue;
DLC.GEO.area_r_tube_empty = area_r_tube_empty;
DLC.GEO.area_r_tube_full = area_r_tube_full;
DLC.GEO.area_inter_tongue_r_tube_empty = area_inter_tongue_r_tube_empty;
DLC.GEO.area_inter_tongue_r_tube_full = area_inter_tongue_r_tube_full;
DLC.GEO.bool_overlaps_all = bool_overlaps_all;
DLC.GEO.bool_tongue_r_tube_empty = bool_tongue_r_tube_empty;
DLC.GEO.bool_tongue_r_tube_full = bool_tongue_r_tube_full;
DLC.GEO.cent_tongue_r_tube_empty_x = cent_tongue_r_tube_empty_x;
DLC.GEO.cent_tongue_r_tube_empty_y = cent_tongue_r_tube_empty_y;
DLC.GEO.cent_tongue_r_tube_full_x = cent_tongue_r_tube_full_x;
DLC.GEO.cent_tongue_r_tube_full_y = cent_tongue_r_tube_full_y;
DLC.GEO.area_l_tube_empty = area_l_tube_empty;
DLC.GEO.area_l_tube_full = area_l_tube_full;
DLC.GEO.area_inter_tongue_l_tube_empty = area_inter_tongue_l_tube_empty;
DLC.GEO.area_inter_tongue_l_tube_full = area_inter_tongue_l_tube_full;
DLC.GEO.bool_tongue_l_tube_empty = bool_tongue_l_tube_empty;
DLC.GEO.bool_tongue_l_tube_full = bool_tongue_l_tube_full;
DLC.GEO.cent_tongue_l_tube_empty_x = cent_tongue_l_tube_empty_x;
DLC.GEO.cent_tongue_l_tube_empty_y = cent_tongue_l_tube_empty_y;
DLC.GEO.cent_tongue_l_tube_full_x = cent_tongue_l_tube_full_x;
DLC.GEO.cent_tongue_l_tube_full_y = cent_tongue_l_tube_full_y;

fprintf(' --> Completed. \n')
end

%% Function: Sort licks
function [DLC, EXPERIMENT_PARAMS] = sort_licks(DLC,EXPERIMENT_PARAMS, params, funcs)
fprintf('Sorting licks ...');

bool_tongue_r_tube_full = DLC.GEO.bool_tongue_r_tube_full;
bool_tongue_l_tube_full = DLC.GEO.bool_tongue_l_tube_full;
area_r_tube = DLC.GEO.area_r_tube_empty + DLC.GEO.area_r_tube_full;
area_l_tube = DLC.GEO.area_l_tube_empty + DLC.GEO.area_l_tube_full;
r_tube_food_perc = DLC.GEO.area_r_tube_full./area_r_tube;
l_tube_food_perc = DLC.GEO.area_l_tube_full./area_l_tube;

tip_tongue_x = DLC.POINTS.tip_tongue_x;
tip_tongue_y = DLC.POINTS.tip_tongue_y;
r_tongue_x = DLC.POINTS.r_tongue_x;
r_tongue_y = DLC.POINTS.r_tongue_y;
l_tongue_x = DLC.POINTS.l_tongue_x;
l_tongue_y = DLC.POINTS.l_tongue_y;

r_nose_x =  DLC.POINTS.r_nose_x;
r_nose_y = DLC.POINTS.r_nose_y;
l_nose_x =  DLC.POINTS.l_nose_x;
l_nose_y = DLC.POINTS.l_nose_y;

r_tube_r_x = DLC.POINTS.r_tube_r_x;
r_tube_r_y = DLC.POINTS.r_tube_r_y;
r_tube_l_x = DLC.POINTS.r_tube_l_x;
r_tube_l_y = DLC.POINTS.r_tube_l_y;
% r_food_y = DLC.POINTS.r_food_y;
% r_tube_food_perc = DLC.data.r_food_y./(DLC.FILE.vid_height - (DLC.data.r_tube_r_y + DLC.data.r_tube_l_y)/2);

l_tube_r_x = DLC.POINTS.l_tube_r_x;
l_tube_r_y = DLC.POINTS.l_tube_r_y;
l_tube_l_x = DLC.POINTS.l_tube_l_x;
l_tube_l_y = DLC.POINTS.l_tube_l_y;
% l_food_y = DLC.POINTS.l_food_y;
% l_tube_food_perc = DLC.data.l_food_y./((DLC.data.l_tube_r_y + DLC.data.l_tube_l_y)/2);

ind_d_lick_max = DLC.IND.ind_d_lick_max;
ind_lick_onset = DLC.IND.ind_lick_onset;
ind_lick_offset = DLC.IND.ind_lick_offset;
time_1K = DLC.TIME.time_1K';

% Classify inner, outer edge, under tube, and grooming
is_r_inner_tube_lick = false(size(ind_lick_onset));
is_l_inner_tube_lick = false(size(ind_lick_onset));
is_r_outer_edge_lick = false(size(ind_lick_onset));
is_l_outer_edge_lick = false(size(ind_lick_onset));
is_r_under_tube_lick = false(size(ind_lick_onset));
is_l_under_tube_lick = false(size(ind_lick_onset));
is_grooming_lick = false(size(ind_lick_onset));

% tongue markers relative to tube markers
tip_crosses_far_edge_r_tube = tip_tongue_x > r_tube_r_x + std(r_tube_r_x);
tip_crosses_far_edge_l_tube = tip_tongue_x > l_tube_r_x + std(l_tube_r_x);
tip_crosses_close_edge_r_tube = tip_tongue_x < r_tube_l_x + std(r_tube_l_x);
tip_crosses_close_edge_l_tube = tip_tongue_x < l_tube_l_x + std(l_tube_l_x);

tip_x_in_r_tube = tip_tongue_x < r_tube_r_x & tip_tongue_x > r_tube_l_x; 
r_x_in_r_tube = r_tongue_x < r_tube_r_x & r_tongue_x > r_tube_l_x;
l_x_in_r_tube = l_tongue_x < r_tube_r_x & l_tongue_x > r_tube_l_x; 
tip_x_in_l_tube = tip_tongue_x < l_tube_r_x & tip_tongue_x > l_tube_l_x; 
r_x_in_l_tube = r_tongue_x < l_tube_r_x & r_tongue_x > l_tube_l_x; 
l_x_in_l_tube = l_tongue_x < l_tube_r_x & l_tongue_x > l_tube_l_x;  
tip_y_in_r_tube = tip_tongue_y > (r_tube_r_y+r_tube_l_y)/2 - 1; 
r_y_in_r_tube = r_tongue_y > (r_tube_r_y+r_tube_l_y)/2 - 1; 
l_y_in_r_tube = l_tongue_y > (r_tube_r_y+r_tube_l_y)/2 - 1; 
tip_y_in_l_tube = tip_tongue_y < (l_tube_r_y+l_tube_l_y)/2 + 1; 
r_y_in_l_tube = r_tongue_y < (l_tube_r_y+l_tube_l_y)/2 + 1; 
l_y_in_l_tube = l_tongue_y < (l_tube_r_y+l_tube_l_y)/2 + 1; 

tip_in_r_tube = tip_x_in_r_tube & tip_y_in_r_tube; 
r_in_r_tube = r_x_in_r_tube & r_y_in_r_tube; 
l_in_r_tube = l_x_in_r_tube & l_y_in_r_tube; 
tip_in_l_tube = tip_x_in_l_tube & tip_y_in_l_tube;
r_in_l_tube = r_x_in_l_tube & r_y_in_l_tube; 
l_in_l_tube = l_x_in_l_tube & l_y_in_l_tube;

for i = 1:length(ind_lick_onset)
    done = false;
    ind = ind_lick_onset(i):ind_lick_offset(i);

%   INNER TUBE LICK DEFINITION
%   - there are frames where the tip, l, r tongue markers are all in 
%     the tube during the lick
%   - when tongue markers are in the y-coordinate range of the tube, 
%     they don't ever cross either of the edges of the tube

    % r inner tube lick
    tongue_in_tube = tip_in_r_tube(ind) & ...
        r_in_r_tube(ind) & ...
        l_in_r_tube(ind);
    tongue_crosses_edge = (tip_y_in_r_tube(ind) & ~tip_x_in_r_tube(ind)) | ...
        (r_y_in_r_tube(ind) & ~r_x_in_r_tube(ind)) | ...
        (l_y_in_r_tube(ind) & ~l_x_in_r_tube(ind));
    if isempty(find(tongue_crosses_edge, 1)) && ~isempty(find(tongue_in_tube, 1))
        is_r_inner_tube_lick(i) = true;
        done = true;
    end

    % l inner tube lick
    tongue_in_tube = tip_in_l_tube(ind) & ...
        r_in_l_tube(ind) & ...
        l_in_l_tube(ind);
    tongue_crosses_edge = (tip_y_in_l_tube(ind) & ~tip_x_in_l_tube(ind)) | ...
        (r_y_in_l_tube(ind) & ~r_x_in_l_tube(ind)) | ...
        (l_y_in_l_tube(ind) & ~l_x_in_l_tube(ind));
    if isempty(find(tongue_crosses_edge, 1)) && ~isempty(find(tongue_in_tube, 1))
        is_l_inner_tube_lick(i) = true;
        done = true;
    end

%   UNDER TUBE LICK DEFINITION
%   - not an inner tube
%   - the tip tongue marker crosses the close edge of the tube when it is
%     in the y_coordinate range of the tube
%   - otherwise, there are frames where the tip tongue marker is in the 
%     y_coordinate range of the tube but never crosses the far edge of the tube
    tongue_crosses_close_edge = tip_crosses_close_edge_r_tube(ind) & tip_y_in_r_tube(ind);
    tongue_crosses_far_edge = tip_crosses_far_edge_r_tube(ind) & tip_y_in_r_tube(ind);
    % r under tube lick
    if ~done && (~isempty(find(tongue_crosses_close_edge, 1))) || ...
            (~done && (isempty(find(tongue_crosses_far_edge, 1))) && ~isempty(find(tip_y_in_r_tube(ind), 1)))
        is_r_under_tube_lick(i) = true;
        done = true;
    end

    tongue_crosses_close_edge = tip_crosses_close_edge_l_tube(ind) & tip_y_in_l_tube(ind);
    tongue_crosses_far_edge = tip_crosses_far_edge_l_tube(ind) & tip_y_in_l_tube(ind);
    % l under tube lick
    if ~done && (~isempty(find(tongue_crosses_close_edge, 1))) || ...
            (~done && (isempty(find(tongue_crosses_far_edge, 1))) && ~isempty(find(tip_y_in_l_tube(ind), 1)))
         is_l_under_tube_lick(i) = true;
        done = true;
    end

%   OUTER EDGE LICK DEFINITION
%   - not an inner tube or under tube lick
%   - there are frames where the tip tongue marker has crossed the far 
%     edge when it is in the y-coordinate range of the tube

    % r outer edge lick
    tongue_crosses_far_edge = tip_crosses_far_edge_r_tube(ind) & tip_y_in_r_tube(ind);
    if ~done && ~isempty(find(tongue_crosses_far_edge, 1))
        is_r_outer_edge_lick(i) = true;
        done = true;
    end

     % l outer edge lick
    tongue_crosses_far_edge = tip_crosses_far_edge_l_tube(ind) & tip_y_in_l_tube(ind);
    if ~done && ~isempty(find(tongue_crosses_far_edge, 1))
        is_l_outer_edge_lick(i) = true;
        done = true;
    end

    % GROOMING LICK DEFINITION
    % not inner, outer edge, or under tube lick
    if ~done
        is_grooming_lick(i) = true;
    end
end

% grooming
ind_grooming_lick = find(is_grooming_lick);
% inner reward/noreward
is_r_reward_inner_tube_lick = bool_tongue_r_tube_full == 1 & is_r_inner_tube_lick;
is_r_noreward_inner_tube_lick = is_r_inner_tube_lick & ~is_r_reward_inner_tube_lick;
is_l_reward_inner_tube_lick = bool_tongue_l_tube_full == 1 & is_l_inner_tube_lick;
is_l_noreward_inner_tube_lick = is_l_inner_tube_lick & ~is_l_reward_inner_tube_lick;
ind_r_reward_inner_tube_lick = find(is_r_reward_inner_tube_lick);
ind_r_noreward_inner_tube_lick = find(is_r_noreward_inner_tube_lick);
ind_l_reward_inner_tube_lick = find(is_l_reward_inner_tube_lick);
ind_l_noreward_inner_tube_lick = find(is_l_noreward_inner_tube_lick);
% outer edge reward/noreward
is_r_reward_outer_edge_lick = bool_tongue_r_tube_full == 1 & r_tube_food_perc >= 0.90 & is_r_outer_edge_lick;
is_r_noreward_outer_edge_lick = is_r_outer_edge_lick & ~is_r_reward_outer_edge_lick;
is_l_reward_outer_edge_lick = bool_tongue_l_tube_full == 1 & l_tube_food_perc >= 0.90 & is_l_outer_edge_lick;
is_l_noreward_outer_edge_lick = is_l_outer_edge_lick & ~is_l_reward_outer_edge_lick;
ind_r_reward_outer_edge_lick = find(is_r_reward_outer_edge_lick);
ind_r_noreward_outer_edge_lick = find(is_r_noreward_outer_edge_lick);
ind_l_reward_outer_edge_lick = find(is_l_reward_outer_edge_lick);
ind_l_noreward_outer_edge_lick = find(is_l_noreward_outer_edge_lick);
% under tube reward/noreward
is_r_reward_under_tube_lick = bool_tongue_r_tube_full == 1 & r_tube_food_perc >= 0.90 & is_r_under_tube_lick;
is_r_noreward_under_tube_lick = is_r_under_tube_lick & ~is_r_reward_under_tube_lick;
is_l_reward_under_tube_lick = bool_tongue_l_tube_full == 1 & l_tube_food_perc >= 0.90 & is_l_under_tube_lick;
is_l_noreward_under_tube_lick = is_l_under_tube_lick & ~is_l_reward_under_tube_lick;
ind_r_reward_under_tube_lick = find(is_r_reward_under_tube_lick);
ind_r_noreward_under_tube_lick = find(is_r_noreward_under_tube_lick);
ind_l_reward_under_tube_lick = find(is_l_reward_under_tube_lick);
ind_l_noreward_under_tube_lick = find(is_l_noreward_under_tube_lick);

% inds for reward-driven licks
ind_lick_onset_r_reward = ind_lick_onset(sort([ind_r_reward_inner_tube_lick;ind_r_reward_outer_edge_lick;ind_r_reward_under_tube_lick; ind_r_noreward_inner_tube_lick;ind_r_noreward_outer_edge_lick;ind_r_noreward_under_tube_lick]));
time_lick_onset_r_reward = time_1K(ind_lick_onset_r_reward);
ind_lick_offset_r_reward = ind_lick_offset(sort([ind_r_reward_inner_tube_lick;ind_r_reward_outer_edge_lick;ind_r_reward_under_tube_lick; ind_r_noreward_inner_tube_lick;ind_r_noreward_outer_edge_lick;ind_r_noreward_under_tube_lick]));
time_lick_offset_r_reward = time_1K(ind_lick_offset_r_reward);
ind_lick_onset_l_reward = ind_lick_onset(sort([ind_l_reward_inner_tube_lick;ind_l_reward_outer_edge_lick;ind_l_reward_under_tube_lick; ind_l_noreward_inner_tube_lick;ind_l_noreward_outer_edge_lick;ind_l_noreward_under_tube_lick]));
time_lick_onset_l_reward = time_1K(ind_lick_onset_l_reward);
ind_lick_offset_l_reward = ind_lick_offset(sort([ind_l_reward_inner_tube_lick;ind_l_reward_outer_edge_lick;ind_l_reward_under_tube_lick; ind_l_noreward_inner_tube_lick;ind_l_noreward_outer_edge_lick;ind_l_noreward_under_tube_lick]));
time_lick_offset_l_reward = time_1K(ind_lick_offset_l_reward);

% add vars to DLC
DLC.CLASS.is_grooming_lick = is_grooming_lick;
DLC.CLASS.is_r_reward_inner_tube_lick = is_r_reward_inner_tube_lick;
DLC.CLASS.is_r_reward_outer_edge_lick = is_r_reward_outer_edge_lick;
DLC.CLASS.is_r_reward_under_tube_lick = is_r_reward_under_tube_lick;
DLC.CLASS.is_r_noreward_inner_tube_lick = is_r_noreward_inner_tube_lick;
DLC.CLASS.is_r_noreward_outer_edge_lick = is_r_noreward_outer_edge_lick;
DLC.CLASS.is_r_noreward_under_tube_lick = is_r_noreward_under_tube_lick;

DLC.CLASS.ind_lick_onset_r_reward = ind_lick_onset_r_reward;
DLC.TIME.time_lick_onset_r_reward = time_lick_onset_r_reward;
DLC.CLASS.ind_lick_offset_r_reward = ind_lick_offset_r_reward;
DLC.TIME.time_lick_offset_r_reward = time_lick_offset_r_reward;

DLC.CLASS.is_l_reward_inner_tube_lick = is_l_reward_inner_tube_lick;
DLC.CLASS.is_l_reward_outer_edge_lick = is_l_reward_outer_edge_lick;
DLC.CLASS.is_l_reward_under_tube_lick = is_l_reward_under_tube_lick;
DLC.CLASS.is_l_noreward_inner_tube_lick = is_l_noreward_inner_tube_lick;
DLC.CLASS.is_l_noreward_outer_edge_lick = is_l_noreward_outer_edge_lick;
DLC.CLASS.is_l_noreward_under_tube_lick = is_l_noreward_under_tube_lick;

DLC.CLASS.ind_lick_onset_l_reward = ind_lick_onset_l_reward;
DLC.TIME.time_lick_onset_l_reward = time_lick_onset_l_reward;
DLC.CLASS.ind_lick_offset_l_reward = ind_lick_offset_l_reward;
DLC.TIME.time_lick_offset_l_reward = time_lick_offset_l_reward;

if EXPERIMENT_PARAMS.flag_figure_debug == 1
    figure;
    hold on;
    % plot(tip_tongue_x(ind_d_lick_max),tip_tongue_y(ind_d_lick_max), '.k');
    plot(tip_tongue_x(ind_d_lick_max(ind_grooming_lick)),tip_tongue_y(ind_d_lick_max(ind_grooming_lick)), 'ob');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_noreward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_noreward_inner_tube_lick)), 'og');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_reward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_reward_inner_tube_lick)), 'or');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_noreward_outer_edge_lick)),tip_tongue_y(ind_d_lick_max(ind_r_noreward_outer_edge_lick)), 'oc');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_reward_outer_edge_lick)),tip_tongue_y(ind_d_lick_max(ind_r_reward_outer_edge_lick)), 'om');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_noreward_under_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_noreward_under_tube_lick)), 'oy');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_reward_under_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_reward_under_tube_lick)), 'o', 'Color', [0.7 0.7 0.7]);
    plot(r_nose_x, r_nose_y,'ok');
    plot(l_nose_x, l_nose_y,'ok');
    plot(r_tube_r_x(ind_d_lick_max),r_tube_r_y(ind_d_lick_max),'sk');
    plot(r_tube_l_x(ind_d_lick_max),r_tube_l_y(ind_d_lick_max),'sk');
     plot(tip_tongue_x(ind_d_lick_max(ind_l_noreward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_noreward_inner_tube_lick)), 'og');
    plot(tip_tongue_x(ind_d_lick_max(ind_l_reward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_reward_inner_tube_lick)), 'or');
    plot(tip_tongue_x(ind_d_lick_max(ind_l_noreward_outer_edge_lick)),tip_tongue_y(ind_d_lick_max(ind_l_noreward_outer_edge_lick)), 'oc');
    plot(tip_tongue_x(ind_d_lick_max(ind_l_reward_outer_edge_lick)),tip_tongue_y(ind_d_lick_max(ind_l_reward_outer_edge_lick)), 'om');
    plot(tip_tongue_x(ind_d_lick_max(ind_l_noreward_under_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_noreward_under_tube_lick)), 'oy');
    plot(tip_tongue_x(ind_d_lick_max(ind_l_reward_under_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_reward_under_tube_lick)), 'o', 'Color', [0.7 0.7 0.7]);
    plot(l_tube_r_x(ind_d_lick_max), l_tube_r_y(ind_d_lick_max),'sk');
    plot(l_tube_l_x(ind_d_lick_max), l_tube_l_y(ind_d_lick_max),'sk');
    xlabel('x position (mm)');
    ylabel('y position (mm)');
    set(gca, 'YDir','reverse')
    xlim([-10 15]);
    ylim([-25 25]);
    title([ EXPERIMENT_PARAMS.file_name ' | Lick sorting summary'], 'interpreter', 'none');
    ESN_Beautify_Plot(gcf, [20 10])
end

fprintf(' --> Completed. \n')
end

%% Function: Quantify food in tube
function [DLC, EXPERIMENT_PARAMS] = quantify_food(DLC,EXPERIMENT_PARAMS, params, funcs)
fprintf('Quantifying food in tube ...');
num_lick = DLC.IND.num_lick;
num_bout = DLC.IND.num_bout;
ind_lick_onset_str_bout = DLC.IND.ind_lick_onset_str_bout;
ind_lick_onset_end_bout = DLC.IND.ind_lick_onset_end_bout ;
ind_lick_onset = DLC.IND.ind_lick_onset;
ind_lick_offset = DLC.IND.ind_lick_offset;

r_tube_r_y = DLC.POINTS.r_tube_r_y;
r_food_x = DLC.POINTS.r_food_x;
r_food_y = DLC.POINTS.r_food_y;
r_food_x(r_food_y < 0) = nan;
r_food_y(r_food_y < 0) = nan;
r_food_max = max(r_food_y);
r_food_min = min(r_food_y);
r_tube_food = (r_food_y - r_food_max)./(r_tube_r_y - r_food_max);


l_tube_r_y = DLC.POINTS.l_tube_r_y;
l_food_x = DLC.POINTS.l_food_x;
l_food_y = DLC.POINTS.l_food_y;
l_food_x(l_food_y > 0) = nan;
l_food_y(l_food_y > 0) = nan;
l_food_max = max(l_food_y);
l_food_min = min(l_food_y);
l_tube_food = (l_food_y - l_food_min)./(l_tube_r_y - l_food_min);

for counter_lick =  1 : 1 : num_lick
    inds_ = ind_lick_onset(counter_lick):ind_lick_offset(counter_lick);
    ind_onset_(counter_lick) = DLC.IND.ind_lick_onset(counter_lick);
    ind_vmax_(counter_lick) = DLC.IND.ind_v_lick_max(counter_lick);
    ind_dmax_(counter_lick) = DLC.IND.ind_d_lick_max(counter_lick);
    ind_vmin_(counter_lick) = DLC.IND.ind_v_lick_min(counter_lick);
    ind_offset_(counter_lick) = DLC.IND.ind_lick_offset(counter_lick);
    for counter_inds_ = 1 : 1: length(inds_)
        r_tube_food_(counter_inds_) = r_tube_food(inds_(counter_inds_));
        r_tube_food_lick_onset(counter_lick,1) = r_tube_food_(1) ;
        r_tube_food_lick_offset(counter_lick,1) = r_tube_food_(end) ;
        l_tube_food_(counter_inds_) = l_tube_food(inds_(counter_inds_));
        l_tube_food_lick_onset(counter_lick,1) = l_tube_food_(1) ;
        l_tube_food_lick_offset(counter_lick,1) = l_tube_food_(end) ;

    end
end

for counter_bout =  1 : 1 : num_bout
    inds_ = ind_lick_onset_str_bout(counter_bout):ind_lick_onset_end_bout(counter_bout);
    for counter_inds_ = 1 : 1: length(inds_)
        r_tube_food_(counter_inds_) = r_tube_food(inds_(counter_inds_));
        r_tube_food_bout_start(counter_bout, 1) = r_tube_food_(1);
        r_tube_food_bout_end(counter_bout, 1) = r_tube_food_(end);
        l_tube_food_(counter_inds_) = l_tube_food(inds_(counter_inds_));
        l_tube_food_bout_start(counter_bout, 1) = l_tube_food_(1);
        l_tube_food_bout_end(counter_bout, 1) = l_tube_food_(end);
    end
end

DLC.FOOD.r_tube_food = r_tube_food;
DLC.FOOD.r_food_x = r_food_x;
DLC.FOOD.r_food_y = r_food_y;
DLC.FOOD.l_tube_food = l_tube_food;
DLC.FOOD.l_food_x = l_food_x;
DLC.FOOD.l_food_y = l_food_y;

DLC.FOOD.r_tube_food_lick_onset = r_tube_food_lick_onset;
DLC.FOOD.r_tube_food_lick_offset = r_tube_food_lick_offset;
DLC.FOOD.l_tube_food_lick_onset = l_tube_food_lick_onset;
DLC.FOOD.l_tube_food_lick_offset = l_tube_food_lick_offset;

DLC.FOOD.r_tube_food_bout_start = r_tube_food_bout_start;
DLC.FOOD.r_tube_food_bout_end = r_tube_food_bout_end;
DLC.FOOD.l_tube_food_bout_start = l_tube_food_bout_start;
DLC.FOOD.l_tube_food_bout_end = l_tube_food_bout_end;
fprintf(' --> Completed. \n')
end

%% Function: Detect harvest str, end, num, times, and duration
function [DLC, EXPERIMENT_PARAMS] = detect_harvest(DLC,EXPERIMENT_PARAMS, params, funcs)
fprintf('Detecting harvest ...');

d_tip = DLC.KINEMATIC.d_tip;
time_1K = DLC.TIME.time_1K';
ind_lick_onset = DLC.IND.ind_lick_onset;
is_grooming_lick = DLC.CLASS.is_grooming_lick;
ind_lick_onset_grooming = (ind_lick_onset(is_grooming_lick));
ind_lick_onset_reward = (ind_lick_onset(~is_grooming_lick));
ind_lick_onset_r_reward = DLC.CLASS.ind_lick_onset_r_reward;
time_lick_onset_r_reward = DLC.TIME.time_lick_onset_r_reward;
ind_lick_onset_l_reward = DLC.CLASS.ind_lick_onset_l_reward;
time_lick_onset_l_reward = DLC.TIME.time_lick_onset_l_reward;
r_tube_food = DLC.FOOD.r_tube_food;
l_tube_food = DLC.FOOD.l_tube_food;

ind_lick_onset_str_bout = DLC.IND.ind_lick_onset_str_bout;
ind_lick_onset_end_bout = DLC.IND.ind_lick_onset_end_bout;

% Check bouts that start or end with grooming lick
% ind_lick_onset_str_bout_grooming = ismember(ind_lick_onset_str_bout, ind_lick_onset_grooming);
% ind_lick_onset_end_bout_grooming = ismember(ind_lick_onset_end_bout, ind_lick_onset_grooming);

for counter_bout = 1 : DLC.IND.num_bout
    % harvest start
    if isempty(find(ind_lick_onset_str_bout(counter_bout) == ind_lick_onset_grooming,1))
        ind_lick_onset_str_harvest(counter_bout,1) = ind_lick_onset_str_bout(counter_bout);
    elseif ~isempty(find(ind_lick_onset_str_bout(counter_bout) == ind_lick_onset_grooming,1))
        shift_ind = 0;
        bool_shift = 1;
        while bool_shift == 1
            if ~isempty(find(ind_lick_onset(find(ind_lick_onset_str_bout(counter_bout) == ind_lick_onset,1) + shift_ind) == ind_lick_onset_grooming, 1)) && ...
                    ind_lick_onset(find(ind_lick_onset_str_bout(counter_bout) == ind_lick_onset,1) + shift_ind) ~= ind_lick_onset(end)
                shift_ind = shift_ind + 1;
            elseif isempty(find(ind_lick_onset(find(ind_lick_onset_str_bout(counter_bout) == ind_lick_onset,1) + shift_ind) == ind_lick_onset_grooming, 1)) || ...
                    ind_lick_onset(find(ind_lick_onset_str_bout(counter_bout) == ind_lick_onset,1) + shift_ind) == ind_lick_onset(end)
                bool_shift = 0;
            end
        end
        ind_lick_onset_str_harvest(counter_bout,1) = ind_lick_onset(find(ind_lick_onset_str_bout(counter_bout) == ind_lick_onset,1) + shift_ind);
    end
    % harvest end
    if isempty(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset_grooming,1))
        ind_lick_onset_end_harvest(counter_bout,1) = ind_lick_onset_end_bout(counter_bout);
    elseif ~isempty(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset_grooming,1))
        shift_ind = 0;
        bool_shift = 1;
        while bool_shift == 1
            if ~isempty(find(ind_lick_onset(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset,1) - shift_ind) == ind_lick_onset_grooming, 1)) && ...
                    ind_lick_onset(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset,1) - shift_ind) ~= ind_lick_onset_str_harvest(counter_bout) && ...
                    ind_lick_onset(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset,1) - shift_ind) ~= ind_lick_onset(1)
                shift_ind = shift_ind + 1;
            elseif isempty(find(ind_lick_onset(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset,1) - shift_ind) == ind_lick_onset_grooming, 1)) || ...
                    ind_lick_onset(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset,1) - shift_ind) == ind_lick_onset_str_harvest(counter_bout) || ...
                    ind_lick_onset(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset,1) - shift_ind) == ind_lick_onset(1)
                bool_shift = 0;
            end
        end
        ind_lick_onset_end_harvest(counter_bout,1) = ind_lick_onset(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset,1) - shift_ind);
    end
end

% Check harcest that start or end with grooming lick
% ind_lick_onset_str_harvest_grooming = ismember(ind_lick_onset_str_harvest, ind_lick_onset_grooming);
% ind_lick_onset_end_harvest_grooming = ismember(ind_lick_onset_end_harvest, ind_lick_onset_grooming);

% Count number of licks in harvest
for counter_harvest = 1 : length(ind_lick_onset_str_harvest)
    inds_ = ind_lick_onset_str_harvest(counter_harvest) : ind_lick_onset_end_harvest(counter_harvest);
    num_lick_harvest(counter_harvest,1) = length(find(inds_ == ind_lick_onset));
end

validity = num_lick_harvest>=3;
num_lick_harvest(~validity) = [];
ind_lick_onset_str_harvest(~validity) = [];
ind_lick_onset_end_harvest(~validity) = [];

time_lick_onset_str_harvest = time_1K(ind_lick_onset_str_harvest);
time_lick_onset_end_harvest = time_1K(ind_lick_onset_end_harvest);
time_harvest_duration = (time_1K(ind_lick_onset_end_harvest) - ...
    time_1K(ind_lick_onset_str_harvest))';

% Determine direction of bouts
num_r = [];
num_l = [];
num_g = [];

for i = 1:length(ind_lick_onset_str_bout)
    inds_lick_onset = ind_lick_onset_str_bout(i):ind_lick_onset_end_bout(i);
    num_r = [num_r; length(find(ismember(inds_lick_onset, ind_lick_onset_r_reward)))];
    num_l = [num_l; length(find(ismember(inds_lick_onset, ind_lick_onset_l_reward)))];
    num_g = [num_g; length(find(ismember(inds_lick_onset, ind_lick_onset_grooming)))];

end
is_bout_r = (num_r > num_l) & (num_r > num_g) & num_r > 2;
is_bout_l = (num_l > num_r) & (num_l > num_g) & num_l > 2;
is_bout_g = (num_g > num_r) & (num_g > num_l) & num_g > 2;

% Determine direction of harvest
num_r = [];
num_l = [];

for i = 1:length(ind_lick_onset_str_harvest)
    inds_lick_onset = ind_lick_onset_str_harvest(i):ind_lick_onset_end_harvest(i);
    num_r = [num_r; length(find(ismember(inds_lick_onset, ind_lick_onset_r_reward)))];
    num_l = [num_l; length(find(ismember(inds_lick_onset, ind_lick_onset_l_reward)))];
end
is_harvest_r = num_r > num_l;
is_harvest_l = num_l > num_r;


DLC.IND.ind_lick_onset_str_harvest = ind_lick_onset_str_harvest;
DLC.IND.ind_lick_onset_end_harvest = ind_lick_onset_end_harvest;
DLC.IND.num_lick_harvest = num_lick_harvest;
DLC.TIME.time_lick_onset_str_harvest = time_lick_onset_str_harvest;
DLC.TIME.time_lick_onset_end_harvest = time_lick_onset_end_harvest;
DLC.TIME.time_harvest_duration = time_harvest_duration;
DLC.CLASS.is_bout_r = is_bout_r;
DLC.CLASS.is_harvest_r = is_harvest_r;
DLC.CLASS.is_bout_l = is_bout_l;
DLC.CLASS.is_harvest_l = is_harvest_l;
DLC.CLASS.is_bout_g = is_bout_g;
if EXPERIMENT_PARAMS.flag_figure_debug == 1
    figure;
    hold on;
    plot(time_1K,d_tip,'-k');
    plot(time_lick_onset_str_harvest,20,'*g');
    plot(time_lick_onset_end_harvest,20,'*r');
    ylim([0 25])
    xlabel('Time (s)');
    ylabel('Displacement (mm)');
    if ~isempty(time_lick_onset_r_reward)
        plot(time_lick_onset_r_reward,19,'.r');
    end
    if ~isempty(time_lick_onset_l_reward)
        plot(time_lick_onset_l_reward,19,'.b');
    end
    yyaxis right
    plot(time_1K,r_tube_food,'.-r');
    %     plot(time_1K,l_tube_food,'.-b');
    ylabel('Reward capacity (0:Empty | 1:Full)')
    title([ EXPERIMENT_PARAMS.file_name ': Harvest | ' num2str(sum(DLC.CLASS.is_bout_l)) ' left bouts | ' num2str(sum(DLC.CLASS.is_bout_r)) ' right bouts'], 'interpreter', 'none');
    ESN_Beautify_Plot(gcf, [20 10])
end

fprintf(' --> Completed. \n')
end


%% Function: PGH_plot_sample_trace
function PGH_plot_sample_trace(DLC)
    % load alignment file for sample diff
    path_to_analyzed_tongue = DLC.FILE.path_to_analyzed;
    align_file = dir([path_to_analyzed_tongue '*aligned.mat']);
    load([path_to_analyzed_tongue align_file.name], 'align_PD');
    sample_diff = align_PD.sample_diff;
    exp_start_time = align_PD.BEHAVE_PD_xcorr_time_1K(1,1);
    
    % load SACS_ALL_DATA for trial data
    path_to_analyzed_eye = [DLC.FILE.path_to_analyzed '..' filesep 'eye' filesep];
    sac_file = dir([path_to_analyzed_eye '*RECAL.mat']);
    load([path_to_analyzed_eye sac_file.name], 'TRIALS_DATA');
    time_end_trial = TRIALS_DATA.time_end';
    time_end_trial = time_end_trial - exp_start_time + sample_diff/1000;
    
    file_name = align_file.name(1:13);
    
    % plot
    figure;
    hold on;
    plot(DLC.TIME.time_1K,DLC.KINEMATIC.d_tip,'-k');
    plot(DLC.TIME.time_d_lick_max_abs,25,'.k', 'MarkerSize', 20);
    plot(time_end_trial,25,'.r', 'MarkerSize', 20)
    ylim([0 26])
    xlabel('Time (s)');
    ylabel('Displacement (mm)');
    set(gca, 'Ycolor','k')
    yyaxis right
    plot(DLC.TIME.time_1K,ESN_smooth(DLC.FOOD.r_tube_food),'-r');
    %     plot(DLC.TIME.time_1K,DLC.TIME.l_tube_food,'.-b');
    ylim([0 1.5])
    ylabel('Reward capacity (0:Empty | 1:Full)')
    set(gca, 'Ycolor','k')
    title(file_name, 'interpreter', 'none');
    ESN_Beautify_Plot(gcf, [20 10])
end
