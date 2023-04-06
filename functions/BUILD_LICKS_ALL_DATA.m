%% MAIN Function: Build LICKS_ALL_DATA
function [LICKS_ALL_DATA, EXPERIMENT_PARAMS] = BUILD_LICKS_ALL_DATA(DLC, EXPERIMENT_PARAMS, params, funcs)
fprintf(['Building LICK_DATA_ALL: ' EXPERIMENT_PARAMS.mat_FileName ' ... ' '\n'])
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs;
%% Build validity
validity = ones(1, DLC.IND.num_lick);
is_filter_outlier = 0;
if is_filter_outlier == 1
    validity(isoutlier(DLC.KINEMATIC.d_lick_max')) = 0;
    validity(isoutlier(DLC.KINEMATIC.v_lick_max')) = 0;
    validity(isoutlier(DLC.KINEMATIC.v_lick_min')) = 0;
    validity(isoutlier(DLC.TIME.time_lick_duration')) = 0;
end
LICKS_ALL_DATA.validity = validity;

%% Build tag and lick_tag_list
lick_tag_list = params.lick.tag_name_list;
lick_tag_bout_list = params.lick.tag_bout_name_list;
lick_tag_harvest_list = params.lick.tag_harvest_name_list;

is_groom= logical(DLC.CLASS.is_grooming_lick(LICKS_ALL_DATA.validity == 1));
is_inner_tube_success =  logical(DLC.CLASS.is_r_reward_inner_tube_lick(LICKS_ALL_DATA.validity == 1) + DLC.CLASS.is_l_reward_inner_tube_lick(LICKS_ALL_DATA.validity == 1));
is_inner_tube_fail = logical(DLC.CLASS.is_r_noreward_inner_tube_lick(LICKS_ALL_DATA.validity == 1) + DLC.CLASS.is_l_noreward_inner_tube_lick(LICKS_ALL_DATA.validity == 1));
is_outer_edge_success = logical(DLC.CLASS.is_r_reward_outer_edge_lick(LICKS_ALL_DATA.validity == 1) + DLC.CLASS.is_l_reward_outer_edge_lick(LICKS_ALL_DATA.validity == 1));
is_outer_edge_fail = logical(DLC.CLASS.is_r_noreward_outer_edge_lick(LICKS_ALL_DATA.validity == 1) + DLC.CLASS.is_l_noreward_outer_edge_lick(LICKS_ALL_DATA.validity == 1));
is_under_tube_success =  logical(DLC.CLASS.is_r_reward_under_tube_lick(LICKS_ALL_DATA.validity == 1) + DLC.CLASS.is_l_reward_under_tube_lick(LICKS_ALL_DATA.validity == 1));
is_under_tube_fail = logical(DLC.CLASS.is_r_noreward_under_tube_lick(LICKS_ALL_DATA.validity == 1) + DLC.CLASS.is_l_noreward_under_tube_lick(LICKS_ALL_DATA.validity == 1));
is_bout_start = ismember(DLC.IND.ind_lick_onset(LICKS_ALL_DATA.validity == 1) ,DLC.IND.ind_lick_onset_str_bout);
is_bout_end = ismember(DLC.IND.ind_lick_onset(LICKS_ALL_DATA.validity == 1) ,DLC.IND.ind_lick_onset_end_bout);
is_harvest_start = ismember(DLC.IND.ind_lick_onset(LICKS_ALL_DATA.validity == 1) ,DLC.IND.ind_lick_onset_str_harvest);
is_harvest_end = ismember(DLC.IND.ind_lick_onset(LICKS_ALL_DATA.validity == 1) ,DLC.IND.ind_lick_onset_end_harvest);

tag = zeros(1,length(DLC.IND.ind_lick_onset(LICKS_ALL_DATA.validity == 1)));
tag_bout = zeros(1,length(DLC.IND.ind_lick_onset(LICKS_ALL_DATA.validity == 1)));
tag_harvest = zeros(1,length(DLC.IND.ind_lick_onset(LICKS_ALL_DATA.validity == 1)));

tag(is_groom) = 1;
tag(is_inner_tube_success) = 2;
tag(is_inner_tube_fail) = 3;
tag(is_outer_edge_success) = 4;
tag(is_outer_edge_fail) = 5;
tag(is_under_tube_success) = 6;
tag(is_under_tube_fail) = 7;

LICKS_ALL_DATA.tag = tag;
EXPERIMENT_PARAMS.lick_tag_list = lick_tag_list;

tag_bout(is_bout_start) = 1;
tag_bout(is_bout_end) = 2;

LICKS_ALL_DATA.tag_bout = tag_bout;
EXPERIMENT_PARAMS.lick_tag_bout_list = lick_tag_bout_list;

tag_harvest(is_harvest_start) = 1;
tag_harvest(is_harvest_end) = 2;

LICKS_ALL_DATA.tag_harvest = tag_harvest;
EXPERIMENT_PARAMS.lick_tag_harvest_list = lick_tag_harvest_list;
%% Build time tags
LICKS_ALL_DATA.time_onset = DLC.TIME.time_lick_onset(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.time_vmax = DLC.TIME.time_v_lick_max_abs(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.time_dmax = DLC.TIME.time_d_lick_max_abs(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.time_vmin = DLC.TIME.time_v_lick_min_abs(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.time_offset = DLC.TIME.time_lick_offset(LICKS_ALL_DATA.validity == 1)';

%% Build kinematics
% max amp, vm+/-, ang
LICKS_ALL_DATA.tongue_dm_max = DLC.KINEMATIC.d_lick_max(LICKS_ALL_DATA.validity == 1)'; % mm
LICKS_ALL_DATA.tongue_vm_max = DLC.KINEMATIC.v_lick_max(LICKS_ALL_DATA.validity == 1)'; % mm/s
LICKS_ALL_DATA.tongue_vm_min = DLC.KINEMATIC.v_lick_min(LICKS_ALL_DATA.validity == 1)'; % mm/s
LICKS_ALL_DATA.tongue_ang_max = DLC.KINEMATIC.angle_lick_max(LICKS_ALL_DATA.validity == 1)';  % deg w.r.t face normal vector

% displacement, velocity, and angle of tongue during lick onset to
% offset
LICKS_ALL_DATA.tongue_dm = DLC.KINEMATIC.d_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_vm = DLC.KINEMATIC.v_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_ang = DLC.KINEMATIC.angle_lick(LICKS_ALL_DATA.validity == 1,:)';

% durations: licks, bout, harvest
LICKS_ALL_DATA.duration_lick = DLC.TIME.time_lick_duration(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.duration_bout = DLC.TIME.time_bout_duration';
LICKS_ALL_DATA.duration_harvest = DLC.TIME.time_harvest_duration';

% displacement, velocity, and angle stream data
LICKS_ALL_DATA.tongue_dm_stream = DLC.KINEMATIC.d_tip;
LICKS_ALL_DATA.tongue_vm_stream = DLC.KINEMATIC.v_tip;
LICKS_ALL_DATA.tongue_ang_stream = DLC.KINEMATIC.angle_midtip;

% time stream data
LICKS_ALL_DATA.time_1K_stream = DLC.TIME.time_1K';


%% Build positions (x, y)
% position: tongue
LICKS_ALL_DATA.tongue_tip_px = DLC.KINEMATIC.tip_tongue_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_tip_py = DLC.KINEMATIC.tip_tongue_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_r_px = DLC.KINEMATIC.r_tongue_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_r_py = DLC.KINEMATIC.r_tongue_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_l_px = DLC.KINEMATIC.l_tongue_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_l_py = DLC.KINEMATIC.l_tongue_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_mid_px = DLC.KINEMATIC.mid_tongue_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_mid_py = DLC.KINEMATIC.mid_tongue_y_lick(LICKS_ALL_DATA.validity == 1,:)';
% position at specific kinematic event (onset, vmax, dmax, vmin, offset)
LICKS_ALL_DATA.tongue_tip_px_onset = DLC.KINEMATIC.tip_tongue_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_py_onset = DLC.KINEMATIC.tip_tongue_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_px_onset = DLC.KINEMATIC.r_tongue_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_py_onset = DLC.KINEMATIC.r_tongue_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_px_onset = DLC.KINEMATIC.l_tongue_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_py_onset = DLC.KINEMATIC.l_tongue_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_px_onset = DLC.KINEMATIC.mid_tongue_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_py_onset = DLC.KINEMATIC.mid_tongue_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_px_vmax = DLC.KINEMATIC.tip_tongue_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_py_vmax = DLC.KINEMATIC.tip_tongue_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_px_vmax = DLC.KINEMATIC.r_tongue_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_py_vmax = DLC.KINEMATIC.r_tongue_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_px_vmax = DLC.KINEMATIC.l_tongue_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_py_vmax = DLC.KINEMATIC.l_tongue_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_px_vmax = DLC.KINEMATIC.mid_tongue_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_py_vmax = DLC.KINEMATIC.mid_tongue_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_px_dmax = DLC.KINEMATIC.tip_tongue_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_py_dmax = DLC.KINEMATIC.tip_tongue_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_px_dmax = DLC.KINEMATIC.r_tongue_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_py_dmax = DLC.KINEMATIC.r_tongue_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_px_dmax = DLC.KINEMATIC.l_tongue_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_py_dmax = DLC.KINEMATIC.l_tongue_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_px_dmax = DLC.KINEMATIC.mid_tongue_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_py_dmax = DLC.KINEMATIC.mid_tongue_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_px_vmin = DLC.KINEMATIC.tip_tongue_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_py_vmin = DLC.KINEMATIC.tip_tongue_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_px_vmin = DLC.KINEMATIC.r_tongue_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_py_vmin = DLC.KINEMATIC.r_tongue_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_px_vmin = DLC.KINEMATIC.l_tongue_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_py_vmin = DLC.KINEMATIC.l_tongue_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_px_vmin = DLC.KINEMATIC.mid_tongue_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_py_vmin = DLC.KINEMATIC.mid_tongue_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_px_offset = DLC.KINEMATIC.tip_tongue_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_py_offset = DLC.KINEMATIC.tip_tongue_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_px_offset = DLC.KINEMATIC.r_tongue_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_py_offset = DLC.KINEMATIC.r_tongue_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_px_offset = DLC.KINEMATIC.l_tongue_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_py_offset = DLC.KINEMATIC.l_tongue_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_px_offset = DLC.KINEMATIC.mid_tongue_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_py_offset = DLC.KINEMATIC.mid_tongue_y_offset(LICKS_ALL_DATA.validity == 1);

% position: nose
LICKS_ALL_DATA.nose_r_px = DLC.KINEMATIC.r_nose_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.nose_r_py = DLC.KINEMATIC.r_nose_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.nose_l_px = DLC.KINEMATIC.l_nose_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.nose_l_py = DLC.KINEMATIC.l_nose_y_lick(LICKS_ALL_DATA.validity == 1,:)';
% position at specific kinematic event (onset, vmax, dmax, vmin, offset)
LICKS_ALL_DATA.nose_r_px_onset = DLC.KINEMATIC.r_nose_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_py_onset = DLC.KINEMATIC.r_nose_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_px_onset = DLC.KINEMATIC.l_nose_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_py_onset = DLC.KINEMATIC.l_nose_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_px_vmax = DLC.KINEMATIC.r_nose_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_py_vmax = DLC.KINEMATIC.r_nose_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_px_vmax = DLC.KINEMATIC.l_nose_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_py_vmax = DLC.KINEMATIC.l_nose_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_px_dmax = DLC.KINEMATIC.r_nose_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_py_dmax = DLC.KINEMATIC.r_nose_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_px_dmax = DLC.KINEMATIC.l_nose_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_py_dmax = DLC.KINEMATIC.l_nose_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_px_vmin = DLC.KINEMATIC.r_nose_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_py_vmin = DLC.KINEMATIC.r_nose_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_px_vmin = DLC.KINEMATIC.l_nose_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_py_vmin = DLC.KINEMATIC.l_nose_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_px_offset = DLC.KINEMATIC.r_nose_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_py_offset = DLC.KINEMATIC.r_nose_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_px_offset = DLC.KINEMATIC.l_nose_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_py_offset = DLC.KINEMATIC.l_nose_y_offset(LICKS_ALL_DATA.validity == 1);

% position: reward
LICKS_ALL_DATA.rew_r_px = DLC.KINEMATIC.r_food_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.rew_r_py = DLC.KINEMATIC.r_food_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.rew_l_px = DLC.KINEMATIC.l_food_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.rew_l_py = DLC.KINEMATIC.l_food_y_lick(LICKS_ALL_DATA.validity == 1,:)';
% position at specific kinematic event (onset, vmax, dmax, vmin, offset)
LICKS_ALL_DATA.rew_r_px_onset = DLC.KINEMATIC.r_food_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_py_onset = DLC.KINEMATIC.r_food_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_px_onset = DLC.KINEMATIC.l_food_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_py_onset = DLC.KINEMATIC.l_food_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_px_vmax = DLC.KINEMATIC.r_food_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_py_vmax = DLC.KINEMATIC.r_food_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_px_vmax = DLC.KINEMATIC.l_food_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_py_vmax = DLC.KINEMATIC.l_food_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_px_dmax = DLC.KINEMATIC.r_food_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_py_dmax = DLC.KINEMATIC.r_food_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_px_dmax = DLC.KINEMATIC.l_food_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_py_dmax = DLC.KINEMATIC.l_food_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_px_vmin = DLC.KINEMATIC.r_food_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_py_vmin = DLC.KINEMATIC.r_food_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_px_vmin = DLC.KINEMATIC.l_food_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_py_vmin = DLC.KINEMATIC.l_food_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_px_offset = DLC.KINEMATIC.r_food_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_py_offset = DLC.KINEMATIC.r_food_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_px_offset = DLC.KINEMATIC.l_food_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_py_offset = DLC.KINEMATIC.l_food_y_offset(LICKS_ALL_DATA.validity == 1);

% position: tubes
LICKS_ALL_DATA.rtube_r_px = DLC.KINEMATIC.r_tube_r_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.rtube_r_py = DLC.KINEMATIC.r_tube_r_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.rtube_l_px = DLC.KINEMATIC.r_tube_l_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.rtube_l_py = DLC.KINEMATIC.r_tube_l_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.ltube_r_px = DLC.KINEMATIC.l_tube_r_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.ltube_r_py = DLC.KINEMATIC.l_tube_r_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.ltube_l_px = DLC.KINEMATIC.l_tube_l_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.ltube_l_py = DLC.KINEMATIC.l_tube_l_y_lick(LICKS_ALL_DATA.validity == 1,:)';
% position at specific kinematic event (onset, vmax, dmax, vmin, offset)
LICKS_ALL_DATA.rtube_r_px_onset = DLC.KINEMATIC.r_tube_r_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_py_onset = DLC.KINEMATIC.r_tube_r_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_px_onset = DLC.KINEMATIC.r_tube_l_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_py_onset = DLC.KINEMATIC.r_tube_l_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_px_vmax = DLC.KINEMATIC.r_tube_r_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_py_vmax = DLC.KINEMATIC.r_tube_r_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_px_vmax = DLC.KINEMATIC.r_tube_l_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_py_vmax = DLC.KINEMATIC.r_tube_l_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_px_dmax = DLC.KINEMATIC.r_tube_r_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_py_dmax = DLC.KINEMATIC.r_tube_r_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_px_dmax = DLC.KINEMATIC.r_tube_l_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_py_dmax = DLC.KINEMATIC.r_tube_l_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_px_vmin = DLC.KINEMATIC.r_tube_r_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_py_vmin = DLC.KINEMATIC.r_tube_r_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_px_vmin = DLC.KINEMATIC.r_tube_l_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_py_vmin = DLC.KINEMATIC.r_tube_l_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_px_offset = DLC.KINEMATIC.r_tube_r_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_py_offset = DLC.KINEMATIC.r_tube_r_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_px_offset = DLC.KINEMATIC.r_tube_l_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_py_offset = DLC.KINEMATIC.r_tube_l_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_px_onset = DLC.KINEMATIC.l_tube_r_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_py_onset = DLC.KINEMATIC.l_tube_r_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_px_onset = DLC.KINEMATIC.l_tube_l_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_py_onset = DLC.KINEMATIC.l_tube_l_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_px_vmax = DLC.KINEMATIC.l_tube_r_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_py_vmax = DLC.KINEMATIC.l_tube_r_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_px_vmax = DLC.KINEMATIC.l_tube_l_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_py_vmax = DLC.KINEMATIC.l_tube_l_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_px_dmax = DLC.KINEMATIC.l_tube_r_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_py_dmax = DLC.KINEMATIC.l_tube_r_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_px_dmax = DLC.KINEMATIC.l_tube_l_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_py_dmax = DLC.KINEMATIC.l_tube_l_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_px_vmin = DLC.KINEMATIC.l_tube_r_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_py_vmin = DLC.KINEMATIC.l_tube_r_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_px_vmin = DLC.KINEMATIC.l_tube_l_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_py_vmin = DLC.KINEMATIC.l_tube_l_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_px_offset = DLC.KINEMATIC.l_tube_r_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_py_offset = DLC.KINEMATIC.l_tube_r_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_px_offset = DLC.KINEMATIC.l_tube_l_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_py_offset = DLC.KINEMATIC.l_tube_l_y_offset(LICKS_ALL_DATA.validity == 1);

%% Build reward-tube capacity
% rew capacity
LICKS_ALL_DATA.rew_capacity_r_lick_onset = DLC.FOOD.r_tube_food_lick_onset(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.rew_capacity_r_lick_offset = DLC.FOOD.r_tube_food_lick_offset(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.rew_capacity_l_lick_onset = DLC.FOOD.l_tube_food_lick_onset(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.rew_capacity_l_lick_offset = DLC.FOOD.l_tube_food_lick_offset(LICKS_ALL_DATA.validity == 1)';

LICKS_ALL_DATA.rew_capacity_r_bout_start = DLC.FOOD.r_tube_food_bout_start';
LICKS_ALL_DATA.rew_capacity_r_bout_end = DLC.FOOD.r_tube_food_bout_end';
LICKS_ALL_DATA.rew_capacity_l_bout_start = DLC.FOOD.l_tube_food_bout_start';
LICKS_ALL_DATA.rew_capacity_l_bout_end = DLC.FOOD.l_tube_food_bout_end';

%% QA
LICKS_ALL_DATA.qa = params.lick.qa;
end
