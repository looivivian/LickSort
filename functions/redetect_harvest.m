function [DLC, EXPERIMENT_PARAMS] = redetect_harvest(DLC,EXPERIMENT_PARAMS, params, funcs, RECLASSIFIED)
fprintf('Detecting harvest ...');

% override manually assigned classes (if lick onset is the
% same)
if ~isempty(RECLASSIFIED.ind_lick_onset)
    for i = 1:length(RECLASSIFIED.ind_lick_onset)
        % find lick num from frame number
        lick_num = [find(DLC.IND.ind_lick_onset == RECLASSIFIED.ind_lick_onset(i), 1)...
        find(DLC.IND.ind_lick_onset == RECLASSIFIED.ind_lick_onset(i) + 1, 1)...
        find(DLC.IND.ind_lick_onset == RECLASSIFIED.ind_lick_onset(i) - 1, 1)];
        if ~isempty(lick_num)
            DLC = setDLCfromClass(DLC, RECLASSIFIED.class(i), lick_num);
        end
    end
end

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

fprintf(' --> Completed. \n')
end

function DLC = setDLCfromClass(DLC, class, lick_num)
    ind_lick_onset = DLC.IND.ind_lick_onset(lick_num);
    ind_lick_offset = DLC.IND.ind_lick_offset(lick_num);
    % reset DLC.CLASS for this lick
    DLC.CLASS.is_r_reward_inner_tube_lick(lick_num) = 0;
    DLC.CLASS.is_r_reward_outer_edge_lick(lick_num) = 0;
    DLC.CLASS.is_r_reward_under_tube_lick(lick_num) = 0;
    DLC.CLASS.is_r_noreward_inner_tube_lick(lick_num) = 0;
    DLC.CLASS.is_r_noreward_outer_edge_lick(lick_num) = 0;
    DLC.CLASS.is_r_noreward_under_tube_lick(lick_num) = 0;
    DLC.CLASS.is_l_reward_inner_tube_lick(lick_num) = 0;
    DLC.CLASS.is_l_reward_outer_edge_lick(lick_num) = 0;
    DLC.CLASS.is_l_reward_under_tube_lick(lick_num) = 0;
    DLC.CLASS.is_l_noreward_inner_tube_lick(lick_num) = 0;
    DLC.CLASS.is_l_noreward_outer_edge_lick(lick_num) = 0;
    DLC.CLASS.is_l_noreward_under_tube_lick(lick_num) = 0;
    DLC.CLASS.is_grooming_lick(lick_num) = 0;

    r_reward_ind = find(DLC.CLASS.ind_lick_onset_r_reward == ind_lick_onset, 1);
    l_reward_ind = find(DLC.CLASS.ind_lick_onset_l_reward == ind_lick_onset, 1);
    if ~isempty(r_reward_ind)
        DLC.CLASS.ind_lick_onset_r_reward(r_reward_ind) = [];
        DLC.CLASS.ind_lick_offset_r_reward(r_reward_ind) = [];
    elseif ~isempty(l_reward_ind)
        DLC.CLASS.ind_lick_onset_l_reward(l_reward_ind) = [];
        DLC.CLASS.ind_lick_offset_l_reward(l_reward_ind) = [];
    end
    
    if class.dir == 'g'
        DLC.CLASS.is_grooming_lick(lick_num) = 1;
    elseif class.dir == 'l'
        if class.reward
            if class.type == 'i'
                DLC.CLASS.is_l_reward_inner_tube_lick(lick_num) = 1;
            elseif class.type == 'u'
                DLC.CLASS.is_l_reward_under_tube_lick(lick_num) = 1;
            elseif class.type == 'o'
                DLC.CLASS.is_l_reward_outer_edge_lick(lick_num) = 1;
            end
            % add to l_reward_ind
            DLC.CLASS.ind_lick_onset_l_reward = sort([DLC.CLASS.ind_lick_onset_l_reward; ind_lick_onset]);
            DLC.CLASS.ind_lick_offset_l_reward = sort([DLC.CLASS.ind_lick_onset_l_reward; ind_lick_offset]);
        else
            if class.type == 'i'
                DLC.CLASS.is_l_noreward_inner_tube_lick(lick_num) = 1;
            elseif class.type == 'u'
                DLC.CLASS.is_l_noreward_under_tube_lick(lick_num) = 1;
            elseif class.type == 'o'
                DLC.CLASS.is_l_noreward_outer_edge_lick(lick_num) = 1;
            end
        end
    elseif class.dir == 'r'
        if class.reward
            if class.type == 'i'
                DLC.CLASS.is_r_reward_inner_tube_lick(lick_num) = 1;
            elseif class.type == 'u'
                DLC.CLASS.is_r_reward_under_tube_lick(lick_num) = 1;
            elseif class.type == 'o'
                DLC.CLASS.is_r_reward_outer_edge_lick(lick_num) = 1;
            end
            % add to r_reward_ind
            DLC.CLASS.ind_lick_onset_r_reward = sort([DLC.CLASS.ind_lick_onset_r_reward; ind_lick_onset]);
            DLC.CLASS.ind_lick_offset_r_reward = sort([DLC.CLASS.ind_lick_onset_r_reward; ind_lick_offset]);
        else
            if class.type == 'i'
                DLC.CLASS.is_r_noreward_inner_tube_lick(lick_num) = 1;
            elseif class.type == 'u'
                DLC.CLASS.is_r_noreward_under_tube_lick(lick_num) = 1;
            elseif class.type == 'o'
                DLC.CLASS.is_r_noreward_outer_edge_lick(lick_num) = 1;
            end
        end
    end
end
