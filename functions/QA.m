%% MAIN Function: Quality assurance
function [DLC,EXPERIMENT_PARAMS] = QA(DLC,EXPERIMENT_PARAMS)
fprintf(['Performing QA: ' EXPERIMENT_PARAMS.file_name])
% load video
%vid_file = dir(strcat(DLC.FILE.path_to_analyzed, '*.mp4'));
%vid_obj = VideoReader(strcat(DLC.FILE.path_to_analyzed, vid_file.name));
vid_obj = VideoReader(strcat(DLC.FILE.path_to_analyzed, EXPERIMENT_PARAMS.file_name(1:17), '.mp4'));
[DLC,EXPERIMENT_PARAMS] = qa_stationary(DLC,EXPERIMENT_PARAMS, vid_obj);

check_lighting_condition = 0;
[DLC,EXPERIMENT_PARAMS] = qa_food_and_tongue(DLC,EXPERIMENT_PARAMS, vid_obj, check_lighting_condition);

save_qa(DLC,EXPERIMENT_PARAMS);
fprintf(' --> Completed. \n')

end


%% Function: QA on stationary markers
function [DLC, EXPERIMENT_PARAMS] = qa_stationary(DLC,EXPERIMENT_PARAMS,vid_obj)
%FPS = EXPERIMENT_PARAMS.FPS;
data_ = table2array(DLC.data);

% video
frame = read(vid_obj, 1);
vid_height = vid_obj.height;
vid_width = vid_obj.width;

% nose markers
fprintf("\nCorrecting nose markers...");
r_nose_x = data_(:,9);
r_nose_y = data_(:,10);
l_nose_x = data_(:,11);
l_nose_y = data_(:,12);
% Eliminate erroneous tracking from nose markers
[r_nose_x_qa, r_nose_y_qa, ~] = find_erroneous_markers(r_nose_x, r_nose_y, [0 0 vid_width vid_height], "r_nose", frame);
[l_nose_x_qa, l_nose_y_qa, ~] = find_erroneous_markers(l_nose_x, l_nose_y, [0 0 vid_width vid_height], "l_nose", frame);
r_nose_x_qa = fillmissing(r_nose_x_qa, 'linear');
r_nose_y_qa = fillmissing(r_nose_y_qa, 'linear');
l_nose_x_qa = fillmissing(l_nose_x_qa, 'linear');
l_nose_y_qa = fillmissing(l_nose_y_qa, 'linear');
DLC.data.r_nose_x = r_nose_x_qa;
DLC.data.r_nose_y = r_nose_y_qa;
DLC.data.l_nose_x = l_nose_x_qa;
DLC.data.l_nose_y = l_nose_y_qa;

% tube markers
fprintf("\nCorrecting tube markers...");
r_tube_r_x = data_(:,17);
r_tube_r_y = data_(:,18);
r_tube_l_x = data_(:,19);
r_tube_l_y = data_(:,20);
l_tube_r_x = data_(:,21);
l_tube_r_y = data_(:,22);
l_tube_l_x = data_(:,23);
l_tube_l_y = data_(:,24);
% check if the tube markers are at a reasonable position
l_l_length = median(l_tube_l_y);
l_r_length = median(l_tube_r_y);
if  abs(l_l_length - l_r_length) >= (l_l_length/4) || abs(l_l_length - l_r_length) >= (l_r_length/4)
    % manual correction
    [l_tube_l_x_qa, l_tube_l_y_qa, ~] = manual_correction(l_tube_l_x, l_tube_l_y, "l\_tube\_l", frame);
    [l_tube_r_x_qa, l_tube_r_y_qa, ~] = manual_correction(l_tube_r_x, l_tube_r_y, "l\_tube\_r", frame);
else
    [l_tube_r_x_qa, l_tube_r_y_qa, ~] = find_erroneous_markers(l_tube_r_x, l_tube_r_y, [0 0 vid_width vid_height/2], "l_tube_r", frame);
    [l_tube_l_x_qa, l_tube_l_y_qa, ~] = find_erroneous_markers(l_tube_l_x, l_tube_l_y, [0 0 vid_width vid_height/2], "l_tube_l", frame);
end
    
r_l_length = vid_height - median(r_tube_l_y);
r_r_length = vid_height - median(r_tube_r_y);
if abs(r_l_length - r_r_length) >= (r_l_length/4) || abs(r_l_length - r_r_length) >= (r_r_length/4)
    % manual correction
    [r_tube_l_x_qa, r_tube_l_y_qa, ~] = manual_correction(r_tube_l_x, r_tube_l_y, "r_tube_l", frame);
    [r_tube_r_x_qa, r_tube_r_y_qa, ~] = manual_correction(r_tube_r_x, r_tube_r_y, "r_tube_r", frame);
else
    [r_tube_r_x_qa, r_tube_r_y_qa, ~] = find_erroneous_markers(r_tube_r_x, r_tube_r_y, [0 vid_height/2 vid_width vid_height], "r_tube_r", frame);
    [r_tube_l_x_qa, r_tube_l_y_qa, ~] = find_erroneous_markers(r_tube_l_x, r_tube_l_y, [0 vid_height/2 vid_width vid_height], "r_tube_l", frame);
end
% Perform linear interpolation
l_tube_r_x_qa = fillmissing(l_tube_r_x_qa, 'linear');
l_tube_r_y_qa = fillmissing(l_tube_r_y_qa, 'linear');
l_tube_l_x_qa = fillmissing(l_tube_l_x_qa, 'linear');
l_tube_l_y_qa = fillmissing(l_tube_l_y_qa, 'linear');
r_tube_r_x_qa = fillmissing(r_tube_r_x_qa, 'linear');
r_tube_r_y_qa = fillmissing(r_tube_r_y_qa, 'linear');
r_tube_l_x_qa = fillmissing(r_tube_l_x_qa, 'linear');
r_tube_l_y_qa = fillmissing(r_tube_l_y_qa, 'linear');

if EXPERIMENT_PARAMS.flag_figure_debug == 1
    figure
    hold on;
    set(gca,'Ydir','reverse')
    title("Stationary markers (Blue: DLC | Green: QA)")
    plot(l_tube_r_x, l_tube_r_y, 'b*')
    plot(l_tube_l_x, l_tube_l_y, 'b*')
    plot(r_tube_r_x, r_tube_r_y, 'b*')
    plot(r_tube_l_x, r_tube_l_y, 'b*')
    plot(r_nose_x, r_nose_y, 'b*')
    plot(l_nose_x, l_nose_y, 'b*')
    plot(l_tube_r_x_qa, l_tube_r_y_qa, 'g*')
    plot(l_tube_l_x_qa, l_tube_l_y_qa, 'g*')
    plot(r_tube_r_x_qa, r_tube_r_y_qa, 'g*')
    plot(r_tube_l_x_qa, r_tube_l_y_qa, 'g*')
    plot(r_nose_x_qa, r_nose_y_qa, 'g*')
    plot(l_nose_x_qa, l_nose_y_qa, 'g*')
end

DLC.data.l_tube_r_x = l_tube_r_x_qa;
DLC.data.l_tube_r_y = l_tube_r_y_qa;
DLC.data.l_tube_l_x = l_tube_l_x_qa;
DLC.data.l_tube_l_y = l_tube_l_y_qa;
DLC.data.r_tube_r_x = r_tube_r_x_qa;
DLC.data.r_tube_r_y = r_tube_r_y_qa;
DLC.data.r_tube_l_x = r_tube_l_x_qa;
DLC.data.r_tube_l_y = r_tube_l_y_qa;
fprintf(' --> Completed. \n')
    
% helper function to detect shift, outlier, and outside expected boundary
% rect = [x_lower, y_lower, x_upper, y_upper]      
function [x, y, ind_DLC_ERR] = find_erroneous_markers(x, y, rect, type, frame)
    mov_x = movmean(x, 1000);
    mov_y = movmean(y, 1000);
    
    if ~isempty(regexp(type, regexptranslate('wildcard', '*nose'))) || (max(mov_x) - min(mov_x) < 5 && max(mov_y) - min(mov_y) < 5)
        ind_DLC_ERR = [];
        % isoutlier(median) method
        x_err = isoutlier(x, "median");
        y_err = isoutlier(y, "median");
        ind_DLC_ERR = unique([find(x_err);find(y_err)]);

        % boundary method
        for i = 1:length(x)
            if x(i) < rect(1) || x(i) > rect(3) || y(i) < rect(2) || y(i) > rect(4)
                ind_DLC_ERR = [ind_DLC_ERR; i];
            end
        end
        ind_DLC_ERR = unique(ind_DLC_ERR);
        x(ind_DLC_ERR) = NaN;
        y(ind_DLC_ERR) = NaN;
    else
        % manual correction
        [x, y, ind_DLC_ERR] = manual_correction(x, y, type, frame);
    end

    % prepare data for interpolation
    if isnan(x(length(x)))
        x(length(x)) = x(find(~isnan(x), 1, "last"));
        y(length(y)) = y(find(~isnan(y), 1, "last"));
    end
    if isnan(x(1))
        x(1) = x(find(~isnan(x), 1, "first"));
        y(1) = y(find(~isnan(y), 1, "first"));
    end

end
function [x, y, ind_DLC_ERR] = manual_correction(x, y, type, frame)
    fig = figure;
    imshow(frame)
    title(type, 'Interpreter', 'none')
    [x_ref, y_ref] = getpts;
    x_ref = x_ref(1);
    y_ref = y_ref(1);
    ind_DLC_ERR = [];
    ind_DLC_ERR = [ind_DLC_ERR; find(abs(x - x_ref)>5)];
    ind_DLC_ERR = [ind_DLC_ERR; find(abs(y - y_ref)>5)];
    ind_DLC_ERR = unique(ind_DLC_ERR);
    if length(ind_DLC_ERR) == length(x)
        x = x_ref*ones(size(x));
        y = y_ref*ones(size(y));
    else
        x(ind_DLC_ERR) = NaN;
        y(ind_DLC_ERR) = NaN;
    end
    close(fig)
end
function [x, y] = kalman_filtering(x, y, FPS)
    %% define main variables for kalman filter
    dt = 1 / FPS;
    u = 0; % acceleration magnitude
    Q = [x(1); y(1); 0; 0]; % initialized state--it has four components: [positionX; positionY; velocityX; velocityY] of the hexbug
    Q_estimate = Q;  % estimate of initial location (what we are updating)
    accel_noise_mag = 2; % process noise: the variability in how fast the tracking is speeding up (stdv of acceleration: meters/sec^2)
    n_x = 10;  % measurement noise in the x direction
    n_y = 10;  % measurement noise in the y direction
    Ez = [n_x 0; 0 n_y];
    Ex = [dt^4/4 0 dt^3/2 0; ...
        0 dt^4/4 0 dt^3/2; ...
        dt^3/2 0 dt^2 0; ...
        0 dt^3/2 0 dt^2].*accel_noise_mag^2; % Ex convert the process noise (stdv) into covariance matrix
    P = Ex; % estimate of initial position variance (covariance matrix)
    %% Define update equations in 2D (Coefficent matrices): A physics based model for where we expect the marker to be [state transition (state + velocity)] + [input control (acceleration)]
    A = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1]; % state update matrix
    B = [(dt^2/2); (dt^2/2); dt; dt]; % input control matrix
    C = [1 0 0 0; 0 1 0 0];  % measurement function C (applied to the state estimate Q to get the expected next new measurement)
    %% Initize result variables
    Q_loc_meas = []; % the path extracted by the tracking algo
    %% Initize estimation variables
    p_estimate = []; % position estimate
    v_estimate = []; % velocity estimate
    predic_state = [];
    predic_var = [];
    for t = 1:size(x)
        % load the given tracking
        Q_loc_meas(:,t) = [x(t); y(t)];
        %% do the kalman filter
        % Predict next state of the with the last state and predicted motion.
        Q_estimate = A * Q_estimate + B * u;
        predic_state = [predic_state; Q_estimate(1)] ;
        % Predict next covariance
        P = A * P * A' + Ex;
        predic_var = [predic_var; P] ;
        % Predicted measurement covariance
        % Kalman Gain
        K = P*C'*inv(C*P*C'+Ez);
        % Update the state estimate.
        if ~isnan(Q_loc_meas(:,t))
            Q_estimate = Q_estimate + K * (Q_loc_meas(:,t) - C * Q_estimate);
        end
        % update covariance estimation.
        P = (eye(4)-K*C)*P;
        %% Store data
        p_estimate = [p_estimate; [Q_estimate(1) Q_estimate(2)]];
        v_estimate = [v_estimate; [Q_estimate(3) Q_estimate(4)]];
    end
    x = p_estimate(:, 1);
    y = p_estimate(:, 2);
end
end

%% Function: QA on food and tongue
function [DLC, EXPERIMENT_PARAMS] = qa_food_and_tongue(DLC,EXPERIMENT_PARAMS,vid_obj,light_condition)
data_ = table2array(DLC.data);
path_to_raw = [EXPERIMENT_PARAMS.mat_PathName '..' filesep '..' filesep '..' filesep 'raw_data'];
% video
vid_height = vid_obj.height;
vid_width = vid_obj.width;

r_nose_x = DLC.data.r_nose_x;
r_nose_y = DLC.data.r_nose_y;
l_nose_x = DLC.data.l_nose_x;
l_nose_y = DLC.data.l_nose_y;
if light_condition == 1
    fprintf("Checking light condition in video ...");

    % 1) Lighting condition (for both tongue and food markers)
    avg_pi = [];
    mean_nose_x = round((mean(r_nose_x) + mean(l_nose_x))/2);
    mean_nose_y = round((mean(r_nose_y) + mean(l_nose_y))/2);

    for i = 1:length(r_nose_x)
        frame = read(vid_obj, i);
        avg_pi = [avg_pi; frame(mean_nose_y, mean_nose_x, :)];
    end
    % Find dark frames
    avg_pi = mean(avg_pi, 3);
    TF = isoutlier(avg_pi,'movmedian',length(avg_pi)/2);
    dark = [];
    for i = 1:length(TF)
        if TF(i)
            if avg_pi(i) < mean(avg_pi)
                dark = [dark; avg_pi(i)];
            else
                TF(i) = 0;
            end
        end
    end

    y0 = (mean(r_nose_y) + mean(l_nose_y ))/2;
    y_lower_lim = 10;
    y_upper_lim = vid_height - 10;
    x_upper_lim = vid_width;
    fprintf(' --> Completed. \n')
else
    TF = [];
    mean_nose_x = round((mean(r_nose_x) + mean(l_nose_x))/2);
    y0 = (mean(r_nose_y) + mean(l_nose_y ))/2;
    y_lower_lim = 10;
    y_upper_lim = vid_height - 10;
    x_upper_lim = vid_width;
end

% QA for food markers
fprintf("Correcting food markers...");
r_food_x = data_(:,13);
r_food_y = data_(:,14);
l_food_x = data_(:,15);
l_food_y = data_(:,16);
r_tube_r_x = DLC.data.r_tube_r_x;
r_tube_r_y = DLC.data.r_tube_r_y;
r_tube_l_x = DLC.data.r_tube_l_x;
r_tube_l_y = DLC.data.r_tube_l_y;
l_tube_r_x = DLC.data.l_tube_r_x;
l_tube_r_y = DLC.data.l_tube_r_y;
l_tube_l_x = DLC.data.l_tube_l_x;
l_tube_l_y = DLC.data.l_tube_l_y;

% 2) Boundary
TF_left = TF;
TF_right = TF;
mean_r_tube_l_x = mean(r_tube_l_x);
mean_r_tube_r_x = mean(r_tube_r_x);
mean_l_tube_l_x = mean(l_tube_l_x);
mean_l_tube_r_x = mean(l_tube_r_x);
for i = 1:length(l_food_y)
    if l_food_y(i) > y0 || l_food_x(i) < mean_l_tube_l_x ||  l_food_x(i) > mean_l_tube_r_x || l_food_y(i) <= y_lower_lim
        TF_left(i) = 1;
    end
    if r_food_y(i) < y0 || r_food_x(i) < mean_r_tube_l_x ||  r_food_x(i) > mean_r_tube_r_x || r_food_y(i) >= y_upper_lim
        TF_right(i) = 1;
    end
end
l_food_x_qa = l_food_x;
l_food_y_qa = l_food_y;
r_food_x_qa = r_food_x;
r_food_y_qa = r_food_y;
% l_food
l_food_incorrect_inds = find(TF_left);
l_food_x_qa(l_food_incorrect_inds) = NaN;
l_food_y_qa(l_food_incorrect_inds) = NaN;
if isnan(l_food_x_qa(1))
    l_food_x_qa(1) = l_food_x_qa(find(~TF_left, 1));
    l_food_y_qa(1) = l_food_y_qa(find(~TF_left, 1));
end
if isnan(l_food_x_qa(length(l_food_x_qa)))
    l_food_x_qa(length(l_food_x_qa)) = l_food_x_qa(find(~TF_left, 1, 'last'));
    l_food_y_qa(length(l_food_x_qa)) = l_food_y_qa(find(~TF_left, 1, 'last'));
end
% r_food
r_food_incorrect_inds = find(TF_right);
r_food_x_qa(r_food_incorrect_inds) = NaN;
r_food_y_qa(r_food_incorrect_inds) = NaN;
if isnan(r_food_x_qa(1))
    r_food_x_qa(1) = r_food_x_qa(find(~TF_right, 1));
    r_food_y_qa(1) = r_food_y_qa(find(~TF_right, 1));
end
if isnan(r_food_x_qa(length(r_food_x_qa)))
    r_food_x_qa(length(r_food_x_qa)) = r_food_x_qa(find(~TF_right, 1, 'last'));
    r_food_y_qa(length(r_food_x_qa)) = r_food_y_qa(find(~TF_right, 1, 'last'));
end
% Fixing all incorrection food markers with makima interpolation
l_food_x_qa = fillmissing(l_food_x_qa, 'makima');
l_food_y_qa = fillmissing(l_food_y_qa, 'makima');
r_food_x_qa = fillmissing(r_food_x_qa, 'makima');
r_food_y_qa = fillmissing(r_food_y_qa, 'makima');

% 3) Detect jumping markers
% l_food
max_diff = mean([median(l_tube_r_y), median(l_tube_l_y)]) / 4;
dy = diff(l_food_y_qa);
jump_inds = zeros(size(l_food_y_qa));
rise = NaN;
for i = 1:length(dy)
    if isnan(rise)
        % find jumping start
        if abs(dy(i)) > max_diff
            if dy(i) > 0
                rise = true;
            else
                rise = false;
            end
                jump_inds(i+1) = 1;
        end
    else
        % find jumping end
        if abs(dy(i)) > max_diff
            if (rise && dy(i) < 0) || (~rise && dy(i) > 0)
                jump_inds(i+1) = -1;
                rise = NaN;
            end
        end
    end
end
jumping = false;
for i = 1:length(jump_inds)
    if jump_inds(i) == 1 && ~isempty(find(jump_inds(i:i+300)==-1, 1))
            jumping = true;
    end
    if jump_inds(i) == -1
        jumping = false;
    end
    if jumping
        l_food_y_qa(i) = NaN;
        l_food_x_qa(i) = NaN;
        l_food_incorrect_inds = [l_food_incorrect_inds; i];
    end
end

% r_food
max_diff = mean([median(r_tube_r_y), median(r_tube_l_y)]) / 4;
dy = diff(r_food_y_qa);
jump_inds = zeros(size(r_food_y_qa));
rise = NaN;
for i = 1:length(dy)
    if isnan(rise)
        % find jumping start
        if abs(dy(i)) > max_diff
            if dy(i) > 0
                rise = true;
            else
                rise = false;
            end
                jump_inds(i+1) = 1;
        end
    else
        % find jumping end
        if abs(dy(i)) > max_diff
            if (rise && dy(i) < 0) || (~rise && dy(i) > 0)
                jump_inds(i+1) = -1;
                rise = NaN;
            end
        end
    end
end
jumping = false;
for i = 1:length(jump_inds)
    if jump_inds(i) == 1 && ~isempty(find(jump_inds(i:i+300)==-1, 1))
            jumping = true;
    end
    if jump_inds(i) == -1
        jumping = false;
    end
    if jumping
        r_food_y_qa(i) = NaN;
        r_food_x_qa(i) = NaN;
        r_food_incorrect_inds = [r_food_incorrect_inds; i];
    end
end

l_food_x_qa = fillmissing(l_food_x_qa, 'makima');
l_food_y_qa = fillmissing(l_food_y_qa, 'makima');
r_food_x_qa = fillmissing(r_food_x_qa, 'makima');
r_food_y_qa = fillmissing(r_food_y_qa, 'makima');

l_food_incorrect_inds = unique(sort(l_food_incorrect_inds));
r_food_incorrect_inds = unique(sort(r_food_incorrect_inds));
if EXPERIMENT_PARAMS.flag_figure_debug == 1
    figure
    hold on
    set(gca,'Ydir','reverse')
    title("Food markers (Blue/Red: QA-left/right | *: DLC)")
    plot(l_food_incorrect_inds, l_food_y(l_food_incorrect_inds), 'b*')
    plot(r_food_incorrect_inds, r_food_y(r_food_incorrect_inds), 'r*')
    plot(l_food_y_qa, 'b-')
    plot(r_food_y_qa, 'r-')
end
% Save variables
DLC.data.l_food_x = l_food_x_qa;
DLC.data.l_food_y = l_food_y_qa;
DLC.data.r_food_x = r_food_x_qa;
DLC.data.r_food_y = r_food_y_qa;
fprintf(' --> Completed. \n')

% QA for tongue markers
fprintf("Correcting tongue markers...");
tip_tongue_x = data_(:,1);
tip_tongue_y = data_(:,2);
r_tongue_x = data_(:,3);
r_tongue_y = data_(:,4);
l_tongue_x = data_(:,5);
l_tongue_y = data_(:,6);
mid_tongue_x = data_(:,7);
mid_tongue_y = data_(:,8);

TF_tip = TF;
TF_l = TF;
TF_r = TF;
TF_mid = TF;
% 1) Boundary
for i = 2:length(tip_tongue_y)
    if tip_tongue_y(i) <= y_lower_lim || tip_tongue_y(i) >= y_upper_lim || abs(angle(i) - angle(i-1)) > 80
        TF_tip(i) = 1;
    end
    if l_tongue_y(i) <= y_lower_lim || l_tongue_y(i) >= y_upper_lim
        TF_l(i) = 1;
    end
    if r_tongue_y(i) <= y_lower_lim || r_tongue_y(i) >= y_upper_lim
        TF_r(i) = 1;
    end
    if mid_tongue_y(i) <= y_lower_lim || mid_tongue_x(i) >= x_upper_lim || mid_tongue_y(i) >= y_upper_lim
        TF_mid(i) = 1;
    end
end

% 2) Abrupt change in distance between the 4 markers
% 3) Abnormal shape (sum of angles not equal to 360)
tl_lengths = zeros(size(tip_tongue_x));
lm_lengths = zeros(size(tip_tongue_x));
mr_lengths = zeros(size(tip_tongue_x));
rt_lengths = zeros(size(tip_tongue_x));
tlm_angles = zeros(size(tip_tongue_x));
lmr_angles = zeros(size(tip_tongue_x));
mrt_angles = zeros(size(tip_tongue_x));
rtl_angles = zeros(size(tip_tongue_x));
for i = 1:length(tip_tongue_x)
    tip = [tip_tongue_x(i), tip_tongue_y(i)];
    mid = [mid_tongue_x(i), mid_tongue_y(i)];
    l = [l_tongue_x(i), l_tongue_y(i)];
    r = [r_tongue_x(i), r_tongue_y(i)];
    
    % populate length
    tl_lengths(i) = norm(tip-l);
    lm_lengths(i) = norm(l-mid);
    mr_lengths(i) = norm(mid-r);
    rt_lengths(i) = norm(r-tip);
    % populate angles
    tlm_angles(i) = find_angle(tip, l, mid);
    lmr_angles(i) = find_angle(l, mid, r);
    mrt_angles(i) = find_angle(mid, r, tip);
    rtl_angles(i) = find_angle(r, tip, l);
end

inds = [];
[~, ind] = findpeaks(diff(tl_lengths), "MinPeakProminence", 30);
inds = [inds; ind];
[~, ind] = findpeaks(diff(lm_lengths), "MinPeakProminence", 30);
inds = [inds; ind];
[~, ind] = findpeaks(diff(mr_lengths), "MinPeakProminence", 30);
inds = [inds; ind];
[~, ind] = findpeaks(diff(rt_lengths), "MinPeakProminence", 30);
inds = [inds; ind];
angles = rtl_angles + mrt_angles + lmr_angles + tlm_angles;
inds = [inds; find(angles < 359 | angles > 361)];

% Fixing all incorrect tongue markers using makima interpolation
tip_tongue_x_qa = tip_tongue_x;
tip_tongue_y_qa = tip_tongue_y;
l_tongue_x_qa = l_tongue_x;
l_tongue_y_qa = l_tongue_y;
mid_tongue_x_qa = mid_tongue_x;
mid_tongue_y_qa = mid_tongue_y;
r_tongue_x_qa = r_tongue_x;
r_tongue_y_qa = r_tongue_y;

tip_incorrect_inds = unique(sort([find(TF_tip); inds]));
l_incorrect_inds = unique(sort([find(TF_l); inds]));
r_incorrect_inds = unique(sort([find(TF_r); inds]));
mid_incorrect_inds = unique(sort([find(TF_mid); inds]));
tip_tongue_x_qa(tip_incorrect_inds) = NaN;
tip_tongue_y_qa(tip_incorrect_inds) = NaN;
l_tongue_x_qa(l_incorrect_inds) = NaN;
l_tongue_y_qa(l_incorrect_inds) = NaN;
r_tongue_x_qa(r_incorrect_inds) = NaN;
r_tongue_y_qa(r_incorrect_inds) = NaN;
mid_tongue_x_qa(mid_incorrect_inds) = NaN;
mid_tongue_y_qa(mid_incorrect_inds) = NaN;
tip_tongue_x_qa = fillmissing(tip_tongue_x_qa, 'makima');
tip_tongue_y_qa = fillmissing(tip_tongue_y_qa, 'makima');
l_tongue_x_qa = fillmissing(l_tongue_x_qa, 'makima');
l_tongue_y_qa = fillmissing(l_tongue_y_qa, 'makima');
r_tongue_x_qa = fillmissing(r_tongue_x_qa, 'makima');
r_tongue_y_qa = fillmissing(r_tongue_y_qa, 'makima');
mid_tongue_x_qa = fillmissing(mid_tongue_x_qa, 'makima');
mid_tongue_y_qa = fillmissing(mid_tongue_y_qa, 'makima');

if EXPERIMENT_PARAMS.flag_figure_debug == 1
    tip_dist = sqrt((tip_tongue_x).^2 + (tip_tongue_y).^2);
    l_dist = sqrt((l_tongue_x).^2 + (l_tongue_y).^2);
    r_dist = sqrt((r_tongue_x).^2 + (r_tongue_y).^2);
    mid_dist = sqrt((mid_tongue_x).^2 + (mid_tongue_y).^2);

    tip_dist_qa = sqrt((tip_tongue_x_qa).^2 + (tip_tongue_y_qa).^2);
    l_dist_qa = sqrt((l_tongue_x_qa).^2 + (l_tongue_y_qa).^2);
    r_dist_qa = sqrt((r_tongue_x_qa).^2 + (r_tongue_y_qa).^2);
    mid_dist_qa = sqrt((mid_tongue_x_qa).^2 + (mid_tongue_y_qa).^2);

    figure
    tip_inds = find(abs(tip_dist(tip_incorrect_inds)-tip_dist_qa(tip_incorrect_inds)) > 5);
    l_inds = find(abs(l_dist(l_incorrect_inds)-l_dist_qa(l_incorrect_inds)) > 5);
    r_inds = find(abs(r_dist(r_incorrect_inds)-r_dist_qa(r_incorrect_inds)) > 5);
    mid_inds = find(abs(mid_dist(mid_incorrect_inds)-mid_dist_qa(mid_incorrect_inds)) > 5);
    tip_incorrect_inds = tip_incorrect_inds(tip_inds);
    l_incorrect_inds = l_incorrect_inds(l_inds);
    r_incorrect_inds = r_incorrect_inds(r_inds);
    mid_incorrect_inds = mid_incorrect_inds(mid_inds);

    subplot(4, 2, 1)
    plot(tip_tongue_x_qa, '-')
    hold on
    plot(tip_incorrect_inds, tip_tongue_x(tip_incorrect_inds), '*')
    title("Tip tongue x QA (*: DLC)")

    subplot(4, 2, 2)
    plot(tip_tongue_y_qa, '-')
    hold on
    plot(tip_incorrect_inds, tip_tongue_y(tip_incorrect_inds), '*')
    title("Tip tongue y QA (*: DLC)")

    subplot(4, 2, 3)
    plot(l_tongue_x_qa, '-')
    hold on
    plot(l_incorrect_inds, l_tongue_x(l_incorrect_inds), '*')
    title("Left tongue x QA (*: DLC)")

    subplot(4, 2, 4)
    plot(l_tongue_y_qa, '-')
    hold on
    plot(l_incorrect_inds, l_tongue_y(l_incorrect_inds), '*')
    title("Left tongue y QA (*: DLC)")

    subplot(4, 2, 5)
    plot(r_tongue_x_qa, '-')
    hold on
    plot(r_incorrect_inds, l_tongue_x(r_incorrect_inds), '*')
    title("Right tongue x QA (*: DLC)")

    subplot(4, 2, 6)
    plot(r_tongue_y_qa, '-')
    hold on
    plot(r_incorrect_inds, r_tongue_y(r_incorrect_inds), '*')
    title("Right tongue y QA (*: DLC)")

    subplot(4, 2, 7)
    plot(mid_tongue_x_qa, '-')
    hold on
    plot(mid_incorrect_inds, mid_tongue_x(mid_incorrect_inds), '*')
    title("Mid tongue x QA (*: DLC)")

    subplot(4, 2, 8)
    plot(mid_tongue_y_qa, '-')
    hold on
    plot(mid_incorrect_inds, mid_tongue_y(mid_incorrect_inds), '*')
    title("Mid tongue y QA (*: DLC)")

%     inds = sort(unique([tip_incorrect_inds; l_incorrect_inds; r_incorrect_inds; mid_incorrect_inds]));
%     for i = 1:length(inds)
%         imwrite(read(vid_obj, inds(i)), strcat(num2str(inds(i)), '.png'))
%     end
end
DLC.data.tip_tongue_x = tip_tongue_x_qa;
DLC.data.tip_tongue_y = tip_tongue_y_qa;
DLC.data.r_tongue_x = r_tongue_x_qa;
DLC.data.r_tongue_y = r_tongue_y_qa;
DLC.data.l_tongue_x = l_tongue_x_qa;
DLC.data.l_tongue_y = l_tongue_y_qa;
DLC.data.mid_tongue_x = mid_tongue_x_qa;
DLC.data.mid_tongue_y = mid_tongue_y_qa;
fprintf(' --> Completed. \n')

function angle = find_angle(a, b, c)
x1 = b(1);
y1 = b(2);
x2 = a(1);
y2 = a(2);
x3 = c(1);
y3 = c(2);
angle = atan2(abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)), ...
            (x2-x1)*(x3-x1)+(y2-y1)*(y3-y1)) * 180/pi;
end
end

%% Function: Update .mat and csv after QA
function save_qa(DLC, EXPERIMENT_PARAMS)
% update csv file
csv_file = dir(strcat(DLC.FILE.path_to_analyzed, '*.csv'));
[csv_table, t, ~] = xlsread(strcat(DLC.FILE.path_to_analyzed, csv_file.name));
csv_table(:,2) = DLC.data.tip_tongue_x;
csv_table(:,3) = DLC.data.tip_tongue_y;
csv_table(:,5) = DLC.data.r_tongue_x;
csv_table(:,6) = DLC.data.r_tongue_y;
csv_table(:,8) = DLC.data.l_tongue_x;
csv_table(:,9) = DLC.data.l_tongue_y;
csv_table(:,11) = DLC.data.mid_tongue_x;
csv_table(:,12) = DLC.data.mid_tongue_y;
csv_table(:,14) = DLC.data.r_nose_x;
csv_table(:,15) = DLC.data.r_nose_y;
csv_table(:,17) = DLC.data.l_nose_x;
csv_table(:,18) = DLC.data.l_nose_y;
csv_table(:,20) = DLC.data.r_food_x;
csv_table(:,21) = DLC.data.r_food_y;
csv_table(:,23) = DLC.data.l_food_x;
csv_table(:,24) = DLC.data.l_food_y;
csv_table(:,26) = DLC.data.r_tube_r_x;
csv_table(:,27) = DLC.data.r_tube_r_y;
csv_table(:,29) = DLC.data.r_tube_l_x;
csv_table(:,30) = DLC.data.r_tube_l_y;
csv_table(:,32) = DLC.data.l_tube_r_x;
csv_table(:,33) = DLC.data.l_tube_r_y;
csv_table(:,35) = DLC.data.l_tube_l_x;
csv_table(:,36) = DLC.data.l_tube_l_y;
writecell(t, strcat(DLC.FILE.path_to_analyzed, csv_file.name));
writematrix(csv_table,strcat(DLC.FILE.path_to_analyzed, csv_file.name), 'WriteMode', 'append');

% update .mat file
data = DLC.data;
save(strcat(DLC.FILE.path_to_analyzed, EXPERIMENT_PARAMS.mat_FileName), 'data')
end
