%% Produces LICKS_ALL_DATA and EXPERIMENT_PARAMS
%  Params
%  mat_file_address:    filepath to a DLC.mat for analysis
%  flag_qa:             bool for whether to perform auto quality assurance
%  flag_figure:         bool for whether to produce figures
%  gui:                 bool for whether to use the GUI
%  params:
%  funcs:
%  
%  Return
%  LICKS_ALL_DATA:      a struct containing analyzed licking data
%  EXPERIMENT_PARAMS:   a struct containing experiment parameters
% 
function [LICKS_ALL_DATA, EXPERIMENT_PARAMS] = PGH_monkey_behavior_lick(mat_file_address, flag_qa, flag_figure, gui, params, funcs)
warning('off','all')
addpath('functions/')
%% Get file_name and file_path

% if there is no inputs, then set pathnames to pwd
if nargin < 1
    [file_name,file_path] = uigetfile(pwd, 'Select DLC mat file');
    flag_qa = 1;
    flag_figure = 1;
    gui = 1;
    params = [];
    funcs = [];
end
if nargin == 1
    [file_path,file_name,file_ext] = fileparts(mat_file_address);
    file_name = [file_name file_ext];
    flag_qa = 1;
    flag_figure = 1;
    gui = false;
    params = [];
    funcs = [];
end
if nargin >= 2
    [file_path,file_name,file_ext] = fileparts(mat_file_address);
    file_name = [file_name file_ext];
end
% add filesep ('/' or '\') to the end of file_path
if ~strcmp(file_path(end), filesep)
    file_path = [file_path filesep];
end

% parts = strsplit(file_path, ["/", "\"]);
% name = erase(parts(length(parts) - 2), '-');
% file_name = char(extractBetween(name, 3, 15));

EXPERIMENT_PARAMS.mat_FileName = file_name;
EXPERIMENT_PARAMS.mat_PathName = file_path;
EXPERIMENT_PARAMS.file_name = file_name;
EXPERIMENT_PARAMS.flag_figure = flag_figure;
EXPERIMENT_PARAMS.flag_figure_debug = 1;

%% Load Data
fprintf('Loading: ')
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_qa flag_figure gui mat_file_address params funcs;
filename = EXPERIMENT_PARAMS.mat_FileName;
path_name = EXPERIMENT_PARAMS.mat_PathName;
load([path_name filename], 'data');
DLC.data = data;
DLC.FILE.path_to_analyzed = path_name;

fprintf([EXPERIMENT_PARAMS.file_name ' --> Completed. \n'])

%% Find video sampling rate (FPS)
% check if the FPS file exists in analyzed figs folder
fprintf('Finding video FPS: ')
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure gui params funcs;
path_name = EXPERIMENT_PARAMS.mat_PathName;
if ~strcmp(path_name(end), filesep);path_name = [path_name filesep];end

datehour = EXPERIMENT_PARAMS.file_name(1:13);

path_to_raw_data = [path_name '..' filesep ...
    '..' filesep '..' filesep 'raw_data' filesep];
path_to_analyzed_figs_tongue = [path_name '..' filesep ...
    '..' filesep '..' filesep 'analyzed_figs' filesep 'behavior_data' filesep 'tongue' filesep];
if isempty(dir(path_to_analyzed_figs_tongue))
    mkdir(path_to_analyzed_figs_tongue);
end

dir_FPS = dir([path_name, '*_video.mat']);
if ~isempty(dir_FPS)
    load([path_name dir_FPS(1).name],'FPS', 'height', 'width', 'duration', 'num_frames')
else
    fprintf('FPS not found, computing now ... \n')
    [LED_FPS, FPS, height, width, duration, num_frames] = PGH_estimate_vid_fps(path_to_raw_data,flag_figure, 30, 150, 1);
    save([path_name datehour '_video.mat'],'FPS','LED_FPS', 'height', 'width', 'duration', 'num_frames');
    saveas(gcf,[path_to_analyzed_figs_tongue  datehour '_FPS'], 'pdf');

end

EXPERIMENT_PARAMS.FPS = FPS;
DLC.FILE.vid_height = height;
DLC.FILE.vid_width = width;
DLC.FILE.duration = duration;
DLC.FILE.num_frames = num_frames;
fprintf([num2str(EXPERIMENT_PARAMS.FPS) ' --> Completed. \n'])

%% Use GUI
if gui 
    app = VL_LickSort;
    app_params(app, EXPERIMENT_PARAMS, DLC, params, funcs);
    while app.done == 0  % polling
      pause(0.05);
    end
    LICKS_ALL_DATA = app.LICKS_ALL_DATA; 
    EXPERIMENT_PARAMS = app.EXPERIMENT_PARAMS;
    app.closeWindow;   % close the app window
else
    %% Analyze DLC data
    fprintf(['Analyzing: ' EXPERIMENT_PARAMS.mat_FileName '\n'])
    clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_qa flag_figure gui mat_file_address params funcs;
    if flag_qa == 1
        [DLC,EXPERIMENT_PARAMS] = QA(DLC,EXPERIMENT_PARAMS);
    end
    [DLC,EXPERIMENT_PARAMS] = ANALYZE(DLC,EXPERIMENT_PARAMS, params, funcs);
    [LICKS_ALL_DATA, EXPERIMENT_PARAMS] = BUILD_LICKS_ALL_DATA(DLC, EXPERIMENT_PARAMS, params, funcs);
    PLOT_LICK_SORTER(LICKS_ALL_DATA, EXPERIMENT_PARAMS, params, funcs, 1);
end
end
