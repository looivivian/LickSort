%% MAIN Function: Plot lick sorter summary
% show_figure = 0: called from gui, shouldn't show gui but save it
% show_figure = 1: called from main script, shouldn't save it
function hFig = PLOT_LICK_SORTER(LICKS_ALL_DATA, EXPERIMENT_PARAMS, params, funcs, show_figure)
fprintf(['Ploting: ' EXPERIMENT_PARAMS.mat_FileName])
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs show_figure;
%% Plot LICKS_ALL_DATA
if show_figure && ~EXPERIMENT_PARAMS.flag_figure
    return;
end
params.cell_name = EXPERIMENT_PARAMS.file_name(1:13);
params.duration = EXPERIMENT_PARAMS.duration_video;
hFig = PGH_plot_lick_sorter(LICKS_ALL_DATA, params, show_figure);
fprintf(' --> Completed. \n')
end