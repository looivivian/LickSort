# LickSort
A data analysis MATLAB Application for marmoset licking experiement recordings

# Starting the Application
To start the application, run the monkey_lick_behavior.m script with the gui parameter set to true. Upon startup, the application will load and analyze the selected .mat file containing the coordinates of tracked points from the experiment recording (converted from DLC.csv produced by DeepLabCut). The application will remain non-interactive until the analysis is completed.

# Visualizing analyzed data
## Visualizing Video Frames
You can select a video frame to be displayed using a lick number (left of the video frame) or a video frame number (below the video frame). The most updated tracking coordinates are overlayed on the displayed frame. If the displayed frame belongs to any lick, the summary information of the lick and the kinematics of the specific frame will be shown.
## Lick Kinematics and Trace
Plots showing the displacement, velocity, trace, and geometrization of the tongue in the displayed frame will be updated if it belongs to any lick. 
## Summary Plots
A Summary plot and an End Position plot are shown on the right panel for an overview of the licking and harvesting behaviors throughout the recording.

# Quality Assurance (QA)
## Auto QA
To run the auto QA algorithm, press the QA button. Analysis will be re-run automatically on the data modified by the QA procedure. Frames with coordinates modified by the QA algorithm will be added to the drop-down list for browsing.
## Manual QA
Fixing individual tracking coordinates
If the tracking from DeepLabCut or the auto QA algorithm is inaccurate, you can fix them by selecting the point you want to modify on the video frame. In a pop-up window, indicate the correct position of the point with a double-click. The coordinates of the selected point will then be updated to the indicated position.
### Fixing lick types
If the analysis algorithm incorrectly classified any lick, you can fix it by indicating the correct lick type using the checkboxes on the left of the video frame. The lick type will be updated once the Re-classify button is pressed.
### Re-analyzing data after Manual QA
To update the lick information and summary plots after performing Manual QA, press the ANALYZE button. The original .mat file updated with data modified in the QA process will be analyzed and displayed.
### Discarding all QA progress
To discard all QA progress and start over with the original .mat data, press the RESTART button. Data from the original .mat file will be analyzed and displayed. 

# Save
To save all the analyzed data and produce the Lick Sorter Plot, press the SAVE button. Analysis will be re-run automatically if no analysis has been done on the data modified by the QA process. If any QA has been done in the analyzed data, an additional \_RECAL.mat file storing all updated coordinates will be saved. 

# Show Plot
After saving all analyzed data, you have the option to view the Lick Sorter Plot by pressing the PLOT button.

# Exiting the Application
You can exit the application by pressing the COMPLETE button. Two structs, [LICKS_ALL_DATA, EXPERIMENT_PARAMS], will be returned to the function call for monkey_lick_behavior.m.


