%% DEVAL PROJECT - MINISCOPE
% CODE TO REPRODUCE FIGURES FROM MANUSCRIPT
%
% Summary created 11/2026 by MMD


%% STEP 0 - CREATION OF MOUSE AND UNIT STRUCTS
% The various functions call on the mouse summary struct and the individual
% unit structs for each day. Since these structs have already been created,
% you have access to the completed forms without having to re-run these
% functions. They are listed here for reference.

mouse = deval_1_createStruct; %creates the mouse struct, which stores information about the mouse identity, sex, reward location, and directories for data
deval_2_add_ratemaps_toStruct(mouse); %initializes the unit structs with unit info and ratemaps, which are stored within the data directories
deval_3_check_place_coding(mouse); %checks whether cells meet a minimum discharge threshold and whether they are spatially stable
deval_4_identifyPlaceFields(mouse); %identifies place fields from ratemaps
deval_5_addSpikeTimes(mouse); %adds the spike times for each unit to the structs

%% NOTES ON OPTIONS IN FIGURE-MAKING FUNCTIONS
% In almost every figure-making functions, there will be the following
% options to set:

% saveOrNot = 1 to save figures, 0 to now
%   Make sure you set the saveDir to the correct directory

% days2use = which days you want to make plots for; the options are:
%   1 = pref1 = first day of preference
%   2 = pref2 = second day of preference
%   3 = preshift = day before the shift day (regardless of preference)
%   4 = downshift = day that the shift occurs
%   5 = postshift = day after shift (regardless of preference)
%   6 = adapt = day mouse showed behavioral adjustment
%   7 = revearsal = day mouse showed a reverse preference
%   8 = ext1 = extinction day 1
%   9 = ext2 = extinction day 2
%       *** The exception is deval_x_remapping_shift, which is written
%       specifically for the shift day. ***

% adpOnly = 1 to only include mice that showed behavioral adjustment

%% FIGURE 1
% Mice develop a preference for large rewards and exhibit behavioral
% adjustment after reward downshift.

% Panels A-D made by Nico
deval_x_frustrativeBehavior(mouse); %Panels E-G; p-values added in illustrator

%% FIGURE 2
% Population code for reward magnitude in dorsal CA1.

% Panels A-E made by Nico
deval_x_preProcess_plotRateMaps(mouse); %Panel F; plots example ratemaps are made for each mouse (both during running and stationary) and some selected for figure
deval_x_normRate_reward(mouse); %Panels G-J; also does the stats for those figures

%% FIGURE 3
% Representation of reward magnitudes is not driven by spatial encoding of
% place cells with place fields at reward locations.

deval_x_preProcess_plot2dPfDetection(mouse); %Panel A; plots example ratemaps are made for each mouse and some selected for figure
deval_x_2dPf_distToRew(mouse); %Panels B-D; also does stats; also saves raw data for Panel D, which Nico used to make the manuscript version of panel D
%   Above function makes other plots as well, including plots when
%   considering fields within 10 cm of the reward instead of 5 and pie
%   charts for each individual mouse.

%% FIGURE 4
% Reward-tethered (RT) place cells provide the cellular substrate for 
% reward-magnitude coding.

deval_x_remapping_shift(mouse); %this function makes all panel figures; see list below for figure names
% A: combinedRewardActiveCells_pie
% B: LargeRewActCells_rateMaps_forNico_Pre_1 and LargeRewActCells_Path_forNico_Pre
% C: SmallRewActCells_rateMaps_forNico_Pre_1 and SmallRewActCells_Path_forNico_Pre
% D: UnstableCells_rateMaps_forNico_Pre_1 and UnstableCells_Path_forNico_Pre
% E-G: rewFr_byRTtype

%% FIGURE 5
% Selective remapping of Large reward-tethered (LRT) place cells following reward magnitude downshift.

deval_x_remapping_shift(mouse); %this function makes all panel figures; see list below for figure names
% A: LargeRewFieldCells_path_forNico_Pre and LargeRewFieldCells_path_forNico_Post and LargeRewFieldCells_rateMaps_forNico_Pre_1 and LargeRewFieldCells_rateMaps_forNico_Post_1
% B: SmallRewFieldCells_path_forNico_Pre and SmallRewFieldCells_path_forNico_Post and SmallRewFieldCells_rateMaps_forNico_Pre_1 and SmallRewFieldCells_rateMaps_forNico_Post_1
% C: LargeRewActCells_path_forNico_Pre and LargeRewActCells_path_forNico_Post and LargeRewActCells_rateMaps_forNico_Pre_1 and LargeRewActCells_rateMaps_forNico_Post_1
% D: SmallRewActCells_path_forNico_Pre and SmallRewActCells_path_forNico_Post and SmallRewActCells_rateMaps_forNico_Pre_1 and SmallRewActCells_rateMaps_forNico_Post_1
% E: UnstableCells_path_forNico_Pre and UnstableCells_path_forNico_Post and UnstableCells_rateMaps_forNico_Pre_1 and UnstableCells_rateMaps_forNico_Post_1
% F: rewardFields_RT_spatCorr
