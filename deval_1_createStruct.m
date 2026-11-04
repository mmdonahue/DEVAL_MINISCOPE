function mouse = deval_1_createStruct(varargin)
% function mouse = deval_1_createStruct(mouse, dataDir)
%
% PURPOSE:
%   Initialize a data struct for the deval project with the basic position
%   data and trial times.
%
% INPUT:
%   mouse = in-progress struct
%   dataDir = directory with raw data
%
% OUTPUT:
%   mouse = data struct
%
% MMD
% 05/2025

if nargin == 0
    mouse = deval_0_hardCode;
    dataDir = 'D:\Devaluation Project Miniscope\dev analysis';
end %set defaults

if nargin == 1
    mouse = varargin{1};
    dataDir = 'D:\Devaluation Project Miniscope\dev analysis';
end %dataDir

if nargin == 2
    mouse = varargin{1};
    dataDir = varargin{2};
end

%% INITIALIZE FOR COORDS

x_arenaSize = 56;
y_arenaSize = 56;
binSz = 4;

%% INITIALIZE FOR READING IN DATA

tNames = {'PreShiftTest', 'Shift', 'postShift', 'Ext'}; %leaving out habituation and pre shift training for now
tFoldNames = {'preshift*', 'Shif*', 'post*', 'Ex*'}; %for indexing folders with slightly different names in each mouse
ttNames = {'Free', 'Forced'};

cd(dataDir)
% mDirs = dir('M*');
% mDirs = mDirs(vertcat(mDirs.isdir)); %delete excel sheets etc.

%% MAKE STRUCT

for m = 1:length(mouse)
    fprintf('Mouse %s\n', mouse(m).name)
    mouse(m).name = mouse(m).name;
    mouse(m).dir = [dataDir '\' mouse(m).name];

    cd(mouse(m).dir)

    for ts = 1:length(tNames)
        fprintf('\t%s\n', tNames{ts})
        mouse(m).task(ts).name = tNames{ts};

        taskDirs = dir(tFoldNames{ts});
        if ts == 1
            taskDirs = flip(taskDirs);
        end %flip pre-shift so we go -3 -2 -1

        if isempty(taskDirs)
            if ts == 4
                continue
            end % some mife don't have extinction (?)
            keyboard
        end %not supposed to happen

        for d = 1:length(taskDirs)
            fprintf('\t\tDay %d/%d\n', d, length(taskDirs))
            % fprintf('\t\t\tTrial')
            cd(taskDirs(d).name)
            mouse(m).task(ts).day(d).name = taskDirs(d).name;
            mouse(m).task(ts).day(d).dir = pwd;

            tmDir = dir('ts*');
            tmConv = readmatrix(tmDir.name);
            posFile = dir('pos*fixed*fixed*');
            if length(posFile) ~= 1
                keyboard
            end %check if there are multiple pos files or none

            tmpPos = readmatrix(posFile.name);
            tmpPos(:,1) = 1:length(tmpPos); %first 6 rows are all 6 (????)
            coords = convert_coords(tmpPos, x_arenaSize, y_arenaSize);
            if m == 1 && ts == 1 && d == 3 || m == 1 && ts == 2 && d == 1 || m == 1 && ts == 3 && d == 1 ...
                    || m == 2 && ts == 1 && d == 3  || m == 2 && ts == 2 && d == 1 || m == 2 && ts == 3 && d == 1 ...
                    || m == 2 && ts == 3 && d == 2 || m ==4 && ts == 1 && d == 3
                coords(:,3) = coords(:,3) * -1 + 56;
            end %corods are flipped on this day for some reason
            mouse(m).task(ts).day(d).coords = coords;

            linCoords = linearize_fig8_coords(coords, binSz);
            mouse(m).task(ts).day(d).linCoords = linCoords;

            events = readtable('EventFile.csv');
            tmpStr = string(events.Var1');

            if ts ~=2
                tInds = find(contains(tmpStr, 'trial'));
            else
                tInds = sort([find(contains(tmpStr, 'Train')) find(contains(tmpStr, 'Test'))]);
                % tIndsInt = sort([find(contains(tmpStr, 'Train')) find(contains(tmpStr, 'Test'))]);
                % checkDs = find(contains(tmpStr, 'Downshift'));
                % tInds = ismember(tIndsInt, checkDs)
                if length(tInds) < 59
                    keyboard
                end
            end %shift days do not include "trial" for the forced/training trials for some reason

            cheatCheck = find(contains(tmpStr, 'Cheat')); %check if any cheating trials needs to be eliminated
            cheatTs = [];
            if ~isempty(cheatCheck)
                tInds(tInds==cheatCheck) = [];

                for ct = 1:length(cheatCheck)
                    tmpT = find(cheatCheck(ct) > tInds, 1, 'Last');
                    cheatTs = [cheatTs tmpT];
                end %find cheating trials
            end %there was cheating

            for t = 1:length(tInds)
                % fprintf('\t%d', t)
                if ismember(t, cheatTs)
                    continue
                end
                if tInds(t) + 4 > height(events)
                    continue
                end
                if startsWith(events.Var1(tInds(t)), 'Test')
                    mouse(m).task(ts).day(d).trial(t).type = 1;
                elseif startsWith(events.Var1(tInds(t)), 'Train')
                    mouse(m).task(ts).day(d).trial(t).type = 2;
                else
                    keyboard
                end %trial type

                tms = zeros(1,5); %initialize - start choice poke ITI ITIend
                tmInds = zeros(1,5); %initialize - start choice poke ITI ITIend

                tms(1) = events.Var2(tInds(t)); %trial start

                if strcmp(events.Var1(tInds(t)+1), 'West') == 0 ...
                        && strcmp(events.Var1(tInds(t)+1), 'East') == 0

                    checkRew = 0; %initialize as bad
                    addCheck = 1;
                    while checkRew == 0
                        addCheck = addCheck + 1;
                        try
                            checkRew = strcmp(events.Var1(tInds(t)+addCheck), 'West') ...
                                || strcmp(events.Var1(tInds(t)+addCheck), 'East');
                        catch; keyboard; end
                    end %while it's still bad
                    mouse(m).task(ts).day(d).trial(t).rewLoc = events.Var1(tInds(t)+addCheck); %saved reward location
                    tms(2) = events.Var2(tInds(t)+addCheck); %choice point door close
                else
                    mouse(m).task(ts).day(d).trial(t).rewLoc = events.Var1(tInds(t)+1); %saved reward location
                    tms(2) = events.Var2(tInds(t)+1); %choice point door close
                end %check that something isn't wrong


                if t ~= length(tInds)
                    pullEvs = events.Var1(tInds(t)+1:tInds(t+1)-1);
                else
                    pullEvs = events.Var1(tInds(t)+1:end);
                end %pull in between events to find first poke
                tmpFunc = @(x)(startsWith(x, 'Poke'));
                findPoke = cellfun(tmpFunc, pullEvs);
                addInds = find(findPoke, 1, 'First');
                if isempty(addInds)
                    tms(3) = NaN; %rat doesn't poke for some reason
                else
                    tms(3) = events.Var2(tInds(t) + addInds);
                end

                % 
                % if ~isnan(tms(3))
                %     if strncmp(mouse(m).task(ts).day(d).trial(t).rewLoc , 'East', 4) && abs(diff([linCoords(tms(3),2) 35])) > 5 ...
                %             || strncmp(mouse(m).task(ts).day(d).trial(t).rewLoc , 'West', 4) && abs(diff([linCoords(tms(3),2) 63])) > 5
                %         keyboard
                %     end
                % end

                % if ~startsWith(events.Var1(tInds(t)+3), 'Poke')
                %
                %     tms(3) = NaN; %first nose poke
                % else
                %     tms(3) = events.Var2(tInds(t)+3); %first nose poke
                % end %check - rat doesn't poke for some reason

                if t ~= length(tInds)
                    if ~startsWith(events.Var1(tInds(t+1)-2), 'ITI starts')
                        pullEvs = events.Var1(tInds(t)+1:tInds(t+1)-1);
                        tmpFunc = @(x)(startsWith(x, 'ITI starts'));
                        findITIs = find(cellfun(tmpFunc, pullEvs));
                        tms(4) = events.Var2(tInds(t) + findITIs); %ITI starts
                        tms(5) = events.Var2(tInds(t) + findITIs + 1); %ITI end
                    end %something's wrong

                    tms(4) = events.Var2(tInds(t+1)-2); %ITI starts
                    tms(5) = events.Var2(tInds(t+1)-1); %ITI end
                else
                    tms(4) = events.Var2(end); %ITI starts
                    tms(5) = NaN; %ITI end
                end %not the very last trial

                try
                    mouse(m).task(ts).day(d).trial(t).tms = transpose(coords(tms,1));
                catch
                    badInds = find(isnan(tms));
                    tmpTms = tms;
                    tmpTms(badInds) = 1; %just temporarily
                    tmpFillTms = transpose(coords(tmpTms,1));
                    tmpFillTms(badInds) = NaN;
                    mouse(m).task(ts).day(d).trial(t).tms = tmpFillTms;
                end
                mouse(m).task(ts).day(d).trial(t).tmInds = tms;

            end %trial
            % fprintf('\n')
            if ~isempty(cheatTs) && length(mouse(m).task(ts).day(d).trial) > cheatTs
                try
                    mouse(m).task(ts).day(d).trial(cheatTs) = [];
                catch; keyboard; end
            end
            
            allTTs = vertcat(mouse(m).task(ts).day(d).trial(:).type);
            CC = bwconncomp(allTTs - 1);
            if ts ~= 2
                mouse(m).task(ts).day(d).phase(1).name = 'free choice';
                mouse(m).task(ts).day(d).phase(1).trials = [1:CC.PixelIdxList{1}(1)-1]';

                mouse(m).task(ts).day(d).phase(2).name = 'forced choice';
                mouse(m).task(ts).day(d).phase(2).trials = [CC.PixelIdxList{1}];
            else
                mouse(m).task(ts).day(d).phase(1).name = 'forced choice pre';
                mouse(m).task(ts).day(d).phase(1).trials = CC.PixelIdxList{1};

                mouse(m).task(ts).day(d).phase(2).name = 'free choice';
                mouse(m).task(ts).day(d).phase(2).trials = [CC.PixelIdxList{1}(end)+1:CC.PixelIdxList{2}(1)-1]';

                mouse(m).task(ts).day(d).phase(3).name = 'forced choice post';
                mouse(m).task(ts).day(d).phase(3).trials = CC.PixelIdxList{2};

            end %set phases

            cd ../
        end %days within the task - shouldn't be up to 10 so alphabetical should work

    end %task
    cd ../
end %mouse


end %function

%% HARD CODE DATA

function mouse = deval_0_hardCode

mouse(1).name = 'M119';
mouse(1).sex = 'F';
mouse(1).largeRew = 'West';
mouse(1).adp = 1;

mouse(2).name = 'M120';
mouse(2).sex = 'F';
mouse(2).largeRew = 'East';
mouse(2).adp = 0;

mouse(3).name = 'M169';
mouse(3).sex = 'M';
mouse(3).largeRew = 'East';
mouse(3).adp = 1;

mouse(4).name = 'M253';
mouse(4).sex = 'M';
mouse(4).largeRew = 'West';
mouse(4).adp = 1;

mouse(5).name = 'M254';
mouse(5).sex = 'M';
mouse(5).largeRew = 'West';
mouse(5).adp = 1;

mouse(6).name = 'M291';
mouse(6).sex = 'M';
mouse(6).largeRew = 'East';
mouse(6).adp = 1;

mouse(7).name = 'M355';
mouse(7).sex = 'F';
mouse(7).largeRew = 'West';
mouse(7).adp = 1;

mouse(8).name = 'M356';
mouse(8).sex = 'F';
mouse(8).largeRew = 'East';
mouse(8).adp = 1;

mouse(9).name = 'M367';
mouse(9).sex = 'F';
mouse(9).largeRew = 'West';
mouse(9).adp = 1;

mouse(10).name = 'M368';
mouse(10).sex = 'M';
mouse(10).largeRew = 'East';
mouse(10).adp = 1;

end %helper function