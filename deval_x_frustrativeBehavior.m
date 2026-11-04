function deval_x_frustrativeBehavior(mouse)
% function deval_x_frustrativeBehavior(mouse)
% 
% PURPOSE:
%   Plot metrics and schematics related to possible frustrative behvaiors
%   displayed by the mice during the deval task.
% 
% INPUTS:
%   mouse = data struct
% 
% OPTIONS (found in function):
%   saveOrNot = 1 to save figures, 0 to not
%   days2use = days you want to create figures for
% 
% MMD
% Ewell Lab
% 02/2026

%% OPTIONS

saveOrNot = 0;

% days2use = [2 4 6 8]; %days I'm interested in; set 1:9 for all
days2use = 4;

adpOnly = 1; %set to 1 to exclude mice that did not adapt

%% INITIALIZE

tNames = {'Preference1', 'Preference2', 'preshift', 'Downshift', 'postshift', 'Adaptation', 'Reversal', 'Extinction1', 'Extinction2'};
cols = {'Blue', 'Red'};
rewNames = {'Large reward', 'Small reward'};

saveDir = 'D:\Devaluation Project Miniscope\Meg analysis\RESULTS\behavior';
% cd(saveDir)

%% GET DATA

for d = days2use
    if d == 4
        continue
    end %shift day is weird and will be redone below

    propPoke = cell(length(mouse), 2); %proportion of trials mice poke to receive reward at
    latTurn = cell(length(mouse), 2);  %latency from ITI end to
    pathLen = cell(length(mouse), 2); %path length from delay end to reward

    for m = 1:length(mouse)
        if adpOnly == 1 && mouse(m).adp == 0
            continue
        end %ignore non-adapted mouse

        if isempty(mouse(m).day(d).dir)
            continue
        end

        events = readtable([mouse(m).day(d).dir '\EventFile.csv']);
        coords = mouse(m).day(d).coords;

        tmpStr = string(events.Var1');
        tInds = sort(find(contains(tmpStr, 'Train')));

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

            if t ~= length(tInds)
                pullEvs = events.Var1(tInds(t)+1:tInds(t+1)-1);
            else
                pullEvs = events.Var1(tInds(t)+1:end);
            end %pull in between events to find first poke

            %large or small reward
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
                rewLoc = events.Var1{tInds(t)+addCheck}; %saved reward location
                tmInd = events.Var2(tInds(t)+addCheck);
            else
                rewLoc = events.Var1{tInds(t)+1}; %saved reward location
                tmInd =  events.Var2(tInds(t)+1);
            end %check that something isn't wrong

            if strcmp(rewLoc, mouse(m).largeRew)
                rewInd = 1;
            else
                rewInd = 2;
            end %large or small reward

            tmpFunc = @(x)(startsWith(x, 'Reward'));
            pokeInd = sum((cellfun(tmpFunc, pullEvs))); %1=poked
            if pokeInd > 1
               pokeInd = 1;
            end
            propPoke{m,rewInd} = [propPoke{m,rewInd}; pokeInd];

            pullCoords = coords(events.Var2(tInds(t)):tmInd,:);
            % plot(pullCoords(:,2), pullCoords(:,3)); xlim([0 56]); ylim([0 56])
            findPass = find(pullCoords(:,2)> 50, 1, 'Last'); %find when the mice pass the 50 cm mark, indicating that they are going into an arm

            latTurn{m,rewInd} = [latTurn{m,rewInd}; coords(events.Var2(tInds(t))+findPass,1) - coords( events.Var2(tInds(t)),1)];

            pullCoords = pullCoords(1:findPass,:);
            xDiff = diff(pullCoords(:,2)); %distance mouse traveled in x dim in between each frame
            yDiff = diff(pullCoords(:,3)); %distance mouse traveled in y dim in between each frame
            tmpSpd = hypot(xDiff, yDiff); %speed in cm/frame
            pathLen{m,rewInd} = [pathLen{m,rewInd}; sum(tmpSpd)];

        end %trial - training aka forced trials only
    end %mouse

    figtitle = ['frustrationBehaviors_' tNames{d}];
    figure('Name', figtitle)

    subplot(1,3,1)
    tmpFunc = @(x)(sum(x) ./ length(x));
    tmpData = cellfun(tmpFunc, propPoke);
    tmpMean = nanmean(tmpData);

    for r = 1:2
        xData = rand(size(tmpData,1),1)*0.5 - 0.25 + r;
        yData = tmpData(:,r);

        hold on
        scatter(xData, yData, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb(cols{r}), 'MarkerFaceAlpha', 0.7)
        line([r-0.2 r+0.2], [tmpMean(r) tmpMean(r)], 'Color', 'k')
    end %reward

    xticks(1:2)
    xticklabels(rewNames)
    ylabel('Proportion rewards sampled')
    ylim([0 1])
    axis square

    tmpData = tmpData(~isnan(tmpData(:,1)),:);
    [p, ~, stats] = signrank(tmpData(:,1), tmpData(:,2));
    title([tNames{d} ' (p = ' num2str(round(p,3)) ')'])



    % [p,tbl,stats,terms] = anovan(statData, {statRew, statPh}, 'model', 'interaction', 'varnames', {'Reward', 'Phase'});



    subplot(1,3,2)
    tmpData = cellfun(@mean, latTurn);
    tmpMean = nanmean(tmpData);

    for r = 1:2
        xData = rand(size(tmpData,1),1)*0.5 - 0.25 + r;
        yData = tmpData(:,r);

        hold on
        scatter(xData, yData, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb(cols{r}), 'MarkerFaceAlpha', 0.7)
        line([r-0.2 r+0.2], [tmpMean(r) tmpMean(r)], 'Color', 'k')
    end %reward

    xticks(1:2)
    xticklabels(rewNames)
    ylabel('Time to reward arm (s)')
    ax = gca;
    ylim([0 ax.YLim(2)])
    axis square

    tmpData = tmpData(~isnan(tmpData(:,1)),:);
    [p, ~, stats] = signrank(tmpData(:,1), tmpData(:,2));
    title([tNames{d} ' (p = ' num2str(round(p,3)) ')'])

    subplot(1,3,3)
    tmpData = cellfun(@mean, pathLen);
    tmpMean = nanmean(tmpData);

    for r = 1:2
        xData = rand(size(tmpData,1),1)*0.5 - 0.25 + r;
        yData = tmpData(:,r);

        hold on
        scatter(xData, yData, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb(cols{r}), 'MarkerFaceAlpha', 0.7)
        line([r-0.2 r+0.2], [tmpMean(r) tmpMean(r)], 'Color', 'k')
    end %reward

    xticks(1:2)
    xticklabels(rewNames)
    ylabel('Path distance to reward arm (cm)')
    ax = gca;
    ylim([0 ax.YLim(2)])
    axis square

    tmpData = tmpData(~isnan(tmpData(:,1)),:);
    [p, ~, stats] = signrank(tmpData(:,1), tmpData(:,2));
    title([tNames{d} ' (p = ' num2str(round(p,3)) ')'])

    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save


end %day

if ismember(4, days2use)
    d = 4;

    phNames = {'Pre-shift', 'Post-shift'};


    propPoke = cell(length(mouse), 2, 2); %proportion of trials mice poke to receive reward at mouse x rew x phase
    latTurn = cell(length(mouse), 2, 2);  %latency from ITI end to
    pathLen = cell(length(mouse), 2, 2); %path length from delay end to reward
    pathsEx = cell(length(mouse), 2, 2); %for plotting example paths

    for m = 1:length(mouse)
        if adpOnly == 1 && mouse(m).adp == 0
            continue
        end %ignore non-adapted mouse

        if isempty(mouse(m).day(d).dir)
            continue
        end

        events = readtable([mouse(m).day(d).dir '\EventFile.csv']);
        coords = mouse(m).day(d).coords;

        tmpStr = string(events.Var1');
        tInds = sort(find(contains(tmpStr, 'Train')));

        cheatCheck = find(contains(tmpStr, 'Cheat')); %check if any cheating trials needs to be eliminated
        cheatTs = [];
        if ~isempty(cheatCheck)
            tInds(tInds==cheatCheck) = [];

            for ct = 1:length(cheatCheck)
                tmpT = find(cheatCheck(ct) > tInds, 1, 'Last');
                cheatTs = [cheatTs tmpT];
            end %find cheating trials
        end %there was cheating

        pullTNames = table2cell(events(tInds,1));
        tmpFunc = @(x)(str2num(x(12:end)));
        tNums = cellfun(tmpFunc,pullTNames);
        tmpBnry = zeros(1,tNums(end));
        tmpBnry(tNums) = 1;

        CC = bwconncomp(tmpBnry);

        for t = 1:length(tInds)
            % fprintf('\t%d', t)
            if ismember(t, cheatTs)
                continue
            end

            if tInds(t) + 4 > height(events)
                continue
            end

            if tNums(t) == CC.PixelIdxList{1}(end)
                checkNext = 0; %initialize as no
                addCheck = 0;
                while checkNext == 0
                    addCheck = addCheck  + 1;
                    checkNext = contains(events.Var1{tInds(t)+addCheck}, 'trial');
                end %while
                pullEvs = events.Var1(tInds(t):tInds(t)+addCheck-1);
            elseif t ~= length(tInds)
                pullEvs = events.Var1(tInds(t)+1:tInds(t+1)-1);
            elseif m == 4 && t == 40
                %something weird happened this day
                checkNext = 0; %initialize as no
                addCheck = 0;
                while checkNext == 0
                    addCheck = addCheck  + 1;
                    checkNext = contains(events.Var1{tInds(t)+addCheck}, 'trial');
                end %while
                pullEvs = events.Var1(tInds(t):tInds(t)+addCheck-1);
            else
                pullEvs = events.Var1(tInds(t)+1:end);
            end %pull in between events to find first poke

            %large or small reward
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
                rewLoc = events.Var1{tInds(t)+addCheck}; %saved reward location
                tmInd = events.Var2(tInds(t)+addCheck);
            else
                rewLoc = events.Var1{tInds(t)+1}; %saved reward location
                tmInd =  events.Var2(tInds(t)+1);
            end %check that something isn't wrong

            if strcmp(rewLoc, mouse(m).largeRew)
                rewInd = 1;
            else
                rewInd = 2;
            end %large or small reward

            tmpFunc = @(x)(ismember(str2num(events.Var1{tInds(t)}(12:end)), x));
            ph = find(cellfun(tmpFunc, CC.PixelIdxList));

            tmpFunc = @(x)(startsWith(x, 'Reward'));
            pokeInd = sum((cellfun(tmpFunc, pullEvs))); %1=poked
            if pokeInd > 1
                keyboard
            end
            propPoke{m,rewInd,ph} = [propPoke{m,rewInd,ph}; pokeInd];

            pullCoords = coords(events.Var2(tInds(t)):tmInd,:);
            % plot(pullCoords(:,2), pullCoords(:,3)); xlim([0 56]); ylim([0 56])
            findPass = find(pullCoords(:,2)> 50, 1, 'Last'); %find when the mice pass the 50 cm mark, indicating that they are going into an arm

            latTurn{m,rewInd,ph} = [latTurn{m,rewInd,ph}; coords(events.Var2(tInds(t))+findPass,1) - coords( events.Var2(tInds(t)),1)];

            pullCoords = pullCoords(1:findPass,:);
            xDiff = diff(pullCoords(:,2)); %distance mouse traveled in x dim in between each frame
            yDiff = diff(pullCoords(:,3)); %distance mouse traveled in y dim in between each frame
            tmpSpd = hypot(xDiff, yDiff); %speed in cm/frame
            pathLen{m,rewInd,ph} = [pathLen{m,rewInd,ph}; sum(tmpSpd)];
            % if sum(tmpSpd) >900
            %     keyboard
            % end

            %flip to align large rewards on the same side always
            if strcmp(mouse(m).largeRew, 'West')
                pullCoords(:,3) = pullCoords(:,3) * -1 + 56;
            end %flip if necessary
            pathsEx{m,rewInd,ph} = cat(2, pathsEx{m,rewInd,ph}, {pullCoords});
        end %trial
    end %mouse

    keyboard

    figtitle = ['frustrationBehaviors_' tNames{d}];
    figure('Name', figtitle)
    spMap = [1 2 3; 4 5 6];

    tmpFunc = @(x)(sum(x) ./ length(x));
    tmpData = cellfun(tmpFunc, propPoke);
    tmpMean = nanmean(tmpData);

    for ph = 1:2
        subplot(2,3,spMap(ph,1))
        for r = 1:2
            xData = rand(size(tmpData,1),1)*0.5 - 0.25 + r;
            yData = tmpData(:,r,ph);

            hold on
            scatter(xData, yData, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb(cols{r}), 'MarkerFaceAlpha', 0.7)
            line([r-0.2 r+0.2], [tmpMean(1,r,ph) tmpMean(1,r,ph)], 'Color', 'k')
        end %reward
        xticks(1:2)
        xticklabels(rewNames)
        ylabel('Proportion rewards sampled')
        ylim([0 1])
        axis square
        title([tNames{d} ' - ' phNames{ph}])

    end %phase

    statData = tmpData(:); %rearrange
    statRew = [ones(size(tmpData,1),1); ones(size(tmpData,1),1)*2; ones(size(tmpData,1),1); ones(size(tmpData,1),1)*2];
    statPh = [ones(size(tmpData,1)*2,1); ones(size(tmpData,1)*2,1)*2];

    statRew = statRew(~isnan(statData));
    statPh = statPh(~isnan(statData));
    statData = statData(~isnan(statData));

    % [p,tbl,stats,terms] = anovan(statData, {statRew, statPh}, 'model', 'interaction', 'varnames', {'Reward', 'Phase'});

    % c = multcompare(stats);



    tmpData = cellfun(@nanmean, latTurn);
    tmpMean = nanmean(tmpData);

    for ph = 1:2
        subplot(2,3,spMap(ph,2))

        for r = 1:2
            xData = rand(size(tmpData,1),1)*0.5 - 0.25 + r;
            yData = tmpData(:,r,ph);

            hold on
            scatter(xData, yData, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb(cols{r}), 'MarkerFaceAlpha', 0.7)
            line([r-0.2 r+0.2], [tmpMean(1,r,ph) tmpMean(1,r,ph)], 'Color', 'k')
        end %reward

        xticks(1:2)
        xticklabels(rewNames)
        ylabel('Time to reward arm (s)')
        ax = gca;
        ylim([0 ax.YLim(2)])
        axis square
        title([tNames{d} ' - ' phNames{ph}])
    end %phase

    tmpData = cellfun(@nanmean, pathLen);
    tmpMean = nanmean(tmpData);

    for ph = 1:2
        subplot(2,3,spMap(ph,3))

        for r = 1:2
            xData = rand(size(tmpData,1),1)*0.5 - 0.25 + r;
            yData = tmpData(:,r,ph);

            hold on
            scatter(xData, yData, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb(cols{r}), 'MarkerFaceAlpha', 0.7)
            line([r-0.2 r+0.2], [tmpMean(1,r,ph) tmpMean(1,r,ph)], 'Color', 'k')
        end %reward

        xticks(1:2)
        xticklabels(rewNames)
        ylabel('Path distance to reward arm (cm)')
        ax = gca;
        ylim([0 ax.YLim(2)])
        axis square
        title([tNames{d} ' - ' phNames{ph}])
    end %phase

    same_axes([2 5]);
    same_axes([3 6]);

    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save

    keyboard

    %% ALT FIG

    figtitle = ['frustrationBehaviors_' tNames{d} '_alt'];
    figure('Name', figtitle)
    tickMap = [1 2; 4 5];

    subplot(1,3,1)
    tmpFunc = @(x)(sum(x) ./ length(x));
    tmpData = cellfun(tmpFunc, propPoke);
    tmpMean = nanmean(tmpData);

    hold on
    for ph = 1:2
        lh = NaN(1,2);
        for r = 1:2
            bgraph = bar(tickMap(ph,r), tmpMean(1,r,ph), 'FaceColor', 'flat');
            lh (r) = bgraph;
            bgraph.CData = rgb(cols{r});

            xData = rand(size(tmpData,1),1)*0.5 - 0.25 + tickMap(ph,r);
            yData = tmpData(:,r,ph);

            % scatter(xData, yData, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb(cols{r}), 'MarkerFaceAlpha', 0.7)
            scatter(xData, yData, 'MarkerEdgeColor', rgb(cols{r}), 'MarkerFaceColor', rgb('White'))

        end %reward
    end %phase

    axis square
    xticks([1.5 4.5])
    xticklabels(phNames)
    ylabel('Proportion rewards sampled')
    ylim([0 1])
    title(tNames{d})
    legend(lh, rewNames, 'Location', 'southoutside', 'Box', 'off')


    statData = tmpData(:); %rearrange
    statRew = [ones(size(tmpData,1),1); ones(size(tmpData,1),1)*2; ones(size(tmpData,1),1); ones(size(tmpData,1),1)*2];
    statPh = [ones(size(tmpData,1)*2,1); ones(size(tmpData,1)*2,1)*2];

    statRew = statRew(~isnan(statData));
    statPh = statPh(~isnan(statData));
    statData = statData(~isnan(statData));

    [p,tbl,stats,terms] = anovan(statData, {statRew, statPh}, 'model', 'interaction', 'varnames', {'Reward', 'Phase'}, 'Display', 'off');
    %
    % c = multcompare(stats); NO MULTCOMPARE, INTERACTION P > 0.05

    subplot(1,3,2)
    tmpData = cellfun(@nanmean, latTurn);
    tmpMean = nanmean(tmpData);

    hold on
    for ph = 1:2
        lh = NaN(1,2);
        for r = 1:2
            bgraph = bar(tickMap(ph,r), tmpMean(1,r,ph), 'FaceColor', 'flat');
            lh (r) = bgraph;
            bgraph.CData = rgb(cols{r});

            xData = rand(size(tmpData,1),1)*0.5 - 0.25 + tickMap(ph,r);
            yData = tmpData(:,r,ph);

            % scatter(xData, yData, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb(cols{r}), 'MarkerFaceAlpha', 0.7)
            scatter(xData, yData, 'MarkerEdgeColor', rgb(cols{r}), 'MarkerFaceColor', rgb('White'))

        end %reward
    end %phase

    axis square
    xticks([1.5 4.5])
    xticklabels(phNames)
    ylabel('Time to reward arm (s)')
    title(tNames{d})
    legend(lh, rewNames, 'Location', 'southoutside', 'Box', 'off')

    combData = vertcat(latTurn(:));
    tmpData = vertcat(combData{:});

    tmpNs = sum(cellfun(@length, latTurn));
    tmpRew = [repmat(1, tmpNs(1,1,1), 1); repmat(2, tmpNs(1,2,1), 1); repmat(1, tmpNs(1,1,2), 1); repmat(2, tmpNs(1,2,2), 1)];
    tmpPh =  [repmat(1, tmpNs(1,1,1), 1); repmat(1, tmpNs(1,2,1), 1); repmat(2, tmpNs(1,1,2), 1); repmat(2, tmpNs(1,2,2), 1)];

    tmpMouse = []; %initialize
    for ph = 1:2
        for r = 1:2
            for m = 1:length(mouse)
                tmpMouse = [tmpMouse; repmat(m, length(latTurn{m,r,ph}), 1)];
            end %mouse
        end %reward
    end %phase

    % ANOVA
    [p,tbl,stats,terms] = anovan(tmpData, {tmpRew, tmpPh}, 'model', 'interaction', 'varnames', {'Reward', 'Phase'}, 'Display', 'off');
    [c,m,h,gnames] = multcompare(stats, 'Dimension', [1 2], 'Display', 'off');
    tbl = array2table(c,"VariableNames", ...
        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl.("Group A")=gnames(tbl.("Group A"));
    tbl.("Group B")=gnames(tbl.("Group B"))

    subplot(1,3,3)
    tmpData = cellfun(@nanmean, pathLen);
    tmpMean = nanmean(tmpData);

    hold on
    for ph = 1:2
        lh = NaN(1,2);
        for r = 1:2
            bgraph = bar(tickMap(ph,r), tmpMean(1,r,ph), 'FaceColor', 'flat');
            lh (r) = bgraph;
            bgraph.CData = rgb(cols{r});

            xData = rand(size(tmpData,1),1)*0.5 - 0.25 + tickMap(ph,r);
            yData = tmpData(:,r,ph);

            % scatter(xData, yData, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb(cols{r}), 'MarkerFaceAlpha', 0.7)
            scatter(xData, yData, 'MarkerEdgeColor', rgb(cols{r}), 'MarkerFaceColor', rgb('White'))

        end %reward
    end %phase

    axis square
    xticks([1.5 4.5])
    xticklabels(phNames)
    ylabel('Path distance to reward arm (cm)')
    title(tNames{d})
    legend(lh, rewNames, 'Location', 'southoutside', 'Box', 'off')

    combData = vertcat(pathLen(:));
    tmpData = vertcat(combData{:});

    tmpNs = sum(cellfun(@length, pathLen));
    tmpRew = [repmat(1, tmpNs(1,1,1), 1); repmat(2, tmpNs(1,2,1), 1); repmat(1, tmpNs(1,1,2), 1); repmat(2, tmpNs(1,2,2), 1)];
    tmpPh =  [repmat(1, tmpNs(1,1,1), 1); repmat(1, tmpNs(1,2,1), 1); repmat(2, tmpNs(1,1,2), 1); repmat(2, tmpNs(1,2,2), 1)];

    tmpMouse = []; %initialize
    for ph = 1:2
        for r = 1:2
            for m = 1:length(mouse)
                tmpMouse = [tmpMouse; repmat(m, length(latTurn{m,r,ph}), 1)];
            end %mouse
        end %reward
    end %phase

    % ANOVA
    [p,tbl,stats,terms] = anovan(tmpData, {tmpRew, tmpPh}, 'model', 'interaction', 'varnames', {'Reward', 'Phase'}, 'Display', 'off');
    [c,m,h,gnames] = multcompare(stats, 'Dimension', [1 2], 'Display', 'off');
    tbl = array2table(c,"VariableNames", ...
        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl.("Group A")=gnames(tbl.("Group A"));
    tbl.("Group B")=gnames(tbl.("Group B"))


    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save

    figtitle = ['ExamplePaths_' tNames{d} '_preShift'];
    figure('Name', figtitle)

    %polygons
    lrPgX = [0.5 2.5 2.5 0.5]*4; %large reward
    lrPgY = [0.5 0.5 3.5 3.5]*4;
    lrPg = polyshape(lrPgX, lrPgY);
    smPgX = [0.5 2.5 2.5 0.5]*4; %small reward
    smPgY = [11.5 11.5 14.5 14.5]*4;
    smPg = polyshape(smPgX, smPgY);
    dPgX = [0.5 2.5 2.5 5.5 5.5 2.5  2.5 0.5]*4; %delay
    dPgY = [3.5 3.5 6.5 6.5 8.5 8.5 11.5 11.5]*4;
    dPg = polyshape(dPgX, dPgY);
    stPgX = [5.5 14.5 14.5 5.5]*4; %stem
    stPgY = [6.5 6.5 8.5 8.5]*4;
    stPg = polyshape(stPgX, stPgY);
    lrArmPgX = [2.5 14.5 14.5 12.5 12.5 2.5]*4; %large reward arm
    lrArmPgY = [0.5 0.5 6.5 6.5 2.5 2.5]*4;
    lrAPg = polyshape(lrArmPgX, lrArmPgY);
    smArmPgX = [2.5 12.5 12.5 14.5 14.5 2.5]*4; %small reward arm
    smArmPgY = [12.5 12.5 8.5 8.5 14.5 14.5]*4;
    smAPg = polyshape(smArmPgX, smArmPgY);

    ph = 1; %when there are differences between the rewards

    cMap = nebula(100);
    cMapInds = linspace(0, 100);

    tInds = [4 NaN 5 3 2 3 1 10 8 2; 5 NaN 8 7 4 4 8 5 1 5]; %to fill in with which example I want to use for each mouse
    for r = 1:2
        for m = 1:length(mouse)
            if isempty(pathsEx{m,r,ph})
                continue
            end
            subplot(2,10,(r-1)*10+m)

            hold on
            plot(lrPg, 'FaceColor', rgb('LightGray'), 'EdgeColor', 'None')
            plot(smPg, 'FaceColor', rgb('LightGray'), 'EdgeColor', 'None')
            plot(dPg, 'FaceColor', rgb('LightGray'), 'EdgeColor', 'None')
            plot(stPg, 'FaceColor', rgb('LightGray'), 'EdgeColor', 'None')
            plot(lrAPg, 'FaceColor', rgb('LightGray'), 'EdgeColor', 'None')
            plot(smAPg, 'FaceColor', rgb('LightGray'), 'EdgeColor', 'None')


            tmpT = tInds(r,m);
            % tmpT = 5;
            tmpCoords = pathsEx{m,r,ph}{tmpT};
            runSpeed = get_runspeed(tmpCoords);
            smRs = smooth_runspeed(runSpeed, 500/1000); %smooth runspeed w 500 ms window
            % plot(tmpCoords(:,2), tmpCoords(:,3));
            % patch([tmpCoords(:,2); NaN], [tmpCoords(:,3); NaN], [cMap(match(smRs(:,2), cMapInds)); NaN], 'EdgeColor','interp');
            % patch([tmpCoords(:,2); NaN], [tmpCoords(:,3); NaN], [match(smRs(:,2), cMapInds); NaN], 'EdgeColor','interp', 'Marker', '.');
            % patch([tmpCoords(:,2); NaN], [tmpCoords(:,3); NaN], [match(smRs(:,2), cMapInds); NaN], 'EdgeColor','interp', 'LineWidth', 1);

            % hold on
            for pt = 1:length(tmpCoords)-1
                % tmp = match(smRs(pt,2), cMapInds)
                plot([tmpCoords(pt,2) tmpCoords(pt+1,2)], [tmpCoords(pt,3) tmpCoords(pt+1,3)], 'Color', cMap(match(smRs(pt,2), cMapInds),:))
            end %pts

            % colormap(nebula)
            % caxis([0 100])
            xlim([0.5 56.5])
            xticks([])
            yticks([])
            ylim([0.5 56.5])
            axis square
            axis off
            title([num2str(round(tmpCoords(end,1)-tmpCoords(1,1))) 's'])


        end %mouse
    end %reward


    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save

    figtitle = ['ExamplePaths_' tNames{d} '_postShift'];
    figure('Name', figtitle)

    ph = 2; %when there are differences between the rewards

    cMap = nebula(100);
    cMapInds = linspace(0, 100);

    tInds = [4 NaN 5 3 2 3 1 10 8 2; 5 NaN 8 7 4 4 8 5 1 5]; %to fill in with which example I want to use for each mouse
    for r = 1:2
        for m = 1:length(mouse)


            hold on
            plot(lrPg, 'FaceColor', rgb('LightGray'), 'EdgeColor', 'None')
            plot(smPg, 'FaceColor', rgb('LightGray'), 'EdgeColor', 'None')
            plot(dPg, 'FaceColor', rgb('LightGray'), 'EdgeColor', 'None')
            plot(stPg, 'FaceColor', rgb('LightGray'), 'EdgeColor', 'None')
            plot(lrAPg, 'FaceColor', rgb('LightGray'), 'EdgeColor', 'None')
            plot(smAPg, 'FaceColor', rgb('LightGray'), 'EdgeColor', 'None')

            if isempty(pathsEx{m,r,ph})
                continue
            end
            subplot(2,10,(r-1)*10+m)
            tmpT = tInds(r,m);
            % tmpT = 5;
            tmpCoords = pathsEx{m,r,ph}{tmpT};
            runSpeed = get_runspeed(tmpCoords);
            smRs = smooth_runspeed(runSpeed, 500/1000); %smooth runspeed w 500 ms window
            % plot(tmpCoords(:,2), tmpCoords(:,3));
            % patch([tmpCoords(:,2); NaN], [tmpCoords(:,3); NaN], [cMap(match(smRs(:,2), cMapInds)); NaN], 'EdgeColor','interp');
            % patch([tmpCoords(:,2); NaN], [tmpCoords(:,3); NaN], [match(smRs(:,2), cMapInds); NaN], 'EdgeColor','interp', 'Marker', '.');
            % patch([tmpCoords(:,2); NaN], [tmpCoords(:,3); NaN], [match(smRs(:,2), cMapInds); NaN], 'EdgeColor','interp', 'LineWidth', 1);

            hold on
            for pt = 1:length(tmpCoords)-1
                % tmp = match(smRs(pt,2), cMapInds)
                plot([tmpCoords(pt,2) tmpCoords(pt+1,2)], [tmpCoords(pt,3) tmpCoords(pt+1,3)], 'Color', cMap(match(smRs(pt,2), cMapInds),:))
            end %pts

            % colormap(nebula)
            % caxis([0 100])
            % caxis([0 100])
            xlim([0.5 56.5])
            xticks([])
            yticks([])
            ylim([0.5 56.5])
            axis square
            axis off
            title([num2str(round(tmpCoords(end,1)-tmpCoords(1,1))) 's'])

        end %mouse
    end %reward


    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save



end %shift day

end %function