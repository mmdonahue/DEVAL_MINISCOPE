function deval_x_normRate_reward(mouse)
% function deval_x_plot_rewardRepresentation(mouse)
%
% PURPOSE:
%   Plot heatmaps showing the firing rates of place cells across the figure
%   8 maze to determine if rewards are over-represented.
%
% MMD
% 05/2025
% Ewell Lab

%% OPTIONS

saveOrNot = 0;

% days2use = [2 4 6 8 9]; %days I'm interested in; set 1:9 for all
days2use = 4; %days I'm interested in; set 1:9 for all

adpOnly = 1; %set to 1 to exclude mice that did not adapt

minDis = 10; %minimum number of "spikes" in a day for a unit to be included

pThresh = 1;

xLim = [-40 10]; %what to set x-axis as

%% INITIALIZE

rewNames = {'Large reward', 'Small reward'};

mapIndsRew = [1:42; 1:14 43:70];
mapIndsRewPlot = cat(1,7:37,[7:14 43:65]);

binAx = ((1:42) - 35)*4; %4 = binSz

ph2use = [2 2 2 3 2 2 2 2 2]; %forced trials always, downshifted forced on shift day
freePh = [1 1 1 2 1 1 1 1 1]; %forced trials always, downshifted forced on shift day

s = RandStream('mt19937ar');

nboot = 1000; %for confidence intervals

rewCols = {'Blue', 'Red'};

tNames = {'Preference1', 'Preference2', 'preshift', 'Downshift', 'postshift', 'Adaptation', 'Reversal', 'Extinction1', 'Extinction2'};

saveDir = 'E:\Devaluation Project Miniscope\Meg analysis\RESULTS\normRate_reward';
cd(saveDir)

%% GET DATA


for d = days2use

    rmsUntracked = [];
    rmsDs = []; %for downsampling
    rmsByMouse = cell(length(mouse),1);
    rmsByMouseByPoke = cell(length(mouse),2,2,length(mouse(end).day(d).phase)); %mouse x reward x poke ind x tt
    rmsByMouseDs = cell(length(mouse),1);
    % rmsByMouseforNico = cell(length(mouse),1);
    statAll = [];
    statAtRew = [];

    if d == 4
        trackLRTs = [];
        trackLRTsMouse = cell(length(mouse),1);
        statAtRewNoLRT = [];
        load("E:\Devaluation Project Miniscope\Meg analysis\RESULTS\LRTcellInfo.mat")
    end %for shift day only for paper

    for m = 1:10
        if adpOnly == 1 && mouse(m).adp == 0
            continue
        end %ignore non-adapted mouse

        if isempty(mouse(m).day(d).dir)
            continue
        end

        load([mouse(m).day(d).dir '\byTrialRatemaps_v3.mat'], 'unit')

        mouseMaps  = [];

        % for down-sampling free choice trials so that they are equal between rewards
        rews = vertcat(mouse(m).day(d).trial(:).rewLoc);
        tmpFunc = @(x)(strncmp(x, mouse(m).largeRew, 4));
        rews = double(cellfun(tmpFunc, rews)); %reward inds
        freeRews = rews(mouse(m).day(d).phase(freePh(d)).trials);
        frSm = length(freeRews) - sum(freeRews); %number of times mouse went to small reward during the free choice

        if frSm < length(freeRews) / 2
            dsTs = zeros(2,frSm);
            dsTs(2,:) = find(freeRews==0)';
            dsTs(1,:) = sort(randsample(s, find(freeRews==1), frSm))';
        elseif frSm > length(freeRews) / 2
            dsTs = zeros(2,length(freeRews)-frSm);
            dsTs(1,:) = find(freeRews==1)';
            dsTs(2,:) = sort(randsample(s, find(freeRews==0), length(freeRews)-frSm))';
        elseif frSm == length(freeRews) / 2
            dsTs = zeros(2,frSm);
            dsTs(1,:) = find(freeRews==1);
            dsTs(2,:) = find(freeRews==0);
        end

        for u = 1:length(unit)
             if unit(u).numDis < minDis
                continue
            end %cell isn't spatially stable within session

          if d == 4
                trackLRTs = [trackLRTs; ismember([m unit(u).ID], LRTinfo, 'rows')];
                trackLRTsMouse{m} = [trackLRTsMouse{m}; ismember([m unit(u).ID], LRTinfo, 'rows')];
            end %shift day - remove LRTs

            uMaps = NaN(1,42,2,length(mouse(m).day(d).phase));
            uMapsDs = NaN(1,42,2,length(mouse(m).day(d).phase));
            uMapsByPoke = NaN(2,42,3,length(mouse(m).day(d).phase));

            maxFr = max([max(unit(u).trajectory(1).rm) max(unit(u).trajectory(2).rm)]);

            rews = vertcat(mouse(m).day(d).trial(:).rewLoc);
            tmpFunc = @(x)(strncmp(x, mouse(m).largeRew, 4));
            rews = double(cellfun(tmpFunc, rews)); %reward inds
            rews(rews==0) = 2; %small reward = 2

            if strncmp(mouse(m).largeRew, 'East', 4)
                mapInds = mapIndsRew;
            else
                mapInds = flip(mapIndsRew);
            end

            for tt = 1:length(mouse(m).day(d).phase)
                % tmpMap = unit(u).trialType(tt).smRm;
                tmpMaps = vertcat(unit(u).trial(mouse(m).day(d).phase(tt).trials).linSmSpkRm);
                tmpRews = rews(mouse(m).day(d).phase(tt).trials);
                tmpPokes = vertcat(mouse(m).day(d).trial(:).tmInds);
                tmpPokes = ~isnan(tmpPokes(:,3));
                tmpPokes = tmpPokes(mouse(m).day(d).phase(tt).trials);
                % maxFr = 0; %initialize

                normMaps = zeros(1,42,2);
                normMapsByPoke = zeros(2,42,2); %poke inds x bins x reward

                dsMaps = zeros(1,42,2);
                for r = 1:2

                    % normMaps(:,:,r) = unit(u).trialType(tt).trajectory(r).rm;
                    % maxFr = max([maxFr max(unit(u).trialType(tt).trajectory(r).rm)]);

                    normMaps(:,:,r) = nanmean(tmpMaps(tmpRews==r,mapInds(r,:)),1);
                    % maxFr = max([maxFr max(normMaps(:,:,r))]);
                    if ~startsWith(mouse(m).day(d).phase(tt).name, 'force')
                        dsMaps(:,:,r) = nanmean(tmpMaps(dsTs(r,:),mapInds(r,:)),1);
                    else
                        dsMaps(:,:,r) = nanmean(tmpMaps(tmpRews==r,mapInds(r,:)),1);
                    end

                    goodPokes = double(tmpRews==r) + double(tmpPokes==1) == 2;
                    badPokes = double(tmpRews==r) + double(tmpPokes==0) == 2;
                    normMapsByPoke(1,:,r) = nanmean(tmpMaps(goodPokes,mapInds(r,:)),1);
                    normMapsByPoke(2,:,r) = nanmean(tmpMaps(badPokes,mapInds(r,:)),1);
                end %reward

                normMaps = normMaps ./ maxFr;
                dsMaps = dsMaps ./ maxFr;

                uMaps(1,:,:,tt) = normMaps;
                uMapsDs(1,:,:,tt) = dsMaps;
                for pi = 1:2
                    uMapsByPoke(pi,:,r,tt) = normMapsByPoke(pi,:,r) ./ maxFr;
                end %poke ind
                % uMapsNico(1,:,r,tt) = normMap(mapIndsRewPlot(r,:));

                % tmpStat = task/day mouse unit reward bin data
                for r = 1:2
                    try
                        tmpData = normMaps(1,25:37,r)';
                        tmpStat = [repmat(d,length(tmpData),1) repmat(m,length(tmpData),1), repmat(u,length(tmpData),1)...
                            repmat(tt,length(tmpData),1), repmat(r,length(tmpData),1), transpose(25:37), tmpData];
                    catch; keyboard; end
                    statAll = [statAll; tmpStat];

                    tmpData = nanmean(normMaps(1,34:36,r));
                    tmpStat = [d m u tt r tmpData];
                    statAtRew = [statAtRew; tmpStat];

                    if d == 4 && ~ismember([m unit(u).ID], LRTinfo, 'rows')
                        statAtRewNoLRT = [statAtRewNoLRT; tmpStat];
                    end %nonLRT
                end %reward for stats

                % normMap = tmpMap ./ maxFr;
                % uMaps = cat(4, uMaps, normMap);
            end %trial type

            rmsUntracked = cat(1, rmsUntracked, uMaps);

            rmsDs = cat(1, rmsDs, uMapsDs);
            rmsByMouse{m} = cat(1, rmsByMouse{m}, uMaps);
           

            % keyboard
            for r = 1:2
                for pi = 1:2
                    for tt = 1:size(uMapsByPoke,4)
                        try
                        rmsByMouseByPoke{m,r,pi,tt} = cat(1, rmsByMouseByPoke{m,r,pi,tt}, uMapsByPoke(pi,:,r,tt));
                        catch; keyboard; end
                    end %tt
                end %poke ind
            end %reward

        end %unit
    end %mouse
    statAll = statAll(~isnan(statAll(:,end)),:);

    statAtRew = statAtRew(~isnan(statAtRew(:,end)),:);
    keyboard
    if d ~=4
        forStatAtRew = statAtRew(statAtRew(:,4)==2,:);
    else
        forStatAtRew = statAtRew(statAtRew(:,4)~=2,:);
    end
    forStatAtRew(:,7) = forStatAtRew(:,6) + 0.0001; %zero pad
    tbl = array2table(forStatAtRew, 'VariableNames', {'day', 'mouse', 'unit', 'trialType', 'reward', 'rate', 'zp_rate'});
    tbl.day = tNames(tbl.day)';
    %         % tbl.day = categorical(tbl.day);
    %         % tbl.day = num2str(tbl.day);
    ttNames = vertcat({mouse(m).day(d).phase(:).name});
    tbl.trialType = ttNames(tbl.trialType)';
    %         % tbl.bin = categorical(tbl.bin);
    % tbl.mouse = categorical(tbl.mouse);
    % tbl.unit = categorical(tbl.unit);
    tbl.mouse = num2str(tbl.mouse);
    tbl.unit = num2str(tbl.unit);
    % tbl.reward = categorical(tbl.reward);
    tbl.reward = rewNames(tbl.reward)';
    if d ~=4
        glme = fitglme(tbl, 'zp_rate ~  reward + (1|mouse) + (1|mouse:unit)', 'DummyVarCoding', 'effects', 'Distribution', 'Gamma'); %reward with random effect of mouse and unit nested in mouse
        emm = emmeans(glme, 'effects');
        H0 = contrasts_wald(glme, emm, [1 -1]);
        fprintf('forced large vs. forced small reward - p = %s\n', num2str(round(H0.pVal,4)))

    else
        glme = fitglme(tbl, 'zp_rate ~ reward*trialType + (1|mouse) + (1|mouse:unit)', 'DummyVarCoding', 'effects', 'Distribution', 'Gamma'); %reward with random effect of mouse and unit nested in mouse
        emm = emmeans(glme, 'effects');
        L = [1 -1 0 0; 1 0 -1 0; 0 0 1 -1; 0 1 0 -1];
        L = [1 -1 0 0];
        H0 = contrasts_wald(glme, emm, L);
        fprintf('forced pre large vs. forced pre small reward - p = %s\n', num2str(round(H0.pVal * 4,4))) %adjust
        L = [0 0 1 -1];
        H0 = contrasts_wald(glme, emm, L);
        fprintf('forced post large vs. forced post small reward - p = %s\n', num2str(round(H0.pVal * 4,4))) %adjust
        L = [1 0 -1 0];
        H0 = contrasts_wald(glme, emm, L);
        fprintf('forced pre large vs. forced post large reward - p = %s\n', num2str(round(H0.pVal * 4,4))) %adjust
        L = [0 1 0 -1];
        H0 = contrasts_wald(glme, emm, L);
        fprintf('forced pre small vs. forced post small - p = %s\n', num2str(round(H0.pVal * 4,4))) %adjust
    end %which day affects the stats

    figtitle = ['Activity_acrossMaze_' tNames{d}];
    figure('Name', figtitle, 'Position', [1 37 748 841])

    rmsByMouseTrim = rmsByMouse(~cellfun(@isempty, rmsByMouse));
    mouseNames = vertcat({mouse(:).name});
    mouseNames = mouseNames(~cellfun(@isempty, rmsByMouse));
    tmpMeans = cellfun(@nanmean,rmsByMouseTrim, 'UniformOutput', false);

    % for tt = 1:size(rmsUntracked,4)
    %     tmpFunc = @(x)(squeeze(mean(x(:,33:36,:,tt),2)));
    %     NicoData = cellfun(tmpFunc, rmsByMouseTrim, 'UniformOutput', false);
    %     for m = 1:length(NicoData)
    %         try
    %         writematrix(NicoData{m}, ['dataForNico_' tNames{d} '_' mouse(end).day(d).phase(tt).name '.xlsx'], 'Sheet', mouseNames{m})
    %         catch; keyboard; end
    %     end %mouse
    % end %tt
 
    for tt = 1:size(rmsUntracked,4)
        % subplot(3,1,tt)
        subplot(3,5,1+((tt-1)*5):4+((tt-1)*5))
        lh = zeros(1,2);
        hold on

        xBox = [-40 -40 -5 -5];
        pgon = polyshape(xBox,[0 0.4 0.4 0]);
        plot(pgon, 'FaceColor', rgb('Yellow'), 'LineStyle', 'None');

        xBox = [-5 -5 5 5];
        pgon = polyshape(xBox,[0 0.4 0.4 0]);
        plot(pgon, 'FaceColor', rgb('Orange'), 'LineStyle', 'None');

        for r = 1:2
            pullMaps = rmsUntracked(:,:,r,tt);
            meanData = nanmean(pullMaps);
            CI = bootci(nboot, {@nanmean, pullMaps}, 'type', 'per');

            hold on;

            lh(r) =  plot_filled_ci(binAx, meanData, CI, rgb(rewCols{r}));
        end %reward type

        ylabel('Normalized Ca^2^+ event rate')
        ylim([0 0.4])

        xlim(xLim)
        xlabel({[], 'Distance to reward (cm)'})
  
        title([mouse(m).day(d).phase(tt).name ' choice trials (n = ' num2str(size(rmsUntracked,1)) ' cells)'])
        leg = legend(lh, rewNames, 'Location', 'northeastoutside', 'AutoUpdate','off');

        text(-30, -0.045, 'Reward approach', 'Color', rgb('GoldenRod'))
        text(-2.5, -0.045, 'Reward', 'Color', rgb('DarkOrange'))

        % [beta,betaNames,statsFixed] = fixedEffects(glme)
        subplot(3,5,tt*5)
        tmpFunc = @(x)reshape((mean(x(1,33:36,:,tt))),1,2);
        pullData = cellfun(tmpFunc, tmpMeans, 'UniformOutput', false);
        pullData = vertcat(pullData{:});

        for r = 1:2
            hold on
            scatter(ones(1,length(pullData))*r, pullData(:,r), 'MarkerFaceColor', rgb(rewCols{r}), 'MarkerEdgeColor', rgb('White'))
        end %reward

        for m = 1:size(pullData,1)
            line([1 2], pullData(m,:), 'Color', 'k')
        end %mouse

        xlim([0.5 2.5])
        xticks(1:2)
        xticklabels(rewNames)
        ylim([0 0.4])
        ylabel('Normalized Ca^2^+ event rate')

        title( [mouse(m).day(d).phase(tt).name ' choice trials (n = ' num2str(size(pullData,1)) ' mice)'])
    end %task


end %day



end %function
