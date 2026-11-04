function deval_x_2dPf_distToRew(mouse)
% function deval_x_2dPf_distToRew(mouse)
%
% PURPOSE:
%   Use 2D place fields to determine if rewards are overrepresented across
%   days in the deval task.
%
% MMD
% Ewell Lab
% 07/2025

%% OPTIONS

saveOrNot = 0;

% days2use = [2 4 6 8]; %days I'm interested in; set 1:9 for all
days2use = 4; %days I'm interested in; set 1:9 for all

ph2use = [2 2 2 3 2 2 2 2 2]; %forced trials always, downshifted forced on shift day

adpOnly = 1; %set to 1 to exclude mice that did not adapt


%% INITIALIZE

minDis = 10; %minimum number of "spikes" in a day for a unit to be included
pThresh = 0.05;

binSz = 4;
binAx = 2:binSz:56;

invBins = flip(1:14);

tNames = {'Preference1', 'Preference2', 'preshift', 'Downshift', 'postshift', 'Adaptation', 'Reversal', 'Extinction1', 'Extinction2'};
rewNames = {'Large reward', 'Small reward', 'Both', 'Neither'};

rewCols = {'SteelBlue', 'Crimson', 'Violet', 'Khaki'}; %large small both neither

rewLocs = [4 6; 4 50]; %east and west
% rewLocInds(:,2) = [12 1; 13 1; 14 1; 12 2; 13 2; 14 2]
rewLocInds(1,:) = sub2ind([14 14], [1 2 3 1 2 3], [1 1 1 2 2 2]);
rewLocInds(2,:) = sub2ind([14 14], [12 13 14 12 13 14], [1 1 1 2 2 2]);


edges = 0:2.5:20;
% edgeAx = edges_to_x_vals((edges - rewInd)*4); %binSz = 4

xLim = [-40 12];

percThresh = 25;
minPfArea = 10; %cm^2

colMap = [rgb('White'); nebula(50)];
% rewCols = {'Blue', 'Red'};

nboot = 1000;

saveDir = 'E:\Devaluation Project Miniscope\Meg analysis\RESULTS\percentage_pf_reward_2d';
cd(saveDir)

%% GET DATA AND PLOT


for d = days2use
    rewCntr5cm = zeros(length(mouse(end).day(d).phase),4);
    rewCntr10cm = zeros(length(mouse(end).day(d).phase),4);

    rewCntr5cmMouse = cell(10,1);
    rewCntr5cmMouse(:) = {zeros(length(mouse(end).day(d).phase),4)};
    rewCntr10cmMouse = cell(10,1);
    rewCntr10cmMouse(:) = {zeros(length(mouse(end).day(d).phase),4)};

    dist2rew = cell(2,3);
    pkDist2rew = cell(2,3);

    allPkLocs = cell(1,3);
    allPkLocs(:) = {zeros(14,14)};

    allPkLocsMouse = cell(1,3);
    allPkLocsMouse(:) = {zeros(14,14,10)};

    distxRate = cell(2,length(mouse(end).day(d).phase));
    distxNormRate = cell(2,length(mouse(end).day(d).phase));
    distxPkRate = cell(2,length(mouse(end).day(d).phase));

    distxSize = cell(2,length(mouse(end).day(d).phase));

    oofFr = cell(2,length(mouse(end).day(d).phase));
    oofFrNorm = cell(2,length(mouse(end).day(d).phase));

    figtitle1 = ['%pf_5cmReward_' tNames{d}];
    f1 = figure('Name', figtitle1);
    figtitle2 = ['%pf_10cmReward_' tNames{d}];
    f2 = figure('Name', figtitle2);

    figtitle3 = ['%pf_5cmReward_' tNames{d} '_byMouse'];
    f3 = figure('Name', figtitle3);
    figtitle4 = ['%pf_10cmReward_' tNames{d} '_byMouse'];
    f4 = figure('Name', figtitle4);
    % figure
    for m = 1:10

        if adpOnly == 1 && mouse(m).adp == 0
            continue
        end %ignore non-adapted mouse

        if isempty(mouse(m).day(d).dir)
            continue
        end %no data for this day for this mouse

        load([mouse(m).day(d).dir '\byTrialRatemaps_v3.mat'], 'unit')

        if strcmp(mouse(m).largeRew, 'West')
            mRewLocs = flip(rewLocs);
            mRewLocInds = flip(rewLocInds);
        else
            mRewLocs = rewLocs;
            mRewLocInds = rewLocInds;
        end %assign for small and large reward

        for tt = 1:length(mouse(m).day(d).phase)
            % allP = [];
            for u = 1:length(unit)
                % allP = [allP unit(u).trialType(tt).pVal];
                % if unit(u).numDis < minDis || unit(u).trialType(ph2use(d)).pVal > pThresh || isnan(unit(u).trialType(ph2use(d)).meanCorr)
                if unit(u).numDis < minDis || min(vertcat(unit(u).trialType(1:length(mouse(m).day(d).phase)).pVal)) >= pThresh...
                        || sum(isnan(vertcat(unit(u).trialType(1:length(mouse(m).day(d).phase)).meanCorr))) == length(mouse(m).day(d).phase)
                    continue
                end %cell isn't spatially stable within session

                ttMap = nanmean(cat(3, unit(u).trial(mouse(m).day(d).phase(tt).trials).smSpkRm),3);
                % pf = get_2d_pfs_v2(ttMap, binSz, percThresh, 0.1*max(max(ttMap)), minPfArea);
                pf = unit(u).pf2d;

                rewPk = zeros(2,2); %initialize as none

                tmpMax = max(max(max(cat(3, unit(u).trialType(1:length(mouse(m).day(d).phase)).smSpkRm)))); %FOR NORMALIZING PURPOSES

                %check if it reached min on this tt
                for p = 1:length(pf)
                    inds = sub2ind(size(ttMap), pf(p).inds(:,2), pf(p).inds(:,1));
                    if max(ttMap(inds)) > 0.2*max(max(unit(u).smSpkRm))

                        %% PLOT CHECK
                        % figure
                        % subplot(1,2,1)
                        % imagesc(unit(u).smSpkRm)
                        % axis xy
                        % axis square
                        % hold on
                        % x = pf(p).inds(:,1);
                        % y = pf(p).inds(:,2);
                        % k = boundary(pf(p).inds(:,1), pf(p).inds(:,2));
                        % if ~isempty(k)
                        %     plot(x(k), y(k), 'Color', 'k', 'LineWidth', 2);
                        % elseif mean(diff(y)) == 0
                        %     line([min(x) max(x)], [y(1) y(1)], 'Color', 'k', 'LineWidth', 2)
                        % elseif mean(diff(x)) == 0
                        %     line([x(1) x(1)], [min(y) max(y)], 'Color', 'k', 'LineWidth', 2)
                        % end
                        % title('Day-averaged map')
                        % colorbar
                        %
                        % subplot(1,2,2)
                        % imagesc(ttMap)
                        % axis xy
                        % axis square
                        % hold on
                        % x = pf(p).inds(:,1);
                        % y = pf(p).inds(:,2);
                        % k = boundary(pf(p).inds(:,1), pf(p).inds(:,2));
                        % if ~isempty(k)
                        %     plot(x(k), y(k), 'Color', 'k', 'LineWidth', 2);
                        % elseif mean(diff(y)) == 0
                        %     line([min(x) max(x)], [y(1) y(1)], 'Color', 'k', 'LineWidth', 2)
                        % elseif mean(diff(x)) == 0
                        %     line([x(1) x(1)], [min(y) max(y)], 'Color', 'k', 'LineWidth', 2)
                        % end
                        %   title([mouse(m).day(d).phase(tt).name ' map'])
                        % colorbar
                        %% END PLOT CHECK


                        [tmpPkRate, tmpPk] = max(ttMap(inds));
                        percFr = unit(u).smSpkRm(inds) ./ max(max(unit(u).smSpkRm));
                        pkBins = [pf(p).inds(tmpPk,1) pf(p).inds(tmpPk,2)]; %NOTE: THIS FLIPS THE ORIENTATION TO MAKE NICO'S FIGURE 2
                        pkLoc = [binAx(pf(p).inds(tmpPk,1)) binAx(pf(p).inds(tmpPk,2))];

                        if strcmp(mouse(m).largeRew, 'West')
                            allPkLocs{tt}(pkBins(1),invBins(pkBins(2))) = allPkLocs{tt}(pkBins(1),invBins(pkBins(2))) + 1;
                            allPkLocsMouse{tt}(pkBins(1),invBins(pkBins(2)),m) = allPkLocsMouse{tt}(pkBins(1),invBins(pkBins(2)),m) + 1;
                        else
                            allPkLocs{tt}(pkBins(1),pkBins(2)) = allPkLocs{tt}(pkBins(1),pkBins(2)) + 1;
                            allPkLocsMouse{tt}(pkBins(1),pkBins(2),m) = allPkLocsMouse{tt}(pkBins(1),pkBins(2),m) + 1;
                        end

                        for r = 1:2
                            tmpPkDist = pdist2(pkLoc, mRewLocs(r,:));
                            tmpAllDist = pdist2([binAx(pf(p).inds(:,1))' binAx(pf(p).inds(:,2))'], mRewLocs(r,:));
                            [~, sortOrd] = sort(tmpAllDist);
                            tmpMinDist = []; %initialize as empty
                            i = 1;
                            while isempty(tmpMinDist)
                                checkFr = ttMap(inds(sortOrd(i)));
                                if checkFr >= max(max(unit(u).smSpkRm))*0.2
                                    tmpMinDist = tmpAllDist(sortOrd(i));
                                else
                                    i = i + 1;
                                end %if firing rate in this bin is at least 20% of the max firing rate from the day-averaged map
                            end %while we haven't found the min dist that meet firing rate min in this tt's ratemap

                            % keyboard

                            % %using peak dist
                            pkDist2rew{r,tt} = [pkDist2rew{r,tt}; tmpPkDist];
                            %  distxRate{r,tt} = [distxRate{r,tt}; tmpPkDist nanmean(ttMap(inds))*60]; %convert to min-1
                            %  distxPkRate{r,tt} = [distxPkRate{r,tt}; tmpPkDist tmpPkRate*60]; %convert to min-1
                            %  distxNormRate{r,tt} = [distxNormRate{r,tt}; tmpPkDist tmpPkRate ./ tmpMax]; %normalize
                            %  rewPk(1,r) = rewPk(1,r) + double(tmpPkDist <= 5);
                            %  rewPk(2,r) = rewPk(2,r) + double(tmpPkDist <= 10);

                            %using min dist
                            dist2rew{r,tt} = [dist2rew{r,tt}; tmpMinDist];
                            distxRate{r,tt} = [distxRate{r,tt}; tmpMinDist nanmean(ttMap(inds))*60]; %convert to min-1
                            distxPkRate{r,tt} = [distxPkRate{r,tt}; tmpMinDist tmpPkRate*60]; %convert to min-1
                            distxNormRate{r,tt} = [distxNormRate{r,tt}; tmpMinDist tmpPkRate ./ tmpMax]; %normalize
                            distxSize{r,tt} = [distxSize{r,tt}; tmpMinDist pf(p).area];
                            rewPk(1,r) = rewPk(1,r) + double(tmpMinDist <= 5);
                            rewPk(2,r) = rewPk(2,r) + double(tmpMinDist <= 10);
                        end %reward
                    end %if it reaches min fr
                end %place field
                rewPk = double(logical(rewPk));

                if sum(rewPk(1,:)) == 0
                    rewCntr5cm(tt,4) = rewCntr5cm(tt,4) + 1;
                    rewCntr5cmMouse{m}(tt,4) = rewCntr5cmMouse{m}(tt,4) + 1;
                elseif sum(rewPk(1,:)) == 1
                    rewCntr5cm(tt,find(rewPk(1,:))) = rewCntr5cm(tt,find(rewPk(1,:))) + 1;
                    rewCntr5cmMouse{m}(tt,find(rewPk(1,:))) = rewCntr5cmMouse{m}(tt,find(rewPk(1,:))) + 1;
                elseif sum(rewPk(1,:)) == 2
                    rewCntr5cm(tt,3) = rewCntr5cm(tt,3) + 1;
                    rewCntr5cmMouse{m}(tt,3) = rewCntr5cmMouse{m}(tt,3) + 1;
                end %what to add

                if sum(rewPk(2,:)) == 0
                    rewCntr10cm(tt,4) = rewCntr10cm(tt,4) + 1;
                    rewCntr10cmMouse{m}(tt,4) = rewCntr10cmMouse{m}(tt,4) + 1;
                elseif sum(rewPk(2,:)) == 1
                    rewCntr10cm(tt,find(rewPk(2,:))) = rewCntr10cm(tt,find(rewPk(2,:))) + 1;
                    rewCntr10cmMouse{m}(tt,find(rewPk(2,:))) = rewCntr10cmMouse{m}(tt,find(rewPk(2,:))) + 1;
                elseif sum(rewPk(2,:)) == 2
                    rewCntr10cm(tt,3) = rewCntr10cm(tt,3) + 1;
                    rewCntr10cmMouse{m}(tt,3) = rewCntr10cmMouse{m}(tt,3) + 1;
                end %what to add

            end %unit
        end %tt

    end %mouse
    keyboard
    figure(f1)
    for tt = 1:size(rewCntr5cm,1)
        subplot(1,3,tt)
        p = pie(rewCntr5cm(tt,:));

        patchHand = findobj(p, 'Type', 'Patch');
        for r = 1:4
            patchHand(r).FaceColor = rgb(rewCols{r}) ;
        end

        lgd = legend(rewNames, 'Location', 'southoutside');
        title({tNames{d}; mouse(m).day(d).phase(tt).name})
    end %tt
    %chi squared test for forced choice pre
    tt = 1;
    n1 = rewCntr5cm(tt,1); N1 = sum(rewCntr5cm(tt,:)); %large
    n2 = rewCntr5cm(tt,2); N2 = sum(rewCntr5cm(tt,:)); %small
    x1 = [repmat('a',N1,1); repmat('b',N2,1)];
    x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
    [tbl,chi2stat,pval] = crosstab(x1,x2);


    if saveOrNot == 1
        saveas(f1, figtitle1, 'fig');
        saveas(f1, figtitle1, 'png');
        saveas(f1, figtitle1, 'epsc');
    end %save

    figure(f2)
    for tt = 1:size(rewCntr10cm,1)
        subplot(1,3,tt)
        p = pie(rewCntr10cm(tt,:));

        patchHand = findobj(p, 'Type', 'Patch');
        for r = 1:4
            patchHand(r).FaceColor = rgb(rewCols{r}) ;
        end

        lgd = legend(rewNames, 'Location', 'southoutside');
        title({tNames{d}; mouse(m).day(d).phase(tt).name})
    end %tt

    if saveOrNot == 1
        saveas(f2,  figtitle2, 'fig');
        saveas(f2, figtitle2, 'png');
        saveas(f2, figtitle2, 'epsc');
    end %save

    figure(f3)
    for m = 1:10
        tmpData = rewCntr5cmMouse{m}(ph2use(d),:);
        if sum(tmpData) == 0
            continue
        end %no data

        subplot(2,5,m)
        p = pie(tmpData);

        patchHand = findobj(p, 'Type', 'Patch');
        for r = 1:4
            patchHand(r).FaceColor = rgb(rewCols{r}) ;
        end
        title([mouse(m).name ' (n = ' num2str(sum(tmpData)) ')'])
    end %mouse

    if saveOrNot == 1
        cd('byMouse')
        saveas(f3, figtitle3, 'fig');
        saveas(f3, figtitle3, 'png');
        saveas(f3, figtitle3, 'epsc');
        cd ../
    end %save

    figure(f4)
    for m = 1:10
        tmpData = rewCntr10cmMouse{m}(ph2use(d),:);
        if sum(tmpData) == 0
            continue
        end %no data

        subplot(2,5,m)
        p = pie(tmpData);

        patchHand = findobj(p, 'Type', 'Patch');
        for r = 1:4
            patchHand(r).FaceColor = rgb(rewCols{r}) ;
        end
        title([mouse(m).name ' (n = ' num2str(sum(tmpData)) ')'])
    end %mouse

    if saveOrNot == 1
        cd('byMouse')
        saveas(f4, figtitle4, 'fig');
        saveas(f4, figtitle4, 'png');
        saveas(f4, figtitle4, 'epsc');
        cd ../
    end %save

    close all
    keyboard

    figtitle6 = ['pfPkLocs' tNames{d}];
    f6 = figure('Name', figtitle6);
    tmpMax = 0;
    for tt = 1:length(mouse(m).day(d).phase)
        subplot(1,3,tt)

        tmpMap = allPkLocs{tt};
        tmpMap = (tmpMap ./ sum(tmpMap(:))) * 100;
        tmpMap(isnan(tmpMap)) = -0.1;
        imagesc(tmpMap)
        axis xy
        axis square
        yticks('')
        xticks([2 13])
        xticklabels(rewNames)
        colormap(colMap)
        caxis([0 4.5])
        title(mouse(m).day(d).phase(tt).name)

        tmpMax = max([tmpMax max(tmpMap(:))]);
        cbr = colorbar('southoutside');
        ylabel(cbr, '% of place fields')
    end %phase

    if saveOrNot == 1
        saveas(f6, figtitle6, 'fig');
        saveas(f6, figtitle6, 'png');
        saveas(f6, figtitle6, 'epsc');
    end %save

    figtitle8 = ['pfDistxRate_' tNames{d}];
    f8 = figure('Name', figtitle8, 'Position', [1 1 669 741]);
    spMap = [1 3 5];

    tmpFunc = @(x)(x(x(:,1)<=5,2));
    pullData = cellfun(tmpFunc, distxRate, 'UniformOutput', false);
    tmpMean = cellfun(@mean, pullData);
    tmpCIFunc =  @(x)(bootci(nboot, {@mean, x}, 'type', 'per'));
    tmpCI = cellfun(tmpCIFunc, pullData, 'UniformOutput', false);
    for tt = 1:length(mouse(m).day(d).phase)
        subplot(3,2,spMap(tt))

        pullCI = horzcat(tmpCI{:,tt});
        bgraph = bar(tmpMean(:,tt), 'FaceColor', 'Flat');
        bgraph.CData = [rgb(rewCols{1}); rgb(rewCols{2})];
        hold on
        errorbar(1:2, tmpMean(:,tt), tmpMean(:,tt)' - pullCI(1,:),  pullCI(end,:)-tmpMean(:,tt)', 'Color', 'k', 'LineStyle', 'None')

        h = dotplot(1:2, pullData(:,tt), 0.1, [rgb(rewCols{1}); rgb(rewCols{2})], [rgb('White'); rgb('White')]);
        xticklabels(rewNames(1:2))
        ylabel('Rate (events/min)')

        % [h,p,ci,stats] = ttest2(pullData{1,tt}, pullData{2,tt});
        [p,h,stats] = ranksum(pullData{1,tt}, pullData{2,tt});

        title({mouse(m).day(d).phase(tt).name, ['p = ' num2str(round(p,3))]})

    end %phase

    save the raw data in an excel file for Nico
    fileName = ['rawRate_' tNames{d} '.xlsx'];
    cd('dataForNico')
    for tt = 1:2
        for r = 1:2
            sheetName = [mouse(end).day(d).phase(tt).name '-' rewNames{r}];
            tmpData = distxRate{r,tt};
            tmpData = tmpData(tmpData(:,1)<=5,2);
            writematrix(tmpData, fileName, 'Sheet', sheetName);
        end %reward
    end %trial type
    cd ../

    tmpFunc = @(x)(x(x(:,1)<=5,2));
    pullData = cellfun(tmpFunc, distxPkRate, 'UniformOutput', false);
    tmpMean = cellfun(@mean, pullData);
    tmpCIFunc =  @(x)(bootci(nboot, {@mean, x}, 'type', 'per'));
    tmpCI = cellfun(tmpCIFunc, pullData, 'UniformOutput', false);
    for tt = 1:length(mouse(m).day(d).phase)
        subplot(3,2,spMap(tt)+1)


        pullCI = horzcat(tmpCI{:,tt});

        bgraph = bar(tmpMean(:,tt), 'FaceColor', 'Flat');
        bgraph.CData = [rgb(rewCols{1}); rgb(rewCols{2})];
        hold on
        errorbar(1:2, tmpMean(:,tt), tmpMean(:,tt)' - pullCI(1,:),  pullCI(end,:)-tmpMean(:,tt)', 'Color', 'k', 'LineStyle', 'None')
        h = dotplot(1:2, pullData(:,tt), 0.1, [rgb(rewCols{1}); rgb(rewCols{2})], [rgb('White'); rgb('White')]);

        xticklabels(rewNames(1:2))
        ylabel('Peak rate (events/min)')
        % title(mouse(m).day(d).phase(tt).name)

        [p,h,stats] = ranksum(pullData{1,tt}, pullData{2,tt});

        title({mouse(m).day(d).phase(tt).name, ['p = ' num2str(round(p,3))]})
    end %phase

    save the peak data in an excel file for Nico
    fileName = ['peakRate_' tNames{d} '.xlsx'];
    cd('dataForNico')
    for tt = 1:2
        for r = 1:2
            sheetName = [mouse(end).day(d).phase(tt).name '-' rewNames{r}];
            tmpData = distxPkRate{r,tt};
            tmpData = tmpData(tmpData(:,1)<=5,2);
            writematrix(tmpData, fileName, 'Sheet', sheetName);
        end %reward
    end %trial type
    cd ../

    if saveOrNot == 1
        saveas(f8, figtitle8, 'fig');
        saveas(f8, figtitle8, 'png');
        saveas(f8, figtitle8, 'epsc');
    end %save

    figtitle9 = ['pfDistxSize_' tNames{d}];
    f9 = figure('Name', figtitle9, 'Position', [1 483 677 243]);

    tmpFunc = @(x)(x(x(:,1)<=5,2));
    pullData = cellfun(tmpFunc, distxSize, 'UniformOutput', false);
    tmpMean = cellfun(@mean, pullData);
    tmpCIFunc =  @(x)(bootci(nboot, {@mean, x}, 'type', 'per'));
    tmpCI = cellfun(tmpCIFunc, pullData, 'UniformOutput', false);
    for tt = 1:length(mouse(m).day(d).phase)
        subplot(1,3,tt)
        % h = dotplot(1:2, pullData(:,tt), 0.1, [rgb(rewCols{1}); rgb(rewCols{2})], [rgb('White'); rgb('White')], 0.2);

        pullCI = horzcat(tmpCI{:,tt});
        bgraph = bar(tmpMean(:,tt), 'FaceColor', 'Flat');
        bgraph.CData = [rgb(rewCols{1}); rgb(rewCols{2})];
        hold on
        errorbar(1:2, tmpMean(:,tt), tmpMean(:,tt)' - pullCI(1,:),  pullCI(end,:)-tmpMean(:,tt)', 'Color', 'k', 'LineStyle', 'None')
        xticklabels(rewNames(1:2))
        ylabel('Place field size (cm^2)')
        title(mouse(m).day(d).phase(tt).name)

        [p,h,stats] = ranksum(pullData{1,tt}, pullData{2,tt});

    end %phase

    %save the size data in an excel file for Nico
    fileName = ['size_' tNames{d} '.xlsx'];
    cd('dataForNico')
    for tt = 1:2
        for r = 1:2
            sheetName = [mouse(end).day(d).phase(tt).name '-' rewNames{r}];
            tmpData = distxSize{r,tt};
            tmpData = tmpData(tmpData(:,1)<=5,2);
            writematrix(tmpData, fileName, 'Sheet', sheetName);
        end %reward
    end %trial type
    cd ../

    if saveOrNot == 1
        saveas(f9, figtitle9, 'fig');
        saveas(f9, figtitle9, 'png');
        saveas(f9, figtitle9, 'epsc');
    end %save

end %day

end %function