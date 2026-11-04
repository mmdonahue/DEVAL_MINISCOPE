function deval_x_preProcess_plot2dPfDetection(mouse)
% function deval_x_preProcess_plot2dPfDetection(mouse)
%
% PURPOSE:
%   Plot ratemaps with detected place cells bounds for a manual check of
%   detection methods.
%
% MMD
% Ewell Lab
% 07/2025


%% OPTIONS/INITIALIZE

d = 4;
minDis = 10;

rs = RandStream('mt19937ar', 'Seed', 1);

colMap = [rgb('White'); jet(100)];

saveDir = 'D:\Devaluation Project Miniscope\Meg analysis\RESULTS\pf_detection\2d';
% cd(saveDir)
% cd(mouse(end).day(d).name)
saveOrNot = 0;

%% GET DATA

for m = 1:length(mouse)
    figtitle1 = [mouse(m).name '_2dpfDetection'];
    f1 = figure('Name', figtitle1, 'Position', [1 4 1680 925]);

    figtitle2 = [mouse(m).name '_2dpfDetection_path'];
    f2 = figure('Name', figtitle2, 'Position', [1 4 1680 925]);

    load([mouse(m).day(d).dir '\byTrialRatemaps_v3.mat'], 'unit')

    keyboard
    allNumDis = vertcat(unit(:).numDis);
    % byTt = vertcat(unit(:).trialType);
    % thisTt = vertcat(byTt(:,tt));
    % allP = vertcat(thisTt(:).pVal);
    % goodUnits = find((double(allNumDis>=minDis) + double(allP<0.05) + double(allP>0) == 3));
    goodUnits = find(allNumDis>=minDis);
    useUnits = sort(randsample(rs, goodUnits, 50));
    useIDs = double(vertcat(unit(useUnits).ID));
    % useUnits = sort(randsample(goodUnits, 10));

    s = ncread([ mouse(m).day(d).dir '\minian_dataset.nc'], 'S');
    unit_id = ncread([ mouse(m).day(d).dir '\minian_dataset.nc'], 'unit_id');
    load([ mouse(m).day(d).dir '\processedData\DriftDetection.mat'])
    goodUnits = ~IsDrifting; %good units are not drifting
    s = s(:,goodUnits);
    unit_id = unit_id(goodUnits);

    tmDir = dir([ mouse(m).day(d).dir '\ts*']);
    tmConv = readmatrix([ mouse(m).day(d).dir '\' tmDir.name]);

    coords = mouse(m).day(d).coords;
    if strcmp(mouse(m).largeRew, 'West')
        coords(:,3) = coords(:,3) * -1 + 56; %flip
    end %if large reward is west, flip it

    runSpeed = get_runspeed(coords);
    smRs = smooth_runspeed(runSpeed, 500/1000); %smooth runspeed w 500 ms window

    uCntr = 0;
    for u = useUnits'
        uCntr = uCntr + 1;
        figure(f1)
        subplot(5,10,uCntr)

        tmpMap = unit(u).smSpkRm;
        tmpMap = tmpMap ./ max(tmpMap(:));
        tmpMap(isnan(tmpMap)) = -0.1;
        pf = unit(u).pf2d;
        pfInds = vertcat({pf(:).inds});

        if strcmp(mouse(m).largeRew, 'West')
            tmpMap = flip(tmpMap);
            flipInds = flip(1:14)';
            tmpFunc = @(x)(cat(2, x(:,1), flipInds(x(:,2))));
            pfInds = cellfun(tmpFunc, vertcat({pf(:).inds}), 'UniformOutput', false);
        end

        imagesc(tmpMap)
        axis xy
        colormap(colMap)
        hold on

        text(5, 12, [num2str(round(max(max(unit(u).smSpkRm)),2)) ' Hz'])
        for p = 1:length(pf)
            % plot(pf(p).inds(:,2), pf(p).inds(:,1))
            %
            % shp = alphaShape(pf(p).inds(:,1), pf(p).inds(:,2))
            % plot(shp)
            % plot(shp,'FaceColor','red','FaceAlpha',0.1)

            % x = pf(p).inds(:,1);
            % y = pf(p).inds(:,2);
            x = pfInds{p}(:,1);
            y = pfInds{p}(:,2);
            k = boundary(x, y);
            if ~isempty(k)
                plot(x(k), y(k), 'Color', 'k', 'LineWidth', 2);
            elseif mean(diff(y)) == 0
                line([min(x) max(x)], [y(1) y(1)], 'Color', 'k', 'LineWidth', 2)
            elseif mean(diff(x)) == 0
                line([x(1) x(1)], [min(y) max(y)], 'Color', 'k', 'LineWidth', 2)
            end
        end %pf

        axis square
        axis off

        figure(f2)
        subplot(5,10,uCntr)

        plot(coords(:,2), coords(:,3), 'Color', rgb('LightGray')) %plot mouse's path in this tt

        [~, spkInds] = findpeaks(s(:,u));

        spkCInds = tmConv(spkInds,2); %convert to coords inds

        checkSpd = smRs(spkCInds,2);
        spkCInds = spkCInds(checkSpd>2);

        hold on
        plot(coords(spkCInds,2), coords(spkCInds,3), 'r.')

        axis square
        axis off

    end %unit

    if saveOrNot == 1
        saveas(f1, figtitle1, 'fig');
        saveas(f1, figtitle1, 'png');
        set(f1,'renderer','Painters');
        saveas(f1, figtitle1, 'epsc');
        saveas(f1, figtitle1, 'svg');

        saveas(f2, figtitle2, 'fig');
        saveas(f2, figtitle2, 'png');
        set(f2,'renderer','Painters');
        saveas(f2, figtitle2, 'epsc');
        saveas(f2, figtitle2, 'svg');
    end %save

    close all

end %mouse




end %function