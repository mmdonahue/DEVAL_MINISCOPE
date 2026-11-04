function deval_x_preProcess_plotRateMaps(mouse)
%
% PURPOSE:
%   Plot ratemaps during running and stationary periods from 20 random
%   cells per mouse from the forced trials from the Pre-1 day.
%
% MMD
% Ewell Lab
% 06/2025

%% OPTIONS/INITIALIZE

d = 2; %pref2
tt = 2; %forced choice trials

minDis = 10;

rs = RandStream('mt19937ar', 'Seed', 1);
colMap = [rgb('White'); jet(100)];

runThresh = 2;

saveDir = 'E:\Devaluation Project Miniscope\Meg analysis\RESULTS\pf_detection\1d';
cd(saveDir)
saveOrNot = 0;

%% BIN CENTERS FOR UN-LINEARIZING

binCtrs(1:14,1) = 1:14;
binCtrs(1:14,2) = 8;
binCtrs(15:21,1) = 14;
binCtrs(15:21,2) = 9:15;
binCtrs(22:35,1) = fliplr(1:14);
binCtrs(22:35,2) = 15;
binCtrs(36:42,1) = 1;
binCtrs(36:42,2) = fliplr(9:15);
binCtrs(43:49,1) = 14;
binCtrs(43:49,2) = fliplr(1:7);
binCtrs(50:63,1) = fliplr(1:14);
binCtrs(50:63,2) = 1;
binCtrs(64:70,1) = 1;
binCtrs(64:70,2) = 1:7;

%% GET DATA

for m = 1:length(mouse)
    figtitle1 = [mouse(m).name '_ratemaps'];
    f1 = figure('Name', figtitle1, 'Position', [1 4 1680 925]);
    figtitle2 = [mouse(m).name '_stationaryRatemaps'];
    f2 = figure('Name', figtitle2, 'Position', [1 4 1680 925]);
    load([mouse(m).day(d).dir '\byTrialRatemaps_v3.mat'], 'unit')
    s = ncread([mouse(m).day(d).dir '\minian_dataset.nc'], 'S');
    load([mouse(m).day(d).dir '\processedData\DriftDetection.mat'])
    goodUnits = ~IsDrifting; %good units are not drifting
    s = s(:,goodUnits);

    allNumDis = vertcat(unit(:).numDis);
    byTt = vertcat(unit(:).trialType);
    thisTt = vertcat(byTt(:,tt));
    allP = vertcat(thisTt(:).pVal);
    goodUnits = find((double(allNumDis>=minDis) + double(allP<0.05) + double(allP>0) == 3));

    coords = mouse(m).day(d).coords;
    linCoords = mouse(m).day(d).linCoords;

    goodTs = find(ismember(vertcat(mouse(m).day(d).trial(:).type), tt));
    tmInds = [mouse(m).day(d).trial(goodTs(1)).tmInds(1) mouse(m).day(d).trial(goodTs(end)).tmInds(end-1)];

    tmDir = dir([mouse(m).day(d).dir '\ts*']);
    tmConv = readmatrix([mouse(m).day(d).dir '\' tmDir.name]);
    trimCoords = coords(tmInds(1):tmInds(end),:);
    trimLinCoords = linCoords(tmInds(1):tmInds(end),:);
    trimConv = tmConv; %same size so we don't have to clear it every loop
    trimConv(:,2) = tmConv(:,2) - tmInds(1);

    useUnits = sort(randsample(rs, goodUnits, 20));
keyboard
    uCntr = 0;
    for u = useUnits'
        uCntr = uCntr + 1;
        rm = unit(u).trialType(tt).smRm;
        normRm = rm ./ max(rm);
        if strncmp(mouse(m).largeRew, 'West', 4)
            normRm = normRm([1:14 43:end 15:42]); %rearrange so it's large reward first
        end %large reward is West

        tmpHm = NaN(15); %initialize 2d heat map
        for bin = 1:70
            tmpHm(binCtrs(bin,2),binCtrs(bin,1)) = normRm(bin);
        end

        pullRms = vertcat(unit(1).trial(goodTs).linSmSpkRm);

        figure(f1)
        subplot(5,4,uCntr)
        tmpHm(isnan(tmpHm)) = -0.1;
        imagesc(tmpHm)
        colormap(colMap)
        axis xy
        axis square
        % camroll(90)
        caxis([-0.1 1])

        xlim([0.5 14.5])
        text(0, 17, 'Small reward')
        text(0, -0.5, 'Large reward')
        text(2, 11, 'Delay')
        text(9,13,num2str(round(max(rm),2)))

    hold on
    line([20/4 20/4], [6 10], 'Color', 'Black', 'LineWidth', 2)
    % pgon = polyshape([0.5 0.75 0.75 0.5],[14.75 14.75 15.25 15.25]);
    % plot(pgon, 'FaceColor', 'k')
    % 
    % pgon = polyshape([0.5 0.75 0.75 0.5],[0.75 0.75 1.25 1.25]);
    % plot(pgon, 'FaceColor', 'k')
      axis off

        title({['Unit' num2str(unit(u).ID)], ''})

        [~, smSpkRm, ~] = stationary_pkRatemap(s(:,u), trimCoords, trimLinCoords, trimConv, runThresh);

        if max(smSpkRm) ~= 0
            normMap = smSpkRm ./ max(smSpkRm);
        else 
            normMap = smSpkRm;
        end
        if strncmp(mouse(m).largeRew, 'West', 4)
            normMap = normMap([1:14 43:end 15:42]); %rearrange so it's large reward first
        end %large reward is West

        tmpHm = NaN(15); %initialize 2d heat map
        for bin = 1:70
            tmpHm(binCtrs(bin,2),binCtrs(bin,1)) = normMap(bin);
        end

        figure(f2)
        subplot(5,4,uCntr)
        tmpHm(isnan(tmpHm)) = -0.1;
        imagesc(tmpHm)
        colormap(colMap)
        axis xy
        axis square
        % camroll(90)
        caxis([-0.1 1])

        xlim([0.5 14.5])
        text(0, 17, 'Small reward')
        text(0, -0.5, 'Large reward')
        % text(2, 9.5, 'Delay')
        text(2, 11, 'Delay')
        text(9,13,num2str(round(max(smSpkRm),2)))


    hold on
    line([20/4 20/4], [6 10], 'Color', 'Black', 'LineWidth', 2)
    % pgon = polyshape([0.5 0.75 0.75 0.5],[14.75 14.75 15.25 15.25]);
    % plot(pgon, 'FaceColor', 'k')
    % 
    % pgon = polyshape([0.5 0.75 0.75 0.5],[0.75 0.75 1.25 1.25]);
    % plot(pgon, 'FaceColor', 'k')
      axis off

        title({['Unit' num2str(unit(u).ID)], ''})

    end %unit

    figure(f1)
    cbr = colorbar('Location', 'southoutside');
    cbr.Position = [ 0.7482    0.07    0.1565    0.015];
    ylabel(cbr, 'Normalized Ca^2^+ event rate')

    figure(f2)
    cbr = colorbar('Location', 'southoutside');
     cbr.Position = [ 0.7482    0.07    0.1565    0.015];
    ylabel(cbr, 'Normalized Ca^2^+ event rate')

    if saveOrNot == 1
        % saveas(f1, figtitle1, 'fig');
        saveas(f1, figtitle1, 'png');
        exportgraphics(f1, [figtitle1 '.eps'], 'ContentType', 'vector')
        saveas(gcf, figtitle1, 'svg');
        % set(f1, 'Renderer', 'Painters')
        % saveas(f1, figtitle1, 'epsc');

        % saveas(f2, figtitle2, 'fig');
        saveas(f2, figtitle2, 'png');
        set(f2, 'Renderer', 'Painters')
        saveas(f2, figtitle2, 'epsc');
    end %saveOrNot

    close all
end %mouse




end %function

%% HELPER FUNCTION - reverse run thresh map
function [spkRm, smSpkRm, tpb] = stationary_pkRatemap(deSig, coords, linCoords, tmConv, runThresh, varargin)
%
% Modified from linearize_fig8_spkRatemap

%% INITIALIZE BINS

numBins = 70;
binSz = 4;

%% INITIALIZE KERNEL FOR SMOOTHING RATEMAP

if nargin == 5
    gWinStd = 8; %cm
    gWinSz = 8;
else
    gWinStd = varargin{2};
    gWinSz = varargin{1};
end %set kernel parameters

gWinSz = gWinSz/binSz; %convert window size to # of bins based on binSz
gWinStd = gWinStd/binSz;

gKrnl = gausskernel(gWinSz, gWinStd);

% This ratemap has to be smoothed in 2 phases: the stem-> E trajectory and
% stem-> W trajectory.

%% RUN SPEED FOR VELOCITY FILTER

runSpeed = get_runspeed(coords);
smRs = smooth_runspeed(runSpeed, 500/1000); %smooth runspeed w 500 ms window
rsBnry = zeros(length(coords),1);
rsBnry(smRs(:,2)<runThresh) = 1;

% This is used for making the time-stamps for when the rat was in each spatial bin more precise
halfFt = mean(diff(coords(:,1)))/2; %half of the frame-rate
%      Mark start time for when the rat was in each bin as the frame-time (mid-point) - half the frame rate
%      and mark the end time as frame-time (mid-point) + half the frame rate

%% GET THE TIMES WHEN THE RAT WAS IN EACH SPATIAL BIN

binTimes = cell(numBins,1);
binInds = cell(numBins,1); %catcher for time indices when rat is in each bin
binTimeSum = zeros(numBins,1);
for x = 1:numBins
    tmpBinInds = linCoords(:,2)==x; %times when mouse was in this bin
    tmpComboBnry = zeros(length(tmpBinInds),1);
    tmpComboBnry(tmpBinInds == 1 & rsBnry == 1) = 1; %times when mouse was in bin AND running at min speed or higher

    timeChunks = bwconncomp(tmpComboBnry,4); %find chunks of time where rat was in this spatial bin

    cntr = 1;
    for c = 1:length(timeChunks.PixelIdxList)
        tmpChunkInds = [timeChunks.PixelIdxList{c}(1) timeChunks.PixelIdxList{c}(end)];
        tmpChunkTimes = coords(tmpChunkInds,1);
        binInds{x}(cntr,:) = tmpChunkInds;
        binTimes{x}(cntr,:) =[tmpChunkTimes(1)-halfFt tmpChunkTimes(2)+halfFt]; %get start and end time for the chunk

        cntr = cntr + 1;
    end %for each chunk of time the rat was in this bin
    try
        binTimeSum(x) = sum(diff(binTimes{x},[],2)); %amt of time rat spent in this spatial bin
    catch; keyboard; end

end %x aka bin

tpb = binTimeSum;

%% GET SPIKE TIMES
% keyboard
pks = find_peaksEH(deSig', 0); %just using stop-gap code bc I can't get the signal processing toolbox while mathworks is down!
spkInds = tmConv(ismember(tmConv(:,1), pks.loc),2); %matching inds for the coords;

%% GET SPIKES FOR EACH SPATIAL BIN

binSpks =  zeros(numBins,1);
for x = 1:numBins
    if ~isempty(binInds{x})
        for t = 1:size(binInds{x},1)
            tmpInds = binInds{x}(t,:);
            % [~, convTmpInds(1)] = min(abs(tmConv(:,2)-tmpInds(1))); %for converting from pos camera to Ca2+ camera times
            % [~, convTmpInds(2)] = min(abs(tmConv(:,2)-tmpInds(2)));

            % if convTmpInds(2) == 1
            %     continue %to next time bin
            % end %Ca2+ camera came on after the position camera

            findSpks = find(spkInds>=tmpInds(1) & spkInds<=tmpInds(2));
            numSpks = length(findSpks);
            binSpks(x) = binSpks(x) + numSpks;
        end %time chunks
    end %there are inds
end %x aka bin

%% CALCULATE RATEMAP

spkRm = binSpks ./ binTimeSum;
% ratemap(isnan(ratemap)) = 0;
spkRm = spkRm';
% tpb = binTimeSum;

%% FIX CORNERS
% mouse can't actually visit the corner bins due to maze edges/reward ports

cornerBins = [22 36 50 64];
intpDelF = interp1(setxor(1:numBins, cornerBins), spkRm(setxor(1:numBins, cornerBins)), cornerBins);
spkRm(cornerBins) = intpDelF;

%% SMOOTH EAST TRAJECTORY

Etraj = spkRm(1:42);
mapLen = length(Etraj);
triEMap = [Etraj Etraj Etraj];
% smEtraj = conv_cir(Etraj, gKrnl)';
smEtraj = nanconv(triEMap, gKrnl');
smEtraj = smEtraj(mapLen+1:2*mapLen);
% smEtraj(isnan(Etraj)) = NaN;

%% SMOOTH WEST TRJAECTORY

Wtraj = spkRm([1:14 43:70]);
mapLen = length(Wtraj);
triWMap = [Wtraj Wtraj Wtraj];
smWtraj = nanconv(triWMap, gKrnl');
smWtraj = smWtraj(mapLen+1:2*mapLen);
% smWtraj = conv_cir(Wtraj, gKrnl)';
% smWtraj(isnan(Wtraj)) = NaN;

%% AVERAGE TOGETHER STEM
% since it got smoothed with both trajectorys
try
    vis = ~isnan(smEtraj(1:14)) + ~isnan(smWtraj(1:14)); %visited once or twice this trial etc.
    smStem = sum([smEtraj(1:14); smWtraj(1:14)], 'omitnan') ./ vis;
catch; keyboard; end

%% PUT IT ALL BACK TOGETHER

smSpkRm = [smStem smEtraj(15:end) smWtraj(15:end)];


end %function
