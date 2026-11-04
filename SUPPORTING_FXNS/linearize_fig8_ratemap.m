function [Caratemap, smCaRm, spkRm, smSpkRm, tpb] = linearize_fig8_ratemap(deSig, coords, linCoords, tmConv, runThresh, varargin)
% function [Caratemap, smCaRm, spkRm, smSpkRm, tpb] = linearize_fig8_ratemap(deSig, coords, linCoords, tmConv, runThresh, gWinSz, gWinStd)
%
% PURPOSE:
%   Create a linearized ratemap for the figure 8 maze. See note below for
%   which bins correspond to which parts of the maze. This assumes a 4 cm
%   bin size (default for linearize_fig8_coords function).
%
% INPUT:
%   deSig = deconvolved signal from minian for one unit
%   coords = (output from convert_coords) array with:
%           coords(:,1) = time in seconds
%           coords(:,2) = x coordinates in cm
%           coords(:,3) = y coordinates in cm
%   linCoords = linearized coordinates for mouse's position.
%       NOTE: use linearize_fig8_coords
%   tmConv = time series for reconciling different sampling rates from
%       position tracker and Ca2+ imaging, with:
%           tmConv(:,1) = timestamps for Ca2+ imaging (~20 fps)
%           tmConv(:,2) = timestamps for position camera (~30 fps)
%   runThresh = minimum speed in cm/smouse must be traveling for activity
%       to be included. Set to 0 for no filter
%
% OUTPUT:
%   Caratemap = unsmooted ratemap with average deltaF/F values for each spatial bin
%   smCaRm = smoothed ratemap with average deltaF/F values for each spatial bin
%   spkRm = unsmoothed ratemap with firing rate for each spatial bin
%       - firing rate calculated from deconvovled signal - time of peak =
%       spike time
%   smSpkRm = smoothed ratemap with firing rate for each spatial bin
%   tpb = mean time mouse spent in each spatial bin
%
% NOTES:
%   - This code is based on the 56 x 56 cm figure 8 maze. This can be
%   changed in the initial lines of this code.
%   - The bins that correspond to various positions are described below:
%       1-14: stem
%       15-21: choice point -> East stem entrance
%       22-35: East stem -> East reward
%       36-42: East reward -> delay zone
%       43-49: choice point -> West stem entrance
%       50-63: West stem -> West reward
%       64-70: West reward -> delay zone
%
% MMD
% 05/2025
% Ewell Lab

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
rsBnry(smRs(:,2)>runThresh) = 1;

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

%% GET DELTAF/F FOR EACH SPATIAL BIN

binDelF =  zeros(numBins,1);
for x = 1:numBins
    if ~isempty(binInds{x})
        for t = 1:size(binInds{x},1)
            tmpInds = binInds{x}(t,:);
            [~, convTmpInds(1)] = min(abs(tmConv(:,2)-tmpInds(1))); %for converting from pos camera to Ca2+ camera times
            [~, convTmpInds(2)] = min(abs(tmConv(:,2)-tmpInds(2)));

            if convTmpInds(2) == 1
                continue %to next time bin
            end %Ca2+ camera came on after the position camera

            CaInds = tmConv(convTmpInds,1);
            CaInds(CaInds==0) = 1; %can't start from bin 0
            try
                tmpDelF = sum(deSig(CaInds(1):CaInds(2)));
            catch; keyboard; end
            binDelF(x) = binDelF(x) + tmpDelF;
        end %time chunks
    end %there are inds
end %x aka bin

%% CALCULATE  CA MAP

Caratemap = binDelF ./ binTimeSum;
% ratemap(isnan(ratemap)) = 0;
Caratemap = Caratemap';
tpb = binTimeSum;

%% FIX CORNERS
% mouse can't actually visit the corner bins due to maze edges/reward ports

cornerBins = [22 36 50 64];
intpDelF = interp1(setxor(1:numBins, cornerBins), Caratemap(setxor(1:numBins, cornerBins)), cornerBins);
Caratemap(cornerBins) = intpDelF;

%% SMOOTH EAST TRAJECTORY

Etraj = Caratemap(1:42);
mapLen = length(Etraj);
triEMap = [Etraj Etraj Etraj];
% smEtraj = conv_cir(Etraj, gKrnl)';
smEtraj = nanconv(triEMap, gKrnl');
smEtraj = smEtraj(mapLen+1:2*mapLen);
% smEtraj(isnan(Etraj)) = NaN;

%% SMOOTH WEST TRJAECTORY

Wtraj = Caratemap([1:14 43:70]);
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

smCaRm = [smStem smEtraj(15:end) smWtraj(15:end)];

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
