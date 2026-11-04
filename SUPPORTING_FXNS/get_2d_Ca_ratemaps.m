function [CaRatemap, spkRatemap, tpb] = get_2d_Ca_ratemaps(deSig, coords, tmConv, binSz, x_arenaSize, y_arenaSize, runThresh, durCrit)
% function [CaRatemap, spkRatemap, tpb] = get_2d_Ca_ratemaps(deSig, coords, tmConv, binSz, x_arenaSize, y_arenaSize, runThresh, durCrit)
%
% PURPOSE:
%   Calculate 2-dimensional ratemaps using the de-convolved deltaF/F Ca2+
%   singal.
%
% INPUT:
%   deSig = deconvolved signal from minian for one unit
%   coords = (output from convert_coords) array with:
%           coords(:,1) = time in seconds
%           coords(:,2) = x coordinates in cm
%           coords(:,3) = y coordinates in cm
%   tmConv = time series for reconciling different sampling rates from
%       position tracker and Ca2+ imaging, with:
%           tmConv(:,1) = timestamps for Ca2+ imaging (~20 fps)
%           tmConv(:,2) = timestamps for position camera (~30 fps)
%   binSz = size of spatial bin for ratemaps
%   x_arenaSize = size of arena in cm in x dimension
%   y_arenaSize = size of arena in cm in y dimension
%   runThresh = minimum speed in cm/smouse must be traveling for activity
%       to be included. Set to 0 for no filter
%   durCrit = minimum amount of time in s that a mouse must occupy a bin
%       for it to be inlcuded in ratemaps
%
% OUTPUT:
%   ratemap = 2d array with average deltaF/F values for each spatial bin
%   tpb = mean time mouse spent in each spatial bin
%
% MMD for Ewell lab based on code by JBT for Colgin Lab
% 05/2025

%get the bins
xBinCtrs = binSz/2:binSz:x_arenaSize;
yBinCtrs = binSz/2:binSz:y_arenaSize;

xBinEdges = 0:binSz:x_arenaSize;
yBinEdges = 0:binSz:y_arenaSize;

%get runspeed and the binary for when mouse is moving over the threshold
runSpeed = get_runspeed(coords);
smRs = smooth_runspeed(runSpeed, 500/1000); %smooth runspeed w 500 ms window
rsBnry = zeros(1, length(coords));
rsBnry(smRs(:,2)>runThresh) = 1;

% This is used for making the time-stamps for when the rat was in each spatial bin more precise
halfFt = mean(diff(coords(:,1)))/2; %half of the frame-rate
%      Mark start time for when the rat was in each bin as the frame-time (mid-point) - half the frame rate
%      and mark the end time as frame-time (mid-point) + half the frame rate

%% GET SPIKES TIMES

pks = find_peaksEH(deSig', 0); %just using stop-gap code bc I can't get the signal processing toolbox while mathworks is down!
spkInds = tmConv(ismember(tmConv(:,1), pks.loc),2); %matching inds for the coords;

%% GET THE TIMES WHEN THE RAT WAS IN EACH SPATIAL BIN

binTimes = cell(length(xBinCtrs), length(yBinCtrs)); %catcher for time-stamps when rat is in each bin
binInds = cell(length(xBinCtrs), length(yBinCtrs)); %catcher for time indices when rat is in each bin
binTimeSum = zeros(length(xBinCtrs), length(yBinCtrs));
for x = 1:length(xBinCtrs)
    xRange = [xBinEdges(x) xBinEdges(x+1)];
    for y = 1:length(yBinCtrs)
        yRange = [yBinEdges(y) yBinEdges(y+1)];

        tmpBinInds = coords(:,2)>=xRange(1) & coords(:,2)<xRange(2) & coords(:,3)>=yRange(1) & coords(:,3)<yRange(2);
        %                   Indices for when the rat is in each spatial bin

        tmpTimeBinary = zeros(1,length(coords(:,1))); %make a binary corresponding to the frames he was in the bin
        tmpTimeBinary(tmpBinInds) = 1;

        tmpComboBnry = zeros(1,length(tmpTimeBinary));
        tmpComboBnry(tmpTimeBinary == 1 & rsBnry == 1) = 1;
        tmpTimeBinary = tmpComboBnry;

        if sum(tmpTimeBinary) > 0 %if rat was in this spatial bin

            try
                timeChunks = bwconncomp(tmpTimeBinary,4); %find chunks of time where rat was in this spatial bin
            catch; keyboard; end
            %NOTE: I'm doing it this 'chunk' way vs. just multiplying # of frames by time-per-frame because it helps in finding spikes below
            %      where I search for all spikes that occurred within each time window

            cntr = 1;
            for c = 1:length(timeChunks.PixelIdxList)
                tmpChunkInds = [timeChunks.PixelIdxList{c}(1) timeChunks.PixelIdxList{c}(end)];
                tmpChunkTimes = coords(tmpChunkInds,1);
                binInds{x,y}(cntr,:) = tmpChunkInds;
                binTimes{x,y}(cntr,:) =[tmpChunkTimes(1)-halfFt tmpChunkTimes(2)+halfFt]; %get start and end time for the chunk

                cntr = cntr + 1;
            end %for each chunk of time the rat was in this bin

            binTimeSum(x,y) = sum(diff(binTimes{x,y},[],2)); %amt of time rat spent in this spatial bin

        end %if the rat was in this spatial bin at all

    end %yRange
end %xRange


%% GET DELTAF/F FOR EACH SPATIAL BIN
binDelF =  zeros(length(xBinCtrs), length(yBinCtrs));
binSpks =  zeros(length(xBinCtrs), length(yBinCtrs));
for x = 1:length(xBinCtrs)
    for y = 1:length(yBinCtrs)
        if ~isempty(binInds{x,y})
            for t = 1:size(binInds{x,y},1)
                tmpInds = binInds{x,y}(t,:);

                [~, convTmpInds(1)] = min(abs(tmConv(:,2)-tmpInds(1))); %for converting from pos camera to Ca2+ camera times
                [~, convTmpInds(2)] = min(abs(tmConv(:,2)-tmpInds(2)));
                if convTmpInds(2) == 1
                    continue %to next time bin
                end %Ca2+ camera came on after the position camera
                try
                    CaInds = tmConv(convTmpInds,1);
                    CaInds(CaInds==0) = 1; %can't start from bin 0
                    tmpDelF = sum(deSig(CaInds(1):CaInds(2)));
                catch; keyboard; end
                binDelF(x,y) = binDelF(x,y) + tmpDelF;

                findSpks = find(spkInds>=tmpInds(1) & spkInds<=tmpInds(2));
                numSpks = length(findSpks);
                binSpks(x,y) = binSpks(x,y) + numSpks;
            end %each chunk of time
        end %if the rat was in this spatial bin
    end %y
end %x

%% FLIP binTimeSum AND binDelF AROUND SO THEY'RE LINED UP WITH COORDS THE RIGHT WAY
%    This happens because coords is referenced to the bottom left corner of a plot but
%    the x & y referencing above starts in the top left (i.e., binTimeSum(1,1) = top left corner).

binTimeSum = flip(rot90(binTimeSum));
binDelF = flip(rot90(binDelF));
binSpks = flip(rot90(binSpks));

%% IF MINIMUM TIME IN SPATIAL BIN CRITERIA BEING ENFORCED...
%   THEN REMOVE SPIKES FROM ANY SPATIAL BIN THAT RAT VISITED < minBinDur ms

binDelF(binTimeSum<(durCrit)) = 0;
binSpks(binTimeSum<(durCrit)) = 0;


%% CALCULATE THE RATE-MAP

CaRatemap = binDelF ./ binTimeSum;
spkRatemap = binSpks ./ binTimeSum;

%% RENAME

tpb = binTimeSum;

end %fxn