function pf = get_circtrack_pfs(rateMap, spatBinSz, varargin)
% function pf = get_circtrack(rateMap, spatBinSz, minMeanFr, minFr, minPfLen)
% 
% PURPOSE:
%   Get circle track place fields.
% 
% INPUTS:
%   rateMap = rateMap across the circle track for this place field,
%       smoothed (in Bieri et al., smoothing w/ Gaussian kernel std = 25
%       cm)
%   spatBinSz = bin size, in cm, used for ratemap
%   minMeanFr = minimum mean firing rate across the ratemap to be
%       considered a potential place cell (in Bieri et al., 2.5 Hz).
%       Default = 0 Hz (not used).
%   minFr (optional) = minimum firing rate for a pin to be considered 
%       part of a potential place field. Default = 1 Hz.
%   minPfLen (optional) = minimum place field length, in cm. Default = 10
%       cm.
% 
% OUTPUT:
%   pf = struct for each place field with fields:
%       inds = ratemap indicies for field
%       radPos = position in degress of place field
%       pkFr = peak firing rate in field
%       pkPos = position in degress of peak firing in field
% 
% MMD
% Colgin Lab - modified for Ewell Lab
% 06/2022 - modified from backward expansion pf code 05/2025

pf = []; %initialize

% rmBinSz = 360/length(rateMap); %ratemap bin size in deg
degBinCtrs = spatBinSz/2:spatBinSz:360; %centers of spatial bins, in deg

minMeanFr = 0;
if nargin > 2
    minMeanFr = varargin{1};
end %check for track diameter

if nanmean(rateMap) <= minMeanFr 
    return
end %doesn't meet minimum

minFr = 1; %Hz
if nargin > 3
    minFr = varargin{2};
end %check 

if max(rateMap) < minFr
    return
end %doesn't meet minimum

minPfLen = 10; %cm
if nargin > 3
    minPfLen = varargin{3};
end %check for track diameter

% minPfLen = 15 * (360/(trackDiameter*pi)); %degrees, apprx. 15 cm in Bieri et al. 2014
contigBins = round(minPfLen/spatBinSz); 

stdMinFr = 0; %Hz, standard minimum firing rate
pkFrCell = max(rateMap);

% minFrUse = max([stdMinFr 0.1*pkFrCell]); %firing rate above the lower of either 1 Hz or 10% of the peak firing rate of the cell

abvPeak = rateMap > minFr;
CC = bwconncomp(abvPeak, 4); %get connected components

inds = [];
for p = 1:CC.NumObjects
    if length(CC.PixelIdxList{1,p}) >= contigBins
        tmpInds = [CC.PixelIdxList{1,p}(1)-1 CC.PixelIdxList{1,p}(end)+1];
        inds = [inds; tmpInds];
    end %if long enough
end %potential place fields

%check if there's anything that wraps around 0 deg
if CC.PixelIdxList{1,1}(1) == 1 && CC.PixelIdxList{1,end}(end) == length(rateMap)
    tmpAllInds = [CC.PixelIdxList{1,end}' CC.PixelIdxList{1,1}'];
    if length(tmpAllInds) >= contigBins
        if ~isempty(inds)
        %delete old potential pfs to replace with more inclusive one IF
        %there are any old pfs
        oldPfs = [];
        tmpOld = [oldPfs find(inds(:,1) == tmpAllInds(1)-1)];
        if ~isempty(tmpOld)
            oldPfs = [oldPfs tmpOld];
        end
        tmpOld =  find(inds(:,2) == tmpAllInds(end)+1);
        if ~isempty(tmpOld)
            oldPfs = [oldPfs tmpOld];
        end
        inds(oldPfs,:) = [];
        end %there are any potential old pfs to delete
        inds = [inds; tmpAllInds(1)-1 tmpAllInds(end)+1];
    end %long enough
end %potential pf around 0

inds(inds==0) = length(rateMap);
inds(inds==length(rateMap)+1) = 1;

for p = 1:size(inds,1)
    if inds(p,2) > inds(p,1)
        pf(p).inds = inds(p,1):inds(p,2);
        pf(p).radPos = degBinCtrs(pf(p).inds);
        
        [pf(p).pkFr, maxInd] = max(rateMap(pf(p).inds));
        pf(p).pkPos = degBinCtrs(pf(p).inds(1) + maxInd - 1);
        pf(p).pkInd = pf(p).inds(1) + maxInd - 1;
    else
        pf(p).inds = [inds(p,1):length(rateMap) 1:inds(p,2)];
        pf(p).radPos = degBinCtrs(pf(p).inds);
        [pf(p).pkFr, maxInd] = max(rateMap(pf(p).inds));
        pf(p).pkPos = pf(p).radPos(maxInd);
        pf(p).pkInd = pf(p).inds(maxInd);
    end %whether cross 0
end %place fields

end %function