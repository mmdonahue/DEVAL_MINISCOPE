function  pf = get_2d_pfs_v2(rateMap, spatBinSz, percThresh, minPkFr, minPfArea)
% function  pf = get_2d_pfs_v2(rateMap, spatBinSz, percThresh, minPkFr, minPfArea)
%
% PURPOSE:
%   To get place field boundaries in x and y dimensions from a 2d rate-map. Based on code by JBT,
%   but simpler.
%
% % INPUT:
%     rateMap = 2d smoothed rate-map for a single unit
%   spatBinSz = the size of bins used to make the rate-map in cm
%               - Employed here to evaluate minimum place field size
%   percThresh = % of max firing rate to be considered part of the place field
%               - If left empty, percThresh = 50%;
%     minPkFr = minimum peak firing rate for a unit to be a place cell
%               - If left empty, minPkFr = 1 Hz;
%   minPfArea = minimum area for a 2d place field, in cm^2
%               - If left empty, minPfSize = 10 cm^2
% OUTPUT:
%       pf = structure array of length #placeFields with subfields:
%            -  pf.inds = the bins, from the linearized rate map, that the place field is in
%            -  pf.area = area, in cm, of the place field
%            - pf.pkFr = peak in-field firing rate (Hz)
%            - pf.meanFr = mean in-field firing rate (Hz)
%
% MMD
% Colgin Lab
% 08/2024

%% CHECK INPUTS AND ASSIGN DEFAULTS

if nargin == 2 || isempty(percThresh)
    percThresh = 50;  %percent
end


if nargin < 4 || isempty(minPkFr)
    minPkFr = 1;  %Hz -- minimum peak FR for a cell to be place cell
end

if nargin < 5 || isempty(minPfArea)
    minPfArea = 10;  %cm^2
end

pf = [];

%% FIND PLACE FIELDS

normRm = rateMap ./ max(max(rateMap));
exBnry = normRm > (percThresh/100);
potPfs = bwconncomp(exBnry, 8); %potential place fields

if potPfs.NumObjects == 0
    return
end

for p = 1:potPfs.NumObjects
    tmpPk = max(rateMap(potPfs.PixelIdxList{p}));

    if tmpPk < minPkFr
        continue
    end %not high enough peak firing rate

    tmpArea = length(potPfs.PixelIdxList{p}) .* spatBinSz;
    if tmpArea < minPfArea
        continue
    end

    pfInd = length(pf) + 1;
    tmpInds = potPfs.PixelIdxList{p};
    [row, col] = ind2sub(size(rateMap), tmpInds);
    pf(pfInd).inds = cat(2, col, row);
    pf(pfInd).area = tmpArea;
    pf(pfInd).pkFr  = tmpPk;
    pf(pfInd).meanFr  = mean(rateMap(potPfs.PixelIdxList{p}));

end %potential pfs


end %funcion