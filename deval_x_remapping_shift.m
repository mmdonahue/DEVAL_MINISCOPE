function deval_x_remapping_shift(mouse)
% function deval_x_remapping_shift(mouse)
%
% PURPOSE:
%   Plot place cell activity across the three phases of the
%   devaluation/shift day.
%
% MMD
% Ewell Lab
% 05/2025

%% OPTIONS

saveOrNot = 0;

adpOnly = 1; %set to 1 to exclude mice that did not adapt

minDis = 10; %minimum number of "spikes" in a day for a unit to be included

pThresh = 0.05;

%% INITIALIZE

rmsForVis = []; %for visualizing ratemaps across three phases

d = 4; %shift day

pathByMouse = cell(10,3);
spkLocs = [];
spkLocsUnstab = [];

spatCorr1d = cell(3,1);
spatCorr1dNoRew = cell(3,1);
rateOverlap1d = cell(3,1);
ifRO1d = cell(3,1); %IN-FIELD rate overlap
lrRO1d = cell(3,1); %LARGE REWARD rate overlap
spatCorr2d = cell(3,1);
rateOverlap2d = cell(3,1);
ifRO2d = cell(3,1); %IN-FIELD rate overlap
lrRO2d = cell(3,1); %LARGE REWARD rate overlap

ifSI2d = cell(3,1); %IN-FIELD rate overlap
lrSI2d = cell(3,1); %LARGE REWARD rate overlap

pkLoc = cell(3,1);

allRewAct = []; %for storing whether unit was active at the large or small reward
allRewField = []; %for storing whether unit had a field at the large or small reward
allRewFr = []; %for storing the firing rate of the unit at each reward
unstabRewFr = []; %for storing the firing rate of the unstable units at each reward
allUnstabRewAct = [];

allMaps1d = cell(3,1);
allMaps2d = cell(3,1); %for storing all the maps for pv corrs
allMaps2dUnstab = cell(3,1); %for storing all the maps for pv corrs

allPfInds = [];

allPfs = [];

lrFrByTrial = cell(3,1); %large reward firing rate by trial

cellInfo = [];
unstabCellInfo = [];

comboNames = {'Forced1 half-maps', 'Forced1/Forced2', 'Forced2 half-maps'};
rewNames = {'Large reward', 'Small reward'};

phCols = {'DeepPink', 'Green', 'SteelBlue'};
phNames = {'Pre', 'Post'};
% colMatrix = []; %initialize
% for c = 1:length(phCols)
%     colMatrix = [colMatrix; rgb(phCols{c})]; %fill in for dotplot
% end %cols

scAx = -1:0.05:1;

colMap = [rgb('White'); jet(100)];
parMap = [rgb('White'); parula(20)];

mapIndsRew = [1:42; 1:14 43:70];
binSz = 4;
binAx = 2:binSz:56;
rewLocs = [4 6; 4 50]; %east and west
rewLocInds(1,:) = sub2ind([14 14], [1 2 3 1 2 3], [1 1 1 2 2 2]);
rewLocInds(2,:) = sub2ind([14 14], [12 13 14 12 13 14], [1 1 1 2 2 2]);

rewCols = {'SteelBlue', 'Crimson', 'Violet', 'Khaki', 'DarkGray'}; %large small both neither
rewNames = {'Large reward', 'Small reward', 'Both', 'Neither'};
nboot = 1000;

saveDir = 'E:\Devaluation Project Miniscope\Meg analysis\RESULTS\shift_day_remap';
cd(saveDir)


%% BIN CENTERS FOR UN-LINEARIZING

binCtrs(1:14,1) = 1:14;
binCtrs(1:14,2) = 8; %stem
binCtrs(15:21,1) = 14;
binCtrs(15:21,2) = fliplr(1:7); %choice point -> East stem entrance
binCtrs(22:35,1) = fliplr(1:14);
binCtrs(22:35,2) = 1; %East stem -> East reward
binCtrs(36:42,1) = 1;
binCtrs(36:42,2) = 1:7; %East reward -> delay zone
binCtrs(43:49,1) = 14;
binCtrs(43:49,2) = 9:15; %choice point -> West stem entrance
binCtrs(50:63,1) = fliplr(1:14);
binCtrs(50:63,2) = 15; %West stem -> West reward
binCtrs(64:70,1) = 1;
binCtrs(64:70,2) = fliplr(9:15); %West reward -> delay zone

%% GET DATA

for m = 1:length(mouse)
    fprintf('Mouse %s\n', mouse(m).name)
    if adpOnly == 1 && mouse(m).adp == 0
        continue
    end %ignore non-adapted mouse

    coords = mouse(m).day(d).coords;
    if strcmp(mouse(m).largeRew, 'West')
        mRewLocs = flip(rewLocs);
        mRewLocInds = flip(rewLocInds);
        mapInds = flip(mapIndsRew);
        coords(:,3) = coords(:,3) * -1 + 56;
    else
        mRewLocs = rewLocs;
        mRewLocInds = rewLocInds;
        mapInds = mapIndsRew;
    end %assign for small and large reward

    %for paths
    runSpeed = smooth_runspeed(get_runspeed(mouse(m).day(d).coords), 500/1000); %smooth runspeed w 500 ms window
    for tt = [1 3]
        startInd = mouse(m).day(d).trial(mouse(m).day(d).phase(tt).trials(1)).tmInds(1);
        endInd = mouse(m).day(d).trial(mouse(m).day(d).phase(tt).trials(end)).tmInds(...
            find(~isnan(mouse(m).day(d).trial(mouse(m).day(d).phase(tt).trials(end)).tmInds), 1, 'Last')  );
        pathByMouse{m,tt} = coords(startInd:endInd,2:3);
    end %tt

    load([mouse(m).day(d).dir '\byTrialRatemaps_v3.mat'], 'unit')
    for u = 1:length(unit)
        if unit(u).numDis < minDis || min(vertcat(unit(u).trialType(1:length(mouse(m).day(d).phase)).pVal)) >= pThresh...
                || sum(isnan(vertcat(unit(u).trialType(1:length(mouse(m).day(d).phase)).meanCorr))) == length(mouse(m).day(d).phase)

            tt = 1; %for FORCED CHOICE PRE
            rews = vertcat(mouse(m).day(d).trial(:).rewLoc);
            tmpFunc = @(x)(strncmp(x, mouse(m).largeRew, 4));
            rews = double(cellfun(tmpFunc, rews)); %reward inds
            rews(rews==0) = 2; %small reward = 2
            tmpMaps = vertcat(unit(u).trial(mouse(m).day(d).phase(tt).trials).linSmSpkRm); %tt = 1 for FORCED CHOICE PRE
            tmpRews = rews(mouse(m).day(d).phase(tt).trials);
            rawMaps = NaN(1,42,2);

            for r = 1:2
                rawMaps(:,:,r) = nanmean(tmpMaps(tmpRews==r,mapInds(r,:)),1);
            end %reward

            storeRewFr = zeros(1,2); %for storing the firing rate
            unStabRewAct = zeros(1,2);
            for r = 1:2
                blFr = rawMaps(1,:,r);

                rewFr = mean(rawMaps(1,34:36,r)); %firing rate the at reward
                storeRewFr(r) = rewFr;
                zFr = (rewFr - nanmean(blFr(:))) ./ nanstd(blFr(:)); %z-score firing rate at the large reward
                unStabRewAct(r) = zFr > 1;
            end %reward

            allUnstabRewAct = [allUnstabRewAct; unStabRewAct];
            unstabRewFr = [unstabRewFr; storeRewFr];

            % do stuff for paths

            tmpLocs = [];
            for tt = [1 3]
                startInd = mouse(m).day(d).trial(mouse(m).day(d).phase(tt).trials(1)).tmInds(1);
                endInd = mouse(m).day(d).trial(mouse(m).day(d).phase(tt).trials(end)).tmInds(...
                    find(~isnan(mouse(m).day(d).trial(mouse(m).day(d).phase(tt).trials(end)).tmInds), 1, 'Last')  );

                ttSpkInds = find(unit(u).spkTms >= mouse(m).day(d).coords(startInd,1) & unit(u).spkTms <= mouse(m).day(d).coords(endInd,1)); %inds for the spikes that are in this tt
                ttSpkSpds = runSpeed(unit(u).spkInds(ttSpkInds),2); %find speed when these spikes occured
                ttSpkInds = ttSpkInds(ttSpkSpds>=2); %only include ones where the mouse was traveling over 2 cm/s
                tmpLocs = cat(1, tmpLocs, {coords(unit(u).spkInds(ttSpkInds),2:3)});
            end %phase
            spkLocsUnstab = cat(1, spkLocsUnstab, {tmpLocs});

            pullMaps2d = cat(3, unit(u).trial(:).smSpkRm);
            rm1 = nanmean(pullMaps2d(:,:,mouse(m).day(d).phase(1).trials),3);
            rm2 = nanmean(pullMaps2d(:,:,mouse(m).day(d).phase(3).trials),3);

            if strcmp(mouse(m).largeRew, 'West')
                allMaps2dUnstab{2} = cat(4, allMaps2dUnstab{2}, cat(3, flip(rm1), flip(rm2)));
            else
                allMaps2dUnstab{2} = cat(4, allMaps2dUnstab{2}, cat(3, rm1, rm2));
            end %whether to flip for east/west

            unstabCellInfo = [unstabCellInfo; m unit(u).ID];
            continue
        end %cell isn't spatially stable within session

        cellInfo = [cellInfo; m unit(u).ID];

        % do stuff for paths

        tmpLocs = [];
        for tt = [1 3]
            startInd = mouse(m).day(d).trial(mouse(m).day(d).phase(tt).trials(1)).tmInds(1);
            endInd = mouse(m).day(d).trial(mouse(m).day(d).phase(tt).trials(end)).tmInds(...
                find(~isnan(mouse(m).day(d).trial(mouse(m).day(d).phase(tt).trials(end)).tmInds), 1, 'Last')  );

            ttSpkInds = find(unit(u).spkTms >= mouse(m).day(d).coords(startInd,1) & unit(u).spkTms <= mouse(m).day(d).coords(endInd,1)); %inds for the spikes that are in this tt
            ttSpkSpds = runSpeed(unit(u).spkInds(ttSpkInds),2); %find speed when these spikes occured
            ttSpkInds = ttSpkInds(ttSpkSpds>=2); %only include ones where the mouse was traveling over 2 cm/s
            tmpLocs = cat(1, tmpLocs, {coords(unit(u).spkInds(ttSpkInds),2:3)});
        end %phase
        spkLocs = cat(1, spkLocs, {tmpLocs});

        %find out if this cell is reward active/has a place field at reward
        tt = 1; %for FORCED CHOICE PRE
        rews = vertcat(mouse(m).day(d).trial(:).rewLoc);
        tmpFunc = @(x)(strncmp(x, mouse(m).largeRew, 4));
        rews = double(cellfun(tmpFunc, rews)); %reward inds
        rews(rews==0) = 2; %small reward = 2
        tmpMaps = vertcat(unit(u).trial(mouse(m).day(d).phase(tt).trials).linSmSpkRm); %tt = 1 for FORCED CHOICE PRE
        tmpRews = rews(mouse(m).day(d).phase(tt).trials);
        rawMaps = NaN(1,42,2);


        for r = 1:2
            rawMaps(:,:,r) = nanmean(tmpMaps(tmpRews==r,mapInds(r,:)),1);
        end %reward

        rewAct = zeros(1,2); %for stating is cell is reward active
        storeRewFr = zeros(1,2); %for storing the firing rate
        for r = 1:2
            blFr = rawMaps(1,:,r);

            rewFr = mean(rawMaps(1,34:36,r)); %firing rate the at reward
            storeRewFr(r) = rewFr;
            zFr = (rewFr - nanmean(blFr(:))) ./ nanstd(blFr(:)); %z-score firing rate at the large reward
            rewAct(r) = zFr > 1;
        end %reward

        allRewFr = [allRewFr; storeRewFr];
        allRewAct = [allRewAct; rewAct];

        %next, check if any of this cells place fields are within 5 cm of the large or small reward
        pf = unit(u).pf2d;
        if ~isempty(pf)
            allInds = vertcat({pf(:).inds}); %inds from every place field
            % keyboard
            if strcmp(mouse(m).largeRew, 'West')
                flipInds = flip(1:14)';
                tmpFunc = @(x)(cat(2, x(:,1), flipInds(x(:,2))));
                pfInds = cellfun(tmpFunc, vertcat({pf(:).inds}), 'UniformOutput', false);
                allPfInds = cat(1, allPfInds, {pfInds});
            else
                allPfInds = cat(1, allPfInds, {horzcat({pf(:).inds})});
            end %wehter needs to flip

            rewField = zeros(1,2);
            for r = 1:2
                tmpFunc = @(x)(min(pdist2([binAx(x(:,1))' binAx(x(:,2))'], mRewLocs(r,:))));
                allMinDist = cellfun(tmpFunc, allInds);
                rewField(r) = min(allMinDist) <=5;
            end %reward
            allRewField = [allRewField; rewField];
        else
            allRewField = [allRewField; 0 0]; %no reward field bc no field at all
        end %if there is a pf

        pullMaps1d = vertcat(unit(u).trial(:).linSmSpkRm);
        pullMaps2d = cat(3, unit(u).trial(:).smSpkRm);
        tmpPfInfo = cat(2, vertcat(pf(:).area), vertcat(pf(:).pkFr), vertcat(pf(:).meanFr));
        allPfs = cat(1, allPfs, {tmpPfInfo});

        %half maps - forced 1

        %1d
        rm1 = nanmean(pullMaps1d(mouse(m).day(d).phase(1).trials(1):round(mouse(m).day(d).phase(1).trials(end)/2),:));
        rm2 = nanmean(pullMaps1d(round(mouse(m).day(d).phase(1).trials(end)/2)+1:mouse(m).day(d).phase(1).trials(end),:));

        if strcmp(mouse(m).largeRew, 'West')
            allMaps1d{1} = cat(4, allMaps1d{1}, cat(3, rm1([1:14 43:70 15:42]), rm2([1:14 43:70 15:42])));
        else
            allMaps1d{1} = cat(4, allMaps1d{1}, cat(3, rm1, rm2));
        end %check large rewrad size

        rm1NR = rm1([1:32 37:60 65:end]);
        rm2NR = rm2([1:32 37:60 65:end]);

        nanInds = isnan(rm1) | isnan(rm2);
        rm1(nanInds) = [];
        rm2(nanInds) = [];

        tmpCorr = corrcoef(rm1, rm2);
        spatCorr1d{1} = [spatCorr1d{1}; tmpCorr(2)];

        % no reward zones
        nanInds = isnan(rm1NR) | isnan(rm2NR);
        rm1NR(nanInds) = [];
        rm2NR(nanInds) = [];

        tmpCorr = corrcoef(rm1NR, rm2NR);
        spatCorr1dNoRew{1} = [spatCorr1dNoRew{1}; tmpCorr(2)];

        %2d
        rm1 = nanmean(pullMaps2d(:,:,mouse(m).day(d).phase(1).trials(1):round(mouse(m).day(d).phase(1).trials(end)/2)),3);
        rm2 = nanmean(pullMaps2d(:,:,round(mouse(m).day(d).phase(1).trials(end)/2)+1:mouse(m).day(d).phase(1).trials(end)),3);

        if strcmp(mouse(m).largeRew, 'West')
            allMaps2d{1} = cat(4, allMaps2d{1}, cat(3, flip(rm1), flip(rm2)));
        else
            allMaps2d{1} = cat(4, allMaps2d{1}, cat(3, rm1, rm2));
        end %whether to flip for east/west

        tmpIfRo = []; %initialize
        tmpSI = [];
        for p = 1:length(pf)
            fr1 = nanmean(rm1(sub2ind([14 14], pf(p).inds(:,2), pf(p).inds(:,1))));
            fr2 = nanmean(rm2(sub2ind([14 14], pf(p).inds(:,2), pf(p).inds(:,1))));

            tmpIfRo = [tmpIfRo; min([fr1 fr2]) ./ max([fr1 fr2])];
            tmpSI = [tmpSI; (fr2 - fr1) ./ (fr2 + fr1)];
        end %pf
        ifRO2d{1} = [ifRO2d{1}; {tmpIfRo}];
        ifSI2d{1} = [ifSI2d{1}; {tmpSI}];

        fr1 = nanmean(rm1(mRewLocInds(1,:))); %large rewaard
        fr2 = nanmean(rm2(mRewLocInds(1,:)));

        lrRO2d{1} = [lrRO2d{1}; min([fr1 fr2]) ./ max([fr1 fr2])];
        lrSI2d{1} = [lrSI2d{1}; (fr2 - fr1) ./ (fr2 + fr1)];

        nanInds = isnan(rm1) | isnan(rm2);
        rm1(nanInds) = [];
        rm2(nanInds) = [];

        tmpCorr = corrcoef(rm1, rm2);
        spatCorr2d{1} = [spatCorr2d{1}; tmpCorr(2)];

        %half maps - forced 2
        rm1 = nanmean(pullMaps1d(mouse(m).day(d).phase(3).trials(1):mouse(m).day(d).phase(3).trials(1)+round(length(mouse(m).day(d).phase(3).trials)/2)-1,:));
        rm2 = nanmean(pullMaps1d(mouse(m).day(d).phase(3).trials(1)+round(length(mouse(m).day(d).phase(3).trials)/2)+1:mouse(m).day(d).phase(3).trials(end),:));

        if strcmp(mouse(m).largeRew, 'West')
            allMaps1d{3} = cat(4, allMaps1d{3}, cat(3, rm1([1:14 43:70 15:42]), rm2([1:14 43:70 15:42])));
        else
            allMaps1d{3} = cat(4, allMaps1d{3}, cat(3, rm1, rm2));
        end %check large rewrad size


        rm1NR = rm1([1:32 37:60 65:end]);
        rm2NR = rm2([1:32 37:60 65:end]);

        nanInds = isnan(rm1) | isnan(rm2);
        rm1(nanInds) = [];
        rm2(nanInds) = [];

        tmpCorr = corrcoef(rm1, rm2);
        spatCorr1d{3} = [spatCorr1d{3}; tmpCorr(2)];

        % no reward zones
        nanInds = isnan(rm1NR) | isnan(rm2NR);
        rm1NR(nanInds) = [];
        rm2NR(nanInds) = [];

        tmpCorr = corrcoef(rm1NR, rm2NR);
        spatCorr1dNoRew{3} = [spatCorr1dNoRew{3}; tmpCorr(2)];

        rm1 = nanmean(pullMaps2d(:,:,mouse(m).day(d).phase(3).trials(1):mouse(m).day(d).phase(3).trials(1)+round(length(mouse(m).day(d).phase(3).trials)/2)-1),3);
        rm2 = nanmean(pullMaps2d(:,:,mouse(m).day(d).phase(3).trials(1)+round(length(mouse(m).day(d).phase(3).trials)/2)+1:mouse(m).day(d).phase(3).trials(end)),3);

        if strcmp(mouse(m).largeRew, 'West')
            allMaps2d{3} = cat(4, allMaps2d{3}, cat(3, flip(rm1), flip(rm2)));
        else
            allMaps2d{3} = cat(4, allMaps2d{3}, cat(3, rm1, rm2));
        end %whether to flip for east/west

        tmpIfRo = []; %initialize
        tmpSI = [];
        for p = 1:length(pf)
            fr1 = nanmean(rm1(sub2ind([14 14], pf(p).inds(:,2), pf(p).inds(:,1))));
            fr2 = nanmean(rm2(sub2ind([14 14], pf(p).inds(:,2), pf(p).inds(:,1))));

            tmpIfRo = [tmpIfRo; min([fr1 fr2]) ./ max([fr1 fr2])];
            tmpSI = [tmpSI; (fr2 - fr1) ./ (fr2 + fr1)];
        end %pf
        ifRO2d{3} = [ifRO2d{3}; {tmpIfRo}];
        ifSI2d{3} = [ifSI2d{3}; {tmpSI}];

        fr1 = nanmean(rm1(mRewLocInds(1,:))); %large rewaard
        fr2 = nanmean(rm2(mRewLocInds(1,:)));

        lrRO2d{3} = [lrRO2d{3}; min([fr1 fr2]) ./ max([fr1 fr2])];
        lrSI2d{3} = [lrSI2d{3}; (fr2 - fr1) ./ (fr2 + fr1)];

        nanInds = isnan(rm1) | isnan(rm2);
        rm1(nanInds) = [];
        rm2(nanInds) = [];

        tmpCorr = corrcoef(rm1, rm2);
        spatCorr2d{3} = [spatCorr2d{3}; tmpCorr(2)];

        %forced 1 vs forced 2
        rm1 = nanmean(pullMaps1d(mouse(m).day(d).phase(1).trials,:));
        rm2 = nanmean(pullMaps1d(mouse(m).day(d).phase(3).trials,:));

        if strcmp(mouse(m).largeRew, 'West')
            allMaps1d{2} = cat(4, allMaps1d{2}, cat(3, rm1([1:14 43:70 15:42]), rm2([1:14 43:70 15:42])));
        else
            allMaps1d{2} = cat(4, allMaps1d{2}, cat(3, rm1, rm2));
        end %check large rewrad size


        rm1NR = rm1([1:32 37:60 65:end]);
        rm2NR = rm2([1:32 37:60 65:end]);

        nanInds = isnan(rm1) | isnan(rm2);
        rm1(nanInds) = [];
        rm2(nanInds) = [];

        tmpCorr = corrcoef(rm1, rm2);
        spatCorr1d{2} = [spatCorr1d{2}; tmpCorr(2)];

        % no reward zones
        nanInds = isnan(rm1NR) | isnan(rm2NR);
        rm1NR(nanInds) = [];
        rm2NR(nanInds) = [];

        tmpCorr = corrcoef(rm1NR, rm2NR);
        spatCorr1dNoRew{2} = [spatCorr1dNoRew{2}; tmpCorr(2)];

        rm1 = nanmean(pullMaps2d(:,:,mouse(m).day(d).phase(1).trials),3);
        rm2 = nanmean(pullMaps2d(:,:,mouse(m).day(d).phase(3).trials),3);

        if strcmp(mouse(m).largeRew, 'West')
            allMaps2d{2} = cat(4, allMaps2d{2}, cat(3, flip(rm1), flip(rm2)));
        else
            allMaps2d{2} = cat(4, allMaps2d{2}, cat(3, rm1, rm2));
        end %whether to flip for east/west

        tmpIfRo = []; %initialize
        tmpSI = [];
        for p = 1:length(pf)
            fr1 = nanmean(rm1(sub2ind([14 14], pf(p).inds(:,2), pf(p).inds(:,1))));
            fr2 = nanmean(rm2(sub2ind([14 14], pf(p).inds(:,2), pf(p).inds(:,1))));

            tmpIfRo = [tmpIfRo; min([fr1 fr2]) ./ max([fr1 fr2])];
            tmpSI = [tmpSI; (fr2 - fr1) ./ (fr2 + fr1)];
        end %pf
        ifRO2d{2} = [ifRO2d{2}; {tmpIfRo}];
        ifSI2d{2} = [ifSI2d{2}; {tmpSI}];

        fr1 = nanmean(rm1(mRewLocInds(1,:))); %large rewaard
        fr2 = nanmean(rm2(mRewLocInds(1,:)));

        lrRO2d{2} = [lrRO2d{2}; min([fr1 fr2]) ./ max([fr1 fr2])];
        lrSI2d{2} = [lrSI2d{2}; (fr2 - fr1) ./ (fr2 + fr1)];

        nanInds = isnan(rm1) | isnan(rm2);
        rm1(nanInds) = [];
        rm2(nanInds) = [];

        tmpCorr = corrcoef(rm1, rm2);
        spatCorr2d{2} = [spatCorr2d{2}; tmpCorr(2)];

        rews = vertcat(mouse(m).day(d).trial(:).rewLoc);
        tmpFunc = @(x)(strncmp(x, mouse(m).largeRew, 4));
        rews = double(cellfun(tmpFunc, rews)); %reward inds
        rews(rews==0) = 2; %small reward = 2
        for tt = [1 3]
            tmpMaps = vertcat(unit(u).trial(mouse(m).day(d).phase(tt).trials).linSmSpkRm);
            tmpRews = rews(mouse(m).day(d).phase(tt).trials);

            lrFr = nanmean(tmpMaps(tmpRews==1,mapInds(1,34:36)),2); %firing rate at the large reward across trials for this phase
            lrFr(lrFr>10) = NaN; %anomolous data
            lrFrByTrial{tt} =  [lrFrByTrial{tt}; lrFr'];
        end %tt

        for tt = 2
            tmpMaps = vertcat(unit(u).trial(mouse(m).day(d).phase(tt).trials).linSmSpkRm);
            tmpRews = rews(mouse(m).day(d).phase(tt).trials);

            lrFr = NaN(20,1);
            lrFr(tmpRews==1) = nanmean(tmpMaps(tmpRews==1,mapInds(1,34:36)),2); %firing rate at the large reward across trials for this phase with NaNs for small reward
            lrFr(lrFr>10) = NaN; %anomolous data
            lrFrByTrial{tt} =  [lrFrByTrial{tt}; lrFr'];
        end %free choice trials

    end %unit

end %mouse


keyboard

%% FIGURE

figtitle = 'spatialCorrelation_shiftDay_1d';
figure('Name', figtitle, 'Position', [1 396 891 396])

subplot(1,2,1)
xLabs = cell(3,1);
for com = 1:3
    c = histc(spatCorr1d{com}, scAx);
    cumCorr = cumsum(c)/sum(c);

    plot(scAx, cumCorr, 'Color', rgb(phCols{com}))
    hold on
end %combo

xlim([min(min(horzcat(spatCorr1d{:}))) 1])
axis square
xlabel({'Spatial correlation coefficient', '\leftarrowLess stable   |   More stable\rightarrow'})
ylabel('Cumulative proportion')
legend(comboNames, 'Location', 'northwest')

subplot(1,2,2)
xLabs = cell(3,1);
delta = 0.1;
for com = 1:3

    [wDist, xAx] = WeightedProportion(spatCorr1d{com}(~isnan(spatCorr1d{com})), scAx(1), scAx(end), delta);

    plot(xAx, wDist, 'Color', rgb(phCols{com}))
    hold on

end %combo

xlim([min(min(horzcat(spatCorr1d{:}))) 1])
axis square
xlabel({'Spatial correlation coefficient', '\leftarrowLess stable   |   More stable\rightarrow'})
ylabel('Proportion')
legend(comboNames, 'Location', 'northwest')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    saveas(gcf, figtitle, 'epsc');
end %save

%% FIGURE

figtitle = 'spatialCorrelation_shiftDay_2d';
figure('Name', figtitle, 'Position', [1 396 891 396])

subplot(1,2,1)
xLabs = cell(3,1);
for com = 1:3
    c = histc(spatCorr2d{com}, scAx);
    cumCorr = cumsum(c)/sum(c);

    plot(scAx, cumCorr, 'Color', rgb(phCols{com}))
    hold on
end %combo

xlim([min(min(horzcat(spatCorr2d{:}))) 1])
axis square
xlabel({'Spatial correlation coefficient', '\leftarrowLess stable   |   More stable\rightarrow'})
ylabel('Cumulative proportion')
legend(comboNames, 'Location', 'northwest')

subplot(1,2,2)
xLabs = cell(3,1);
delta = 0.1;
for com = 1:3

    [wDist, xAx] = WeightedProportion(spatCorr2d{com}(~isnan(spatCorr2d{com})), scAx(1), scAx(end), delta);

    plot(xAx, wDist, 'Color', rgb(phCols{com}))
    hold on

end %combo

xlim([min(min(horzcat(spatCorr2d{:}))) 1])
axis square
xlabel({'Spatial correlation coefficient', '\leftarrowLess stable   |   More stable\rightarrow'})
ylabel('Proportion')
legend(comboNames, 'Location', 'south')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    saveas(gcf, figtitle, 'epsc');
end %save


%% PV CORR HEATMAP

figtitle = 'pvCorr_heatmap_shift_2d';
figure('Name', figtitle, 'Position', [1 396 891 396])

for com = 1:3
    subplot(1,3,com)
    tmpHm = NaN(14,14);

    for xb = 1:14
        for yb = 1:14
            pullBins1 = allMaps2d{com}(xb,yb,1,:);
            pullBins2 = allMaps2d{com}(xb,yb,2,:);

            nanInds = isnan(pullBins1) | isnan(pullBins2);
            pullBins1(nanInds) = [];
            pullBins2(nanInds) = [];

            tmpCorr = corrcoef(pullBins1, pullBins2);

            tmpHm(xb,yb) = tmpCorr(2);
        end %y bin
    end %x bin

    imagesc(tmpHm)
    colormap(colMap)
    axis xy
    axis square
    title(comboNames{com})
    axis off
    clim([-0.1 0.7])
end %combo

cbr = colorbar('Location', 'southoutside');
cbr.Position = [    0.6912    0.2    0.2134    0.02];
ylabel(cbr, 'Correlation')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    saveas(gcf, figtitle, 'epsc');
end %save

figtitle = 'pvCorr_heatmap_shift_1d';
figure('Name', figtitle, 'Position', [1 396 891 396])

for com = 1:3
    subplot(1,3,com)
    tmpHm = NaN(14,14);

    for xb = 1:70
        pullBins1 = allMaps1d{com}(1,xb,1,:);
        pullBins2 = allMaps1d{com}(1,xb,2,:);

        nanInds = isnan(pullBins1) | isnan(pullBins2);
        pullBins1(nanInds) = [];
        pullBins2(nanInds) = [];

        tmpCorr = corrcoef(pullBins1, pullBins2);

        tmpHm(binCtrs(xb,2),binCtrs(xb,1)) = tmpCorr(2);
    end %x bin

    imagesc(tmpHm)
    colormap(colMap)
    axis xy
    axis square
    title(comboNames{com})
    axis off
    clim([-0.1 0.85])
end %combo

cbr = colorbar('Location', 'southoutside');
cbr.Position = [    0.6912    0.2    0.2134    0.02];
ylabel(cbr, 'Correlation')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    saveas(gcf, figtitle, 'epsc');
end %save

%% PV CORR MAP BY RT TYPE

figtitle = 'pvCorr_heatmap_shift_RTcombos_1d';
figure('Name', figtitle, 'Position', [0 396 891 396])

com = 2;
subplot(1,3,com)
tmpHm = NaN(14,14);

pullInds = cell(3,1);
pullInds{1} = ~allRewField(:,1) & ~allRewField(:,2) & allRewAct(:,1) & ~allRewAct(:,2); %cells that are active at teh large reward and not small
pullInds{2} = ~allRewField(:,1) & ~allRewField(:,2) & ~allRewAct(:,1) & allRewAct(:,2); %cells that are active at teh small reward and not large
pullInds{3} = ~allRewField(:,1) & ~allRewField(:,2) & ~allRewAct(:,1) & ~allRewAct(:,2); %cells that are active at neither
tmpTtls = {'LRT', 'SRT', 'nonRT'};

for sp = 1:3
    for xb = 1:70
        pullBins1 = allMaps1d{com}(1,xb,1,pullInds{sp});
        pullBins2 = allMaps1d{com}(1,xb,2,pullInds{sp});

        nanInds = isnan(pullBins1) | isnan(pullBins2);
        pullBins1(nanInds) = [];
        pullBins2(nanInds) = [];

        tmpCorr = corrcoef(pullBins1, pullBins2);

        tmpHm(binCtrs(xb,2),binCtrs(xb,1)) = tmpCorr(2);
    end %x bin

    subplot(1,3,sp)
    imagesc(tmpHm)
    colormap(colMap)
    axis xy
    axis square
    title(tmpTtls{sp})
    axis off
    clim([-0.15 1])

    cbr = colorbar('Location', 'southoutside');
    ylabel(cbr, 'Correlation')
end %subplot

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    saveas(gcf, figtitle, 'epsc');
end %save

%% FOR PAPER - ONLY YELLOW CELLS F1/F2 COMPARISON

figtitle = 'LRT_nonLRT_F1F2remap';
figure('Name', figtitle)

subplot(1,2,1)
hold on
com = 2;
stat1 = [];
stat2 = [];

tmpCols = {'Gold', 'OrangeRed'};

r = 4; %cells with fields at neither rewards
pullInds = ~allRewField(:,1) & ~allRewField(:,2) & allRewAct(:,1);
data1 = spatCorr1d{com}(pullInds);
bgraph = bar(1, nanmean(spatCorr1d{com}(pullInds)), 'FaceColor', 'Flat');
bgraph.CData = rgb(tmpCols{1});
tmpCI = bootci(nboot, {@nanmean, spatCorr1d{com}(pullInds)}, 'type', 'per');
errorbar(1, nanmean(spatCorr1d{com}(pullInds)),  nanmean(spatCorr1d{com}(pullInds)) - min(tmpCI),...
    max(tmpCI) -  nanmean(spatCorr1d{com}(pullInds)), 'Color', 'k', 'LineStyle', 'None');
% scatter((r-0.1)+(rand(1,sum(pullInds))*0.2), spatCorr1d{com}(pullInds), 'MarkerFaceColor', rgb(rewCols{r}),...
%     'MarkerEdgeColor', rgb('White'), 'MarkerFaceAlpha', 0.25)

stat1 = spatCorr1d{com}(pullInds);

% statData = [statData; spatCorr1d{com}(pullInds)];
% statCom = [statCom; repmat(com, size(spatCorr1d{com}(pullInds)))];
% statRew = [statRew; repmat(r, size(spatCorr1d{com}(pullInds)))];
% statZ = [statZ; repmat(1, size(spatCorr1d{com}(pullInds)))];

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & ~allRewAct(:,1);
data2 = spatCorr1d{com}(pullInds);
bgraph = bar(2, nanmean(spatCorr1d{com}(pullInds)), 'FaceColor', 'Flat');
bgraph.CData = rgb(tmpCols{2});
bgraph.FaceAlpha = 0.5;
tmpCI = bootci(nboot, {@nanmean, spatCorr1d{com}(pullInds)}, 'type', 'per');
errorbar(2, nanmean(spatCorr1d{com}(pullInds)),  nanmean(spatCorr1d{com}(pullInds)) - min(tmpCI),...
    max(tmpCI) -  nanmean(spatCorr1d{com}(pullInds)), 'Color', 'k', 'LineStyle', 'None');

stat2 = spatCorr1d{com}(pullInds);

xticks(1:2)
xticklabels({'LRT', 'non-LRT'})
ylabel('Spatial correlation')
yticks(0:0.1:0.5)

[p,h,stats] = ranksum(stat1, stat2);
title(['Z = ' num2str(round(stats.zval,2)) ' p = ' num2str(p)])

% add shuffle - 5 reps of same n as LRT
subplot(1,2,2)
hold on
com = 2;
stat1 = [];
stat2 = [];

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & allRewAct(:,1);
data1 = spatCorr1d{com}(pullInds);
bgraph = bar(1, nanmean(spatCorr1d{com}(pullInds)), 'FaceColor', 'Flat');
bgraph.CData = rgb(tmpCols{1});
tmpCI = bootci(nboot, {@nanmean, spatCorr1d{com}(pullInds)}, 'type', 'per');
errorbar(1, nanmean(spatCorr1d{com}(pullInds)),  nanmean(spatCorr1d{com}(pullInds)) - min(tmpCI),...
    max(tmpCI) -  nanmean(spatCorr1d{com}(pullInds)), 'Color', 'k', 'LineStyle', 'None');
% scatter((r-0.1)+(rand(1,sum(pullInds))*0.2), spatCorr1d{com}(pullInds), 'MarkerFaceColor', rgb(rewCols{r}),...
%     'MarkerEdgeColor', rgb('White'), 'MarkerFaceAlpha', 0.25)

stat1 = spatCorr1d{com}(pullInds);

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & ~allRewAct(:,1);

for rep = 1:5
    rs = RandStream('mt19937ar', 'Seed', rep);

    useUnitsInds = sort(randsample(rs, sum(pullInds), 142));
    useUnits = find(pullInds); useUnits = useUnits(useUnitsInds);

    bgraph = bar(rep+1, nanmean(spatCorr1d{com}(useUnits)), 'FaceColor', 'Flat');
    bgraph.CData = rgb(tmpCols{2});
    tmpCI = bootci(nboot, {@nanmean, spatCorr1d{com}(useUnits)}, 'type', 'per');
    errorbar(rep+1, nanmean(spatCorr1d{com}(useUnits)),  nanmean(spatCorr1d{com}(useUnits)) - min(tmpCI),...
        max(tmpCI) -  nanmean(spatCorr1d{com}(useUnits)), 'Color', 'k', 'LineStyle', 'None');

    stat2 = spatCorr1d{com}(useUnits);

    [p,h,stats] = ranksum(stat1, stat2);
    text(rep+.75, 0.6, ['p = ' num2str(round(p,2))])
end %repetition

xticks(1:6)
xticklabels({'LRT', 'non-LRT rep 1', 'non-LRT rep 2', 'non-LRT rep 3', 'non-LRT rep 4', 'non-LRT rep 5'})
ylabel('Spatial correlation')
yticks(0:0.1:0.5)


if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    saveas(gcf, figtitle, 'epsc');
end %save

%% PAPER FOLLOW UP - REMOVE REWARD ZONES

figtitle = 'LRT_nonLRT_F1F2remap_noRewZones';
figure('Name', figtitle)

subplot(1,2,1)
hold on
com = 2;
stat1 = [];
stat2 = [];

tmpCols = {'Gold', 'OrangeRed'};

r = 4; %cells with fields at neither rewards
pullInds = ~allRewField(:,1) & ~allRewField(:,2) & allRewAct(:,1);
data1 = spatCorr1dNoRew{com}(pullInds);
bgraph = bar(1, nanmean(spatCorr1dNoRew{com}(pullInds)), 'FaceColor', 'Flat');
bgraph.CData = rgb(tmpCols{1});
tmpCI = bootci(nboot, {@nanmean, spatCorr1dNoRew{com}(pullInds)}, 'type', 'per');
errorbar(1, nanmean(spatCorr1dNoRew{com}(pullInds)),  nanmean(spatCorr1dNoRew{com}(pullInds)) - min(tmpCI),...
    max(tmpCI) -  nanmean(spatCorr1dNoRew{com}(pullInds)), 'Color', 'k', 'LineStyle', 'None');
% scatter((r-0.1)+(rand(1,sum(pullInds))*0.2), spatCorr1d{com}(pullInds), 'MarkerFaceColor', rgb(rewCols{r}),...
%     'MarkerEdgeColor', rgb('White'), 'MarkerFaceAlpha', 0.25)

stat1 = spatCorr1dNoRew{com}(pullInds);

% statData = [statData; spatCorr1d{com}(pullInds)];
% statCom = [statCom; repmat(com, size(spatCorr1d{com}(pullInds)))];
% statRew = [statRew; repmat(r, size(spatCorr1d{com}(pullInds)))];
% statZ = [statZ; repmat(1, size(spatCorr1d{com}(pullInds)))];

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & ~allRewAct(:,1);
data2 = spatCorr1dNoRew{com}(pullInds);
bgraph = bar(2, nanmean(spatCorr1dNoRew{com}(pullInds)), 'FaceColor', 'Flat');
bgraph.CData = rgb(tmpCols{2});
bgraph.FaceAlpha = 0.5;
tmpCI = bootci(nboot, {@nanmean, spatCorr1dNoRew{com}(pullInds)}, 'type', 'per');
errorbar(2, nanmean(spatCorr1dNoRew{com}(pullInds)),  nanmean(spatCorr1dNoRew{com}(pullInds)) - min(tmpCI),...
    max(tmpCI) -  nanmean(spatCorr1dNoRew{com}(pullInds)), 'Color', 'k', 'LineStyle', 'None');

stat2 = spatCorr1dNoRew{com}(pullInds);

xticks(1:2)
xticklabels({'LRT', 'non-LRT'})
ylabel('Spatial correlation')
yticks(0:0.1:0.5)

[p,h,stats] = ranksum(stat1, stat2);
title(['Z = ' num2str(round(stats.zval,2)) ' p = ' num2str(p)])

% add shuffle - 5 reps of same n as LRT
subplot(1,2,2)
hold on
com = 2;
stat1 = [];
stat2 = [];

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & allRewAct(:,1);
data1 = spatCorr1dNoRew{com}(pullInds);
bgraph = bar(1, nanmean(spatCorr1dNoRew{com}(pullInds)), 'FaceColor', 'Flat');
bgraph.CData = rgb(tmpCols{1});
tmpCI = bootci(nboot, {@nanmean, spatCorr1dNoRew{com}(pullInds)}, 'type', 'per');
errorbar(1, nanmean(spatCorr1dNoRew{com}(pullInds)),  nanmean(spatCorr1dNoRew{com}(pullInds)) - min(tmpCI),...
    max(tmpCI) -  nanmean(spatCorr1dNoRew{com}(pullInds)), 'Color', 'k', 'LineStyle', 'None');
% scatter((r-0.1)+(rand(1,sum(pullInds))*0.2), spatCorr1d{com}(pullInds), 'MarkerFaceColor', rgb(rewCols{r}),...
%     'MarkerEdgeColor', rgb('White'), 'MarkerFaceAlpha', 0.25)

stat1 = spatCorr1dNoRew{com}(pullInds);

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & ~allRewAct(:,1);

for rep = 1:5
    rs = RandStream('mt19937ar', 'Seed', rep);

    useUnitsInds = sort(randsample(rs, sum(pullInds), 142));
    useUnits = find(pullInds); useUnits = useUnits(useUnitsInds);

    bgraph = bar(rep+1, nanmean(spatCorr1dNoRew{com}(useUnits)), 'FaceColor', 'Flat');
    bgraph.CData = rgb(tmpCols{2});
    tmpCI = bootci(nboot, {@nanmean, spatCorr1dNoRew{com}(useUnits)}, 'type', 'per');
    errorbar(rep+1, nanmean(spatCorr1dNoRew{com}(useUnits)),  nanmean(spatCorr1dNoRew{com}(useUnits)) - min(tmpCI),...
        max(tmpCI) -  nanmean(spatCorr1dNoRew{com}(useUnits)), 'Color', 'k', 'LineStyle', 'None');

    stat2 = spatCorr1dNoRew{com}(useUnits);

    [p,h,stats] = ranksum(stat1, stat2);
    text(rep+.75, 0.6, ['p = ' num2str(round(p,2))])
end %repetition

xticks(1:6)
xticklabels({'LRT', 'non-LRT rep 1', 'non-LRT rep 2', 'non-LRT rep 3', 'non-LRT rep 4', 'non-LRT rep 5'})
ylabel('Spatial correlation')
yticks(0:0.1:0.5)

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf, 'Renderer', 'Painters')
    saveas(gcf, figtitle, 'epsc');
end %save


%% EXAMPLE MAPS FOR NICO - LRT

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & allRewAct(:,1); %cells that are active at the large reward
findInds = find(pullInds);

rs = RandStream('mt19937ar', 'Seed', 1);
useUnits = sort(randsample(rs, length(findInds), 40));

for ph = 1:2
    figtitle = ['LargeRewActCells_rateMaps_forNico_' phNames{ph} '_0'];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))
        tmpMap = allMaps2d{2}(:,:,ph,findInds(uInd));
        normMap = tmpMap ./ max(tmpMap(:));
        normMap(isnan(normMap)) = -0.1;
        imagesc(normMap)
        axis xy

        hold on
        tmpPfInds = allPfInds{findInds(uInd)};
        for p = 1:length(tmpPfInds)
            x = tmpPfInds{p}(:,1);
            y = tmpPfInds{p}(:,2);
            k = boundary(x, y);
            if ~isempty(k)
                plot(x(k), y(k), 'Color', 'k', 'LineWidth', 2);
            elseif mean(diff(y)) == 0
                line([min(x) max(x)], [y(1) y(1)], 'Color', 'k', 'LineWidth', 2)
            elseif mean(diff(x)) == 0
                line([x(1) x(1)], [min(y) max(y)], 'Color', 'k', 'LineWidth', 2)
            end
        end %place field

        axis square
        axis off
        colormap(colMap)
        caxis([-0.1 1])

        % title({['Peak = ' num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'], ['Rew = ' num2str(round(nanmean(tmpMap(rewLocInds(1,:)))*60,1)) ' min^-^1']})
        title([mouse(cellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(findInds(uInd),2))])

        text(3, 4, [num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'])

    end %unit
    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save
end %phase

for ph = 1:2
    figtitle = ['LargeRewActCells_rateMaps_forNico_' phNames{ph} '_1'];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))
        tmpMap = allMaps2d{2}(:,:,ph,findInds(uInd));
        normMap = tmpMap ./ max(max(max(allMaps2d{2}(:,:,:,findInds(uInd)))));
        normMap(isnan(normMap)) = -0.1;
        imagesc(normMap)
        axis xy

        hold on
        tmpPfInds = allPfInds{findInds(uInd)};
        for p = 1:length(tmpPfInds)
            x = tmpPfInds{p}(:,1);
            y = tmpPfInds{p}(:,2);
            k = boundary(x, y);
            if ~isempty(k)
                plot(x(k), y(k), 'Color', 'k', 'LineWidth', 2);
            elseif mean(diff(y)) == 0
                line([min(x) max(x)], [y(1) y(1)], 'Color', 'k', 'LineWidth', 2)
            elseif mean(diff(x)) == 0
                line([x(1) x(1)], [min(y) max(y)], 'Color', 'k', 'LineWidth', 2)
            end
        end %place field

        axis square
        axis off
        colormap(colMap)
        caxis([-0.1 1])

        % title({['Peak = ' num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'], ['Rew = ' num2str(round(nanmean(tmpMap(rewLocInds(1,:)))*60,1)) ' min^-^1']})
        title([mouse(cellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(findInds(uInd),2))])

        text(3, 4, [num2str(round(max(max(max(allMaps2d{2}(:,:,:,findInds(uInd)))))*60,1)) ' min^-^1'])

    end %unit

    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save
end %phase

pathByMouse = pathByMouse(:,[1,3]);

for ph = 1:2
    figtitle = ['LargeRewActCells_path_forNico_' phNames{ph}];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))

        plot(pathByMouse{cellInfo(findInds(uInd),1),ph}(:,1), pathByMouse{cellInfo(findInds(uInd),1),ph}(:,2), 'Color', rgb('LightGray'))

        hold on
        tmpLocs = spkLocs{findInds(uInd)}{ph};
        plot(tmpLocs(:,1), tmpLocs(:,2), 'r.')

        xlim([0 56])
        ylim([0 56])
        axis square
        axis off

        title([mouse(cellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(findInds(uInd),2))])

    end %unit


    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save

end %phase

%% EXAMPLE MAPS FOR NICO - SRT

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & ~allRewAct(:,1) & allRewAct(:,2); %cells that are active at the SMALL reward
findInds = find(pullInds);

rs = RandStream('mt19937ar', 'Seed', 1);
useUnits = sort(randsample(rs, length(findInds), 40));

for ph = 1:2
    figtitle = ['SmallRewActCells_rateMaps_forNico_' phNames{ph} '_0'];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))
        tmpMap = allMaps2d{2}(:,:,ph,findInds(uInd));
        normMap = tmpMap ./ max(tmpMap(:));
        normMap(isnan(normMap)) = -0.1;
        imagesc(normMap)
        axis xy

        hold on
        tmpPfInds = allPfInds{findInds(uInd)};
        for p = 1:length(tmpPfInds)
            x = tmpPfInds{p}(:,1);
            y = tmpPfInds{p}(:,2);
            k = boundary(x, y);
            if ~isempty(k)
                plot(x(k), y(k), 'Color', 'k', 'LineWidth', 2);
            elseif mean(diff(y)) == 0
                line([min(x) max(x)], [y(1) y(1)], 'Color', 'k', 'LineWidth', 2)
            elseif mean(diff(x)) == 0
                line([x(1) x(1)], [min(y) max(y)], 'Color', 'k', 'LineWidth', 2)
            end
        end %place field

        axis square
        axis off
        colormap(colMap)
        caxis([-0.1 1])

        % title({['Peak = ' num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'], ['Rew = ' num2str(round(nanmean(tmpMap(rewLocInds(1,:)))*60,1)) ' min^-^1']})
        title([mouse(cellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(findInds(uInd),2))])

        text(3, 4, [num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'])

    end %unit
    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save
end %phase

for ph = 1:2
    figtitle = ['SmallRewActCells_rateMaps_forNico_' phNames{ph} '_1'];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))
        tmpMap = allMaps2d{2}(:,:,ph,findInds(uInd));
        normMap = tmpMap ./ max(max(max(allMaps2d{2}(:,:,:,findInds(uInd)))));
        normMap(isnan(normMap)) = -0.1;
        imagesc(normMap)
        axis xy

        hold on
        tmpPfInds = allPfInds{findInds(uInd)};
        for p = 1:length(tmpPfInds)
            x = tmpPfInds{p}(:,1);
            y = tmpPfInds{p}(:,2);
            k = boundary(x, y);
            if ~isempty(k)
                plot(x(k), y(k), 'Color', 'k', 'LineWidth', 2);
            elseif mean(diff(y)) == 0
                line([min(x) max(x)], [y(1) y(1)], 'Color', 'k', 'LineWidth', 2)
            elseif mean(diff(x)) == 0
                line([x(1) x(1)], [min(y) max(y)], 'Color', 'k', 'LineWidth', 2)
            end
        end %place field

        axis square
        axis off
        colormap(colMap)
        caxis([-0.1 1])

        % title({['Peak = ' num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'], ['Rew = ' num2str(round(nanmean(tmpMap(rewLocInds(1,:)))*60,1)) ' min^-^1']})
        title([mouse(cellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(findInds(uInd),2))])

        text(3, 4, [num2str(round(max(max(max(allMaps2d{2}(:,:,:,findInds(uInd)))))*60,1)) ' min^-^1'])

    end %unit

    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save
end %phase

if size(pathByMouse,2) == 3
    pathByMouse = pathByMouse(:,[1,3]);
end

for ph = 1:2
    figtitle = ['SmallRewActCells_path_forNico_' phNames{ph}];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))

        plot(pathByMouse{cellInfo(findInds(uInd),1),ph}(:,1), pathByMouse{cellInfo(findInds(uInd),1),ph}(:,2), 'Color', rgb('LightGray'))

        hold on
        tmpLocs = spkLocs{findInds(uInd)}{ph};
        plot(tmpLocs(:,1), tmpLocs(:,2), 'r.')

        xlim([0 56])
        ylim([0 56])
        axis square
        axis off

        title([mouse(cellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(findInds(uInd),2))])

    end %unit


    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save

end %phase

%% EXAMPLE MAPS FOR NICO - LARGE REWARD PF

pullInds = allRewField(:,1) & ~allRewField(:,2); %cells that have pf at large
findInds = find(pullInds);

rs = RandStream('mt19937ar', 'Seed', 1);
useUnits = sort(randsample(rs, length(findInds), 40));

for ph = 1:2
    figtitle = ['LargeRewFieldCells_rateMaps_forNico_' phNames{ph} '_0'];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))
        tmpMap = allMaps2d{2}(:,:,ph,findInds(uInd));
        normMap = tmpMap ./ max(tmpMap(:));
        normMap(isnan(normMap)) = -0.1;
        imagesc(normMap)
        axis xy

        hold on
        tmpPfInds = allPfInds{findInds(uInd)};
        for p = 1:length(tmpPfInds)
            x = tmpPfInds{p}(:,1);
            y = tmpPfInds{p}(:,2);
            k = boundary(x, y);
            if ~isempty(k)
                plot(x(k), y(k), 'Color', 'k', 'LineWidth', 2);
            elseif mean(diff(y)) == 0
                line([min(x) max(x)], [y(1) y(1)], 'Color', 'k', 'LineWidth', 2)
            elseif mean(diff(x)) == 0
                line([x(1) x(1)], [min(y) max(y)], 'Color', 'k', 'LineWidth', 2)
            end
        end %place field

        axis square
        axis off
        colormap(colMap)
        caxis([-0.1 1])

        % title({['Peak = ' num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'], ['Rew = ' num2str(round(nanmean(tmpMap(rewLocInds(1,:)))*60,1)) ' min^-^1']})
        title([mouse(cellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(findInds(uInd),2))])

        text(3, 4, [num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'])

    end %unit
    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save
end %phase

for ph = 1:2
    figtitle = ['LargeRewFieldCells_rateMaps_forNico_' phNames{ph} '_1'];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))
        tmpMap = allMaps2d{2}(:,:,ph,findInds(uInd));
        normMap = tmpMap ./ max(max(max(allMaps2d{2}(:,:,:,findInds(uInd)))));
        normMap(isnan(normMap)) = -0.1;
        imagesc(normMap)
        axis xy

        hold on
        tmpPfInds = allPfInds{findInds(uInd)};
        for p = 1:length(tmpPfInds)
            x = tmpPfInds{p}(:,1);
            y = tmpPfInds{p}(:,2);
            k = boundary(x, y);
            if ~isempty(k)
                plot(x(k), y(k), 'Color', 'k', 'LineWidth', 2);
            elseif mean(diff(y)) == 0
                line([min(x) max(x)], [y(1) y(1)], 'Color', 'k', 'LineWidth', 2)
            elseif mean(diff(x)) == 0
                line([x(1) x(1)], [min(y) max(y)], 'Color', 'k', 'LineWidth', 2)
            end
        end %place field

        axis square
        axis off
        colormap(colMap)
        caxis([-0.1 1])

        % title({['Peak = ' num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'], ['Rew = ' num2str(round(nanmean(tmpMap(rewLocInds(1,:)))*60,1)) ' min^-^1']})
        title([mouse(cellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(findInds(uInd),2))])

        text(3, 4, [num2str(round(max(max(max(allMaps2d{2}(:,:,:,findInds(uInd)))))*60,1)) ' min^-^1'])

    end %unit

    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save
end %phase

if size(pathByMouse,2) == 3
    pathByMouse = pathByMouse(:,[1,3]);
end

for ph = 1:2
    figtitle = ['LargeRewFieldCells_path_forNico_' phNames{ph}];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))

        plot(pathByMouse{cellInfo(findInds(uInd),1),ph}(:,1), pathByMouse{cellInfo(findInds(uInd),1),ph}(:,2), 'Color', rgb('LightGray'))

        hold on
        tmpLocs = spkLocs{findInds(uInd)}{ph};
        plot(tmpLocs(:,1), tmpLocs(:,2), 'r.')

        xlim([0 56])
        ylim([0 56])
        axis square
        axis off

        title([mouse(cellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(findInds(uInd),2))])

    end %unit


    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save

end %phase

%% EXAMPLE MAPS FOR NICO - SMALL REWARD PF

pullInds = ~allRewField(:,1) & allRewField(:,2); %cells that have pf at small
findInds = find(pullInds);

rs = RandStream('mt19937ar', 'Seed', 1);
useUnits = sort(randsample(rs, length(findInds), 40));

for ph = 1:2
    figtitle = ['SmallRewFieldCells_rateMaps_forNico_' phNames{ph} '_0'];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))
        tmpMap = allMaps2d{2}(:,:,ph,findInds(uInd));
        normMap = tmpMap ./ max(tmpMap(:));
        normMap(isnan(normMap)) = -0.1;
        imagesc(normMap)
        axis xy

        hold on
        tmpPfInds = allPfInds{findInds(uInd)};
        for p = 1:length(tmpPfInds)
            x = tmpPfInds{p}(:,1);
            y = tmpPfInds{p}(:,2);
            k = boundary(x, y);
            if ~isempty(k)
                plot(x(k), y(k), 'Color', 'k', 'LineWidth', 2);
            elseif mean(diff(y)) == 0
                line([min(x) max(x)], [y(1) y(1)], 'Color', 'k', 'LineWidth', 2)
            elseif mean(diff(x)) == 0
                line([x(1) x(1)], [min(y) max(y)], 'Color', 'k', 'LineWidth', 2)
            end
        end %place field

        axis square
        axis off
        colormap(colMap)
        caxis([-0.1 1])

        % title({['Peak = ' num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'], ['Rew = ' num2str(round(nanmean(tmpMap(rewLocInds(1,:)))*60,1)) ' min^-^1']})
        title([mouse(cellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(findInds(uInd),2))])

        text(3, 4, [num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'])

    end %unit
    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save
end %phase

for ph = 1:2
    figtitle = ['SmallRewFieldCells_rateMaps_forNico_' phNames{ph} '_1'];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))
        tmpMap = allMaps2d{2}(:,:,ph,findInds(uInd));
        normMap = tmpMap ./ max(max(max(allMaps2d{2}(:,:,:,findInds(uInd)))));
        normMap(isnan(normMap)) = -0.1;
        imagesc(normMap)
        axis xy

        hold on
        tmpPfInds = allPfInds{findInds(uInd)};
        for p = 1:length(tmpPfInds)
            x = tmpPfInds{p}(:,1);
            y = tmpPfInds{p}(:,2);
            k = boundary(x, y);
            if ~isempty(k)
                plot(x(k), y(k), 'Color', 'k', 'LineWidth', 2);
            elseif mean(diff(y)) == 0
                line([min(x) max(x)], [y(1) y(1)], 'Color', 'k', 'LineWidth', 2)
            elseif mean(diff(x)) == 0
                line([x(1) x(1)], [min(y) max(y)], 'Color', 'k', 'LineWidth', 2)
            end
        end %place field

        axis square
        axis off
        colormap(colMap)
        caxis([-0.1 1])

        % title({['Peak = ' num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'], ['Rew = ' num2str(round(nanmean(tmpMap(rewLocInds(1,:)))*60,1)) ' min^-^1']})
        title([mouse(cellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(findInds(uInd),2))])

        text(3, 4, [num2str(round(max(max(max(allMaps2d{2}(:,:,:,findInds(uInd)))))*60,1)) ' min^-^1'])

    end %unit

    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save
end %phase

if size(pathByMouse,2) == 3
    pathByMouse = pathByMouse(:,[1,3]);
end

for ph = 1:2
    figtitle = ['SmallRewFieldCells_path_forNico_' phNames{ph}];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))

        plot(pathByMouse{cellInfo(findInds(uInd),1),ph}(:,1), pathByMouse{cellInfo(findInds(uInd),1),ph}(:,2), 'Color', rgb('LightGray'))

        hold on
        tmpLocs = spkLocs{findInds(uInd)}{ph};
        plot(tmpLocs(:,1), tmpLocs(:,2), 'r.')

        xlim([0 56])
        ylim([0 56])
        axis square
        axis off

        title([mouse(cellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(findInds(uInd),2))])

    end %unit


    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');

    end %save

end %phase

%% EXAMPLE MAPS FOR NICO - UNSTABLE

pullInds = allUnstabRewAct(:,1)  | allUnstabRewAct(:,2);
findInds = find(pullInds);

rs = RandStream('mt19937ar', 'Seed', 1);
useUnits = sort(randsample(rs, length(findInds), 40));

for ph = 1:2
    figtitle = ['UnstableCells_rateMaps_forNico_' phNames{ph} '_0'];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))
        tmpMap = allMaps2dUnstab{2}(:,:,ph,findInds(uInd));
        normMap = tmpMap ./ max(tmpMap(:));
        normMap(isnan(normMap)) = -0.1;
        imagesc(normMap)
        axis xy

        axis square
        axis off
        colormap(colMap)
        caxis([-0.1 1])

        % title({['Peak = ' num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'], ['Rew = ' num2str(round(nanmean(tmpMap(rewLocInds(1,:)))*60,1)) ' min^-^1']})
        title([mouse(unstabCellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(unstabCellInfo(uInd),2))])

        text(3, 4, [num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'])

    end %unit
    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save
end %phase

for ph = 1:2
    figtitle = ['UnstableCells_rateMaps_forNico_' phNames{ph} '_1'];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))
        tmpMap = allMaps2dUnstab{2}(:,:,ph,findInds(uInd));
        normMap = tmpMap ./ max(max(max(allMaps2dUnstab{2}(:,:,:,findInds(uInd)))));
        normMap(isnan(normMap)) = -0.1;
        imagesc(normMap)
        axis xy

        axis square
        axis off
        colormap(colMap)
        caxis([-0.1 1])

        % title({['Peak = ' num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'], ['Rew = ' num2str(round(nanmean(tmpMap(rewLocInds(1,:)))*60,1)) ' min^-^1']})
        title([mouse(unstabCellInfo(findInds(uInd),1)).name ' U' num2str(unstabCellInfo(findInds(uInd),2))])

        text(3, 4, [num2str(round(max(max(max(allMaps2dUnstab{2}(:,:,:,findInds(uInd)))))*60,1)) ' min^-^1'])

    end %unit

    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save
end %phase

if size(pathByMouse,2) == 3
    pathByMouse = pathByMouse(:,[1,3]);
end

for ph = 1:2
    figtitle = ['UnstableCells_path_forNico_' phNames{ph}];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))

        plot(pathByMouse{unstabCellInfo(findInds(uInd),1),ph}(:,1), pathByMouse{unstabCellInfo(findInds(uInd),1),ph}(:,2), 'Color', rgb('LightGray'))

        hold on
        tmpLocs = spkLocsUnstab{findInds(uInd)}{ph};
        plot(tmpLocs(:,1), tmpLocs(:,2), 'r.')

        xlim([0 56])
        ylim([0 56])
        axis square
        axis off

        title([mouse(unstabCellInfo(findInds(uInd),1)).name ' U' num2str(unstabCellInfo(findInds(uInd),2))])

    end %unit


    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save

end %phase

%% EXAMPLE MAPS FOR NICO - NON LRT


pullInds = ~allRewField(:,1) & ~allRewField(:,2) & ~allRewAct(:,1); %cells that are NOT active at the large reward
findInds = find(pullInds);

rs = RandStream('mt19937ar', 'Seed', 1);
useUnits = sort(randsample(rs, length(findInds), 40));


for ph = 1:2
    figtitle = ['nonLRTCells_rateMaps_forNico_' phNames{ph} '_1'];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))
        tmpMap = allMaps2d{2}(:,:,ph,findInds(uInd));
        normMap = tmpMap ./ max(max(max(allMaps2d{2}(:,:,:,findInds(uInd)))));
        normMap(isnan(normMap)) = -0.1;
        imagesc(normMap)
        axis xy

        hold on
        tmpPfInds = allPfInds{findInds(uInd)};
        for p = 1:length(tmpPfInds)
            x = tmpPfInds{p}(:,1);
            y = tmpPfInds{p}(:,2);
            k = boundary(x, y);
            if ~isempty(k)
                plot(x(k), y(k), 'Color', 'k', 'LineWidth', 2);
            elseif mean(diff(y)) == 0
                line([min(x) max(x)], [y(1) y(1)], 'Color', 'k', 'LineWidth', 2)
            elseif mean(diff(x)) == 0
                line([x(1) x(1)], [min(y) max(y)], 'Color', 'k', 'LineWidth', 2)
            end
        end %place field

        axis square
        axis off
        colormap(colMap)
        caxis([-0.1 1])

        % title({['Peak = ' num2str(round(max(tmpMap(:))*60,1)) ' min^-^1'], ['Rew = ' num2str(round(nanmean(tmpMap(rewLocInds(1,:)))*60,1)) ' min^-^1']})
        title([mouse(cellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(findInds(uInd),2))])

        text(3, 4, [num2str(round(max(max(max(allMaps2d{2}(:,:,:,findInds(uInd)))))*60,1)) ' min^-^1'])

    end %unit

    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save
end %phase

pathByMouse = pathByMouse(:,[1,3]);

for ph = 1:2
    figtitle = ['nonLRTCells_rateMaps_forNico_' phNames{ph}];
    figure('Name', figtitle)

    for uInd = useUnits'
        subplot(5,8,match(uInd,useUnits))

        plot(pathByMouse{cellInfo(findInds(uInd),1),ph}(:,1), pathByMouse{cellInfo(findInds(uInd),1),ph}(:,2), 'Color', rgb('LightGray'))

        hold on
        tmpLocs = spkLocs{findInds(uInd)}{ph};
        plot(tmpLocs(:,1), tmpLocs(:,2), 'r.')

        xlim([0 56])
        ylim([0 56])
        axis square
        axis off

        title([mouse(cellInfo(findInds(uInd),1)).name ' U' num2str(cellInfo(findInds(uInd),2))])

    end %unit


    if saveOrNot == 1
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        set(gcf, 'Renderer', 'Painters')
        saveas(gcf, figtitle, 'epsc');
    end %save

end %phase

%% COMPARE SPAT CORR ACROSS CELL TYPES

figtitle = 'rewardFields_RT_spatCorr';
figure('Name', figtitle)
com = 2;

tmpCols = {'Blue', 'Red', 'Purple', 'SteelBlue', 'Crimson', 'Violet', 'Yellow'};
tmpNames = {'large reward pf', 'small reward pf', 'both reward pf', 'LRT', 'SRT', 'BRT', 'non-reward active'};

hold on
pullInds = allRewField(:,1) & ~allRewField(:,2); %cells with fields at large and not small
data1 = spatCorr1d{com}(pullInds);
bgraph = bar(1, nanmean(spatCorr1d{com}(pullInds)), 'FaceColor', 'Flat');
bgraph.CData = rgb(tmpCols{1});
tmpCI = bootci(nboot, {@nanmean, spatCorr1d{com}(pullInds)}, 'type', 'per');
errorbar(1, nanmean(spatCorr1d{com}(pullInds)),  nanmean(spatCorr1d{com}(pullInds)) - min(tmpCI),...
    max(tmpCI) -  nanmean(spatCorr1d{com}(pullInds)), 'Color', 'k', 'LineStyle', 'None');

pullInds = ~allRewField(:,1) & allRewField(:,2); %cells with fields at small and not large
data2 = spatCorr1d{com}(pullInds);
bgraph = bar(2, nanmean(spatCorr1d{com}(pullInds)), 'FaceColor', 'Flat');
bgraph.CData = rgb(tmpCols{2});
tmpCI = bootci(nboot, {@nanmean, spatCorr1d{com}(pullInds)}, 'type', 'per');
errorbar(2, nanmean(spatCorr1d{com}(pullInds)),  nanmean(spatCorr1d{com}(pullInds)) - min(tmpCI),...
    max(tmpCI) -  nanmean(spatCorr1d{com}(pullInds)), 'Color', 'k', 'LineStyle', 'None');

pullInds = allRewField(:,1) & allRewField(:,2); %cells with fields at both
data3 = spatCorr1d{com}(pullInds);
bgraph = bar(3, nanmean(spatCorr1d{com}(pullInds)), 'FaceColor', 'Flat');
bgraph.CData = rgb(tmpCols{3});
tmpCI = bootci(nboot, {@nanmean, spatCorr1d{com}(pullInds)}, 'type', 'per');
errorbar(3, nanmean(spatCorr1d{com}(pullInds)),  nanmean(spatCorr1d{com}(pullInds)) - min(tmpCI),...
    max(tmpCI) -  nanmean(spatCorr1d{com}(pullInds)), 'Color', 'k', 'LineStyle', 'None');

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & allRewAct(:,1) & ~allRewAct(:,2);
data4 = spatCorr1d{com}(pullInds);
bgraph = bar(4, nanmean(spatCorr1d{com}(pullInds)), 'FaceColor', 'Flat');
bgraph.CData = rgb(tmpCols{4});
tmpCI = bootci(nboot, {@nanmean, spatCorr1d{com}(pullInds)}, 'type', 'per');
errorbar(4, nanmean(spatCorr1d{com}(pullInds)),  nanmean(spatCorr1d{com}(pullInds)) - min(tmpCI),...
    max(tmpCI) -  nanmean(spatCorr1d{com}(pullInds)), 'Color', 'k', 'LineStyle', 'None');

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & ~allRewAct(:,1) & allRewAct(:,2);
data5 = spatCorr1d{com}(pullInds);
bgraph = bar(5, nanmean(spatCorr1d{com}(pullInds)), 'FaceColor', 'Flat');
bgraph.CData = rgb(tmpCols{5});
tmpCI = bootci(nboot, {@nanmean, spatCorr1d{com}(pullInds)}, 'type', 'per');
errorbar(5, nanmean(spatCorr1d{com}(pullInds)),  nanmean(spatCorr1d{com}(pullInds)) - min(tmpCI),...
    max(tmpCI) -  nanmean(spatCorr1d{com}(pullInds)), 'Color', 'k', 'LineStyle', 'None');

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & allRewAct(:,1) & allRewAct(:,2);
data6 = spatCorr1d{com}(pullInds);
bgraph = bar(6, nanmean(spatCorr1d{com}(pullInds)), 'FaceColor', 'Flat');
bgraph.CData = rgb(tmpCols{6});
tmpCI = bootci(nboot, {@nanmean, spatCorr1d{com}(pullInds)}, 'type', 'per');
errorbar(6, nanmean(spatCorr1d{com}(pullInds)),  nanmean(spatCorr1d{com}(pullInds)) - min(tmpCI),...
    max(tmpCI) -  nanmean(spatCorr1d{com}(pullInds)), 'Color', 'k', 'LineStyle', 'None');

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & ~allRewAct(:,1) & ~allRewAct(:,2);
data7 = spatCorr1d{com}(pullInds);
bgraph = bar(7, nanmean(spatCorr1d{com}(pullInds)), 'FaceColor', 'Flat');
bgraph.CData = rgb(tmpCols{7});
tmpCI = bootci(nboot, {@nanmean, spatCorr1d{com}(pullInds)}, 'type', 'per');
errorbar(7, nanmean(spatCorr1d{com}(pullInds)),  nanmean(spatCorr1d{com}(pullInds)) - min(tmpCI),...
    max(tmpCI) -  nanmean(spatCorr1d{com}(pullInds)), 'Color', 'k', 'LineStyle', 'None');

xticks(1:7)
xticklabels(tmpNames)
ylabel('Spatial correlation')

statData = [data1; data2; data3; data4; data5; data6; data7];
statCell = [ones(length(data1),1); ones(length(data2),1)*2; ones(length(data3),1)*3; ones(length(data4),1)*4; ones(length(data5),1)*5; ones(length(data6),1)*6; ones(length(data7),1)*7];

[p,tbl,stats,terms] = anovan(statData, {statCell}, 'model', 'full',...
    'varnames', {'CellType'}, 'Display', 'off');
if p < 0.05
    [c,~,h,gnames] = multcompare(stats, 'Display', 'off');

    sigInds = find(c(:,end)<0.05);
    yLine = 0.55:0.025:0.55+(0.025*length(sigInds));
    for sI = 1:length(sigInds)
        line([c(sigInds(sI),1) c(sigInds(sI),2)], [yLine(sI) yLine(sI)], 'Color', 'k')
    end %sigInds
end %p

% alt stats for what Nico used in the paper
statData = [data1; data2; data4; data5; data7];
statCell = [ones(length(data1),1); ones(length(data2),1)*2;  ones(length(data4),1)*4; ones(length(data5),1)*5; ones(length(data7),1)*7];
[p,tbl,stats,terms] = anovan(statData, {statCell}, 'model', 'full',...
    'varnames', {'CellType'}, 'Display', 'off');
if p < 0.05
    [c,~,h,gnames] = multcompare(stats, 'Display', 'off');

    sigInds = find(c(:,end)<0.05);
    yLine = 0.55:0.025:0.55+(0.025*length(sigInds));
    for sI = 1:length(sigInds)
        line([c(sigInds(sI),1) c(sigInds(sI),2)], [yLine(sI) yLine(sI)], 'Color', 'k')
    end %sigInds
end %p

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    saveas(gcf, figtitle, 'epsc');
end %save

%% FIG - REW FRS FOR RTS

rewFrByCell = cell(4,1); %initialize

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & allRewAct(:,1) & ~allRewAct(:,2);
rewFrByCell{1} = allRewFr(pullInds,:)*60;

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & ~allRewAct(:,1) & allRewAct(:,2);
rewFrByCell{2} = allRewFr(pullInds,:)*60;

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & sum(allRewAct,2);
rewFrByCell{3} = allRewFr(pullInds,:)*60;

pullInds = sum(unstabRewFr,2) > 0;
rewFrByCell{4} = unstabRewFr(pullInds,:)*60;

tmpMean = cellfun(@mean, rewFrByCell, 'UniformOutput', false);
tmpMean = vertcat(tmpMean{:});
tmpFunc = @(x)(bootci(nboot, {@mean, x}, 'type', 'per'));
tmpCI = cellfun(tmpFunc, rewFrByCell, 'UniformOutput', false)';

figtitle = 'rewFr_byRTtype';
figure('Name', figtitle)

tmpNames = {'LRT', 'SRT', 'all RT', 'no place field'};

bgraph = bar(tmpMean, 'FaceColor', 'Flat');
hold on

for k = 1:numel(bgraph)
    % Get the precise x-coordinates of the current bar series
    xCoordinates = bgraph(k).XEndPoints;

    %get relevant CI
    tmpFunc = @(x)(x(:,k));
    pullCI = cellfun(tmpFunc, tmpCI, 'UniformOutput', false);
    pullCI = horzcat(pullCI{:});
    % Plot the error bars without connecting lines
    errorbar(xCoordinates, tmpMean(:,k)', tmpMean(:,k)' - pullCI(1,:), pullCI(2,:) - tmpMean(:,k)', 'LineStyle', 'none', 'Color', 'k');

    bgraph(k).CData = rgb(rewCols{k});
end

xticks(1:4)
xticklabels(tmpNames)
ylabel('Firing rate (min^-^1)')
legend(rewNames(1:2))

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    saveas(gcf, figtitle, 'epsc');
end %save

[p,h,stats] = signrank(rewFrByCell{3}(:,1), rewFrByCell{3}(:,2));
[p,h,stats] = signrank(rewFrByCell{4}(:,1), rewFrByCell{4}(:,2));

statData = [rewFrByCell{1}(:,1); rewFrByCell{1}(:,2); rewFrByCell{2}(:,1); rewFrByCell{2}(:,2)];
statCell = [ones(size(rewFrByCell{1},1)*2,1); ones(size(rewFrByCell{2},1)*2,1)*2];
statRew = [ones(size(rewFrByCell{1},1),1); ones(size(rewFrByCell{1},1),1)*2;...
    ones(size(rewFrByCell{2},1),1); ones(size(rewFrByCell{2},1),1)*2];

[p,tbl,stats,terms] = anovan(statData, {statCell, statRew}, 'model', 'interaction', 'varnames', {'CellType', 'Reward'}, 'Display', 'off');
[c,~,h,gnames] = multcompare(stats, 'Dimension', [1 2], 'Display', 'off');


%% ALL RT PIE CHART

rewActInds = allRewAct(:,1) | allRewAct(:,2); %all cells that are active at either reward
actFields = allRewField(rewActInds,:); %pull out reward field info for cells that are active at reward

lRewOnly = sum(actFields(:,1) & ~actFields(:,2)); %cells that have a field at large
sRewOnly = sum(~actFields(:,1) & actFields(:,2)); %cells that have a field at small
bRew = sum(actFields(:,1) & actFields(:,2)); %cells that have a field at both
nRew = sum(~actFields(:,1) & ~actFields(:,2)); %cells that have a field at neither

unStabRew = sum(allUnstabRewAct(:,1) | allUnstabRewAct(:,2));

figtitle = 'combinedRewardActiveCells_pie';
figure('Name', figtitle, 'Position', [383 299 669 594])
tmpCols = {'SteelBlue', 'Crimson', 'Violet', 'Khaki', 'DarkGray'}; %large small both neither
p = pie([ lRewOnly sRewOnly bRew nRew unStabRew]);
patchHand = findobj(p, 'Type', 'Patch');
for t = 1:5
    patchHand(t).FaceColor = rgb(tmpCols{t});
end

lgd = legend({'Large reward place field', 'Small reward place field', 'Both rewards place field', 'Other location place field', 'No place field (spatially unstable)'}, 'Location', 'southoutside');
n = sum([ lRewOnly sRewOnly bRew nRew unStabRew]);

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    saveas(gcf, figtitle, 'epsc');
end %save


%% GET NS for email

nLRA = sum(~allRewField(:,1) & ~allRewField(:,2) & allRewFr(:,1)>0 & allRewFr(:,2)==0); %large reward active
nSRA = sum(~allRewField(:,1) & ~allRewField(:,2) & allRewFr(:,1)==0 & allRewFr(:,2)>0); %small reward active
nBRA = sum(~allRewField(:,1) & ~allRewField(:,2) & allRewFr(:,1)>0 & allRewFr(:,2)>0); %both reward active

%% FOR FUTURE REFERENCE

pullInds = ~allRewField(:,1) & ~allRewField(:,2) & allRewAct(:,1); %cells that are active at the large reward

LRTinfo = cellInfo(pullInds,:);



end %function