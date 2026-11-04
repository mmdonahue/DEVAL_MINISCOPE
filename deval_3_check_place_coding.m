function deval_3_check_place_coding(mouse)
% function deval_3_check_place_coding(mouse)
%
% PURPOSE:
%   Check the properties to help determine if cells are place cells. This
%   code will check and save the number of discharges in a session and
%   determine if spatial stability measures pass a shuffled dataset as
%   Gergely did in his miniscope dataset.
%
% MMD
% 05/2025
% Ewell Lab


%% INITIALIZE FOR RMS

binSz = 4; %cm
runThresh = 2; %cm/s
durCrit = 150/1000; %150 ms

x_arenaSize = 56;
y_arenaSize = 56;

ttInds = {1, 2, [1 2]};
ttNames = {'free', 'forced', 'all'};
shiftTtNames = {'forced1', 'free', 'forced2', 'all'};

numShuff = 500;

%% ADD TO STRUCT

for m = 1:length(mouse)
    fprintf('Mouse %s\n', mouse(m).name)
  
    for ts = 1:length(mouse(m).task)
      
        fprintf('\t%s\n', mouse(m).task(ts).name)
        for d = 1:length(mouse(m).task(ts).day)
            fprintf('\t\tDay %d/%d\n', d, length(mouse(m).task(ts).day))

            load([mouse(m).task(ts).day(d).dir '\byTrialRatemaps_v3.mat'])
            if ~isfield(unit, 'numDis') || ~isfield(unit, 'trialType')
                s = ncread([mouse(m).task(ts).day(d).dir '\minian_dataset.nc'], 'S');
                load([mouse(m).task(ts).day(d).dir '\processedData\DriftDetection.mat'])
                goodUnits = ~IsDrifting; %good units are not drifting
                s = s(:,goodUnits);
            end %load

            if ~isfield(unit, 'numDis')
                fprintf('\t\t\tGetting number of spikes per day\n...')
                % w = waitbar(0, 'Running through units...', 'Name', 'units');
                for u = 1:length(unit)
                    % waitbar(u/sum(goodUnits), w, sprintf('Unit %d of %d', u, sum(goodUnits)));
                    pks = find_peaksEH(s(:,u), 0); %just using stop-gap code bc I can't get the signal processing toolbox while mathworks is down!
                    unit(u).numDis = length(pks.loc);
                end %unit

                save([mouse(m).task(ts).day(d).dir '\byTrialRatemaps_v3.mat'], 'unit') %save it real quick
                % close(w)
            end %number of discharges per unit is not yet a field

            if ~isfield(unit, 'trialType')
                fprintf('\t\t\tGetting spatial stability significance scores...\n')
                coords = mouse(m).task(ts).day(d).coords;
                linCoords = mouse(m).task(ts).day(d).linCoords;
                tmDir = dir([mouse(m).task(ts).day(d).dir '\ts*']);
                tmConv = readmatrix([mouse(m).task(ts).day(d).dir '\' tmDir.name]);

                % shuffle the signal for all units at once
                shuffInd = randsample(20*60+1:length(s)-20*60, numShuff); %shuffle inds that aren't within 1 min of original time points (20 fps * 60 s)
                shMaps = cell(length(unit),numShuff);
                shMaps(:) = {zeros(length(mouse(m).task(ts).day(d).trial),70)};

                %get and save shuffled maps
                w1 = waitbar(0, 'Running through shuffles...', 'Name', 'shuffles');
                for ff = 1:numShuff
                    waitbar(ff/numShuff, w1, sprintf('Shuffle %d', ff));
                    shuffSig = circshift(s,shuffInd(ff));

                    for t = 1:length(mouse(m).task(ts).day(d).trial)
                        tmInds = mouse(m).task(ts).day(d).trial(t).tmInds;
                        tmInds = [tmInds(1) tmInds(5)]; %trial start - ITI end
                        if isnan(tmInds(end))
                            tmInds(end) = length(coords);
                        end %just go to the end of the recording

                        trimCoords = coords(tmInds(1):tmInds(end),:);
                        trimLinCoords = linCoords(tmInds(1):tmInds(end),:);
                        trimConv = tmConv; %same size so we don't have to clear it every loop
                        trimConv(:,2) = tmConv(:,2) - tmInds(1);

                        for u = 1:length(unit)
                            [~, shMaps{u,ff}(t,:), ~] = linearize_fig8_spkRatemap(shuffSig(:,u), trimCoords, trimLinCoords, trimConv, runThresh);

                        end %unit
                    end %trial
                end %shuff
                close(w1)

                w = waitbar(0, 'Running through units...', 'Name', 'units');
                for u = 1:length(unit)

                    waitbar(u/sum(goodUnits), w, sprintf('Unit %d of %d', u, sum(goodUnits)));

                    %now calculated the correlation coefficients and p-values for each trial type
                    for tt = 1:length(mouse(m).task(ts).day(d).phase)
                        % keyboard
                        % tts = vertcat(mouse(m).task(ts).day(d).trial(:).type); %trial types
                        % goodTs = ismember(tts, ttInds{tt});
                        goodTs = mouse(m).task(ts).day(d).phase(tt).trials;
                        unit(u).trialType(tt).name = mouse(m).task(ts).day(d).phase(tt).name;

                        pullMaps = vertcat(unit(u).trial(:).linSmSpkRm); %pull ratemaps for trials that match desired trial type
                        pullMaps = pullMaps(goodTs,:);
                        meanMap = nanmean(pullMaps);
                        unit(u).trialType(tt).smRm = meanMap;

                        %now compare maps from 1st half and 2nd half of session
                        fMaps = nanmean(pullMaps(1:floor(size(pullMaps,1)/2),:));
                        sMaps = nanmean(pullMaps(floor(size(pullMaps,1)/2)+1:end,:));
                        nanInds = isnan(fMaps) | isnan(sMaps);
                        fMaps(nanInds) = [];
                        sMaps(nanInds) = [];
                        fsCorr = corrcoef(fMaps, sMaps);
                        unit(u).trialType(tt).fsCorr = fsCorr(2);

                        %now compare even and odd maps
                        eMaps = nanmean(pullMaps(2:2:size(pullMaps,1),:));
                        oMaps = nanmean(pullMaps(1:2:end,:));
                        nanInds = isnan(eMaps) | isnan(oMaps);
                        eMaps(nanInds) = [];
                        oMaps(nanInds) = [];
                        eoCorr = corrcoef(eMaps, oMaps);
                        unit(u).trialType(tt).eoCorr = eoCorr(2);

                        unit(u).trialType(tt).meanCorr = nanmean([unit(u).trialType(tt).fsCorr unit(u).trialType(tt).eoCorr]);

                        %now pull corr coeffs from the the shuffled dataset
                        shCorrs = NaN(1,numShuff);
                        for ff = 1:numShuff

                            %now compare maps from 1st half and 2nd half of session
                            fMaps = nanmean(shMaps{u,ff}(1:floor(size(pullMaps,1)/2),:));
                            sMaps = nanmean(shMaps{u,ff}(floor(size(pullMaps,1)/2)+1:end,:));
                            nanInds = isnan(fMaps) | isnan(sMaps);
                            fMaps(nanInds) = [];
                            sMaps(nanInds) = [];
                            fsCorr = corrcoef(fMaps, sMaps);

                            %now compare even and odd maps
                            eMaps = nanmean(shMaps{u,ff}(2:2:end,:));
                            oMaps = nanmean(shMaps{u,ff}(1:2:end,:));
                            nanInds = isnan(eMaps) | isnan(oMaps);
                            eMaps(nanInds) = [];
                            oMaps(nanInds) = [];
                            eoCorr = corrcoef(eMaps, oMaps);

                            shCorrs(ff) = mean([fsCorr(2) eoCorr(2)]);
                        end %shuffles

                        pVal = length(find(atanh(shCorrs)>atanh(unit(u).trialType(tt).meanCorr))) ./ length(shCorrs); %Fisher's transformation of corr coef
                        unit(u).trialType(tt).pVal = pVal;
                    end %trial type
                end %unit
                close(w)
            end %trial type - go through and get maps and stability estimates
            save([mouse(m).task(ts).day(d).dir '\byTrialRatemaps_v3.mat'], 'unit') %save it real quick

            % keyboard
            
        end %day
    end %task
end %mouse


end %function