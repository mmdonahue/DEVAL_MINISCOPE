function deval_2_add_ratemaps_toStruct(mouse)
% function deval_2_add_ratemaps_toStruct(mouse)
%
% PURPOSE:
%   Creates "unit" struct saved in raw data folder with unit information
%   and ratemaps for each cell.
%
% INPUT:
%   mouse = output from deval_1_...
%
% OUTPUT:
%   None in workspace.
%
% MMD
% 05/2025
 
% NOTE: For description of "v2" etc. in the filename, see notes at the very
% bottom of the function. 

%% INITIALIZE FOR RMS

binSz = 4; %cm
runThresh = 2; %cm/s
durCrit = 150/1000; %150 ms

x_arenaSize = 56;
y_arenaSize = 56;

%% ADD TO STRUCT

for m = 1:length(mouse)
    fprintf('Mouse %s\n', mouse(m).name)
    for ts = 1:length(mouse(m).task)

        fprintf('\t%s\n', mouse(m).task(ts).name)
        for d = 1:length(mouse(m).task(ts).day)
            fprintf('\t\tDay %d/%d\n', d, length(mouse(m).task(ts).day))


            if isfile([mouse(m).task(ts).day(d).dir '\byTrialRatemaps_v3.mat'])
                % load([mouse(m).task(ts).day(d).dir '\byTrialRatemaps_v2.mat'])
            % else
                 continue   
            end
            unit = struct(); %for saving individual trial ratemaps for this day
                % Load the calcium and spike data
                s = ncread([mouse(m).task(ts).day(d).dir '\minian_dataset.nc'], 'S');
                unit_id = ncread([mouse(m).task(ts).day(d).dir '\minian_dataset.nc'], 'unit_id');
                % a = ncread([mouse(m).task(ts).day(d).dir '\minian_dataset.nc'], 'a');
                load([mouse(m).task(ts).day(d).dir '\processedData\DriftDetection.mat'])
                goodUnits = ~IsDrifting; %good units are not drifting

                tmDir = dir([mouse(m).task(ts).day(d).dir '\ts*']);
                tmConv = readmatrix([mouse(m).task(ts).day(d).dir '\' tmDir.name]);

                s = s(:,goodUnits);
                unit_id = unit_id(goodUnits);

                coords = mouse(m).task(ts).day(d).coords;
                linCoords = mouse(m).task(ts).day(d).linCoords;

                w = waitbar(0, 'Running through units...', 'Name', 'units');
            
                for u = 1:sum(goodUnits)
                    waitbar(u/sum(goodUnits), w, sprintf('Unit %d of %d', u, sum(goodUnits)));

                    unit(u).ID = unit_id(u);

                    % Go and create trial-by-trial ratemaps and save in folder so we aren't making our data struct enormous
                    for t = 1:length(mouse(m).task(ts).day(d).trial)

                        tmInds = mouse(m).task(ts).day(d).trial(t).tmInds;
                        tmInds = [tmInds(1) tmInds(5)]; %trial start - ITI end
                        if isnan(tmInds(end))
                            tmInds(end) = length(coords);
                        end %just go to the end of the recording

                        trimCoords = coords(tmInds(1):tmInds(end),:);
                        trimLinCoords = linCoords(tmInds(1):tmInds(end),:);

                        [ratemap, tpb] = get_2d_Ca_ratemaps(s(:,u), trimCoords, tmConv, binSz, x_arenaSize, y_arenaSize,  runThresh, durCrit);
                        smRm = smooth_2d_Ca_ratemap(ratemap, binSz);
                        unit(u).trial(t).ratemap = ratemap;
                        unit(u).trial(t).smRm = smRm;
        
                        trimConv = tmConv; %same size so we don't have to clear it every loop
                        trimConv(:,2) = tmConv(:,2) - tmInds(1); 

                        [linCaRm, linCaSmRm, linSpkRm, linSmSpkRm, linTpb] = linearize_fig8_ratemap(s(:,u), trimCoords, trimLinCoords, trimConv, runThresh);
                        % if nansum(linCaSmRm) ~=0
                        % keyboard
                        % end
                        unit(u).trial(t).linCaRm = linCaRm;
                        unit(u).trial(t).linCaSmRm = linCaSmRm;
                        unit(u).trial(t).linSpkRm = linSpkRm;
                        unit(u).trial(t).linSmSpkRm = linSmSpkRm;

                    end %trial
                end %unit

                close(w)
                save([mouse(m).task(ts).day(d).dir '\byTrialRatemaps_v3.mat'], 'unit')
        end %day
    end %task
end %mouse

end %function


 %% VERSION HISTORY

 % v1 = initial attempt
 % v2 = fixed smoothing and corners of linearized ratemap, added firing rate map with peak of deconvolved signal = spike time
 % v3 = fixed some flipped coords and one mouse has wrong coords on one day