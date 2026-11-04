function deval_5_addSpikeTimes(mouse)
% function deval_5_addSpikeTimes(mouse)
%
% PURPOSE:
%   For some reason, I never added the "spike times" to the unit struct. I
%   am adding them now.
% 
% MMD



%% ADD TO STRUCT

for m = 1:length(mouse)
    fprintf('Mouse %s\n', mouse(m).name)
    for d = 1:length(mouse(m).day)
        if isempty(mouse(m).day(d).dir)
            continue
        end

        fprintf('\t\t%s\n', mouse(m).day(d).name)
        load([mouse(m).day(d).dir '\byTrialRatemaps_v3.mat'])

        s = ncread([ mouse(m).day(d).dir '\minian_dataset.nc'], 'S');
        unit_id = ncread([ mouse(m).day(d).dir '\minian_dataset.nc'], 'unit_id');

        load([ mouse(m).day(d).dir '\processedData\DriftDetection.mat'])
        goodUnits = ~IsDrifting; %good units are not drifting
        if length(unit) ~= sum(goodUnits)
            keyboard %PROBLEM!
        end

        tmDir = dir([ mouse(m).day(d).dir '\ts*']);
        tmConv = readmatrix([ mouse(m).day(d).dir '\' tmDir.name]);

        w = waitbar(0, 'Running through units...', 'Name', 'units');

        for u = 1:length(unit)
            matchID = find(unit_id==unit(u).ID);
            waitbar(u/sum(goodUnits), w, sprintf('Unit %d of %d', u, sum(goodUnits)));

            pks = find_peaksEH(s(:,matchID)', 0); %just using stop-gap code bc I can't get the signal processing toolbox while mathworks is down!
            spkInds = tmConv(match(pks.loc, tmConv(:,1)),2); %matching inds for the coords;
            spkInds(spkInds==0) = []; %sampling before the tracking comes online
            unit(u).spkInds = spkInds;
            unit(u).spkTms = mouse(m).day(d).coords(spkInds,1);

        end %unit

        close(w)
        save([ mouse(m).day(d).dir '\byTrialRatemaps_v3.mat'], 'unit')


    end %day

end %mouse






end %function