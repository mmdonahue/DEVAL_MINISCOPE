function deval_4_identifyPlaceFields(mouse)
% function deval_4_identifyPlaceFields(mouse)
%
% PURPOSE:
%   Determine if units had an identifiable place field.
%
% MMD
% Ewell Lab
% 05/2025

%% INITIALIZE FOR RMS

binSz = 4; %cm
minSize = 10; %cm

eastInds = 1:42;
westInds = [1:14 43:70];
rInds = [eastInds; westInds];
badInds = [36:42 64:70];

trNames = {'Large reward', 'Small reward'};


minMeanFr = 0;
minPkFr = 0;

%% ADD TO STRUCT

for m = 1:length(mouse)
    fprintf('Mouse %s\n', mouse(m).name)
    % for ts = 1:length(mouse(m).task)
    for ts = 1:4
        fprintf('\t%s\n', mouse(m).task(ts).name)
        for d = 1:length(mouse(m).task(ts).day)
            fprintf('\t\tDay %d/%d\n', d, length(mouse(m).task(ts).day))
            load([mouse(m).task(ts).day(d).dir '\byTrialRatemaps_v3.mat'])
            % if isfield(unit, 'pf')
            %     continue
            % end %already done
            % if ~isfield(unit, 'trialType')
            %     continue
            % end

            rews = vertcat(mouse(m).task(ts).day(d).trial(:).rewLoc);
            tmpFunc = @(x)(strncmp(x, mouse(m).largeRew, 4));
            rews = double(cellfun(tmpFunc, rews)); %reward inds
            rews(rews==0) = 2; %small reward = 2

            for u = 1:length(unit)
           
                for tt = 1:length(mouse(m).task(ts).day(d).phase)
                    tmpMax = 0;

                    for r = 1:2
                        unit(u).trialType(tt).trajectory(r).name = trNames{r};

                        goodTs = mouse(m).task(ts).day(d).phase(tt).trials;

                        rewMaps = vertcat(unit(u).trial(ismember(1:length(rews),goodTs)' + double(rews==r) == 2).linSmSpkRm);
                        meanMap = nanmean(rewMaps,1);
                        if isempty(meanMap)
                            meanMap = NaN(1,70);
                        end %never went to this reward in this trialtype

                        if strncmp(mouse(m).largeRew, 'West', 4)
                            meanMap = meanMap([1:14 43:end 15:42]); %rearrange so it's large reward first
                        end %large reward is West
                        unit(u).trialType(tt).trajectory(r).rm = meanMap(rInds(r,:));

                        tmpMax = max([tmpMax max(unit(u).trialType(tt).trajectory(r).rm)]);

                    end %reward

                    for r = 1:2
                        if sum(~isnan(unit(u).trialType(tt).trajectory(r).rm)) == 0
                            unit(u).trialType(tt).trajectory(r).pf = [];
                            continue %to next reward type
                        end %to turns in this firection for this tt
                        try
                            pf =  get_circtrack_pfs(unit(u).trialType(tt).trajectory(r).rm, binSz, minMeanFr, tmpMax*0.1, minSize);
                        catch; keyboard; end
                        if ~isempty(pf)
                            pf = rmfield(pf, 'radPos');
                            pf = rmfield(pf, 'pkPos');
                            unit(u).trialType(tt).trajectory(r).pf = pf;
                        else
                            unit(u).trialType(tt).trajectory(r).pf = [];
                        end

                    end %large and small trajectory
                end %trial type

            end %unit
            save([mouse(m).task(ts).day(d).dir '\byTrialRatemaps_v3.mat'], 'unit') %save it real quick

        end %day
    end %task
end %mouse





end %function