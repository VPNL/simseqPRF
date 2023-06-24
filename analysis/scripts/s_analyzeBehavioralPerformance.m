%% s_analyzeBehavioralPerformance.m

% Define params
subjnrs     = [1,2,3,7,8,9,10,11,12,13];
projectDir  = fullfile(simseqRootPath);
saveFigs    = true; % Save figures or not

%% SIMSEQ DATA

% Preallocate space
nrRuns  = 8;
hitRate = NaN(length(subjnrs),nrRuns);
HTS     = hitRate;
MSS     = hitRate;
FA      = hitRate;
faRate  = hitRate;

% Loop over subjects
for s = 1:length(subjnrs)
    fprintf('\nLoading subject S%02d..',subjnrs(s))
    pths = getSubjectPaths(projectDir,subjnrs(s));
    saveFigDir = fullfile(projectDir, 'results','group', 'figures');
    
    nrStimFiles = length(pths.stimFiles);
    
    % Load stim files with behavior
    for ii = 1:nrStimFiles
        fprintf('.')
        p = load(pths.stimFiles{ii});
        
        hitRate(s,ii) = p.perf.hitRate;
        HTS(s,ii)     = length(p.perf.hitTr);
        MSS(s,ii)     = length(p.perf.missTr);
        FA(s,ii)      = length(p.perf.falseFlip);
        faRate(s,ii)  = p.perf.falseRate;
    end
    fprintf('Finished!\n')
end

% Replace first run of S8, when buttonbox didn't work with a NaN
hitRate(hitRate==0) = NaN;

fprintf('Average performance = %2.0f %% correct (with SD = %2.1f)\n', ...
    100*mean(mean(hitRate,2,'omitnan')),100*std(mean(hitRate,2,'omitnan')))
fprintf('Average false alarm rate = %02d %% (%02d trials)', ...
    mean(mean(faRate,2,'omitnan')),mean(mean(FA,2,'omitnan')))

%% RETINOTOPY DATA

% Preallocate space
nrRuns  = 4;
hitRate = NaN(length(subjnrs),nrRuns);
HTS     = hitRate;
MSS     = hitRate;
FA      = hitRate;
faRate  = hitRate;

% Loop over subjects
for s = 1:length(subjnrs)
    fprintf('\nLoading subject S%02d..',subjnrs(s))
    pths = getSubjectPaths(projectDir,subjnrs(s));
    saveFigDir = fullfile(projectDir, 'results','group', 'figures');
    
    toonDir = fullfile(pths.dataDirToon, 'behavior');
    responseFiles = dir([toonDir '/*.mat']);
    nrStimFiles =  length(responseFiles);
    
    % Load stim files with behavior
    for ii = 1:nrStimFiles
        fprintf('.')
        p = load(fullfile(responseFiles(ii).folder,responseFiles(ii).name));
        dt = length(p.response.keyCode)/p.params.scanDuration;
        hitR = []; falseA = [];
        hits = find(p.response.keyCode);
        hitsClean = [];
        for jj = 1:length(hits)
            if jj==length(hits)
                hitsClean = [hitsClean, hits(jj)];
            elseif diff([hits(jj),hits(jj+1)])==1
                % do nothing
            else
                hitsClean = [hitsClean, hits(jj)];
            end
        end
        fixChanges = find(diff(p.stimulus.fixSeq));
        timeOfButtonPress = p.response.secs(hitsClean);
        rspWindow = (120/dt)+1; %floor((p.params.fix.responseTime*1000)/8);
        
        % for now, this only deals with hits, no false alarms...
        for hh = hitsClean
            
            targetRange = (hh - rspWindow(1)):hh;
            targetRange = targetRange(targetRange<length(p.response.keyCode));
            if any(ismember(fixChanges',targetRange))
                hitR = [hitR hh]; % this is a hit
            else
                falseA = [falseA hh];
            end
        end
        fprintf('Subject %s\n',pths.toon)
        fprintf('All buttonpresses: %s\n',num2str(timeOfButtonPress))
        fprintf('All hits: %s\n',num2str(hitR))
        fprintf('Hit rate: %2.1f%%',length(hitR)./length(fixChanges)*100)
%         hts_sec = find(p.response.secs);        
        
        fixChangesAll{s,ii} = fixChanges;
        buttonPressAll{s,ii} = hitsClean;
        hitRAll{s,ii} = hitR;
        falseAAll{s,ii} = falseA;
        hitRate(s,ii) = length(hitR)./length(fixChanges);
        HTS(s,ii)     = length(hitR);
        MSS(s,ii)     = length(hits)-length(hitR);
        FA(s,ii)      = length(falseA);
    end
    fprintf('Finished!\n')
end

hitRate(hitRate==0)=NaN;

fprintf('Average performance = %2.0f %% correct (with SD = %2.1f)\n', ...
    100*mean(mean(hitRate,2,'omitnan'),'omitnan'),...
    100*std(mean(hitRate,2,'omitnan'),[],'omitnan'));
fprintf('Nr false alarm = %2.2f\n', mean(mean(FA,2,'omitnan'),'omitnan'))
