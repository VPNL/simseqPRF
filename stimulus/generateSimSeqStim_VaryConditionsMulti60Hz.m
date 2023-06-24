function generateSimSeqStim_VaryConditionsMulti60Hz()
% Function to generate variations of the sequential and simultaneous
% presentation of square stimuli in the upper right visual field, as in
% Kastner et al '98 (Science) and exp 1 in '01 (J Neurophys). Varations in-
% clude: longer durations, more transients in sequential presentation and
% more stimuli covering the field of view.
%
% Vary stim area and duration for both SIM & SEQ:
%       (6)  SEQ 16x 1 square presented 3456 ms (with 34 ms ISI).
%       (12) SIM 1x 16 squares presented 3456 ms, followed by blank.
%
% This function will save the image sequence, images themselves and movie.
% Images can be used for simulator and the second version of the simseq
% fMRI experiment.
%
% Written by Eline Kupers -- Stanford University 2020

%% 1 Set general params
clear;
showMovie   = false;   % whether or not to play a movie of the stimulus sequence
saveFigs    = true;    % whether or not to save figures
saveData    = true;    % whether or not to save sequence and images
createVideo = false;    % whether or not to save movie as .avi
verbose     = true;    % whether or not to plot figures

trialType = 'varyStimDur_x_Area4_2x2design';
version = 99;

if saveFigs || saveData
    saveFolder  = fullfile(simseqRootPath, 'stimulus', 'stim_simulation','varyConditions', [trialType '_v' num2str(version)]);
    if ~exist(saveFolder,'dir'); mkdir(saveFolder); end
    makeprettyfigures;
end

labels = {'SEQ1-4x200ms-4deg2','SEQ2-4x1000ms-4deg2',...
    'SEQ3-4x200ms-16deg2','SEQ4-4x1000ms-16deg2',...
    'SIM1-4x200ms-4deg2','SIM2-4x1000ms-4deg2',...
    'SIM3-4x200ms-16deg2','SIM4-4x1000ms-16deg2'};

%% 2 Set spatial screen/stimulus params
screenHeightPx  = 101;  % number of pixels of screen vertically
screenWidthPx   = 101;  % number of pixels of screen horizontally
screenHeightDeg = 24;   % screen diameter degrees of visual angle
pxPerDeg        = screenHeightPx/screenHeightDeg; %assumed pixels per degree:
ifi             = 1/60; % refresh rate for 60 hz monitor

%max eccentricity to use. So for each version we will try numLocations between -maxEcc and +maxEcc
maxEcc = screenHeightDeg/2; % deg visual angle

% Define grid, points are square centers
gridSpaceDeg = 1; %sqrt(1); % space in degrees between grid points, we use sqrt(1.9) for pilot 2 and sqrt(2) for pilot 2 to have eccen of 2
nrGridPoints = 15; %11;
hemifield    = 'bothf'; % choose 'urhf' (upper right, default if empty),
                        % 'llhf' (lower left) or 'bothf' (both urhf llhf)

%% 3 Set temporal screen/stimulus params
TR = 1000;    % ms

params = getStimParams(trialType, 'gridSpaceDeg',gridSpaceDeg, ...
    'nrGridPoints',nrGridPoints, 'hemifield',hemifield,'verbose',verbose);
params.durations.gap               = 2*ifi;     % (sec) min of 2 frames for blank between stim frames
params.durations.initialBlankDurS  = 0;       % s --> Add blank later to combine runs
params.durations.blankDurS         = 0;       % s
params.gridSpaceDeg                = gridSpaceDeg;
params.nrGridPoints                = nrGridPoints;
params.hemifield                   = hemifield;

%% Define run sequence
params.initialBlankDurS  = 12; % s
params.finalBlankDurS    = 12; % s
params.blankDurS         = 12; % s
framesOfInitialBlank     = params.initialBlankDurS/ifi; % nr frames
defaultFramesOfInterBlockBlank  = params.blankDurS/ifi; % nr frames
framesOfFinalBlank       = params.finalBlankDurS/ifi;  % nr frames

totalNrUniqueRuns = 2; %3;
nrHemis           = 2;
dummyTimePoints   = 35000; % arbitrary nr of time frames, just make sure its long enough to capture a single run
seqRun            = zeros(totalNrUniqueRuns,nrHemis,dummyTimePoints);
blockStartStop    = cell(1,totalNrUniqueRuns);
blockStartStopWithBlank = blockStartStop;

simNr = [6,12;6,12];

%% Loop over variations
count = 1;
for nn = 1:length(params.locations)
    for ll = 1:length(params.durations.simSingleFrameDurS)
        for conditionOrder = [1,2] % 1 = seq, 2 = sim, 3 = sim, same nr of transients as seq
            clear sequence images images_seq;
            p = params;
            p.locations   = squeeze(params.locations{nn}(nn,:,:))'; % stim size, [x,y], nr
            p.locationIdx = params.locationIdx{nn};
            p.durations.simSingleFrameDurS = params.durations.simSingleFrameDurS(ll);
            p.durations.seqSingleFrameDurS = params.durations.seqSingleFrameDurS(ll);
            p.nrSeqRepetitions    = params.nrSeqRepetitions(ll);
            p.totalTrialsPerBlock = params.totalTrialsPerBlock{nn}(ll);
            p.nrOfBlocks = params.nrOfBlocks{nn}(conditionOrder,ll);
            for bb = 1:p.nrOfBlocks
                sequence(bb,:,:) = ...
                    generateSequence(p.locations, ifi, p.durations, ...
                    p.nrSeqRepetitions, conditionOrder, p.totalTrialsPerBlock, hemifield);
            end
            sequenceAll{conditionOrder,count} = sequence +((nn-1)*6); clear sequence; % we add six because seq with 4 squares has blank+5 unique images
            conditionOrderAll{conditionOrder,count} = count*ones(1,p.nrOfBlocks);   % 1 = Seq, 2 = Sim
        end
        count = count+1;
    end
end

figure(21); clf; ncols = 2; 
nrows = size(sequenceAll,2); 
for ii = 1:size(sequenceAll,2)
    subplot(nrows,ncols,ii); hold all;
    plot(squeeze(sequenceAll{1,ii}(1,:,:))'); title(sprintf('SEQ %i', ii)); drawnow;
    xlim([0 size(sequenceAll{1,ii}(1,:,:),3)])
    
    subplot(nrows,ncols,ii+4); hold all;
    plot(squeeze(sequenceAll{2,ii}(1,:,:))'); title(sprintf('SIM %i', ii)); drawnow;
    xlim([0 size(sequenceAll{1,ii}(1,:,:),3)])
end

conditionOrderSelect = conditionOrderAll;

%% Concatenate SIM and SEQ conditions such that 1-6 SEQ / 7-10 SIM
allConditionBlocks = [];
for jj = 1:size(conditionOrderSelect,1)
    for ii = 1:size(conditionOrderSelect,2)
        allConditionBlocks = [allConditionBlocks, conditionOrderSelect{jj,ii} + (jj-1)*size(conditionOrderSelect,2)];
    end
end

% Get nr of repetition blocks
nrBlocksPerCond = histc(allConditionBlocks,1:max(unique(allConditionBlocks)));

% Shuffle all blocks
conditionOrderAllShuffled = Shuffle(allConditionBlocks(:));
% conditionOrderAllShuffled = allConditionBlocks(:);

% Concatenate sequences so they match 1-12 order
sequenceAll2 = {sequenceAll{1,:,:},sequenceAll{2,:,:}};
% sequenceAll2 = {sequenceAll{1,:},sequenceAll{2,:},sequenceAll{3,:}};

% Add derivatives to params struct so we can save it with images
params.conditionOrderAllShuffled = conditionOrderAllShuffled;
params.sequenceAll = sequenceAll2;
params.nrBlocksPerCond = nrBlocksPerCond;

params.sequenceAllTruncated = sequenceAll2;
for ii = 1:length(sequenceAll2), len(ii) = (size(sequenceAll2{ii},3)); end

% When to split giant run?
splitPoint = (sum((len+defaultFramesOfInterBlockBlank).*nrBlocksPerCond) - 2*(params.initialBlankDurS/ifi))/totalNrUniqueRuns;

%% Build unique run sequence

% reset counters
timeCounter = 1;
blockCounter = 1;
runNr = 1;
buildSequence = 1;
framesPerTR = (TR/1000)/ifi;
seqBlock = zeros(totalNrUniqueRuns,dummyTimePoints);
while buildSequence    
    
    for ii = 1:length(conditionOrderAllShuffled)
        cond = conditionOrderAllShuffled(ii);
        condIdx = sum(conditionOrderAllShuffled(1:ii)==cond);
        
        % Add stimulus
        stimTimeFrames = timeCounter:(timeCounter+(size(sequenceAll2{cond},3))-1);
        seqRun(runNr,:,stimTimeFrames) = sequenceAll2{cond}(condIdx,:,:);
        seqBlock(runNr,stimTimeFrames) = cond*ones(1,length(stimTimeFrames));
        
        blockStartStop{runNr}(blockCounter,:) = [cond,ifi*stimTimeFrames(1),ifi*(stimTimeFrames(end)+1)];
        
        % update counter
        timeCounter = timeCounter+size(sequenceAll2{cond},3);
        
        % Check if counter is at TR integer
        timeFromfullTR = mod(timeCounter,framesPerTR);
        if timeFromfullTR>0
            if timeFromfullTR>(0.75*framesPerTR)
                framesToAdd = framesPerTR-timeFromfullTR;
                framesOfInterBlockBlank = defaultFramesOfInterBlockBlank+framesToAdd;
            elseif timeFromfullTR<=(0.75*framesPerTR)
                framesOfInterBlockBlank = defaultFramesOfInterBlockBlank-timeFromfullTR;
            end
        else
            framesOfInterBlockBlank = defaultFramesOfInterBlockBlank;
        end
        
        % add 8 s of blank after block presentation
        blankTimePoints = timeCounter:(timeCounter+framesOfInterBlockBlank-1);
        seqRun(runNr,:,blankTimePoints) = ones(1,2,framesOfInterBlockBlank);
        seqBlock(runNr,blankTimePoints) = zeros(1,length(blankTimePoints));
        
        blockStartStopWithBlank{runNr}(blockCounter,:) = [cond,stimTimeFrames(1),blankTimePoints(end)+1];

        % update counter
        timeCounter = timeCounter+framesOfInterBlockBlank;
        blockCounter = blockCounter+1;
        if ii<length(conditionOrderAllShuffled)-1
            nextTrialLen = len(conditionOrderAllShuffled(ii+1));
        end
        if ii<length(conditionOrderAllShuffled)
            if timeCounter>splitPoint
                runNr = runNr+1;
                blockCounter = 1;
                timeCounter = 1;
            elseif (nextTrialLen > defaultFramesOfInterBlockBlank) && (timeCounter+nextTrialLen) > splitPoint
                runNr = runNr+1;
                blockCounter = 1;
                timeCounter = 1;
            end
        end
    end
    if runNr > totalNrUniqueRuns
        warning('More than %d runs! Will shuffle conditions and try again!\n',totalNrUniqueRuns);
        conditionOrderAllShuffled = Shuffle(conditionOrderAllShuffled');
        runNr = 1; timeCounter = 1;
        seqRun  = zeros(totalNrUniqueRuns,nrHemis,dummyTimePoints);
        seqBlock = zeros(totalNrUniqueRuns,dummyTimePoints);
        blockStartStop = cell(1,totalNrUniqueRuns);
        blockStartStopWithBlank = blockStartStop;
        blockCounter = 1;
        buildSequence = buildSequence+1;
    elseif runNr == totalNrUniqueRuns
        buildSequence = 0;
    end
    % Quit if we can't find a solution
    if buildSequence > 15
        error('Yikes, can''t find a solution! Check splitpoint for runs!\n')
    end
end

params.blockStartStop = blockStartStop;
params.blockStartStopWithBlank = blockStartStopWithBlank;
clear nextTrialLen timeCounter blockCounter blankTimePoints timeFromfullTR stimTimePoints

% set blank 7 to blank 1
seqRun(seqRun==7) = 1;
seqRun(seqRun==0) = 1;

% Remove all possible blanks at the end, while preserving the last stim
% block across runs
for r = 1:size(seqRun,1)
    postBlankFrames(r) = length(seqRun(r,1,:))-(blockStartStop{r}(end,3)/ifi);
end
[truncVal,truncIdx] = min(postBlankFrames);
startTrunc = size(seqRun(truncIdx,1,:),3)-truncVal;
seqRun(:,:,startTrunc:end) = [];
seqBlock(:,startTrunc:end) = []; clear startTrunc truncVal truncIdx

% Now add pre/post run blanks
seqRun = cat(3, ones(size(seqRun,1),size(seqRun,2),framesOfInitialBlank), seqRun, ones(size(seqRun,1),size(seqRun,2),framesOfFinalBlank));
seqBlock = cat(2, zeros(size(seqBlock,1),framesOfInitialBlank), seqBlock, zeros(size(seqBlock,1),framesOfFinalBlank));

% round stimulus to a full second
if mod(size(seqRun,3)/TR,1)>0
    seqRun = cat(3, seqRun, ones(size(seqRun,1),size(seqRun,2),(framesPerTR-mod(size(seqRun,3),framesPerTR))));
    seqBlock = cat(2, seqBlock, zeros(size(seqBlock,1),(framesPerTR-mod(size(seqBlock,2),framesPerTR))));
end

params.seqBlock = seqBlock;

fprintf('Single run is %5.0f frames, equal to %3.0f seconds, or ~ %1.2f min\n', size(seqRun,3), size(seqRun,3)/framesPerTR, size(seqRun,3)/(60*framesPerTR))
for r = 1:length(blockStartStop)
    blockStartStop{r}(:,2) = blockStartStop{r}(:,2)+framesOfInitialBlank*ifi;
    blockStartStop{r}(:,3) = blockStartStop{r}(:,3)+framesOfInitialBlank*ifi;
end

%% Debug figure; check if runs are roughly the same length
figure(97); clf; set(gcf,'Position',[73,148,1773,814]);
clf;
for rx = 1:size(seqRun,1)
    subplot(size(seqRun,1)+1,1,1);
    plot(squeeze(seqRun(rx,1,:))); hold on;
    title('All runs','FontSize',10);     ylabel('condition nr'); box off;
    set(gca, 'TickLength', [0.005,0.005]); set(gca, 'FontSize', 10)
    xlim([0 size(seqRun,3)])
    
    subplot(size(seqRun,1)+1,1,1+rx);
    plot(squeeze(seqRun(rx,1,:))); hold on;
    plot(blockStartStop{rx}(:,2)/ifi,ones(1,length(blockStartStop{rx}(:,1))),'gx');
    plot(blockStartStop{rx}(:,3)/ifi,ones(1,length(blockStartStop{rx}(:,1))),'rx');
    title(sprintf('Run %d',rx), 'FontSize',10); box off;
    set(gca, 'TickLength', [0.005,0.005])
    ylabel('condition nr'); set(gca, 'FontSize', 10)
    xlim([0 size(seqRun,3)])
end
xlabel('time (s)'); set(gca, 'FontSize', 10)
% Save debug figure
if saveFigs
    fname = sprintf('runOverview_simseq_%s_1deg2_%s_v%d', ...
        trialType, hemifield, version);
    print(fullfile(saveFolder, fname), '-dpng');
    savefig(fullfile(saveFolder, [fname '.fig']))
end

%% Debug figure; Check if all conditions are represented with approx equal TRs
allFrames = size(seqRun,1)*size(seqRun,3);
totalStimSeconds = zeros(size(seqRun,1),length(sequenceAll2));
totalStimSecondsWithBlanks = totalStimSeconds;

for ii = 1:size(seqRun,1)
    dur = blockStartStop{ii}(:,3)-blockStartStop{ii}(:,2);
    theseConds = blockStartStop{ii}(:,1);
    for dd = 1:length(dur)
        totalStimSeconds(ii,theseConds(dd)) = totalStimSeconds(ii,theseConds(dd))+dur(dd);
    end
    
    dur = blockStartStopWithBlank{ii}(:,3)-blockStartStopWithBlank{ii}(:,2);
    theseConds = blockStartStopWithBlank{ii}(:,1);
    for dd = 1:length(dur)
        totalStimSecondsWithBlanks(ii,theseConds(dd)) = totalStimSecondsWithBlanks(ii,theseConds(dd))+dur(dd);
    end
end

edges1 = cumsum(sum(totalStimSecondsWithBlanks,1));
edges2 = cumsum(sum(totalStimSeconds,1));
figure(96); clf;  set(gcf, 'Position', [194,263,1520,699]);
cmap = turbo(length(sequenceAll2));
for ed = length(edges1):-1:1
    ax1 = subplot(1,2,1); hold on;
    bar(1,edges1(ed),'Stacked', 'FaceColor', cmap(ed,:), 'barWidth',0.75)
    ax2 = subplot(1,2,2); hold on;
    bar(1,edges2(ed),'Stacked', 'FaceColor', cmap(ed,:), 'barWidth',0.75)
end
subplot(121); xlim([0.5 1.5])
plot([0.5,1.5],[1,1].*allFrames,'k--')
title('Condition distribution with blanks')
subplot(122); ylim([0 1400]); xlim([0.5 1.5])
title('Condition distribution only stim frames')

h = findobj(gca);
legend(h([2:1:length(sequenceAll2)+1]),labels, 'Location','EastOutside', 'FontSize',9); legend boxoff;

% Save debug figure
if saveFigs
    fname = sprintf('stimCondTRs_simseq_%s_1deg2_%s_v%d', ...
        trialType, hemifield, version);
    print(fullfile(saveFolder, fname), '-dpng');
    savefig(fullfile(saveFolder, [fname '.fig']))
end
clear edges1 edges2 totalStimFramesWithBlanks totalStimFrames theseConds dur zeroIdx allFrames

%% Create stimulus images

for nn = 1:length(params.locations)
    p = params;
    p.locations   = params.locations{nn};
    p.locationIdx = params.locationIdx{nn};
    
    %centerX and y of the screen, in pixels
    ctrX = screenWidthPx/2; % px
    ctrY = ctrX; % px
    
    nLocs = size(p.locations);
    if length(nLocs)==2, nLocs(3) = 1; end
    
    % Store images at each position (+ all stim + blank)
    if strcmp(hemifield,'bothf') % 8 + 3 = 2* 4 indiv stim + 2 with all + 1 blank
        images{nn} = zeros(screenHeightPx, screenWidthPx, (nLocs(1)*nLocs(3))+3);
    else % 4+2 = 4 individual stim, 1 with all stim, 1 blank
        images{nn} = zeros(screenHeightPx, screenWidthPx, (nLocs(1)*nLocs(3))+2);
    end
    
    
    count = 0;
    for ll = 1:length(params.stimAreaDeg2)
        count = count + 2;
        %stimulus dimensions in pixels
        stimWidthPx  = params.stimWidthDeg(ll)*pxPerDeg; % px
        stimHeightPx = params.stimHeightDeg(ll)*pxPerDeg; %px
        rect = [0 0 stimHeightPx stimWidthPx];
        
        
        %         for mi = 1:nLocs(2)
        for ni = 1:nLocs(3)
            %loop through the number of simulations, each with a different 
            % number of unique stimus locations
            xyLoc = squeeze(p.locations(ll,:,ni));
            
            %horizontal and vertical center positions, in pixels, relative to left
            %edge of screen
            curCtrX = ctrX + (xyLoc(1)*pxPerDeg);
            curCtrY = ctrY - (xyLoc(2)*pxPerDeg);
            
            % Select pixels for stimulus
            newRect = CenterRectOnPoint(rect,curCtrX,curCtrY);
            
            startX = floor(newRect(1));
            endX   = floor(newRect(3));
            
            startY = floor(newRect(2));
            endY   = floor(newRect(4));
            
            lY = (startY:endY);
            lX = (startX:endX);
                
            if length(lY)>length(lX)
                lY = lY(1:length(lX));
            elseif length(lX)>length(lY)
                lX = lX(1:length(lY));
            end
                        
            % Insert in image array
            images{nn}(lY, lX, count) = 1;
            
            % concatenate in stim in 6th column for all images sequentially
            images{nn}(lY, lX, simNr(nn,ll)) = 1;
            
            count = count+1;
        end
        
    end
    %     end
    
    if verbose %% debug: plot figure
        figure; clf; set(gca, 'Position', [144 37 1166 925])
        nrows = 5;
        ncols = ceil(size(images{nn},3)/nrows);
        for im = 1:size(images{nn},3)
            subplot(nrows, ncols,im);
            imagesc(images{nn}(:,:,im)); hold on; title(im)
            axis image; axis off; colormap gray;
            set(gca, 'FontSize', 15);
        end
        if saveFigs
            fname = sprintf('stimUniqueFrames%d_simseq_%s_%d_%ddeg2_%s_v%d', ...
                nn, trialType, params.stimAreaDeg2(1),params.stimAreaDeg2(2), hemifield, version);
            print(fullfile(saveFolder, fname), '-dpng');
            savefig(fullfile(saveFolder, [fname '.fig']))
        end
    end
end

%% Create sequence
allIm(1,:,:,:) = images{1}(:,:,1:12);
allIm(2,:,:,:) = images{2}(:,:,1:12);
clear images;

% loop through stimuli
for run = 1:size(seqRun,1)
    sequence = squeeze(seqRun(run,:,:));
    images_seq = NaN(size(allIm,2),size(allIm,3),size(sequence,2));
    
    for si = 1:length(sequence(1,:))
        idx = sequence(1,si);
        images_seq(:,:,si) = squeeze(allIm(1,:,:,idx));
        
        if strcmp(hemifield,'bothf') && (size(sequence,1) > 1)
            idx2 = sequence(2,si);
            images_seq(:,:,si) = ...
                [images_seq(:,:,si) + squeeze(allIm(2,:,:,idx2))];
        end
    end
    
    % Convert stim values to logical
    images = logical(images_seq); clear images_seq
    
    %% Plot sequence
    t = 0:1/ifi:((length(sequence)-1)/ifi);
    figure(102); clf; set(gcf, 'Color', 'w', 'Position', [107,233,1696,729])
    subplot(211)
    plot(t,sequence(1,:)-1,'k', 'LineWidth',2); hold all;
    if size(sequence,1) >1
        plot(t,sequence(2,:)-1,'r:', 'LineWidth',2);
        legend('Upper Right', 'Lower Left'); legend boxoff
    end
    xlabel('time (s)')
    ylabel('condition nr')
    set(gca, 'FontSize', 20, 'TickDir','out', 'yTick',[0:1:max(unique(sequence))]); box off;
    xlim([0,max(t)]);
    
    subplot(212); cla;
    plot(t,squeeze(any(any(images,1),2)),'k', 'LineWidth',2);
    xlabel('time (s)')
    ylabel('pixel contrast')
    set(gca, 'FontSize', 20, 'TickDir','out', 'YTick', [0 1]); box off;
    xlim([0,max(t)])
    
    
    % Save debug figures, images and sequence
    fname = sprintf('stim_simseq_run%d_v%d',run, version);
    
    if saveFigs
        print(fullfile(saveFolder, fname), '-dpng');
        savefig(fullfile(saveFolder, [fname '.fig']))
    end
    
    if saveData
        save(fullfile(saveFolder, [fname '.mat']), 'images', 'sequence', 'params', '-v7.3');
    end
    
    %% Write par file
    condNames = labels;
    colors = turbo(length(sequenceAll2)); % white + rainbow
    par = struct();
    % Define initial blank period
    par.onset(1) = 0;
    par.code(1)  = 0;
    par.cond{1}  = 'blank';
    par.color{1} = [1 1 1];
    count = 2;
    for ii = 1:size(blockStartStop{run},1)
        if  ii==1
            par.onset(count) = blockStartStop{run}(ii,2)+params.initialBlankDurS-2; % onset in seconds
        else
            par.onset(count) = blockStartStop{run}(ii,2)+par.onset(2);
        end
        par.onset(count) = blockStartStop{run}(ii,2); % onset in seconds
        par.code(count)  = blockStartStop{run}(ii,1); % condition code (blank=0, seq4 = 1-3, seq16 = 4-6, sim4=7,9, sim16=10,12)
        par.cond{count}  = condNames{par.code(count)};
        par.color{count} = colors(par.code(count),:);
        count = count+1;
        
        par.onset(count) = blockStartStop{run}(ii,3); % onset in seconds
        par.code(count)  = 0; % condition code (blank=0, seq4 = 1-3, seq16 = 4-6, sim4=7,9, sim16=10,12)
        par.cond{count}  = 'blank';
        par.color{count} = [1 1 1];
        count = count+1;
    end
    
    outFile = sprintf('simseq_%s_run%d_v%d.par', trialType, run, version);
    
    fidout = fopen(fullfile(saveFolder,outFile),'w');
    for pp=1:length(par.onset)
        fprintf(fidout,'%d \t %d \t', par.onset(pp), par.code(pp));
        fprintf(fidout,'%s \t', par.cond{pp});
        fprintf(fidout,'%1.3f %1.3f %1.3f \n', par.color{pp});
    end
    fclose(fidout);
    fclose('all');
    fprintf('Wrote parfile %s successfully.\n', fullfile(saveFolder,outFile));
end


%% debug: show move and plot sequence
if showMovie
    for fr=1:10:size(images,3)
        figure(1); clf; set(gcf,'Position',[740   486   598   476])
        imagesc(images(:,:,fr)); title(sprintf('Frame %d - %3.3fs',fr, fr*ifi));
        colormap gray; axis square;
        pause(0.0001);
    end
end
if createVideo
    % Prepare the new file.
    vidObj = VideoWriter(fullfile(saveFolder,[fname '.mp4']));
    open(vidObj);
    % Create an animation.
    figure(101); clf;
    axis equal tight off
    for k = 1:size(images,3)
        imagesc(images(:,:,k)); title(sprintf('Frame %d - %3.3fs',k, k*ifi));
        axis square; colormap gray; drawnow;
        %             pause(0.01)
        % Write each frame to the file.
        currFrame = getframe(101);
        writeVideo(vidObj,currFrame);
    end
    % Close the file.
    close(vidObj);
end

end




