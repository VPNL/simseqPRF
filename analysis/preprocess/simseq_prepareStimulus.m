function simseq_prepareStimulus(sessionDir)
% sessionDir = vistasoft session in simseq project folder.
% for example, '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal/experiments/simseq/data/subj01/session1'

% Get number of runs/files
d = dir(fullfile(sessionDir,'simseqPRF*.mat'));

% Define params
fs = 1000; % ms resolution
screenHeightPx  = 101;  % number of pixels of fit screen vertically
screenWidthPx   = 101;  % number of pixels of fit screen horizontally
screenHeightDeg = 24;   % screen diameter degrees of visual angle
pixPerDeg       = screenHeightPx/screenHeightDeg;

minTimeGap = 0.033; % sec (force a 33ms minimum time gap between stimuli, 
                 % needed for temporal model to see them as separate stim.

% Loop over files
for run = 1:length(d)
    
    load(fullfile(d(run).folder, d(run).name));
    
    % Extract params used for grid and stimulus locations
    stimWidthPx     = params.squareSizeDeg*pixPerDeg;
    stimHeightPx    = stimWidthPx;
    
    % Create grid from experiment
    gr = params.gridSpaceDeg*([1:params.grid]-ceil(params.grid/2));
    [X,Y]=(meshgrid(gr));X=X';Y=Y';
%     if scan.coil == 16 % flip left/right (put this in params next time!).
% %         X = -X;
%     end
    XY = [X(:), Y(:)];

    % Get stim locations
    locations = unique(cell2mat(flip.pos)); % degrees from central fixation (deg)
    
    % Store images at each position
    images = zeros(screenHeightPx, screenWidthPx, length(locations)+1); % add one with all stimuli and one blank
    
    % Loop through the number of unique stimulus locations, skip blank
    for ni = 2:length(locations)
        
        stimX = XY(locations(ni),1)*pixPerDeg;
        stimY = XY(locations(ni),2)*pixPerDeg;
        
        curCtrX = ceil(screenHeightPx/2) + stimX;
        curCtrY = ceil(screenWidthPx/2) + stimY;
        
        % Select pixels for stimulus
        startX = round(curCtrX - floor(stimWidthPx/2));
        endX   = round(curCtrX + ceil(stimWidthPx/2) - 1);
        startY = round(curCtrY - floor(stimHeightPx/2));
        endY   = round(curCtrY + ceil(stimHeightPx/2) -1);
        
        % Insert in image array
        images(startY:endY, startX:endX, ni) = 1;
        
        % concatenate in last column for all images sequentially
        images(startY:endY, startX:endX, length(locations)+1) = 1;
    end
    
    % convert to single integer array
    images = single(images);
    
    % debug figure to plot binarized stimuli
    figure(1); clf; set(gca, 'Position', [543,608,1235,322])
    for ii = 1:size(images,3)
        subplot(1, size(images,3), ii);
        imagesc(images(:,:,ii)); hold on; title(ii)
        axis square; colormap gray;
        set(gca, 'FontSize', 15)
        if ii == 1
            ylabel('Pixels (#)')
        elseif ii == 3
            xlabel('Pixels (#)')
        end
    end
    
    
    %% Create sequence
    sequence = zeros(length(trial),trial(1).dur*fs);
    cutoff = ((trial(1).dur-minTimeGap)*fs);
    for ii = 1:length(trial)
        if (ii == 5) || (ii == 0)
            sequence(ii,:) = repmat(trial(ii).cond, [1 trial(ii).dur*fs]);
        elseif (ii >= 1) || (ii<=4)
            sequence(ii,1:cutoff) = repmat(trial(ii).cond, [1 (trial(ii).dur-minTimeGap)*fs]);
            sequence(ii,(cutoff+1):size(sequence,2)) = zeros(1,(minTimeGap*fs));
        end
    end
    
    % rectify and add 1 for blank (1), single stim (2:5), all four stim (6)
    sequenceFinal = sequence';
    sequenceFinal = sequenceFinal(:)+1;
    
    sequence = sequenceFinal;
    save(fullfile(sessionDir, 'Stimuli',sprintf('images_and_sequenceWith33msGap_run%d.mat',run)), 'images', 'sequence')
       
    % make truncated sequence for testing
    sequence = sequenceFinal(1:44000);
    save(fullfile(sessionDir, 'Stimuli',sprintf('images_and_sequenceTruncatedWith33msGap_run%d.mat',run)), 'images', 'sequence')
    
%     % Make separate sequences for two conditions
%     sequenceSeq = sequenceFinal;
%     sequenceSim = sequenceFinal;
%     sequenceSeq(sequence==6) = 1;
%     sequenceSim(((sequence>1)&sequence<6)) = 1;
%     
%     sequence = sequenceSim;
%     save(fullfile(sessionDir, 'Stimuli',sprintf('images_and_sequenceSimOnly_run%d.mat',run)), 'images', 'sequence')
%     
%     sequence = sequenceSeq;
%     save(fullfile(sessionDir, 'Stimuli',sprintf('images_and_sequenceSeqOnly_run%d.mat',run)), 'images', 'sequence')

end

