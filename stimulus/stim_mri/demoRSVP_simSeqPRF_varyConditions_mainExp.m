function demoRSVP_simSeqPRF_varyConditions_mainExp(subject)
% Demo of 2 blocks of single experimental run of the
% simultaneous-sequential spatiotemporal pRF experiment
%
% INPUTS
% [subject]       : (string) ID of subject (default to 'test')

%%% at the scanner:
% set trigger box to these options (per JG):
% change modes / manual configure / HHSC 1X5D / USB / HID NAR NO5

if ~exist('subject','var')
    scan.subj = 'test'; else, scan.subj  = subject;
end

scan.coil = 16;
scan.runNum = 1;
scan.versionNr = 5;
scan.atScanner = 0;
scan.useEyelink = 0;
scan.offset = [0,0];
scan.loadTrials = [];

flipLR = 0;
flipUD = 0;

input('Hit enter to proceed.');

%%% PTB preferences
Screen('Preference', 'SkipSyncTests', 1);
Screen('Preference','TextRenderer', 0);

%%%% Screen
params.backgroundColor = [128 128 128];    % color RGB
params.textColor       = [255 255 255];    % color RGB
scan.screenWidth       = 38; %104;              % in cm; % laptop=27.5, office=43, CNI = 104cm at both 16ch and 32ch
scan.screenHeight      = 20; %104;              % in cm; % laptop=27.5, office=43, CNI = 104cm at both 16ch and 32ch                                     % in cm; CNI = 270-273cm or 272 at 32ch, 265 at 16ch;
scan.viewingDist       = 25;                 % OLD: 265 cm;


%%%% open screen for refresh rate and window size
screen = max(Screen('Screens'));
[win, rect] = Screen('OpenWindow',screen,params.backgroundColor);

%%% move xc and yc if needed
xc = rect(3)/2 + scan.offset(1); % rect and center, with the flixibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2 + scan.offset(2);
params.xc = xc; params.yc = yc;

[w,h] = RectSize(Screen('TextBounds',win,'Loading textures.. Starting soon!'));
DrawFormattedText(win, 'Loading textures.. Starting soon!', xc-w/2,yc-h/2, params.textColor(1,:),[], flipLR, flipUD);
Screen(win, 'Flip', 0);

%%%% keyboard checks
[keyboardIndices, productNames, ~] = GetKeyboardIndices; %#ok<ASGLU>
laptop = max(keyboardIndices); % don't differentiate scanner/laptop
scanner = laptop;

%%%% set-up rand
rand('twister', sum(100*clock));
scan.rand = rand;

%%%% files and things
scan.root = pwd;
scan.date = datestr(now,30);

%%%% images used
im2deg = load(fullfile(simseqRootPath,'stimulus', 'stim_mri','cartoon','cartoonChopped_pilot3','cartoonChopped_2deg_FOV38by29cm_rms0.1.mat'), 'images');
im4deg = load(fullfile(simseqRootPath,'stimulus', 'stim_mri','cartoon','cartoonChopped_pilot3','cartoonChopped_4deg_FOV38by29cm_rms0.1.mat'), 'images');
trig = load(fullfile(simseqRootPath,'stimulus', 'stim_mri',sprintf('stim_simseq_varyStimDur_x_Area4_2x2design_1deg2_bothf_mid_fullrun%d_v%d.mat',scan.runNum, scan.versionNr)));
trialSetsDir = fullfile(simseqRootPath,'stimulus', 'stim_mri','trialSets_VaryConditions_mainExp');
if ~exist(trialSetsDir, 'dir'), mkdir(trialSetsDir); end

%%%% load calibration file
load(fullfile(simseqRootPath,'stimulus', 'stim_mri','calibration','gamma.mat'), 'gamma');
scan.gammaTable = gamma; clear gamma;

%%%% image order
ims.numIDs = size(im2deg.images,4);
ims.order = Shuffle(1:ims.numIDs);

%%%% timing
params.countDown      = 5;%10;                 % countdown at the beginning of the scan - 10
params.TRlength       = 1;                 % in seconds

%%%% layout
params.grid           = 15; %11;           % number of grid elements
params.gridSpaceDeg   = 1; %sqrt(1.9);     % center to center spacing of the grid
params.squareSizeDeg  = [2,4]; %1;         % size of each square (width/height)
params.fixRadDeg      = .2;                % in degrees, the size of the biggest white dot in the fixation
params.cuePix         = 2;                 % in pixels, the thickness of the cue ring. now draws *inside* fixRadDeg, so that the cue is not made bigger by increasing this param

%%%% scale the stims for the screen
scan.ppd = pi* rect(4) / (atan((scan.screenHeight*0.5)/scan.viewingDist)) / 360;
scan.squareSize  = (params.squareSizeDeg*scan.ppd);                 % in degrees, the size of our squares
scan.fixRad      = (params.fixRadDeg*scan.ppd);
scan.littleFix   = (scan.fixRad*.25);

%%% define grid of stimulus center locations
gridDeg = params.gridSpaceDeg*([1:params.grid]-ceil(params.grid/2));
[X,Y] = meshgrid(gridDeg);
X = X'; Y = Y';
if flipLR, X = -X; end
params.locationsDeg = [X(:) Y(:)]; clear X Y

%%% sample centers
params.locationsIdx   = [43, 183]; % upper right [5,-5], lower left [5,-5]
params.mainCentersDeg = params.locationsDeg(params.locationsIdx,:); % LH 5 deg, RH 5 deg
gapFromCenter         = 0.41; % in deg -- note: same as for first pilot

scan.rects = cell(1,length(params.squareSizeDeg));
for sz = 1:length(params.squareSizeDeg)
    halfSquareSize = params.squareSizeDeg(sz)/2;
    for ll = 1:length(params.locationsIdx)
        % coords: 1 x [x1,y1,x2,y2] 
        upperright = CenterRectOnPoint([0 0 params.squareSizeDeg(sz) params.squareSizeDeg(sz)], ...
            params.mainCentersDeg(ll,1) + (gapFromCenter + halfSquareSize), ...
            params.mainCentersDeg(ll,2) - (gapFromCenter + halfSquareSize));
        
        upperleft = CenterRectOnPoint([0 0 params.squareSizeDeg(sz) params.squareSizeDeg(sz)], ...
            params.mainCentersDeg(ll,1) - (gapFromCenter + halfSquareSize), ... x
            params.mainCentersDeg(ll,2) - (gapFromCenter + halfSquareSize));  % y
        
        lowerright = CenterRectOnPoint([0 0 params.squareSizeDeg(sz) params.squareSizeDeg(sz)], ...
            params.mainCentersDeg(ll,1) + (gapFromCenter + halfSquareSize), ... x
            params.mainCentersDeg(ll,2) + (gapFromCenter + halfSquareSize));  % y
        
        lowerleft =  CenterRectOnPoint([0 0 params.squareSizeDeg(sz) params.squareSizeDeg(sz)], ...
            params.mainCentersDeg(ll,1) - (gapFromCenter + halfSquareSize), ... x
            params.mainCentersDeg(ll,2) + (gapFromCenter + halfSquareSize));  % y
        
        % Rects field has the following dimensions:
        % 4 locations, 4 coords [x1,y1,x2,y2], 2 hemis 
        locPix = cat(2,lowerright',upperright',lowerleft',upperleft') * scan.ppd;
        locPixCentered =  locPix + repmat([xc;yc;xc;yc],1,4);
        scan.rects{sz}(:,:,ll) = locPixCentered;
    end
end
params.locationsPix = params.locationsDeg*scan.ppd;
scan.centers = [params.locationsPix(:,1)+xc params.locationsPix(:,2)+yc];

if flipUD
    scan.centers = repmat([rect(3) rect(4)],length(scan.centers),1)-scan.centers;
    yc = rect(4)-yc; xc = rect(3)-xc;
    for sz = 1:length(params.squareSizeDeg)
        for ll = 1:length(params.locationsIdx)
            scan.rects{sz}(:,:,ll)
        end
    end
end

%%%% tasks
params.refreshRateHz = (1/60); %Screen('GetFlipInterval',win); % 60 Hz
params.framesPerSec  = (1/params.refreshRateHz);
params.taskBack      = 1;                        % 1-back
params.taskFreq      = 40;                       % every 40th frame, so new letter every (40/60) seconds (frequency of new RSVP letter)
params.task          = ['Please fixate.' num2str(params.taskBack) '-Back on Letters: Press button when letter repeats.']; % Text shown at the beginning of the run
params.taskColor     = [0 0 0; 1 1 1].*255;      % black & white
params.taskProb      = (1/9);                    % proportion of trials with nback ID targets
params.respWindow    = 50;   % number of flips in which we'll count a correct response
params.str           = ['A' 'S' 'D' 'F' 'G' 'H' 'J' 'K', 'B', 'P'];% alphabet([],1);%           % entire alphabet, capitalized
params.imageRateHz   = params.refreshRateHz;  % Hz

%%%%%%%%%%%%%%%
% trial setup %
%%%%%%%%%%%%%%%
subsetOfFrames           = 200:3500; %400:2400;

params.numConds          = 8; % 2 sizes, 2 timings, sim/seq (prev pilot 2 had 12; % 2 sizes, 3 timings, sim & seq conditions)
params.condNames         = {'SEQ-4deg2-200ms','SEQ-4deg2-1000ms', ...
                            'SEQ-16deg2-200ms','SEQ-16deg2-1000ms', ...
                            'SIM-4deg2-200ms','SIM-4deg2-1000ms', ...
                            'SIM-16deg2-200ms','SIM-16deg2-1000ms'};
params.blockNr           = trig.params.seqBlock(scan.runNum,:);
params.seq               = trig.sequence(:,subsetOfFrames);
params.seqTiming         = 0:(params.refreshRateHz):(params.refreshRateHz*length(params.seq(1,:))); % seconds
params.seqTiming         = params.seqTiming(1:end-1);

params.pos               = {trig.params.locationIdx{1},...% LH
                            trig.params.locationIdx{2}}; % RH

params.loc               = scan.rects; 

% Define response keys
responseKeys = zeros(1,256);
responseKeys(KbName({'1!';'2@';'3#';'4$';'5%';'6^'}))=1;

% Check if frames per images is a round nr, to avoid display delays
params.framesPerImage    = params.refreshRateHz/params.imageRateHz;
if mod(params.framesPerImage,1)~=0
    warning('Frames per image is not a multiple of screen refreshrate! This can cause inaccurate timing')
end

%%%% timing optimization
scan.screenRefreshRate  = params.refreshRateHz;
scan.imageRate          = params.imageRateHz;
slack                   = params.refreshRateHz/2;
scan.time               = params.countDown+params.seqTiming(end);   % time in seconds
scan.allFlips           = 0:params.refreshRateHz:scan.time; %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RSVP targ 1 = repeat (1-back)
% RSVP targs happen on blank trials as well
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(scan.loadTrials)
    load(fullfile(trialSetsDir, sprintf('%s.mat',scan.loadTrials)),'trial');
else
    trial = createTrialSet_varyArea(params, ims, scan, trialSetsDir);
end

%% indexing stimulus properties for each screenflip - now incl countdown
if params.countDown > 0
      flip.trials      = [repmat(-1*Expand(fliplr(1:params.countDown),params.framesPerSec,1),2,1) params.seq];
else
    flip.trials      = params.seq;
end
flip.rsvptarg  = zeros(1, length(flip.trials));
flip.rsvpID    = flip.rsvptarg;
flip.rsvpColor = zeros(length(flip.trials),3);
flip.IDs       = cell(1, length(flip.trials));
flip.pos       = flip.IDs;

expStart = params.countDown*params.framesPerSec;

for n = 1:length(trial)
    flip.rsvptarg(n+expStart)    = trial(n).rsvptarg;
    flip.IDs{n+expStart}         = trial(n).IDs;
    flip.rsvpID(n+expStart)      = trial(n).rsvpID;
    flip.pos{n+expStart}         = trial(n).pos;
    flip.loc{n+expStart}         = trial(n).location;
    flip.rsvpColor(n+expStart,:) = trial(n).rsvpColor;
    flip.block(n+expStart,:)     = trial(n).block;
end
%%
HideCursor;
Priority(9);
Screen(win, 'TextSize', 18);

Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%%%%% flipping images upside down for the 16ch headcoil
squares2deg = cell(1,size(im2deg.images,4));
squares4deg = cell(1,size(im4deg.images,4));
for n = 1:size(im2deg.images,4)
    
    if flipUD == 1 && flipLR == 0
        squares2deg{n} = flipud(im2deg.images(:,:,:,n));
        squares4deg{n} = flipud(im4deg.images(:,:,:,n));
    end
    
    if flipLR == 1 && flipUD == 0
        squares2deg{n} = fliplr(im2deg.images(:,:,:,n));
        squares4deg{n} = fliplr(im4deg.images(:,:,:,n));
    end
    
    if flipLR == 1 && flipUD == 1
        squares2deg{n} = flipud(fliplr(im2deg.images(:,:,:,n))); %#ok<FLUDLR>
        squares4deg{n} = flipud(fliplr(im4deg.images(:,:,:,n))); %#ok<FLUDLR>
    end
    
    if ((flipLR == 0) && (flipUD == 0))
        squares2deg{n} = im2deg.images(:,:,:,n);
        squares4deg{n} = im4deg.images(:,:,:,n);
    end
end

clear im2deg im4deg; clear trig;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       experiment                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% start recording the response
scan.pressFlip = [];
scan.pressedTime = [];
scan.pressedKeyCode = [];

% Get block onsets to print in edf file
blockOnsets = find(diff(flip.block)>0)+1;

%%% Make the textures
fprintf('\nMaking textures.')
fTex = cell(1,length(scan.allFlips)-1);
for n = 1:length(scan.allFlips)
    if any(flip.IDs{n} > 0) % draw multiple squares
        trialImages = flip.IDs{n};
        for ii = 1:length(trialImages)
            if any(ismember(flip.trials(:,n),(2:6)))
                f = squares2deg{trialImages(ii)};
                fTex{n}(ii) = Screen('MakeTexture',win,f);
            elseif any(ismember(flip.trials(:,n),(8:12)))
                f = squares4deg{trialImages(ii)};
                fTex{n}(ii) = Screen('MakeTexture',win,f);
            end
        end
    end
    if mod(n,length(scan.allFlips)/10)==0
        fprintf('.')
    end
end
fprintf('\nDone!\n')

% Save some memory
clear trialImages squares2deg squares4deg f

%%% initial window - wait for backtick
startText = 'Waiting for trigger, or press "t" key to start scan.';
Screen('DrawText',win, startText, 10,10,params.textColor);

[w,h] = RectSize(Screen('TextBounds',win,params.task));
DrawFormattedText(win, params.task, xc-(w/2),yc-(h/2), params.taskColor(1,:),[], flipLR, flipUD);
Screen(win, 'Flip', 0);

%%%% This will record keypresses
KbQueueCreate(scanner,responseKeys); % set to scanner
KbQueueStart(); % 

% Load Gamma table
Screen('LoadNormalizedGammaTable', win, scan.gammaTable); % 3rd input = loadOnNextFlip flag

% Release cpu
WaitSecs(0.01);
% KbQueueFlush(scanner)

%%%%%%% START drawing textures before task flipping
for n = 1:length(scan.allFlips)
    if any(flip.IDs{n} > 0) % draw multiple squares
        Screen('DrawTextures', win, fTex{n}',[],flip.loc{n});
        Screen('Close', fTex{n});
    end
    
    Screen('FillOval', win,[255 255 255],[xc-1 yc-1 xc+1 yc+1]); % small fixation dot
    tr = flip.trials(1,n);
    if flip.trials(:,n)<0 % countdown
        [w,h] = RectSize(Screen('TextBounds',win,num2str(-1*tr)));
        DrawFormattedText(win,num2str(-1*tr),xc-(w/2), yc-(h/2)+15, [0 0 0],[],flipLR,flipUD); % yc+4
        
    elseif flip.rsvpID(n) > 0  % draw the rsvp task
        % RSVP letter/digit
        [w,h] = RectSize(Screen('TextBounds',win,params.str(flip.rsvpID(n))));
        DrawFormattedText(win,params.str(flip.rsvpID(n)), xc-(w/2), yc-(h/2)+15, flip.rsvpColor(n,:),[],flipLR,flipUD); %xc-6, yc+4,
        
        %     elseif flip.rsvpID(n) == 0
        %         Screen('FillOval', win,[180 0 0 255*.5], [xc-3 yc-3 xc+3 yc+3]); % small fixation dot
        %         Screen('FillOval', win,[255 255 255], [xc-1 yc-1 xc+1 yc+1]); % small fixation dot
    end

    %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if n == 1
        KbTriggerWait(KbName('t'), laptop); % check laptop to get trigger. %checkKeyPress('t');%
        
        [VBLT, scan.startRun, endFlipTime, missed] = Screen(win, 'Flip', 0);
        scan.stimOnsetTime(n) = scan.startRun;
    else
        [VBLT, scan.stimOnsetTime(n), endFlipTime, missed] = Screen(win, 'Flip', scan.startRun + scan.allFlips(n) - slack);
    end
    
    % Record screen timing into 'response' struct
    response.VBLTimestamp(n) = VBLT;
    response.stimOnsetTime(n) = scan.stimOnsetTime(n);
    response.endFlipTime(n) = endFlipTime;
    response.missed(n) = missed;
    
    % listen for response  - correct if you respond to previous 3 letters
%     [pressed,pressedTime,pressedKeyCode] = KbCheck();
    [pressed, pressedKeyCode, pressedTime] = KbQueueCheck(scanner);
    if pressed == 1
        scan.pressFlip = [scan.pressFlip n];
        scan.pressedTime = [scan.pressedTime pressedTime(pressedTime>0)];
        scan.pressedKeyCode = [scan.pressedKeyCode find(pressedKeyCode)];
        KbQueueFlush(scanner);
    end
    KbQueueFlush(scanner); 
%     KbQueueStop(scanner);   
end

%%%% to show the very last flip screen for its 250ms
[response.VBLTimestamp(n+1), response.stimOnsetTime(n+1), ...
    response.endFlipTime(n+1), response.missed(n+1)] = ...
    Screen(win, 'Flip', scan.startRun + scan.allFlips(length(scan.allFlips)) - slack);

%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%

scan.runTime = GetSecs - scan.startRun;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute performance
perf = struct('task',params.task,'hitTr',[],'falseFlip',[]);
% for now, this only deals with hits, no false alarms...
for p = scan.pressFlip((diff(scan.pressFlip(2:end))>1))+1 % scan.flip(2:end) if KbQueueCheck
    targetRange = p-params.respWindow:p-1;
    if sum(flip.rsvptarg(targetRange)>0) % we only care about task 1 in this version
        perf.hitTr = [perf.hitTr targetRange(find(flip.rsvptarg(targetRange)))]; % this is a hit
    else
        perf.falseFlip = [perf.falseFlip p];
    end
end

% Put behavior to struct perf
perf.targTrials = find(flip.rsvptarg);
perf.hitRate    = length(perf.hitTr)/length(perf.targTrials);
perf.falseTr    = length(unique(perf.falseFlip));
perf.missTr     = setdiff(perf.targTrials,perf.hitTr);
perf.falseRate  = length(perf.falseTr)/length(perf.targTrials);

% Print performance on screen
feedback = sampleFeedbackText(perf.hitRate);
perfText = sprintf('%s Hit rate is %2.1f%%, nr of false alarms: %2.0f', feedback, perf.hitRate*100,perf.falseTr);
[w,h] = RectSize(Screen('TextBounds',win,perfText));
DrawFormattedText(win, perfText, xc-w/2,yc-h/2, [0 0 0],[], flipLR, flipUD);
Screen(win, 'Flip', 0);

pause(5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save data
curPath = pwd;
cd(fullfile(simseqRootPath));
saveDir = fullfile(simseqRootPath, 'data',['S' scan.subj]);
if ~exist(fullfile(saveDir),'dir'); mkdir(saveDir); end
save(fullfile(saveDir, ['simseqPRF_' scan.subj '_run' num2str(scan.runNum) '_v' num2str(versionNr) '_' scan.date 'demo.mat']), 'params', 'scan', 'trial', 'flip', 'response','perf', '-v7.3');

%% Clean up
KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');

% Print performance
fprintf('%s: hit rate %.2f%%\n',params.task,100*perf.hitRate);
fprintf('False alarm presses: %2.0f\n',perf.falseTr);

fclose all;
cd(curPath);
