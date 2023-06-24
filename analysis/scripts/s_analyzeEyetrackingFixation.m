%% s_analyzeEyetrackingFixation.m
% EDF files need to be convert to ascii files, using the EDF2ASC software
% provided by SR-Research for Eyelink1000:
% https://www.sr-research.com/support/thread-13.html
% You need to make an account to access this support page!
%
% This analysis requires the I2MC toolbox: https://github.com/royhessels/I2MC
% Hessels, R.S., Niehorster, D.C., Kemner, C., & Hooge, I.T.C. (2017). 
% Noise-robust fixation detection in eye-movement data - Identification by 2-means clustering (I2MC). 
% Behavior Research Methods, 49(5): 1802--1823. doi: 10.3758/s13428-016-0822-1
%
% This script will preprocess eye gaze data in the following steps:
% * Remove blinks (automatically labeled by Eyelink algorithm)
% * 


%% Define params
subjnrs     = [1,2,3,7,12]; % S10, S11 and S13 do not have eyetracking data
% Data from S8 and S9 are too noisy to analyze
projectDir  = '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal/';
plotFigures = false;
saveFigs    = true;   % Save figures or not
saveData    = true;   % Save preprocessed data or not

% Add necessary paths:
eyeCodeFolder = fullfile(projectDir,'experiments','simseq','code','eyetracking');
addpath(genpath(fullfile(eyeCodeFolder)));
addpath(genpath(fullfile('~/matlab/git/toolboxes/I2MC')));

% Screen params
params.pixScreen = [1920 1080]; % pixels
params.cmScreen  = [38 20];  % cm
params.viewingDistance = 25; % cm;
ppd = pi* params.pixScreen(2) / (atan((params.cmScreen(2)*0.5)/params.viewingDistance)) / 360;

% Eyetracking analysis params
params.blinkWinSec  = [0.1,0.1]; % sec % Blink window is pretty conservative right now. (Could also be less conservative [0.2 0.35])
params.fs           = 1000; % hz
params.blankPeriod  = 12*params.fs; % ms
params.preExpPeriod = 30*params.fs; % ms
params.k            = 3;     % k-means
params.velthresh    = 400;   % deg/s (more than this will be seen as noise)
params.absthresh    = 10;    % deg (more than 10 deg radius will fall beyond display and considered noise)
params.stdthresh    = 2.5;     % final std threshold for outliers
params.vThres       = 6;     % deg/s velocity threshold for microsaccade detection
params.msMinDur     = 6;     % ms - min duration for microsaccade detection
params.pdd          = ppd;

% To compute velocity
if params.fs >= 500
    params.velocityType = 3;
else
    params.velocityType = 2;
end

% I2CM clustering algorithm options
opt.xres          = params.pixScreen(1); % maximum value of horizontal resolution in pixels
opt.yres          = params.pixScreen(2); % maximum value of vertical resolution in pixels
opt.missingx      = NaN; % missing value for horizontal position in eye-tracking data (example data uses -xres). used throughout functions as signal for data loss
opt.missingy      = NaN; % missing value for vertical position in eye-tracking data (example data uses -yres). used throughout functions as signal for data loss
opt.freq          = params.fs; % sampling frequency of data (check that this value matches with values actually obtained from measurement!)
opt.scrSz         = params.cmScreen; % screen size in cm
opt.disttoscreen  = params.viewingDistance; % distance to screen in cm.
opt.pixperdeg     = ppd;
params.I2MC       = opt;

% Colormaps and conditions
darkred    = [0.48,0.016, 0.011];
darkblue   = [0.267, 0.33, 0.767];

% Allocate space to store preprocessed data
allSubjectData = struct('xyDataDeg',[],'tData',[],'timeSeqIdx',[],...
    'timeSimIdx',[],'blockTimeStart',[],'blockTimeStop',[],...
    'blockCondition',[],'params',[], 'raw',[]);

makeprettyfigures;

if saveFigs
    ses =  getSessionNrMainExp(subjnrs(1));
    pths = getSubjectPaths(projectDir,subjnrs(1),ses);
    saveFigDir = fullfile(fileparts(fileparts(pths.figureDir)), 'average', 'figures', ...
        'main_eyetrackingPerformance_SimSeq');
    if ~exist(saveFigDir, 'dir'); mkdir(saveFigDir); end
end

%% Loop over subjects
for s = 1:length(subjnrs)
    fprintf('\nLoading subject S%02d..',subjnrs(s))
    ses =  getSessionNrMainExp(subjnrs(s));
    pths = getSubjectPaths(projectDir,subjnrs(s),ses);
    edffile = dir(fullfile(pths.dataDirSimSeq,pths.subjID,pths.session, 'eye','*.asc'));
    
    for r = 1:length(edffile)
        if plotFigures
            fH1 = figure(109+r); clf; set(gcf,'Position', [0  0 1980 780],'Color', 'w');
        end
        
        % ===============================================
        % ================ Load run data ================
        % ===============================================
        % Columns are [Time by X gaze by Y gaze by pupil size by X velocity
        % by Y Velocity by X res by Y Res by input (?)]
        % more info at https://www.sr-research.com/support/thread-7675.html
        asc     = read_eyelink_asc(fullfile(edffile(r).folder,edffile(r).name));
        asc.dat = asc.dat';
        
        if sum(asc.dat(:,2)~=0)>1
            tDataRaw            = asc.dat(:,1);
            xyPosRaw            = [asc.dat(:,2),asc.dat(:,3)];
            xyVelRaw            = [asc.dat(:,5),asc.dat(:,6)];
            blinkData.startTime = table2array(asc.eblink(:,2));
            blinkData.endTime   = table2array(asc.eblink(:,3));
            saccEL              = table2array(asc.ssacc(:,2)); % NB: EyeLink saccades are usually not very robust, at least not for microsaccades
            
            % Reset time to 0
            tData               = tDataRaw - tDataRaw(1);
            blinkData.startTime = blinkData.startTime - tDataRaw(1);
            blinkData.endTime   = blinkData.endTime - tDataRaw(1);
            saccEL              = saccEL - tDataRaw(1);
            
            % Recompute velocity for microsaccade analysis
            xyVel = vecvel(xyPosRaw,params.fs,params.velocityType); % pixels / s
            
            % Remove first 30s before experiment onset:
            xyPos = xyPosRaw;
            xyPos(1:params.preExpPeriod,:) = NaN;
            xyVel(1:params.preExpPeriod,:) = NaN;
            blinkData.startTime = blinkData.startTime(blinkData.startTime>params.preExpPeriod);
            blinkData.endTime = blinkData.endTime(blinkData.endTime>params.preExpPeriod);
            saccEL = saccEL(saccEL>params.preExpPeriod);
            
            % Convert to deg visual angle
            nSamp  = length(~isnan(xyPos));
            xyDataDeg  = xyPos./ppd;
            xyVelDeg   = xyVel./ppd;
            
            % ===============================================
            % ================ Get meta data ================
            % ===============================================
            
            % Define the start and end time of experiment, as well as start of
            % experiment blocks
            jj=1; blockTimeStart = []; blockCondition = [];
            for ii = 1:size(asc.msg,1)
                if regexp(asc.msg{ii,3},'.\RECORD*')
                    startIdx =  ii;
                    startTime = str2num(asc.msg{startIdx,2});
                end
                
                if regexp(asc.msg{ii,3},'BLOCK_ID*')
                    blockID =  ii;
                    blockTimeStart(jj) = str2num(asc.msg{blockID,2});
                    blockCondition(jj) = str2num(asc.msg{blockID,3}(end-2:end-1)) - 48;
                    jj = jj +1;
                end
                
                if regexp(asc.msg{ii,3},'TRIAL_END')
                    endIdx =  ii;
                    endTime = str2num(asc.msg{endIdx,2});
                    break
                end
            end
            blockTimeStop  = [blockTimeStart(2:end)-params.blankPeriod, blockTimeStart(end)+params.blankPeriod];
            blockTimeStart = blockTimeStart - tDataRaw(1);
            blockTimeStop  = blockTimeStop - tDataRaw(1);
            
            % Remove blinks and no data periods
            xyPosDegNoBlink = xyDataDeg;
            for iBlink = 1:length(blinkData.startTime)
                blinkInd = tData >= blinkData.startTime(iBlink)-params.blinkWinSec(1)*params.fs & ...
                    tData <= blinkData.endTime(iBlink)+params.blinkWinSec(2)*params.fs;
                xyPosDegNoBlink(blinkInd,:) = NaN;
            end
            
            %% I2MC
            data = [];
            data.time   = tData;
            data.left.X = xyPosDegNoBlink(:,1);
            data.left.Y = xyPosDegNoBlink(:,2);
            [fix,ndata] = I2MCfunc(data,opt);
            
            % PLOT RESULTS
            if plotFigures
                if saveFigs
                    fName = sprintf('S%d_I2CM_output_run%d', s,r);
                    plotResults(data,fix,fullfile(saveFigDir,fName),[opt.xres opt.yres]);
                else
                    plotResults(data,fix,[],[opt.xres opt.yres]);
                end
            end
            
            % Fixation centering
            clear mnsX mnsY
            tmpx = xyPosDegNoBlink(:,1); tmpy = xyPosDegNoBlink(:,2);
            if size(xyPosDegNoBlink,1)>10 && length(fix.xpos)>params.k % need at least 10 datapoints
                
                % compute the means for found fixation periods (2 or 3)
                mnsXidx = kmeans(fix.xpos',params.k);
                mnsYidx = kmeans(fix.ypos',params.k);
                for mm = 1:params.k
                    mnsX(mm) = mean(fix.xpos(mnsXidx==mm));
                    mnsY(mm) = mean(fix.ypos(mnsYidx==mm));
                    
                    % Get timepoints for each mean and subtract
                    curFixPeriodsX = find(mnsXidx==mm);
                    for tt = curFixPeriodsX'
                        fixTimePoints = fix.startT(tt):fix.endT(tt);
                        fixTimePoints(fixTimePoints>length(tmpx))=[];
                        tmpx(fixTimePoints) =  tmpx(fixTimePoints) - mnsX(mm);
                    end
                    
                    curFixPeriodsY = find(mnsYidx==mm);
                    for tt = curFixPeriodsY'
                        fixTimePoints = fix.startT(tt):fix.endT(tt);
                        fixTimePoints(fixTimePoints>length(tmpy))=[];
                        tmpy(fixTimePoints)  =  tmpy(fixTimePoints) - mnsY(mm);
                    end
                end
                
                xyDataDegCentered = [tmpx,tmpy];
                
                %% Remove noise outliers
                outlrsVel = find(abs(xyVelDeg(:,1))>params.velthresh | abs(xyVelDeg(:,2))>params.velthresh);
                outlrsVel = outlrsVel+[-2:2];
                
                outlrsAbs = find(abs(xyDataDegCentered(:,1))>params.absthresh | abs(xyDataDegCentered(:,2))>params.absthresh);
                outlrsAbs = outlrsAbs+[-2:2];
                
                outlrsAll = unique([outlrsVel(:);outlrsAbs(:)]);
                outlrsAll = sort(outlrsAll);
                outlrsAll(outlrsAll>length(tData))=[];
                
                xyPosDegFiltered = xyDataDegCentered;
                xyPosDegFiltered(outlrsAll,:) = NaN;
                
                %% Final STD thresh
                stdThresh = params.stdthresh*mean(std(xyPosDegFiltered,[],1,'omitnan'));
                outlrsSTD = find(std(xyPosDegFiltered,[],2,'omitnan') > stdThresh);
                outlrsSTD = outlrsSTD+[-2:2];
                outlrsSTD = unique(outlrsSTD(:));
                outlrsSTD(outlrsSTD>length(tData))=[];
                xyPosDegFiltered(outlrsSTD,:) = NaN;
                outliers = cat(1,outlrsAll,outlrsSTD);
                
                %% =================================================
                %  =============== Detect saccades =================
                %  =================================================
                %
                % Eyelink also detects saccades, but the algorithm isn't that reliable.
                % This is only if you actually want to do analyses with the saccades (e.g.,
                % look at how microsaccades are modulated by attention). If just for
                % monitoring goodness of fixation, you can probably just rely on the
                % Eyelink detected saccades (eyd.saccades)
                % Outputs are:
                % * onset of saccade,
                % * end of saccade,
                % * peak velocity of saccade (vpeak)
                % * horizontal component     (dx)
                % * vertical component       (dy)
                % * horizontal amplitude     (dX)
                % * vertical amplitude       (dY)
                %         [sacRaw1,radius1] = microsacc(xyPos,xyVel,vThres,msMinDur);
                %         [sacRaw2,radius2] = microsacc(xyPosDegNoBlink,xyVelRaw,vThres,msMinDur);
                
                [sacRaw,radius] = microsacc(xyPosDegFiltered,xyVelRaw,params.vThres,params.msMinDur);
                % remove the ones that occurr closely together (overshoot)
                numSacs = size(sacRaw,1);
                minInterSamples = ceil(0.01*1000);
                interSac = sacRaw(2:end,1)- sacRaw(1:end-1,2);
                sac = sacRaw([1; find(interSac > minInterSamples)+1],:);
                fprintf('%d rejected for close spacing\n', numSacs - size(sac,1));
                fprintf('%d saccades detected\n', size(sac,1));
                fprintf('number of saccades detected: %d, detectRadius: [%0.3f %0.3f]\n', ...
                    size(sac,1), radius(1), radius(2));
                
                mcrsaccIndx = [];
                for sc = 1:size(sac,1)
                    mcrsaccIndx = cat(2, mcrsaccIndx, sac(sc,1):sac(sc,2));
                end
                params.mcrsacc.sacRaw = sacRaw;
                params.mcrsacc.radius = radius;
                params.mcrsacc.numSacs = numSacs;
                params.mcrsacc.interSac = interSac;
                params.mcrsacc.sac = sac;
                params.mcrsacc.mcrsaccIndx = mcrsaccIndx;
                
                % invert y so 0 is bottom left
                params.invertY = true;
                xyPosDegFiltered(:,2) = -xyPosDegFiltered(:,2);
                
                % Get epochs for SIM & SEQ blocks
                blockCmap = [];
                for bb = 1:length(blockCondition)
                    if blockCondition(bb)>5
                        cmap = darkblue;
                    else
                        cmap = darkred;
                    end
                    blockCmap = cat(1,blockCmap,cmap);
                end
                
                timeSimIdx = zeros(size(tData)); timeSeqIdx = zeros(size(tData));
                for jj = 1:length(blockTimeStart)
                    tempIdx = ((tData > blockTimeStart(jj)) & tData < blockTimeStop(jj));
                    if blockCondition(jj) < 5
                        timeSeqIdx(tempIdx) = 1;
                    elseif blockCondition(jj) >= 5
                        timeSimIdx(tempIdx) = 1;
                    end
                end
                timeSeqIdx = logical(timeSeqIdx);
                timeSimIdx = logical(timeSimIdx);
                
                %% Store data
                allSubjectData(s,r).raw = asc;
                allSubjectData(s,r).params = params;
                allSubjectData(s,r).xyDataDeg = xyPosDegFiltered;
                allSubjectData(s,r).tData = tData;
                allSubjectData(s,r).timeSeqIdx = timeSeqIdx;
                allSubjectData(s,r).timeSimIdx = timeSimIdx;
                allSubjectData(s,r).blockTimeStart = blockTimeStart;
                allSubjectData(s,r).blockTimeStop = blockTimeStop;
                allSubjectData(s,r).blockCondition = blockCondition;
                
                
                if plotFigures
                    %% Plot X COORDINATES
                    figure(fH0); clf
                    sgtitle(sprintf('Eyetracking Subject S%02d, run %d',subjnrs(s), r))
                    
                    % RAW data
                    subplot(3,4,[1:3]); cla
                    mx = max(xyDataDeg(:));
                    mn = min(xyDataDeg(:));
                    plot(tData/1000,xyDataDeg(:,1), 'Color', [.7 .7 .7]); hold all;
                    plot(tData/1000,xyDataDeg(:,2), 'Color', [0 0 0]); hold all;
                    plot(mcrsaccIndx/1000, mean([mn,mx]).*ones(1,length(mcrsaccIndx)),'cx')
                    
                    % Plot per stimulus condition
                    for bb = 1:length(blockCondition)
                        plot([blockTimeStart(bb), blockTimeStart(bb)]./1000,[mn mx], ':','color',blockCmap(bb,:));
                        plot([blockTimeStop(bb), blockTimeStop(bb)]./1000,[mn mx], ':','color','k');
                    end
                    grid on;
                    ylim([mn-10, mx+10]);
                    set(gca, 'FontSize', 15); box off;
                    ylabel('Raw Position (deg)')
                    
                    % Recentered data
                    subplot(3,4,[5:7]); cla;
                    mx = max(xyDataDegCentered(:));
                    mn = min(xyDataDegCentered(:));
                    plot(tData/1000,xyDataDegCentered(:,1), 'Color', [.7 .7 .7]); hold all;
                    plot(tData/1000,xyDataDegCentered(:,2), 'Color', [0 0 0]);
                    plot(tData(outliers)/1000,xyDataDegCentered(outliers,1), 'rx');
                    plot(tData(outliers)/1000,xyDataDegCentered(outliers,2), 'rx');
                    
                    % Plot per stimulus condition
                    for bb = 1:length(blockCondition)
                        plot([blockTimeStart(bb), blockTimeStart(bb)]./1000,[mn mx], ':','color',blockCmap(bb,:));
                        plot([blockTimeStop(bb), blockTimeStop(bb)]./1000,[mn mx], ':','color','k');
                    end
                    grid on;
                    ylim(1.1.*[mn, mx]);
                    ylabel('Center corrected position (deg)')
                    set(gca, 'FontSize', 15); box off;
                    
                    subplot(3,4,[9:11]); cla;
                    mx = max(xyPosDegFiltered(:));
                    mn = min(xyPosDegFiltered(:));
                    plot(tData/1000,xyPosDegFiltered(:,1), 'Color', [0.7 0.7 0.7]); hold all;
                    plot(tData/1000,xyPosDegFiltered(:,2), 'Color', [0 0 0]);
                    
                    % Plot per stimulus condition
                    for bb = 1:length(blockCondition)
                        plot([blockTimeStart(bb), blockTimeStart(bb)]./1000,[mn mx], ':','color',blockCmap(bb,:));
                        plot([blockTimeStop(bb), blockTimeStop(bb)]./1000,[mn mx], ':','color','k');
                    end
                    grid on;
                    ylim([mn-1, mx+1]);
                    xlabel('Time (s)');
                    ylabel('Outlier corrected position')
                    set(gca, 'FontSize', 15); box off;
                    
                    % Plot XY ON GRID
                    subplot(3,4,8);
                    plot(xyDataDegCentered(:,1),xyDataDegCentered(:,2),'o-','Color',[.7 .7 .7]); hold on;
                    scatter(mean(xyDataDegCentered(:,1),'omitnan'),mean(xyDataDegCentered(:,2),'omitnan'),'r+')
                    axis square; grid on;
                    xlabel('X (pixels)');  ylabel('Y (pixels)');
                    set(gca, 'FontSize', 15);
                    xlim([-20 20]); ylim([-20 20])
                    
                    subplot(3,4,12); cla;
                    plot(xyPosDegFiltered(:,1),xyPosDegFiltered(:,2),'o-','Color',[.7 .7 .7]);   hold all;
                    scatter(mean(xyPosDegFiltered(:,1),'omitnan'),mean(xyPosDegFiltered(:,2),'omitnan'),'r+');
                    scatter(xyPosDegFiltered(timeSeqIdx,1),xyPosDegFiltered(timeSeqIdx,2),10,darkred, 'filled');
                    scatter(xyPosDegFiltered(timeSimIdx,1),xyPosDegFiltered(timeSimIdx,2),10,darkblue, 'filled');
                    axis square; grid on;
                    xlabel('X (deg)');  ylabel('Y (deg)');
                    set(gca, 'FontSize', 15);
                    xlim([-10 10]); ylim([-10 10])
                    
                    if saveFigs
                        fName = sprintf('S%d_eyetracking_run%d_preproc',subjnrs(s),r);
                        saveas(fH0,fullfile(saveFigDir,[fName '.png']))
                        print(fH0,fullfile(saveFigDir,fName),'-depsc');
                    end
                    %% PLOT CLEAN DATA
                    figure(fH1)
                    subplot(1,4,[1:3]);
                    mx = max(xyPosDegFiltered(:));
                    mn = min(xyPosDegFiltered(:));
                    plot(tData/1000,xyPosDegFiltered(:,1), 'Color', [0.7 0.7 0.7]); hold all;
                    plot(tData/1000,xyPosDegFiltered(:,2), 'Color', [0 0 0]);
                    
                    % Plot per stimulus condition
                    for bb = 1:length(blockCondition)
                        plot([blockTimeStart(bb), blockTimeStart(bb)]./1000,[-10 10], ':','color',blockCmap(bb,:));
                        plot([blockTimeStop(bb), blockTimeStop(bb)]./1000,[-10 10], ':','color','k');
                    end
                    grid on;
                    ylim(1.1*[mn mx]);
                    xlabel('Time (s)'); ylabel('X or Y (deg)');
                    title(sprintf('Subject S%02d: Fixation traces clean, run %d',subjnrs(s),r))
                    set(gca, 'FontSize', 15); box off;
                    
                    % Plot XY ON GRID
                    subplot(1,4,4); cla;  hold all;
                    p1 = patch([0.59, 0.59, 4.59, 4.59],[0.59, 4.59, 4.59, 0.59],'k'); hold all;
                    p1.FaceAlpha = 0.1;
                    p2 = patch([2.59, 2.59, 4.59, 4.59],[2.59, 4.59, 4.59, 2.59],'k');
                    p2.FaceAlpha = 0.2;
                    p3 = patch(-1.*[0.59, 0.59, 4.59, 4.59],-1.*[0.59, 4.59, 4.59, 0.59],'k');
                    p3.FaceAlpha = 0.1;
                    p4 = patch(-1.*[2.59, 2.59, 4.59, 4.59],-1.*[2.59, 4.59, 4.59, 2.59],'k');
                    p4.FaceAlpha = 0.2;
                    
                    plot(xyPosDegFiltered(:,1),xyPosDegFiltered(:,2),'o-','Color',[.7 .7 .7]);
                    scatter(mean(xyPosDegFiltered(:,1),'omitnan'),mean(xyPosDegFiltered(:,2),'omitnan'),'r+');
                    scatter(xyPosDegFiltered(timeSeqIdx,1),xyPosDegFiltered(timeSeqIdx,2),10,darkred, 'filled');
                    scatter(xyPosDegFiltered(timeSimIdx,1),xyPosDegFiltered(timeSimIdx,2),10,darkblue, 'filled');
                    
                    axis square; grid on;
                    xlabel('X (deg)');  ylabel('Y (deg)');
                    set(gca, 'FontSize', 15);
                    xlim([-5 5]); ylim([-5 5])
                    
                    if saveFigs
                        fName = sprintf('S%d_eyetracking_run%d_clean',subjnrs(s),r);
                        saveas(fH1,fullfile(saveFigDir,[fName '.png']))
                        print(fH1,fullfile(saveFigDir,fName),'-depsc');
                    end
                end
            end
        end
    end
end

if saveData
    saveDataDir = fullfile(fileparts(fileparts(pths.figureDir)), 'eyetracking');
    if ~exist(saveDataDir, 'dir'); mkdir(saveDataDir); end
    save(fullfile(saveDataDir,[sprintf('S%i_',subjnrs) 'eyetrackingDataClean']),'allSubjectData')
end


return



