function [] = concatParfilesUniqueRuns(projectDir, subjnr, sesNr)
% script to concatenate unique runs in pilot 3. 

% subjnr     = 6; %9;
% sesNr      = 1;
% projectDir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal/';
% projectDir = '/share/kalanit/biac2/kgs/projects/simseqPRF/';

cd(fullfile(projectDir,'experiments/simseq/'))

pilotNr    = 3;
pths       = getSubjectPaths(projectDir,subjnr,sesNr);
version    = pths.expversionNr;

dataFolder = fullfile(pths.dataDirSimSeq);
saveFolder = fullfile(dataFolder,sprintf('subj%02d',subjnr), pths.session,'Stimuli','parfiles');
scanFolder = fullfile(dataFolder,sprintf('subj%02d',subjnr), pths.session, 'Gray', 'MotionComp_RefScan1');

% tmp = load(fullfile(scanFolder,'TSeries','Scan1','tSeries1.mat'));
% TRsPerRun = size(tmp.tSeries,1);

if subjnr == 3
    StartRun2inSeconds = 325;
elseif ismember(subjnr,[2,11,12])
    StartRun2inSeconds = 323; 
else
    StartRun2inSeconds = 324; % i.e. next run begins at 325 TRs (or 324 seconds counting from 0)
end


% Get conditions, colors
if pilotNr == 2
    runs  = pths.runOrder(1:3);
    stimFileNameBase = 'stim_simseq_varyTransients_x_Area_x_StimDur_2loc_3timings_1deg2_bothf_mid_fullrun';
    condNames = {'SEQ1-4xorig','SEQ2-4xtrans','SEQ3-4xdur','SEQ4-16xorig','SEQ5-16xtrans','SEQ6-16xdur',...
        'SIM1-4xorig','','SIM3-4xdur','SIM4-16xorig','','SIM6-16xdur'};
    outFile = sprintf('simseq_varyConditions_run123_v%d.par', version);
    runAddTRs = 4; % we remove 6 TRs (6 s) after scan start to avoid inhomogeities

elseif pilotNr == 3
    stimFileNameBase = 'stim_simseq_varyStimDur_x_Area4_2x2design_1deg2_bothf_mid_fullrun';
    runs  = pths.runOrder(1:2); 
    condNames = {'SEQ1-4x200ms-4deg2','SEQ2-4x1000ms-4deg2',...
                    'SEQ3-4x200ms-16deg2','SEQ4-4x1000ms-16deg2',...
                    'SIM1-4x200ms-4deg2','SIM2-4x1000ms-4deg2',...
                    'SIM3-4x200ms-16deg2','SIM4-4x1000ms-16deg2'};
    outFile = sprintf('simseq_varyStimDur_x_Area_run12_v%d_corrected2.par', version);
    runRemoveTRs = 2; % we remove 2x 8 s (first 6 s count down + 2s blank,
    nrBlankInit  = 12; % TRs = sec
    nrBlankEnd   = 12; % TRs = sec
    countDown    = 6;  % TRs = sec
                   % then another 8s blank) after scan start to avoid inhomogeities
%     deleteEndTR = 1; % 0 for s3 (stimulus sequence happened to be longer)
%     deleteFrontTR = 8;
end
colors = turbo(length(condNames)); % rainbow

% Create par struct, loop over runs to fill in
par = struct();
count = 1;
for r = 1:length(runs)
    run = runs(r);
    stimFileName = sprintf('%s%d_v%d.mat',...
        stimFileNameBase, run, version);
    load(fullfile(dataFolder,sprintf('subj%02d', subjnr),pths.session,'Stimuli',...
        stimFileName));
    runLength(r) = size(images,3)/60;
    
    blockStartStop = params.blockStartStop;
    condNum = blockStartStop{r}(:,1); % integers
    onsets  = blockStartStop{r}(:,2); % in seconds
    offsets = blockStartStop{r}(:,3); % in seconds
    onsets(1) = 0; % Reset to 0 as start (not a fraction)
    
    onsets = onsets + nrBlankInit - runRemoveTRs;
    offsets = offsets + nrBlankInit - runRemoveTRs;
    
    % Define initial blank period
    if run == runs(1)
        newRunStart = 0;
        par.onset(count) = 0; 
        par.code(count)  = 0;
        par.cond{count}  = 'blank';
        par.color{count} = [1 1 1];
    else
        newRunStart = StartRun2inSeconds; 
        par.onset(count) = newRunStart; % last and first blank are 12 TRs, but we remove first 6 TRs (6 s) after scan start to avoid inhomogeities
        par.code(count)  = 0;
        par.cond{count}  = 'blank';
        par.color{count} = [1 1 1];
    end
     count = count+1;
    
    
    for ii = 1:length(onsets)
%         if ii==1
%             par.onset(count) = newRunStart;
%         else
            par.onset(count) = newRunStart + onsets(ii);
%         end
        par.code(count)  = condNum(ii); % condition code (blank=0, seq4 = 1-3, seq16 = 4-6, sim4=7,9, sim16=10,12)
        par.cond{count}  = condNames{par.code(count)};
        par.color{count} = colors(par.code(count),:);
        count = count+1;
        
        % Define blank at offset of block
        par.onset(count) = newRunStart + offsets(ii); % onset in seconds
        par.code(count)  = 0; % condition code (blank=0, seq4 = 1-3, seq16 = 4-6, sim4=7,9, sim16=10,12)
        par.cond{count}  = 'blank';
        par.color{count} = [1 1 1];
        count = count+1;
    end
end


fidout = fopen(fullfile(saveFolder,outFile),'w');
for pp=1:length(par.onset)
    fprintf(fidout,'%d \t %d \t', par.onset(pp), par.code(pp));
    fprintf(fidout,'%s \t', par.cond{pp});
    fprintf(fidout,'%1.3f %3.3f %1.3f \n', par.color{pp});
end
fclose(fidout);
fclose('all');
fprintf('Wrote parfile %s successfully.\n', fullfile(saveFolder,outFile));
