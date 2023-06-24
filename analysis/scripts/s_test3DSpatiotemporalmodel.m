%% s_test3DSpatiotemporalmodel.m
% Script to test different stages of 3ch spatio temporal model
tic
projectDir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal/';
% projectDir = '/Volumes/kalanit/biac2/kgs/projects/spatiotemporal/';
subjnr  = 1;
pilotNr = 8;

if ismember(pilotNr,[6,7,8,9])
    ses    = 1;
else
    if subjnr ~= 2
        ses    = 2;
    else ses = 3;
    end
end
pths   = getSubjectPaths(projectDir, subjnr,ses);
simSeqSessionDir = fullfile(pths.dataDirSimSeq, pths.subjID, pths.session);


% Define params
upsampleStimType    = 'upsampleFrom1Hz';
stimRun             = 1;
hemi                = 'both';%'both';
dt                  = 'Averages';
roi                 = 'artificial';%'VO1';%'artificial';
roiType             = 'stimcorner4_area4sq_eccen5'; %'stimcorner';
roifile             = [roi '_toon_simseq' roiType];
veThresh            = 0.1;
idx_voxel           = [1:4]; %[1:10]; % for artifical pRFs: idx 1 and 2 are small pRFs; 3 and 4 are big pRFs
yl                  = [-0.1,0.4];
cyany              = [0,188,227]./255;
if strcmp(roi,'artificial') || isempty(regexp(roiType,'stimcorner','ONCE'))
    roiSubFolder = [];
else
    roiSubFolder = 'simseqROIs';
end

params.analysis.spatialModel          = 'onegaussianFit'; % Choose: 'onegaussianFit' or 'cssFit'
params.analysis.temporalModel         = '3ch-stLN'; %'1ch-glm';%;%'2ch-exp-sig'; %'3ch-stLN';%'1ch-glm';%'3ch-stLN';%'3ch-stLN';       % Choose: '3ch-stLN', '2ch-exp-sig', '1ch-glm', '1ch-dcts'
params.analysis.normNeuralChan        = true;
params.stim.framePeriod               = 1; % TR (seconds)
params.saveDataFlag                   = true;
params.stim.sparsifyFlag              = false; % leave out zero response in conv(prf,stim)
params.recomputePredictionsFlag       = true;
params.analysis.predFile              = 'tmp.mat';
params.analysis.zeroPadPredNeuralFlag = true;
params.analysis.reluFlag              = true;
params.analysis.useMedianROIexponent  = true;


%% Generate stimulus
clear msStim;
switch pilotNr
    case -1 % use full field on-off
        gap = 0.033; %s
        idx_on_off = [2:(0.25+gap):3]*1000;
        msStim = zeros(101,101,30000);
        msStim(:,:,idx_on_off(1):(idx_on_off(1)+249)) = 1;
        msStim(:,:,idx_on_off(2):(idx_on_off(2)+249)) = 1;
        msStim(:,:,idx_on_off(3):(idx_on_off(3)+249)) = 1;
        msStim(:,:,idx_on_off(4):(idx_on_off(4)+249)) = 1;
        msStim(:,:,15000:15250) = 1;
        msStim = reshape(msStim, 101*101,[]);
        stim.sequence = zeros(1,30000);
        stim.sequence(idx_on_off(1):(idx_on_off(1)+249)) = 1;
        stim.sequence(idx_on_off(2):(idx_on_off(2)+249)) = 1;
        stim.sequence(idx_on_off(3):(idx_on_off(3)+249)) = 1;
        stim.sequence(idx_on_off(4):(idx_on_off(4)+249)) = 1;

        stim.sequence(15000:16000)= 1;
        
    case 1  % small squares
        stimfile = fullfile(simSeqSessionDir, 'Stimuli', 'stim_simseq_varyStimDur_x_Area4_2x2design_1deg2_bothf_mid_fullrun1_v3.mat');
        stim = load(stimfile);
        
    case 2  % big squares
        stimfile = fullfile(simSeqSessionDir, 'Stimuli', 'stim_simseq_varyStimDur_x_Area4_2x2design_1deg2_bothf_mid_fullrun1_v4.mat');
        stim = load(stimfile);
        
    case 3  % small and big squares
        stimfile = fullfile(simSeqSessionDir, 'Stimuli', 'stim_simseq_varyStimDur_x_Area4_2x2design_1deg2_bothf_mid_fullrun1_v2.mat');
        stim = load(stimfile);
        
    case 4  % small, big squares and full field
        stimfile = fullfile(simSeqSessionDir, 'Stimuli', 'stim_simseq_varyStimDur_x_Area4_2x2design_1deg2_bothf_mid_fullrun1_v2.mat');
        stim = load(stimfile);
        
    case 5 % Pilot 3: Vary stim dur (100,200,1000ms) and stim area (big/small squares)
        stimfile = fullfile(simSeqSessionDir, 'Stimuli', sprintf('stim_simseq_varyStimDur_x_Area4_2x2design_1deg2_bothf_mid_fullrun1_v5.mat'));
        stim = load(stimfile);
        fullfieldCond = max(unique(stim.sequence))+1;
        
    case 6 % Pilot 1: SEQ block followed by SIM block repeated many times
        stimfile = fullfile(simSeqSessionDir, 'Stimuli', sprintf('images_and_sequenceWith33msGap_run%d.mat',stimRun));
        stim = load(stimfile);
        xl = [0 50]; % only one block
        
    case 7 % Insub's stRetinotopy experiment
        stRetDir = fullfile(projectDir, 'experiments','stRet','data','subj01','session1');
        stimfile = fullfile(stRetDir, 'Stimuli', 'images_and_params_run01.mat');
        load(stimfile)
        stim.images = images;
        stim.sequence = sequence;
        
     case 8 % Pilot 1: 1 SEQ trial followed by SIM trial
        stimfile = fullfile(simSeqSessionDir, 'Stimuli', sprintf('images_and_sequenceWith33msGap_run%d.mat',stimRun));
        stim = load(stimfile);
        shortSequence = [stim.sequence(10000:13000)', ones(1,10000), stim.sequence(26000:29000)', ones(1,10000)];
        stim.sequence = shortSequence;
        xl = [0 26]; % only one block
        
    case 9 % Pilot 1: 1 SEQ trial followed by SIM trial 1000 ms
        stimfile = fullfile(simSeqSessionDir, 'Stimuli', sprintf('images_and_sequenceWith33msGap_run%d.mat',stimRun));
        stim = load(stimfile);
        shortSequence = ones(1,35000);
        shortSequence(2000:2200) = 6;
        
        shortSequence(16000:16200) = 4;
        shortSequence(16233:16433) = 2;
        shortSequence(16467:16667) = 3;
        shortSequence(16700:16900) = 5;
        
%         shortSequence(2000:3000) = 4;
%         shortSequence(3033:4033) = 2;
%         shortSequence(4066:5066) = 3;
%         shortSequence(5099:6099) = 5;
%         shortSequence(15000:16000) = 6;
        stim.sequence = shortSequence;
        xl = [0 26]; % only one block
end

if pilotNr == -1
    
%     % Get upsampled stimulus sequence
%     trLength = size(msStim,2)/60;
%     xq = 0:0.001:(trLength-(1/60));
%     t  = 0:(1/60):((size(msStim,2)-1)/60);
%     for jj = 1:size(msStim,2)
%         newSeq(:,jj) = interp1(t,stim.sequence(:,jj),xq, 'previous');
%     end
%     [Bnew,Nnew] = RunLength_M(newSeq(:,1));
%     [Bold,Nold] = RunLength_M(stim.sequence(:,1));
%     assert(isequal(Bnew,Bold));
%     
%     repeats = (Nnew-Nold);
%     msStim2 = zeros(101,101,length(newSeq),'single');
%     count = 1;
%     for mm = 1:length(Bold)
%         origIm = find(stim.sequence(:,1)==Bold(mm));
%         origIm = origIm(1);
%         reps = Nnew(mm);
%         msStim2(:,:,count:(count+reps-1)) = repmat(msStim(:,:,origIm),[1 1 reps]);
%         count = count+reps;
%     end
%     
%     msStim = reshape(msStim2,101*101,[]);
    
elseif ismember(pilotNr,[1:5])
    newSeq = [];
    if size(stim.sequence,2) > size(stim.sequence,1)
        stim.sequence = stim.sequence';
    end
    % Get upsampled stimulus sequence
    trLength = length(stim.sequence)/60;
    xq = 0:0.001:(trLength-(1/60));
    t  = 0:1/60:((length(stim.sequence)-1)/60);
    for jj = 1:size(stim.sequence,2)
        newSeq(:,jj) = interp1(t,stim.sequence(:,jj),xq, 'previous');
    end
    
    [Bnew,Nnew] = RunLength_M(newSeq(:,1));
    [Bold,Nold] = RunLength_M(stim.sequence(:,1));
    assert(isequal(Bnew,Bold));
    
    repeats = (Nnew-Nold);
    msStim = zeros(size(stim.images,1),size(stim.images,2),length(newSeq),'single');
    count = 1;
    for mm = 1:length(Bold)
        origIm = find(stim.sequence(:,1)==Bold(mm));
        origIm = origIm(1);
        reps = Nnew(mm);
        msStim(:,:,count:(count+reps-1)) = repmat(stim.images(:,:,origIm),[1 1 reps]);
        count = count+reps;
    end
    
    msStim = reshape(msStim,101*101,[]);
%     if pilotNr==5
%         msStim = msStim(:,1:end);
%     end
elseif ismember(pilotNr,[6,7,8,9])
    % Get upsampled stimulus sequence
     if pilotNr==6
        msStim = zeros(size(stim.images,1),size(stim.images,2),30000);   
     else
        msStim = zeros(size(stim.images,1),size(stim.images,2),length(stim.sequence));
     end
    for ii = 1:length(stim.sequence)
        if stim.sequence(ii)==0
            msStim(:,:,ii) = zeros(101,101);
        else
            msStim(:,:,ii) = stim.images(:,:,stim.sequence(ii));
        end
    end
    msStim = reshape(msStim, size(msStim,1)*size(msStim,2),[]);
    if pilotNr==6
        msStim = msStim(:,1:50000);
    end
end

if pilotNr == 4
    % add fullfield stim at the end
    gap = 0.034; %s
    idx_on_off = [2:(0.1+gap):3]*1000;
    fullf = zeros(101*101,30000);
    for ii = 1:4
        fullf(:,idx_on_off(ii)+[0:(gap*1000)]) = 1;
    end
    fullf(:,15000:19000) = 1;
    msStim = cat(2,msStim, fullf);
end

%% Apply the model

params = getSpatialParams(params,1);
params = getTemporalParams(params);

if strcmp(params.analysis.temporalModel, '3ch-stLN')
    params.analysis.combineNeuralChan     = [1, 2, 2];
else
    params.analysis.combineNeuralChan     = [1:params.analysis.temporal.num_channels];
end

% Get spatial prf params
if strcmp(roi,'artificial')
    pRFParams = simseq_simulatePRFParams(hemi, params.analysis.spatialModel, veThresh);
else
    pRFParams = simseq_loadPRFParams(pths, dt, hemi, roifile, params.analysis.spatialModel, veThresh, [], roiSubFolder);
end

% Add extra params and replace spatial field
if strcmp(hemi,'lh')
	params.analysis.spatial.lh = pRFParams.lh;
elseif strcmp(hemi,'rh')
    params.analysis.spatial.rh = pRFParams.rh;
elseif strcmp(hemi,'both')
    params.analysis.spatial.lh = pRFParams.lh;
    params.analysis.spatial.rh = pRFParams.rh;
end

clear prfParams;
fNames = {'lh','rh'};
if strcmp(params.analysis.spatialModel,'cssFit')
    if strcmp(roi,'artificial')
        tmp = simseq_simulatePRFParams(hemi, 'cssFit', veThresh);
    else
        tmp = simseq_loadPRFParams(pths, dt, hemi, roifile, 'cssFit', veThresh, [],'simseqROIs');
    end
    
    for fn = 1:length(fNames)
        params.analysis.spatial.(fNames{fn}).exponent = tmp.(fNames{fn}).exponent;
    end
    
elseif strcmp(params.analysis.temporalModel,'3ch-stLN')
    if strcmp(roi,'artificial')
        tmp = simseq_simulatePRFParams(hemi, 'cssFit', veThresh);
    else
        tmp = simseq_loadPRFParams(pths, dt, hemi, roifile, 'cssFit', veThresh, [],'simseqROIs');
    end
    params.analysis.temporal.param = rmfield(params.analysis.temporal.param,'exponent');
    for fn = 1:length(fNames)
        params.analysis.temporal.param.exponent.(fNames{fn}) = tmp.(fNames{fn}).exponent;
    end
end

% if params.analysis.useMedianROIexponent
%     if strcmp(params.analysis.temporalModel,'3ch-stLN')
%         % Use average ROI spatial exponent as static temporal exponent
%         medianExp = median([params.analysis.temporal.param.exponent.lh, ...
%             params.analysis.temporal.param.exponent.rh],'omitnan');
%         params.analysis.temporal.param.exponent = medianExp;
%     end
%     if strcmp(params.analysis.spatialModel,'cssFit')
%         % Use average ROI spatial exponent
%         medianExp = median([params.analysis.spatial.exponent.lh, ...
%             params.analysis.spatial.exponent.rh],'omitnan');
%         params.analysis.spatial.exponent.lh = medianExp;
%         params.analysis.spatial.exponent.rh = medianExp;
%     end
% end

%% Run it!
for fn = 1:length(fNames)
    params.analysis.spatial.(fNames{fn}).x0            = params.analysis.spatial.(fNames{fn}).x0(idx_voxel);
    params.analysis.spatial.(fNames{fn}).y0            = params.analysis.spatial.(fNames{fn}).y0(idx_voxel);
    params.analysis.spatial.(fNames{fn}).sigmaMajor    = params.analysis.spatial.(fNames{fn}).sigmaMajor(idx_voxel);
    params.analysis.spatial.(fNames{fn}).sigmaMinor    = params.analysis.spatial.(fNames{fn}).sigmaMinor(idx_voxel);
    params.analysis.spatial.(fNames{fn}).theta         = params.analysis.spatial.(fNames{fn}).theta(idx_voxel);
    params.analysis.spatial.(fNames{fn}).effectiveSize = params.analysis.spatial.(fNames{fn}).effectiveSize(idx_voxel);
    params.analysis.spatial.(fNames{fn}).varexpl       = params.analysis.spatial.(fNames{fn}).varexpl(idx_voxel);
    params.analysis.spatial.(fNames{fn}).veMask        = params.analysis.spatial.(fNames{fn}).veMask(idx_voxel);
    params.analysis.spatial.(fNames{fn}).exponent      = params.analysis.spatial.(fNames{fn}).exponent(idx_voxel);
    if strcmp(params.analysis.temporalModel,'3ch-stLN')
        params.analysis.temporal.(fNames{fn}).param.exponent = params.analysis.temporal.param.exponent.(fNames{fn});
    end
end
if strcmp(params.analysis.temporalModel,'3ch-stLN')
    params.analysis.temporal.param.exponent = [params.analysis.temporal.lh.param.exponent, params.analysis.temporal.rh.param.exponent];
end
params.analysis.normAcrossRuns = true;
params.analysis.hrf.func = [];
params.analysis.hrf.type = 'spm';
params.useGPU = false;

predictions = stPredictBOLDFromStim(params, logical(msStim));

% Build spatial, temporal, spatiotemporal pRF models
[linearPRFModel, params] = get3DSpatiotemporalpRFs(params);

toc
%% Make plots

makeprettyfigures;

% Plot filters
figure(1); clf;
imagesc(reshape(predictions.prfs(:,idx_voxel(1)),101,101));
axis square; colorbar; title(sprintf('%s - Spatial pRF filter',params.analysis.spatialModel))

figure(10); clf;
imagesc(reshape(predictions.prfs(:,idx_voxel(4)),101,101));
axis square; colorbar; title(sprintf('%s - Spatial pRF filter',params.analysis.spatialModel))

% 
% figure(2); clf;
% title(sprintf('%s - Temporal pRF filter',params.analysis.temporalModel),'Interpreter', 'none')
% plot(linearPRFModel.temporal(:,1),'b', 'lineWidth',2); hold on;
% plot(1:size(linearPRFModel.temporal(:,1),1),zeros(size(linearPRFModel.temporal(:,1))),'k', 'lineWidth', 0.5)
% 
% if strcmp(params.analysis.temporalModel,'3ch-stLN')
%     hold on; plot(linearPRFModel.temporal(:,2), 'r', 'lineWidth',2);
%     plot(linearPRFModel.temporal(:,3), 'g', 'lineWidth',2);
%     legend('Sustained IRF', '','Transient (odd)','Transient (even)')
% end
% xlabel('Time (ms)'); ylabel('Amplitude'); box off;

%% Plot spatial, spatiotemporal linear prf responses
stim1D = sum(msStim,1);
t_s = 0:params.analysis.temporal.tr:(size(predictions.predBOLD,1)-params.analysis.temporal.tr);

dt = (1/params.analysis.temporal.fs);
t_stim = dt:dt:(dt*length(stim1D));
t_ms = 0:dt:((dt*length(predictions.predNeural))-dt);
xl = [0 max(t_ms)];

if size(predictions.prfResponse,3)>1
    ncols = 5;
else
    ncols = 4;
end
for idx = 1:length(idx_voxel)
    figure(28+idx); clf; set(gcf,'Position',[681    55   930   922]);
    
    % Stimulus
%     subplot(ncols,1,1); plot(t_ms,stim1D./max(stim1D), 'k'); title('Stimulus'); box off;
%     ylabel({'Fraction pixels', 'with contrast'})
%     xlim(xl); ylim([-0.2 1.2]); 
    
    % Spatial pRF response
%     subplot(ncols,1,2); plot(t_ms,spatial_prfResponse', 'k');
%     title(sprintf('Spatial only %s PRF model * stim', params.analysis.spatialModel),'Interpreter', 'none');
%     ylabel({'Pred neural','resp (a.u.)'})
%     xlim(xl); box off;
%     Spatiotemporal pRF response (sustained/linear)
%     
    subplot(ncols,1,1); hold on; plot(t_stim,predictions.prfResponse(:,idx_voxel(idx),1), 'k');
    title(sprintf('PRF * stim', params.analysis.temporalModel, linearPRFModel.names{1}),'Interpreter', 'none');
    ylabel({'Pred neural','resp (a.u.)'})
    xlim(xl); ylim([-0.2 1.2]); box off;
    
    subplot(ncols,1,2); hold on; 
    plot(t_ms,predictions.predNeural(:,idx_voxel(idx),1), 'k');
    title(sprintf('Predicted Neural response', params.analysis.temporalModel, linearPRFModel.names{1}),'Interpreter', 'none');
    ylabel({'Pred neural','resp (a.u.)'})
    xlim(xl); ylim([-0.2 1.2]); box off;
    
    subplot(ncols,1,3); 
    plot(t_s,zeros(1,size(predictions.predBOLD,1)),'k'); hold on;
    plot(t_s,predictions.predBOLD(:,idx_voxel(idx),1), 'k','lineWidth',3);
        title(sprintf('Predicted BOLD response', params.analysis.temporalModel, linearPRFModel.names{1}),'Interpreter', 'none');
        ylabel({'Pred BOLD','(% signal change)'})
        xlim(xl); ylim([-0.04 0.25]); box off;
    
    if size(predictions.predNeural,3)>1
    % Spatiotemporal pRF response (transient odd)
        subplot(ncols,1,2); hold on;
        plot(t_ms,predictions.predNeural(:,idx_voxel(idx),2), 'color',cyany);
%         title(sprintf('Spatiotemporal %s %s PRF model * stim', params.analysis.temporalModel, linearPRFModel.names{2}),'Interpreter', 'none');
%         ylabel({'Pred neural','resp (a.u.)'})
%         xlim(xl); ylim([-0.2 1.2]); box off;
        
        subplot(ncols,1,3);
%         plot(t_s,zeros(1,size(predictions.predBOLD,1)),'k'); hold on;
        plot(t_s,predictions.predBOLD(:,idx_voxel(idx),2), 'color',cyany,'lineWidth',3);
%         title(sprintf('Predicted Neural and BOLD response', params.analysis.temporalModel, linearPRFModel.names{1}),'Interpreter', 'none');
%         ylabel({'Pred BOLD','(% signal change)'})
%         xlim(xl); ylim([-0.04 0.3]); box off;
    
        subplot(ncols,1,4); hold on;
        %     plot(t_ms,predictions.prfResponse(:,idx_voxel(idx),1), 'k');
        hold on; plot(t_s,zeros(1,size(predictions.predBOLD,1)),'k')
%         plot(t_ms,predictions.predNeural(:,idx_voxel(idx),2), 'color',cyany);
        plot(t_s,predictions.predBOLD(:,idx_voxel(idx),1)+predictions.predBOLD(:,idx_voxel(idx),2), 'k','lineWidth',3);
        title(sprintf('Pred BOLD combined', params.analysis.temporalModel, linearPRFModel.names{1}),'Interpreter', 'none');
        ylabel({'Pred BOLD','(% signal change)'})
        xlim(xl); ylim([-0.04 0.25]); box off;
        
    end
%     if size(predictions.predNeural,3)>1
%         % Spatiotemporal pRF response (transient odd)
%         subplot(ncols,1,4); hold on; plot(t_s,predictions.predBOLD(:,idx_voxel(idx),2), 'r');
%         title(sprintf('Spatiotemporal %s %s PRF model * stim', params.analysis.temporalModel, linearPRFModel.names{2}),'Interpreter', 'none');
%         xlim(xl); ylim([-0.04 0.075]); box off;
%         % Spatiotemporal pRF response (transient even)
%         subplot(ncols,1,4); hold on; plot(t_s,predictions.predBOLD(:,idx_voxel(idx),2), 'g:');
%         title(sprintf('Spatiotemporal %s %s PRF model * stim', params.analysis.temporalModel, linearPRFModel.names{3}),'Interpreter', 'none');
%         xlabel('Time (s)');
%         ylabel({'Pred BOLD','(% signal change)'})
% %         xlim(xl); ylim([-0.04 0.075]); box off;
%     end
    subSaveFolder = '~/Desktop';
    printName = sprintf('%s%s%d_200ms_simFirst',params.analysis.spatialModel,params.analysis.temporalModel,idx_voxel(idx));
    print(gcf,'-depsc2','-painters','-r300','-loose',fullfile(subSaveFolder,printName));
end

return
%%
if size(predictions.predNeural,3)>1
    ncols = 6;
else
    ncols = 4;
end

yl = 1.1.*[min(predictions.predBOLD(:,idx_voxel)), max(predictions.predBOLD(:,idx_voxel))];

figure(29); clf; set(gcf,'Position',[681    55   930   922]);
% subplot(ncols,1,1);
% plot(t_ms,relu_prfResponse{1}, 'k');
% xlim(xl);  box off;
% % ylim(yl);
% title('Linear pRF response + relu')

makeprettyfigures

t_ms = dt:dt:(size(predictions.predNeural,1)./1000);
for idx = 1:length(idx_voxel)
    
    for tt = 1:length(t_ms)
    	subplot(ncols,1,2); cla; hold on;
        
        plot(t_ms(1:tt),predictions.predNeural(1:tt,idx_voxel(idx),1), 'k');
        xlim(xl); drawnow; pause(0.001);%ylim(yl); box off
    end
    if size(predictions.predBOLD,3)==1
        title(sprintf('Non-linear: %s PRF response + relu', predictions.params.analysis.nonlinearity),'Interpreter', 'none'); box off;
        
        subplot(ncols,1,3); cla;
        hold on; hold on; plot(t_s,zeros(size(predictions.predBOLD)),'k')
        plot(t_s,predictions.predBOLD(:,idx_voxel(idx),1),  'k','lineWidth',4); hold on
%         hold on; plot(t_ms,predictions.predNeural(:,idx_voxel(idx),1), 'k', 'lineWidth',2); hold on
        xlim(xl); ylim(yl);
        title(sprintf('BOLD response: %s PRF + relu', predictions.params.analysis.nonlinearity),'Interpreter', 'none'); box off;
        
    elseif size(predictions.predBOLD,3)>1
        subplot(ncols,1,1); hold on;
        hold on; plot(t_s,predictions.predBOLD(:,idx_voxel(idx),1), 'lineWidth',2); hold on
        title(sprintf('Non-linear: %s %s PRF response + relu', linearPRFModel.names{1}, predictions.params.analysis.nonlinearity),'Interpreter', 'none'); box off;
        
        subplot(ncols,1,2); hold on;
        plot(t_ms,predictions.predNeural(:,idx_voxel(idx),2));
        title(sprintf('Non-linear: %s %s PRF response + relu', linearPRFModel.names{2}, predictions.params.analysis.nonlinearity),'Interpreter', 'none');
        hold on; plot(t_s,predictions.predBOLD(:,idx_voxel(idx),2), 'lineWidth',2);
        box off; xlim(xl); %ylim(yl);
        
        %     subplot(ncols,1,3); plot(t_ms,predictions.predNeural(:,idx_voxel,3), 'k');
        %     title(sprintf('Non-linear: %s %s PRF response + relu', linearPRFModel.names{3}, predictions.params.analysis.nonlinearity),'Interpreter', 'none'); box off;
        %     hold on; plot(t_s,predictions.predBOLD(:,idx_voxel,3), 'g', 'lineWidth',2);
        %     xlim(xl); %ylim(yl);
        
        subplot(ncols,1,4); hold on;
        plot(t_s,predictions.predBOLD(:,idx_voxel(idx),1), 'lineWidth',2); hold on;
        plot(t_s,predictions.predBOLD(:,idx_voxel(idx),2), 'lineWidth',2);
        %     plot(t_s,predictions.predBOLD(:,idx_voxel,3), 'g:', 'lineWidth',3);
        xlim(xl); %ylim(yl);
        title(sprintf('BOLD responses: %s PRF + relu', predictions.params.analysis.nonlinearity),'Interpreter', 'none'); box off;
        legend({'Sustained','Transient'}, 'Location','Best'); legend boxoff
        
        subplot(ncols,1,5); hold on;
        plot(t_s,sum(predictions.predBOLD(:,idx_voxel(idx),:),3), 'lineWidth',2); hold on;
        xlim(xl); %ylim(yl);
        title(sprintf('Sum of all BOLD responses: %s PRF + relu', predictions.params.analysis.nonlinearity),'Interpreter', 'none'); box off;
    end
    %%
    if size(predictions.predBOLD,3)>1
        yl2 = [-0.1,1];
        t_ms = dt:dt:(size(predictions.predNeural,1)./1000);

        scaleFactor = max(predictions.predNeural(:,idx_voxel(idx),1))/max(predictions.predBOLD(:,idx_voxel(idx),1));
        figure(30);  set(gcf,'Position',[681    55   930   922]); clf;
        subplot(411); hold on;
        plot(t_ms,predictions.predNeural(:,idx_voxel(idx),1), 'k');
        plot(t_s,predictions.predBOLD(:,idx_voxel(idx),1), 'r','lineWidth',3); title('Sustained BOLD response'); box off;
        xlim([0 max(t_s)]);  ylim(yl2)
        subplot(412); hold on;
        plot(t_ms,predictions.predNeural(:,idx_voxel(idx),2),'k');
        plot(t_s,scaleFactor*predictions.predBOLD(:,idx_voxel(idx),2), 'g', 'lineWidth',3); title('SCALED Transient odd BOLD response'); box off;
        xlim([0 max(t_s)]);  ylim(yl2)
%         subplot(413);hold on;
%         plot(t_ms,predictions.predNeural(:,idx_voxel(idx),3));
%         plot(t_s,scaleFactor*predictions.predBOLD(:,idx_voxel(idx),3), 'g', 'lineWidth',3);  title('SCALED Transient even BOLD response'); box off;
%         xlim([0 max(t_s)]); ylim(yl2)
        
        sumBothChannels = (1/2)*predictions.predBOLD(:,idx_voxel(idx),1) + (1/2)*scaleFactor*predictions.predBOLD(:,idx_voxel(idx),2);
        subplot(414); hold on;
        plot(t_s, zeros(size(t_s)),'lineWidth',0.5); hold on
        plot(t_s,sumBothChannels, 'b', 'lineWidth',3); title('50/50 S/T weighted sum of all channel BOLD responses'); box off;
        xlim([0 max(t_s)]); ylim(yl2)
    end
end


%
% % Spatial
% spatial_prfResponse = msStim'*linearPRFModel.spatial.prfs; % time x voxels
%
% % Then temporal
% linST_prfResponse = NaN(size(spatial_prfResponse,1),size(spatial_prfResponse,2),size(linearPRFModel.temporal,2));
% for n = 1:length(linearPRFModel.names)
%     linST_prfResponse(:,:,n) = convCut2(spatial_prfResponse, linearPRFModel.temporal(:,n), size(spatial_prfResponse,1));
% end

%% Apply ReLU
% [relu_prfResponse, params] = applyReLU(linST_prfResponse,params);
%
% % Apply nonlinearity (spatial compression, temporal compression, or both)
% [nonlin_prfResponse, params] = applyNonlinearity(relu_prfResponse,params);
%
% for ii = 1:size(nonlin_prfResponse,2)
%    nonlin_prfResponse(:,:,ii) = normMax(nonlin_prfResponse(:,:,ii));
% end
%
% if length(nonlin_prfResponse)>1
%     transientSum = nonlin_prfResponse(:,:,2) + nonlin_prfResponse(:,:,3);
%     nonlin_prfResponse = cat(3,nonlin_prfResponse,transientSum);
% end
% % Define HRF and convolve with neural prediction
% hrf = canonical_hrf(1 / params.analysis.temporal.fs, [5 14 28]);
%
% % Convolve neural response with HRF per channel, and downsample to TR
% predBOLD = getPredictedBOLDResponse(params, nonlin_prfResponse, hrf);
