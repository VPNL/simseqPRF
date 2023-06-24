function fH = makeSupplFigure2_SimulateModelPredictions(saveFigs)
% Simulate neural and BOLD predictions by pRF models
projectDir = simseqRootPath;
subjnr  = 1;
pths    = getSubjectPaths(projectDir, subjnr);
simSeqSessionDir = fullfile(pths.dataDirSimSeq);

% Define params
stimRun             = 1;
hemi                = 'both';
dt                  = 'Averages';
roi                 = 'artificial';
roiType             = 'stimcorner4_area4sq_eccen5'; 
roifile             = [roi '_toon_simseq' roiType];
veThresh            = 0.1;
idx_voxel           = [1:4]; % four artifical pRFs: idx 1 and 2 are small pRFs; 3 and 4 are big pRFs
colors              = getColormapPRFModels(0);

if strcmp(roi,'artificial') || isempty(regexp(roiType,'stimcorner','ONCE'))
    roiSubFolder = [];
else
    roiSubFolder = 'simseqROIs';
end

params.analysis.spatialModel          = 'onegaussianFit'; % Choose: 'onegaussianFit' or 'cssFit'
params.analysis.temporalModel         = '3ch-stLN';       % Choose: '1ch-glm' or '3ch-stLN'
params.analysis.normNeuralChan        = true;
params.stim.framePeriod               = 1; % TR (seconds)
params.saveDataFlag                   = true;
params.stim.sparsifyFlag              = false; % leave out zero response in conv(prf,stim)
params.recomputePredictionsFlag       = true;
params.analysis.predFile              = 'tmp.mat';
params.analysis.zeroPadPredNeuralFlag = true;
params.analysis.reluFlag              = true;
params.analysis.useMedianROIexponent  = true;

% Set colors
if strcmp(params.analysis.spatialModel,'onegaussianFit') && ...
        strcmp(params.analysis.temporalModel,'3ch-stLN')
    colorBOLD = colors(3,:);
elseif strcmp(params.analysis.spatialModel,'cssFit') && ...
        strcmp(params.analysis.temporalModel,'1ch-glm')
     colorBOLD = colors(2,:);
elseif strcmp(params.analysis.spatialModel,'onegaussianFit') && ...
        strcmp(params.analysis.temporalModel,'1ch-glm')
     colorBOLD = colors(1,:);
end
%% Generate stimulus
        
% 1 SEQ trial followed by SIM trial
stimfile = fullfile(simSeqSessionDir, sprintf('images_and_sequenceWith33msGap_run%d.mat',stimRun));

if ~exist(stimfile,'dir') 
    simseq_prepareStimulus(fullfile(simseqRootPath,'data','simseq','behavior','subj01'))
end
stim = load(stimfile);
shortSequence = [stim.sequence(10000:13000)', ones(1,10000), stim.sequence(26000:29000)', ones(1,10000)];
stim.sequence = shortSequence;

% Get upsampled stimulus sequence
for ii = 1:length(stim.sequence)
    if stim.sequence(ii)==0
        msStim(:,:,ii) = zeros(101,101);
    else
        msStim(:,:,ii) = stim.images(:,:,stim.sequence(ii));
    end
end
msStim = reshape(msStim, size(msStim,1)*size(msStim,2),[]);

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
    tmp = simseq_simulatePRFParams(hemi, 'cssFit', veThresh);
    for fn = 1:length(fNames)
        params.analysis.spatial.(fNames{fn}).exponent = tmp.(fNames{fn}).exponent;
    end
    
elseif strcmp(params.analysis.temporalModel,'3ch-stLN')
    tmp = simseq_simulatePRFParams(hemi, 'cssFit', veThresh);
    params.analysis.temporal.param = rmfield(params.analysis.temporal.param,'exponent');
    for fn = 1:length(fNames)
        params.analysis.temporal.param.exponent.(fNames{fn}) = tmp.(fNames{fn}).exponent;
    end
end

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

figure(2); clf;
title(sprintf('%s - Temporal pRF filter',params.analysis.temporalModel),'Interpreter', 'none')
plot(linearPRFModel.temporal(:,1),'b', 'lineWidth',2); hold on;
plot(1:size(linearPRFModel.temporal(:,1),1),zeros(size(linearPRFModel.temporal(:,1))),'k', 'lineWidth', 0.5)

if strcmp(params.analysis.temporalModel,'3ch-stLN')
    hold on; plot(linearPRFModel.temporal(:,2), 'r', 'lineWidth',2);
    plot(linearPRFModel.temporal(:,3), 'g', 'lineWidth',2);
    legend('Sustained IRF', '','Transient (odd)','Transient (even)')
end
xlabel('Time (ms)'); ylabel('Amplitude'); box off;

%% Plot spatial, spatiotemporal linear prf responses
stim1D = sum(msStim,1);
t_s = 0:params.analysis.temporal.tr:(size(predictions.predBOLD,1)-params.analysis.temporal.tr);

dt = (1/params.analysis.temporal.fs);
t_stim = dt:dt:(dt*length(stim1D));
t_ms = 0:dt:((dt*length(predictions.predNeural))-dt);
xl = [0 max(t_ms)];

if size(predictions.prfResponse,3)>1
    ncols = 6;
else
    ncols = 3;
end
for idx = 1:length(idx_voxel)
    figure(28+idx); clf; set(gcf,'Position',[681    55   930   922]);
    
    % Stimulus
    subplot(ncols,1,1); plot(t_stim,stim1D./max(stim1D), 'k'); title('Stimulus'); box off;
    ylabel({'Fraction pixels', 'with contrast'})
    xlim(xl); ylim([-0.2 1.2]); 
    
    % Neural pRF response
    subplot(ncols,1,2); hold on; 
    plot(t_ms,predictions.predNeural(:,idx_voxel(idx),1), 'k');
    title(sprintf('Predicted Neural response %s %s', params.analysis.temporalModel, linearPRFModel.names{1}),'Interpreter', 'none');
    ylabel({'Pred neural','resp (a.u.)'})
    xlim(xl); ylim([-0.2 1.2]); box off;
    
    % BOLD pRF response
    subplot(ncols,1,3); 
    plot(t_s,zeros(1,size(predictions.predBOLD,1)),'k'); hold on;
    plot(t_s,predictions.predBOLD(:,idx_voxel(idx),1), 'k--','lineWidth',3);
        title(sprintf('Predicted BOLD response %s %s', params.analysis.temporalModel, linearPRFModel.names{1}),'Interpreter', 'none');
        ylabel({'Pred BOLD','(% signal change)'})
        xlim(xl); ylim([-0.04 0.25]); box off;
    
    if size(predictions.predNeural,3)>1
    % Spatiotemporal pRF response (transient odd)
        subplot(ncols,1,4); hold on;
        plot(t_ms,predictions.predNeural(:,idx_voxel(idx),2), 'k');
        title(sprintf('Predicted Neural response %s %s', params.analysis.temporalModel, linearPRFModel.names{2}),'Interpreter', 'none');
        ylabel({'Pred neural','resp (a.u.)'})
        xlim(xl); ylim([-0.2 1.2]); box off;
        
        subplot(ncols,1,5);
        plot(t_s,predictions.predBOLD(:,idx_voxel(idx),2), 'color',colorBOLD,'lineWidth',3, 'lineStyle',':');
        title(sprintf('Predicted BOLD response %s %s', params.analysis.temporalModel, linearPRFModel.names{1}),'Interpreter', 'none');
        ylabel({'Pred BOLD','(% signal change)'})
        xlim(xl); ylim([-0.04 0.3]); box off;
    
        subplot(ncols,1,6); hold on;
        hold on; plot(t_s,zeros(1,size(predictions.predBOLD,1)),'k')
        plot(t_s,predictions.predBOLD(:,idx_voxel(idx),1)+predictions.predBOLD(:,idx_voxel(idx),2), 'color',colorBOLD,'lineWidth',3, 'lineStyle','-');
        title(sprintf('Pred BOLD combined %s %s', params.analysis.temporalModel, linearPRFModel.names{1}),'Interpreter', 'none');
        ylabel({'Pred BOLD','(% signal change)'})
        xlim(xl); ylim([-0.04 0.25]); box off;
        
    end
    
    if saveFigs
        saveFigDir = fullfile(simseqRootPath, 'results','group');
        subDir = 'supplFig2';
        if ~exist(fullfile(saveFigDir,subDir),'dir')
            mkdir(fullfile(saveFigDir,subDir)); 
        end
        printName = sprintf('%s%s%d_200ms_SimSeq',params.analysis.spatialModel,params.analysis.temporalModel,idx_voxel(idx));
        print(gcf,'-depsc2','-painters','-r300','-loose',fullfile(saveFigDir,subDir,printName));
    end
end

return

