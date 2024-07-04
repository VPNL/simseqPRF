function fH = makeSupplFigure3_SimulateModelPredictions(saveFigs,saveFigDir)
% Function to reproduce supplementary figure 3:
% panel B: Simulate neural and BOLD predictions by pRF models
%
% From the paper:
% Title:   Rethinking simultaneous suppression in visual cortex via
%          compressive spatiotemporal population receptive fields.
% Authors: Kupers, Kim, Grill-Spector (2024).
% Journal: Nature Communications
% DOI:     XXX
%
% Requires getting MRI data from OSF (see downloadDataTableFromOSF.m)
%
% Code written by E.R. Kupers (2024) Stanford University
%
% INPUTS (required):
% - saveFigs        : save figures or not?
% - saveFigDir      : folder to save figures
%
% OUTPUTS:
% - fH         : figure handle
%
%%
projectDir = simseqRootPath;
subjnr     = 1;
pths       = getSubjectPaths(projectDir, subjnr);

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
spatialModels       = {'onegaussianFit','cssFit','onegaussianFit'};
temporalModels      = {'1ch-glm','1ch-glm','3ch-stLN'};
modellabel          = {'LSS','CSS','CST'};
if strcmp(roi,'artificial') || isempty(regexp(roiType,'stimcorner','ONCE'))
    roiSubFolder = [];
else
    roiSubFolder = 'simseqROIs';
end

%% Generate stimulus

% 1 SEQ trial followed by SIM trial
stimfile = fullfile(fullfile(simseqRootPath,'data','stimuli','simulation'), ...
    sprintf('images_and_sequenceWith33msGap_run%d.mat',stimRun));

if ~exist(stimfile,'file')
    simseq_prepareStimulus(fullfile(simseqRootPath,'data','simseq','behavior','subj01'))
end
stim = load(stimfile);
shortSequence = zeros(1,35000);

shortSequence(2000:2200) = 4;
shortSequence(2233:2433) = 2;
shortSequence(2467:2667) = 3;
shortSequence(2700:2900) = 5;
shortSequence(15000:15200) = 6;

uniqueConds = unique(shortSequence);
stimToSimulate = zeros(101,101,length(uniqueConds));
for uc = 2:length(uniqueConds)
    origIm = find(stim.sequence==uniqueConds(uc));
    origIm = origIm(1,1);
    stimToSimulate(:,:,uc) = origIm;
end

for ii = 1:length(shortSequence)
    if shortSequence(ii)==0
        msStim(:,:,ii) = zeros(101,101);
    else
        msStim(:,:,ii) = stimToSimulate(:,:,shortSequence(ii));
    end
    
end

msStim = reshape(msStim, size(msStim,1)*size(msStim,2),[]);
stim.sequence = shortSequence;
stim1D = sum(msStim,1);
%% Set parameters

for nn = 1:length(spatialModels)
    params.analysis.spatialModel          = spatialModels{nn}; % Choose: 'onegaussianFit' or 'cssFit'
    params.analysis.temporalModel         = temporalModels{nn};       % Choose: '1ch-glm' or '3ch-stLN'
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
            params.analysis.temporal.(fNames{fn}).param.tau_s  = params.analysis.temporal.param.tau_s.*ones(size(params.analysis.spatial.(fNames{fn}).x0));
            params.analysis.temporal.(fNames{fn}).param.tau_t  = params.analysis.temporal.(fNames{fn}).param.tau_s;
        end
    end
    if strcmp(params.analysis.temporalModel,'3ch-stLN')
        params.analysis.temporal.param.exponent = [params.analysis.temporal.lh.param.exponent, params.analysis.temporal.rh.param.exponent];
        params.analysis.temporal.param.tau_s = [params.analysis.temporal.lh.param.tau_s, params.analysis.temporal.rh.param.tau_s];
        params.analysis.temporal.param.tau_t = [params.analysis.temporal.lh.param.tau_t, params.analysis.temporal.rh.param.tau_t];
    end
    params.analysis.normAcrossRuns = true;
    params.analysis.hrf.func = [];
    params.analysis.hrf.type = 'spm';
    params.useGPU = false;
    
    predictions = stPredictBOLDFromStim(params, logical(msStim));
    
    % Build spatial, temporal, spatiotemporal pRF models
    % [linearPRFModel, params] = get3DSpatiotemporalpRFs(params);
    % % Plot spatial and temporal components
    % figure(1); clf;
    % imagesc(reshape(predictions.prfs(:,idx_voxel(1)),101,101));
    % axis square; colorbar; title(sprintf('%s - Spatial pRF filter',params.analysis.spatialModel))
    %
    % figure(10); clf;
    % imagesc(reshape(predictions.prfs(:,idx_voxel(4)),101,101));
    % axis square; colorbar; title(sprintf('%s - Spatial pRF filter',params.analysis.spatialModel))
    %
    % figure(2); clf;
    % title(sprintf('%s - Temporal pRF filter',params.analysis.temporalModel),'Interpreter', 'none')
    % plot(linearPRFModel.temporal(:,1,idx_voxel(1)),'b', 'lineWidth',2); hold on;
    % plot(1:size(linearPRFModel.temporal(:,1,idx_voxel(1)),1),zeros(size(linearPRFModel.temporal(:,1,idx_voxel(1)))),'k', 'lineWidth', 0.5)
    %
    % if strcmp(params.analysis.temporalModel,'3ch-stLN')
    %     hold on; plot(linearPRFModel.temporal(:,2,idx_voxel(1)), 'r', 'lineWidth',2);
    %     plot(linearPRFModel.temporal(:,3,idx_voxel(1)), 'g', 'lineWidth',2);
    %     legend('Sustained IRF', '','Transient (odd)','Transient (even)')
    % end
    % xlabel('Time (ms)'); ylabel('Amplitude'); box off;
    % xlim([0 200]);
    toc
    
    %% Plot spatial, spatiotemporal linear prf responses
    
    t_s = 0:params.analysis.temporal.tr:(size(predictions.predBOLD,1)-params.analysis.temporal.tr);
    
    dt = (1/params.analysis.temporal.fs);
    t_stim = dt:dt:(dt*length(stim1D));
    t_ms = 0:dt:((dt*length(predictions.predNeural))-dt);
    xl = [0 max(t_ms)];
    
    ncols = 3;
    
    for idx = 4
        fH = figure; clf; set(gcf,'Position',[681    55   930   922]);
        % Stimulus
        subplot(ncols,1,1); plot(t_stim,stim1D./max(stim1D), 'k'); title('Stimulus'); box off;
        ylabel({'Fraction pixels', 'with contrast'})
        xlim(xl); ylim([-0.2 1.2]);
        
        % Neural pRF response
        subplot(ncols,1,2); hold on;
        yyaxis left
        plot(t_ms,predictions.predNeural(:,idx_voxel(idx),1), 'k');
        xlim(xl); ylim([-0.2 1.2]); box off;
        ylabel({'Pred neural','resp (a.u.)'})
        % ADD BOLD pRF response
        yyaxis right
        plot(t_s,predictions.predBOLD(:,idx_voxel(idx),1), '--','lineWidth',3, 'color',colorBOLD);
        xlim(xl); ylim([-0.2 1.2]); box off;
        ylabel({'Pred BOLD','(% signal change)'})
        title(sprintf('Sustained channel Predicted Neural and BOLD response %s', modellabel{nn}),'Interpreter', 'none');
        
        if size(predictions.predNeural,3)>1
            % Spatiotemporal pRF response (transient odd)
            subplot(ncols,1,3); hold on;
            yyaxis left
            plot(t_ms,predictions.predNeural(:,idx_voxel(idx),2), 'k');
            xlim(xl); ylim([-0.2 1.2]); box off;
            ylabel({'Pred neural','resp (a.u.)'})

            yyaxis right
            plot(t_s,predictions.predBOLD(:,idx_voxel(idx),2), 'color',colorBOLD,'lineWidth',3, 'lineStyle','-');
            xlim(xl); ylim([-0.2 1.2]); box off;
            ylabel({'Pred BOLD','(% signal change)'})
            title(sprintf('Combined Transient channel Predicted Neural & BOLD response %s', modellabel{nn}),'Interpreter', 'none');
        end
    end

    if saveFigs
        subDir = 'supplFig3';
        if ~exist(fullfile(saveFigDir,subDir),'dir')
            mkdir(fullfile(saveFigDir,subDir));
        end
        fName = sprintf('supplFig3_%s%s%d_200ms_SimSeq_%s',params.analysis.spatialModel,params.analysis.temporalModel,idx_voxel(idx),modellabel{nn});
        saveas(gcf, fullfile(fullfile(saveFigDir,subDir), [fName '.png']))
        %         print(gcf,'-depsc2','-painters','-r300','-loose',fullfile(saveFigDir,subDir,fName));
    end
end

return

