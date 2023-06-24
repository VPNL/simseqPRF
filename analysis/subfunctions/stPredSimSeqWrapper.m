function stPred = stPredSimSeqWrapper(params)

% preallocate space
stPred = struct();
stPred.predHRF      = [];
stPred.predNeural   = [];
stPred.params       = [];
stPred.prfResponse  = [];
stPred.rf           = [];
stPred.hrf          = [];

% Use average ROI spatial exponent as static exponent
if params.analysis.useMedianROIexponent
    if strcmp(params.analysis.temporalModel, '3ch-stLN')
        catExp = [params.analysis.temporal.param.exponent.lh,params.analysis.temporal.param.exponent.rh];
        medianTemporalExp = median(catExp,'omitnan').*ones(size(catExp));
    end
    if strcmp(params.analysis.spatialModel,'cssFit')
        medianSpatialExp = median([params.analysis.spatial.lh.exponent,params.analysis.spatial.rh.exponent],'omitnan');
    end
end

% Get temporal and spatial params
params = getTemporalParams(params);
params = getSpatialParams(params,1);
params.analysis.temporal.fs = 1000/params.analysis.subsampleRateMs;

% Get number of voxels across hemisphere
if isfield(params.analysis.spatial,'lh')
    numVoxLeft = length(params.analysis.spatial.lh.x0); 
else, numVoxLeft = 0; end
if isfield(params.analysis.spatial,'rh')
    numVoxRight = length(params.analysis.spatial.rh.x0); 
else, numVoxRight = 0; end
numVoxels = numVoxLeft+numVoxRight;

% If fixed exponent or use median exponent across ROIs and ues ST model: 
if strcmp(params.analysis.temporalModel, '3ch-stLN')
    if params.analysis.useMedianROIexponent
        params.analysis.temporal.param.exponent = medianTemporalExp;
    elseif ~isnan(params.analysis.useFixedExponent)
        params.analysis.temporal.param.exponent = ones(1,numVoxels).*params.analysis.useFixedExponent;
    elseif isempty(params.analysis.temporal.param.exponent) && ~params.analysis.useMedianROIexponent && isempty(params.analysis.useFixedExponent)
    	params.analysis.temporal.param.exponent = ones(1,numVoxels);
    elseif ~isempty(params.analysis.temporal.exponent)
        params.analysis.temporal.param.exponent = params.analysis.temporal.exponent;
    end
end

% if median exponent is requested and CSS model is used
if params.analysis.useMedianROIexponent &&  strcmp(params.analysis.spatialModel,'cssFit')
    if isfield(params.analysis.spatial,'lh')
        params.analysis.spatial.lh.exponent = medianSpatialExp.*ones(size(params.analysis.spatial.lh.exponent));
    end
    if isfield(params.analysis.spatial,'rh')
        params.analysis.spatial.rh.exponent = medianSpatialExp.*ones(size(params.analysis.spatial.rh.exponent));
    end
end


% Loop over stimulus files
for numStim = 1:size(params.stim,1)
    
    % Get stim
    stim3Dflattened = double(params.stim(numStim).images_unconvolved);
    fprintf('[%s]: Generating pRF time courses for %s model...%d voxels (stimulus = %d) \n', ...
        mfilename, params.analysis.temporalModel,numVoxels, numStim);
    
    % Generate neural prediction with spatiotemporal model (non-)linearities
    params.analysis.spatial.pRFModelType  = 'unitVolume';
    params.analysis.spatial.keepPixels    = [];
    params.saveDataFlag                   = false; % we'll do it later
    params.analysis.zeroPadPredNeuralFlag = true;
    params.analysis.normNeuralChan        = false;
    params.analysis.reluFlag              = true;
    params.analysis.hrf.type              = 'spm';
    params.analysis.normAcrossRuns        = false;
    predictions = stPredictBOLDFromStim(params, stim3Dflattened);
    
    prf2D = reshape(predictions.prfs,sqrt(size(predictions.prfs,1)),sqrt(size(predictions.prfs,1)),[]);
    stPred(numStim).predHRF      = predictions.predBOLD;
    stPred(numStim).predNeural   = predictions.predNeural;
    stPred(numStim).params       = predictions.params;
    stPred(numStim).prfResponse  = predictions.prfResponse;
    stPred(numStim).rf           = prf2D;
    stPred(numStim).hrf          = predictions.hrf;
end

%%
if params.verbose
    
    timeSnip = [0 size(stim3Dflattened,2)]; timeSnipZoom = [10000 size(stim3Dflattened,2)/4];
%     timeSnip = [10000 45500]; timeSnipZoom = [10000 20000];
    green = [0 0.7 0]; red = [0.7 0 0];
    figure(2); clf; set(gcf, 'Color', 'w', 'Position',[619 1 1359 960]);
    for ll = 1:size(prf2D,3)
        clf;
        % STIMULUS PIXELS
        subplot(4,3,[1 2]); cla; 
        plot(any(stim3Dflattened), 'k', 'lineWidth',2); 
        title('Stimulus time course (matrix flattened across pixels)')
        xlabel('Time (ms)');  ylabel({'Pixel contrast'});
        xlim(timeSnip)
        box off; set(gca,'TickDir','out','FontSize',12)
        
        % 2D PRF 
        subplot(4,3,3); cla; imagesc(squeeze(prf2D(:,:,ll))); axis square; colorbar; %set(gca,'CLim',[0 1])
        title(sprintf('pRF %d',ll))
        hold on; plot([50 50],[1 101],'w-',[1 101],[50 50],'w-')
        set(gca, 'XTick',[1 50 101],'XTickLabel',{'-12' '0' '12'}, ...
            'YTick',[1 50 101],'YTickLabel',{'-12' '0' '12'}, ...
            'TickDir','out','FontSize',12);
        
        % PRF STIMULUS TIME COURSE
        subplot(4,3,[4 5]); plot(predictions.prfResponse(:,ll,1), 'k', 'lineWidth',2);
        title('pRF*stimulus');
        box off; set(gca,'TickDir','out','FontSize',12)
        xlabel('Time (ms)');  xlim(timeSnip); ylabel({'pRF response'; 'not norm. (a.u.)'})
        
        % ZOOM OF PRF STIMULUS TIME COURSE
        subplot(4,3,[6]); plot(predictions.prfResponse(:,ll,1), 'k', 'lineWidth',2);
        box off; set(gca,'TickDir','out','FontSize',12)
        xlabel('Time (ms)');  xlim(timeSnipZoom); title('Zoom: pRF*stim');
        
        % PREDICTED NEURAL TIME COURSE (by st model)
        subplot(4,3,[7 8]); cla; plot(predictions.predNeural(:,ll,1), 'k', 'lineWidth',2); 
        title('Predicted neural timecourse ')
        box off; set(gca,'TickDir','out','FontSize',12)
        xlabel('Time (ms)'); xlim(timeSnip); ylabel({'Neural response'; 'scaled to sust. (a.u.)'})
        
        % ZOOM PREDICTED NEURAL TIME COURSE (by st model)
        subplot(4,3,[9]); cla; plot(predictions.predNeural(:,ll,1), 'k', 'lineWidth',2); 
        box off; set(gca,'TickDir','out','FontSize',12)
        xlabel('Time (ms)'); xlim(timeSnipZoom); box off;
        title('Zoom: neural sustained'); 
        
        % PREDICTED BOLD RESPONSE
        subplot(4,3,[10 11]);plot(predictions.predBOLD(:,ll,1), 'k', 'lineWidth',2); 
        title('Predicted BOLD time course');
        xlabel('Time (s)'); xlim(timeSnip./1000)
        box off; set(gca,'TickDir','out','FontSize',12)
        ylabel({'BOLD response'; 'scaled to sust. (% change)'}) 

        if size(predictions.predBOLD,3)>1
            % PREDICTED NEURAL TIME COURSE
            subplot(4,3,[7 8]); 
            if size(predictions.predNeural,3)==2
                hold on; plot(predictions.predNeural(:,ll,2), ':','color',red, 'lineWidth',1);
                legend({'Sustained','Transient'}, 'Location','NorthEast'); legend boxoff;
            elseif size(predictions.predNeural,3)==3
                hold on; plot(predictions.predNeural(:,ll,2), ':','color',red, 'lineWidth',1);
                hold on; plot(predictions.predNeural(:,ll,3), ':','color',green, 'lineWidth',1);
                legend({'Sustained','Transient On', 'Transient Off'}, 'Location','NorthEast'); legend boxoff;
            end
            % ZOOM NEURAL TIME COURSE
            subplot(4,3,[12]); 
            if size(predictions.predBOLD,3)==2
                hold on; plot(predictions.predNeural(:,ll,2), '-','color',red, 'lineWidth',2);
            elseif size(predictions.predBOLD,3)==3
                hold on; plot(predictions.predNeural(:,ll,2), ':','color',red, 'lineWidth',1);
                hold on; plot(predictions.predNeural(:,ll,3), ':','color',green, 'lineWidth',1);
            end
            title('Zoom: neural transient'); xlabel('Time (ms)'); xlim(timeSnipZoom)
            box off; set(gca,'TickDir','out','FontSize',12); 
            % PREDICTED BOLD RESPONSE
            subplot(4,3,[10 11]);
            if size(predictions.predBOLD,3)==2
                hold on; plot(predictions.predBOLD(:,ll,2), '-', 'color',red,'lineWidth',2);
                lbls = {'Sustained','Transient'};
            elseif size(predictions.predBOLD,3)==3
                hold on; plot(predictions.predBOLD(:,ll,2), '-', 'color',red,'lineWidth',2);                
                hold on; plot(predictions.predBOLD(:,ll,3), '-', 'color',green,'lineWidth',2);
                lbls = {'Sustained','Transient On', 'Transient Off'};
            end
            hold on; plot(sum(predictions.predBOLD(:,ll,:),3)./size(predictions.predBOLD,3), '-', 'color','b','lineWidth',2);
            legend({lbls{:}, 'Sum'}, 'Location','NorthWest'); legend boxoff; 
        end
        % save fig
        if params.saveFigs
            [folderName1,folderName2] = fileparts(params.analysis.predFile);
            if ~exist(fullfile(folderName1,folderName2),'dir'); mkdir(fullfile(folderName1,folderName2)); end
            saveas(2, fullfile(folderName1,folderName2,sprintf('Stim_and_modelprediction_voxel%d.png',ll)));
        end
    end % voxels
end % verbose

return

