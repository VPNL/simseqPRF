function [params] = simseq_loadstRetPRFParams(pths, hemi, roi, ...
    spatialModel, temporalModel, veThresh, hvol, subFolder)
% Function to select pRF params from toonotopy experiment within ROI

% cd to toon data
cd(fullfile(pths.dataDirToon, pths.toon))

if ~exist('veThresh', 'var') || isempty(veThresh)
    veThresh = 0;
end

if ~exist('subFolder', 'var')
    subFolder = '';
end

if ~exist('hvol','var') || isempty(hvol)
    %% Load hidden Gray/Volume
    cd(fullfile(pths.dataDirSimSeq, pths.subjID, pths.session))
    setVAnatomyPath('./3DAnatomy/t1.nii.gz');
    hvol = initHiddenGray;
    hvol = viewSet(hvol, 'curdt', 'MotionComp_RefScan1');
end

if strcmp(hemi,'both'), hemis = {'lh','rh'}; else, hemis = {hemi}; end

%% Load results from pRF model

% %%%%% StRet data, within new ROI %%%%%
% rm_stRet = load(fullfile(pths.stRetResultsDir,'df_cv0_defaultHRF.mat'));
tmodelfull = strrep(temporalModel,'-','_');
tmodel = strsplit(temporalModel,'-');
tmodel = tmodel{2};

a = load(fullfile(pths.stRetResultsDir,sprintf('%s_stRet_DT_CSTopt_DNST_matchingVoxels_simseqSession.mat',pths.subjID)));
DT_simseq = a.DT_simseq{ismember(a.DT_simseq_names,temporalModel)};

for hi = 1:length(hemis)
    h = hemis{hi};
    roiName = strsplit(roi,'_');
    roiName = [h '_' roiName{1}];
    
    simseqRoiFile  = [h '_' roi '.mat'];
    
    hvol = loadROI(hvol,fullfile('./3DAnatomy','ROIs',subFolder,simseqRoiFile),[],[],1,0);
    if hi == 2 && length(hvol.ROIs)==1, hi = 1; end
    if ~isempty(hvol.ROIs) 
        if isfield(hvol.ROIs,'coords')
            % Load roi indices on cortex
            [coordsIdx, coords3D] = roiIndices(hvol, hvol.ROIs(hi).coords,1);
            fprintf('\n[%s] loading ROI: %s\n',mfilename, roiName)
            assert(isequal(coords3D,hvol.ROIs(hi).coords))
        
            % Select subject, model, and ROI pRF data
            idx = ismember(DT_simseq.subjID,pths.subjID) & ...
            ismember(DT_simseq.model,'st') & ...
            ismember(DT_simseq.tmodel, tmodelfull) & ...
            ismember(DT_simseq.roiName, roiName);
        
            roiData = DT_simseq(idx,:);
        
            % double check if ROI coordinates and st ret param coordinates 
            % are from the same anatomical surface:
            % * hvol.coords = cortical surface coords
            % * hvol.ROIs(hi).coords = roi coords on surface (drawn from Toonotopy data)
            % * coords3D = st Ret DataTable coords (spatiotemporal retinotopy experiment)
            [~, indROI, indGray] = intersectCols(roiData.("3DCoords")',hvol.coords);
            [~, indROI2, indGray2] = intersectCols(hvol.ROIs(hi).coords,hvol.coords);
            [~, indROI3, indGray3] = intersectCols(coords3D,hvol.coords);
            assert(isequal(indGray,indGray2))
            assert(isequal(indGray,indGray3))
        
        if ~isempty(roiData)
            
            % Get ROI and gray indices
            %         [~, indROI, indGray] = intersectCols(hvol.ROIs(hi).coords,hvol.coords);
            params.(h).indROI        = indROI3;
            params.(h).indGray       = indGray3;
            params.(h).hvolROI       = hvol.ROIs(hi);
            params.(h).hvolcoords    = hvol.coords;
            
            % Get x,y,sigma, etc params
            params.(h).x0         = roiData.x0(indROI3)'; 
            params.(h).y0         = roiData.y0(indROI3)'; 
            params.(h).sigmaMajor = roiData.sigma(indROI3)'; 
            params.(h).sigmaMinor = roiData.sigma(indROI3)';
            params.(h).theta      = zeros(size(roiData.x0)); 
            params.(h).beta       = roiData.beta(indROI3)'; 
            params.(h).exponent   = roiData.exponent(indROI3)'; % should be 1 for glm and dcts
            params.(h).varexpl    = roiData.cv_varexp(indROI3)';
            
            % Get variance explained and its ROI mask
            params.(h).veMask      = params.(h).varexpl>veThresh;
            params.(h).varexpl_all = params.(h).varexpl;
            params.(h).varexpl     = params.(h).varexpl(params.(h).veMask);
            
            % Remove prf below variance explained threshold
            params.(h).x0          = params.(h).x0(params.(h).veMask);
            params.(h).y0          = params.(h).y0(params.(h).veMask);
            params.(h).sigmaMajor  = params.(h).sigmaMajor(params.(h).veMask);
            params.(h).sigmaMinor  = params.(h).sigmaMinor(params.(h).veMask);
            params.(h).theta       = params.(h).theta(params.(h).veMask);
            params.(h).exponent    = params.(h).exponent(params.(h).veMask);
            params.(h).beta        = params.(h).beta(params.(h).veMask);
            
            if strcmp(temporalModel, '1ch-dcts')
                tparams = cell2mat(roiData.temporal);
                params.(h).temporal.param.tau1      = tparams(indROI3,1)'; 
                params.(h).temporal.param.weight    = tparams(indROI3,2)'; 
                params.(h).temporal.param.tau2      = tparams(indROI3,3)';
                params.(h).temporal.param.n         = tparams(indROI3,4)';
                params.(h).temporal.param.sigma     = tparams(indROI3,5)';
                params.(h).effectiveSize            = params.(h).sigmaMajor./sqrt(params.(h).temporal.param.n);
                
                params.(h).temporal.param.tau1      = params.(h).temporal.param.tau1(params.(h).veMask);
                params.(h).temporal.param.weight    = params.(h).temporal.param.weight(params.(h).veMask);
                params.(h).temporal.param.tau2      = params.(h).temporal.param.tau2(params.(h).veMask);
                params.(h).temporal.param.n         = params.(h).temporal.param.n(params.(h).veMask);
                params.(h).temporal.param.sigma     = params.(h).temporal.param.sigma(params.(h).veMask);
                
            elseif strcmp(temporalModel, '3ch-stLN')
                tparams = cell2mat(roiData.temporal);
                params.(h).temporal.param.exponent   = roiData.exponent(indROI3)';
                params.(h).temporal.param.tau_s      = tparams(indROI3,1)'; % IRF expects a 10 Hz sampled tau parameter
                params.(h).temporal.param.tau_t      = tparams(indROI3,2)'; % IRF expects a 10 Hz sampled tau parameter
                params.(h).effectiveSize             = params.(h).sigmaMajor./sqrt(params.(h).temporal.param.exponent);
                
                params.(h).temporal.param.exponent   = params.(h).temporal.param.exponent(params.(h).veMask);
                params.(h).temporal.param.tau_s      = params.(h).temporal.param.tau_s(params.(h).veMask);
                params.(h).temporal.param.tau_t      = params.(h).temporal.param.tau_t(params.(h).veMask);
            else
                params.(h).temporal.exponent   = [];
            end
            
            if strcmp(spatialModel, 'oneGaussianFit') && strcmp(temporalModel, '1ch-glm')
                params.(h).spatial.exponent = ones(size(params.(h).x0)); % set to one
                params.(h).effectiveSize = params.(h).sigmaMajor;
            end
            
            params.(h).temporal.param.all    = roiData.temporal(params.(h).veMask)';
            params.(h).effectiveSize = params.(h).effectiveSize(params.(h).veMask);
            
            
            params.(h).indROI        = params.(h).indROI(params.(h).veMask);
            params.(h).indGray       = params.(h).indGray(params.(h).veMask);
            params.(h).hvolcoords    = params.(h).hvolcoords(:,params.(h).veMask);
            
            %         if strcmp(spatialModel, 'differenceOfGaussiansFit')
            %             params.(h).sigmaSurround = params.(h).sigmaSurround(params.(h).veMask);
            %             params.(h).betaSurround  = params.(h).betaSurround(params.(h).veMask);
            %         end
        else
            params.(h) = [];
            hvol.ROIs(hi).name = [];
            params.(h).hvolROI.coords = [];
            params.(h).varexpl     = [];
            params.(h).veMask      = [];
            params.(h).x0          = [];
            params.(h).y0          = [];
            params.(h).sigmaMajor  = [];
            params.(h).sigmaMinor  = [];
            params.(h).theta       = [];
            params.(h).exponent    = [];
            params.(h).beta        = [];
            params.(h).effectiveSize = [];
        end
    else
        params.(h) = [];
        hvol.ROIs(hi).name = [];
        params.(h).hvolROI.coords = [];
        params.(h).varexpl     = [];
        params.(h).veMask      = [];
        params.(h).x0          = [];
        params.(h).y0          = [];
        params.(h).sigmaMajor  = [];
        params.(h).sigmaMinor  = [];
        params.(h).theta       = [];
        params.(h).exponent    = [];
        params.(h).beta        = [];
        params.(h).effectiveSize = [];
        end 
    end
    
    % Add XY fov grid
    params.(h).nSamples = 50; %hvol.rm.retinotopyParams.analysis.numberStimulusGridPoints; % 50 corresponds to 101 x 101
    params.(h).sampleRate = 101; %hvol.rm.retinotopyParams.analysis.sampleRate; % deg visual angle
    params.(h).fieldSize  = 12; %hvol.rm.retinotopyParams.analysis.fieldSize; % deg visual angle
    
    x = single( linspace(- params.(h).fieldSize,  params.(h).fieldSize, 1+ params.(h).nSamples*2) );
    y = x;
    [X,Y] = meshgrid(x,y);
    
    % Update the sampling grid to reflect the sample points used.
    params.(h).X = X(:);
    params.(h).Y = -1*Y(:);

end

return