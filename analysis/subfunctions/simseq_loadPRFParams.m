function [params, hvol] = simseq_loadPRFParams(pths, dt, hemi, roi, spatialModel, veThresh, hvol, subFolder)
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
    setVAnatomyPath('./3DAnatomy/t1.nii.gz');
    hvol = initHiddenGray;
    hvol = viewSet(hvol, 'curdt', dt);
end

%% Load results from pRF model
rm_model = dir(fullfile('./Gray', dt, sprintf('*%s*fFit.mat',spatialModel)));

% Use search fit for DoG model, at least for subj01 (need to check if
% necessary for other subjects as well)
if strcmp(spatialModel, 'differenceOfGaussiansFit')
    rm_model = dir(fullfile('./Gray', dt, sprintf('*%s*sFit.mat',spatialModel)));
end

hvol = rmSelect(hvol,1,fullfile(rm_model.folder,rm_model.name));

if strcmp(hemi,'both'), hemis = {'lh','rh'}; else, hemis = {hemi}; end

for hi = 1:length(hemis)
    h = hemis{hi};
    roifile  = [h '_' roi '.mat'];
    
    hvol = loadROI(hvol,fullfile('./3DAnatomy','ROIs',subFolder,roifile),[],[],1,0);
    fprintf('\n[%s] loading ROI: %s\n',mfilename, roifile)
    
    if (~isempty(hvol.ROIs) && (hi==1)) || ((hi==2) && size(hvol.ROIs,2)>1)
        
        assert(isequal(hvol.ROIs(hi).name,[h '_' roi]));
        
        % Get ROI and gray indices
        [~, indROI, indGray] = intersectCols(hvol.ROIs(hi).coords,hvol.coords);
        params.(h).indROI        = indROI;
        params.(h).indGray       = indGray;
        params.(h).hvolROI       = hvol.ROIs(hi);
        params.(h).hvolcoords    = hvol.coords;
        
        % Get x,y,sigma, etc params
        params.(h).x0         = hvol.rm.retinotopyModels{1}.x0(indGray);
        params.(h).y0         = hvol.rm.retinotopyModels{1}.y0(indGray);
        params.(h).sigmaMajor = hvol.rm.retinotopyModels{1}.sigma.major(indGray);
        params.(h).sigmaMinor = hvol.rm.retinotopyModels{1}.sigma.minor(indGray);
        params.(h).theta      = zeros(size(params.(h).x0));
        params.(h).beta       = hvol.rm.retinotopyModels{1}.beta(:,:,1);
        if strcmp(spatialModel, 'cssFit')
            params.(h).exponent   = hvol.rm.retinotopyModels{1}.exponent(indGray);
            params.(h).exponent(params.(h).exponent > 1) = 1; % cap exp at 1
            params.(h).effectiveSize = params.(h).sigmaMajor ./ sqrt(params.(h).exponent);
        else
            params.(h).exponent = ones(size(params.(h).x0)); % set to one
            params.(h).effectiveSize = params.(h).sigmaMajor;
        end
        
        
        % Get sigma and beta of (negative) surround Gaussian, when using DoG model
        if strcmp(spatialModel, 'differenceOfGaussiansFit')
            params.(h).sigmaSurround = hvol.rm.retinotopyModels{1}.sigma2.major(indGray);
            params.(h).betaSurround  = hvol.rm.retinotopyModels{1}.beta(:,:,2);
            params.(h).varexplPos  = rmGet(hvol.rm.retinotopyModels{1}, 'varexppos');
            params.(h).varexplNeg  = rmGet(hvol.rm.retinotopyModels{1}, 'varexpneg');
        end
        
        % Get variance explained and its ROI mask
        params.(h).varexpl     = rmGet(hvol.rm.retinotopyModels{1}, 'varexplained');
        
        params.(h).varexpl  = params.(h).varexpl(indGray);
        params.(h).veMask   = params.(h).varexpl>veThresh;
        
        % Remove prf below variance explained threshold
        params.(h).x0          = params.(h).x0(params.(h).veMask);
        params.(h).y0          = params.(h).y0(params.(h).veMask);
        params.(h).sigmaMajor  = params.(h).sigmaMajor(params.(h).veMask);
        params.(h).sigmaMinor  = params.(h).sigmaMinor(params.(h).veMask);
        params.(h).theta       = params.(h).theta(params.(h).veMask);
        params.(h).exponent    = params.(h).exponent(params.(h).veMask);
        params.(h).beta        = params.(h).beta(params.(h).veMask);
        
        params.(h).indROI        = indROI(params.(h).veMask);
        params.(h).indGray       = indGray(params.(h).veMask);
        params.(h).hvolROI       = hvol.ROIs(hi);
        params.(h).hvolcoords    = hvol.coords(:,params.(h).veMask);
        
        if strcmp(spatialModel, 'differenceOfGaussiansFit')
            params.(h).sigmaSurround = params.(h).sigmaSurround(params.(h).veMask);
            params.(h).betaSurround  = params.(h).betaSurround(params.(h).veMask);
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
    
    % Add XY fov grid
    
    params.(h).nSamples = hvol.rm.retinotopyParams.analysis.numberStimulusGridPoints; % 50 corresponds to 101 x 101
    params.(h).sampleRate = hvol.rm.retinotopyParams.analysis.sampleRate; % deg visual angle
    params.(h).fieldSize  = hvol.rm.retinotopyParams.analysis.fieldSize; % deg visual angle
    
    x = single( linspace(- params.(h).fieldSize,  params.(h).fieldSize, 1+ params.(h).nSamples*2) );
    y = x;
    [X,Y] = meshgrid(x,y);
    
    % Update the sampling grid to reflect the sample points used.
    params.(h).X = X(:);
    params.(h).Y = -1*Y(:);
end

return
