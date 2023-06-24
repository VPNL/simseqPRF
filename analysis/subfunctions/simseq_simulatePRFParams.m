function params = simseq_simulatePRFParams(hemi, spatialModel, veThresh)
% Function to simulate artificual pRF params within ROI

if ~exist('veThresh', 'var') || isempty(veThresh)
    veThresh = 0;
end

if strcmp(hemi,'both'), hemis = {'lh','rh'}; else, hemis = {hemi}; end

sigmas = [0.5, 2];
x0     = [4.23, 5.8];
y0     = [4.23, 5.8];

for hi = 1:length(hemis)
    h = hemis{hi};
    if strcmp(h,'rh')
        a = -1;
    else
        a = 1;
    end
    % Get x,y,sigma, etc params
    params.(h).x0         = a*repmat(x0,1,length(sigmas));
    params.(h).y0         = a*repmat(y0,1,length(sigmas));
    params.(h).sigmaMajor = repmat(sigmas,length(x0),1);
    params.(h).sigmaMajor = params.(h).sigmaMajor(:)';
    
    params.(h).sigmaMinor = params.(h).sigmaMajor;
    params.(h).theta      = zeros(size(params.(h).x0));
    params.(h).beta       = ones(size(params.(h).x0));
    if strcmp(spatialModel, 'cssFit')
        params.(h).exponent = [0.8 0.8 0.4 0.4];
        params.(h).effectiveSize = params.(h).sigmaMajor ./ sqrt(params.(h).exponent);
    else
        params.(h).exponent = ones(size(params.(h).x0));
        params.(h).effectiveSize = params.(h).sigmaMajor;
    end
    
%     % Get sigma and beta of (negative) surround Gaussian, when using DoG model
%     if strcmp(spatialModel, 'differenceOfGaussiansFit')
%         params.(h).sigmaSurround = hvol.rm.retinotopyModels{1}.sigma2.major(indGray);
%         params.(h).betaSurround  = hvol.rm.retinotopyModels{1}.beta(:,:,2);
%         params.(h).varexplPos  = rmGet(hvol.rm.retinotopyModels{1}, 'varexppos');
%         params.(h).varexplNeg  = rmGet(hvol.rm.retinotopyModels{1}, 'varexpneg');
%     end
    
    % Get variance explained and its ROI mask
    params.(h).varexpl     = ones(size(params.(h).x0));
    params.(h).veMask      = params.(h).varexpl>veThresh;
    
    % Remove prf below variance explained threshold
    params.(h).x0          = params.(h).x0(params.(h).veMask);
    params.(h).y0          = params.(h).y0(params.(h).veMask);
    params.(h).sigmaMajor  = params.(h).sigmaMajor(params.(h).veMask);
    params.(h).sigmaMinor  = params.(h).sigmaMinor(params.(h).veMask);
    params.(h).theta       = params.(h).theta(params.(h).veMask);
    params.(h).exponent    = params.(h).exponent(params.(h).veMask);
    params.(h).beta        = params.(h).beta(params.(h).veMask);

%     if strcmp(spatialModel, 'differenceOfGaussiansFit')
%         params.(h).sigmaSurround = params.(h).sigmaSurround(params.(h).veMask);
%         params.(h).betaSurround  = params.(h).betaSurround(params.(h).veMask);
%     end
%     
    % Add XY fov grid
    params.(h).nSamples   = 50; % 50 corresponds to 101 x 101
    params.(h).fieldSize  = 12; % deg visual angle
    params.(h).sampleRate = params.(h).fieldSize/params.(h).nSamples; % deg visual angle
    
    x = single( linspace(- params.(h).fieldSize,  params.(h).fieldSize, 1+ params.(h).nSamples*2) );
    y = x;
    [X,Y] = meshgrid(x,y);
    
    % Update the sampling grid to reflect the sample points used.
    params.(h).X = X(:);
    params.(h).Y = -1*Y(:);
end
