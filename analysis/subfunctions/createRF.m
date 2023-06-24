function [rf2D, rf1D] = createRF(params)


for ii = 1:length(params.x0)
    % Build 2D Gaussian pRFs
    % We use build-in vistasoft function because that's what we use to
    % solve retinotopy data.
    rf(:,ii) = rfGaussian2d(params.X, params.Y,...
            params.sigmaMajor(ii), ... 
            params.sigmaMinor(ii), ...
            params.theta(ii), ...
            params.x0(ii), ...
            params.y0(ii));
        
%     This 2D Gaussian pRF model converts to unit volume and truncates at 5 SD   
%   rf(:,ii) = pmGaussian2d(params.X, params.Y,...
%             params.sigmaMajor(ii), ... 
%             params.sigmaMinor(ii), ...
%             params.theta(ii), ...
%             params.x0(ii), ...
%             params.y0(ii));

end

% Check if there are any negative pRFs if using DoG model, build them and
% subtract from center pRF
if strcmp(params.spatialModel,'differenceOfGaussiansFit')
        
    closeToZero = 0.001;
    DoG_idx = params.sigmaSurround-params.sigmaMajor>closeToZero;
    
    if sum(DoG_idx)>0
        for ii = 1:length(DoG_idx)
            rf2  = rfNormGaussian2d(params.X, -1*params.Y,...
            params.sigmaSurround(DoG_idx(mii)), ... 
            params.sigmaSurround(DoG_idx(ii)), ...
            params.theta(DoG_idx(ii)), ...
            params.x0(DoG_idx(ii)), ...
            params.y0(DoG_idx(ii)));
            
            rf(:,DoG_idx(ii)) = rf(:,DoG_idx(ii)) - rf2;
        end
    end
end

% Mask RF at 5 sd to remove trailing edge (~99% of volume)
% if params.trimRFFlag
    for ii = 1:length(params.x0)
        xx = sqrt(size(rf(:,ii),1));
        pRF = reshape(rf(:,ii),xx,xx);
        r = (params.sigmaMajor(ii)*2.5); % diameter of 5 SD pRF
        ctrPRF = [params.y0(ii);params.x0(ii)]./params.sampleRate; % get center 
        cutoff = r./params.sampleRate; % define eccentricity of cutoff
        
        cx = find((min(params.X):params.sampleRate:max(params.X))==0);
        cy = find((min(params.Y):params.sampleRate:max(params.Y))==0);
        mask = makecircleimage(xx,cutoff,[],[],[],0,[cy cx]+ctrPRF);
        truncPRF  =  pRF.*mask;
        truncPRF_sum = truncPRF ./ sum(truncPRF(:));
        rf(:,ii) = truncPRF_sum(:);
    end
% end



rf2D = reshape(rf,[sqrt(size(rf,1)) sqrt(size(rf,1)) size(rf,2)]);
rf2D = flipud(rf2D); % our y-axis is flipped
rf1D = reshape(rf2D, [size(rf,1), size(rf,2)]); 
% end

return
        