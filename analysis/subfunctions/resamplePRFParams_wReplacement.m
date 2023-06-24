function [median_resampledPRFSz,median_resampledCSSExp,median_resampledCSTExp, ...
          mean_resampledBetavalSust, mean_resampledBetavalTrans, ...
          resampledCVR2, diff_resampledCVR2, mean_diffCVR2, maxNC] = ...
            resamplePRFParams_wReplacement(ds)

%% Accumulate and resample effective pRF size, CSS exponent, CST exponent, cv-R2

% Define params 
nboot       = 1000; % nr of bootstraps for resampling
subjnrs     = double(unique(ds.Subject)');
roisToPlot  = unique(ds.ROI,'stable');
newROIOrder = [1,2,3,4,5,8,9,6,7];
roisToPlot  = roisToPlot(newROIOrder);

% Sample median
medianPRFSz     = NaN(length(subjnrs),length(roisToPlot));
medianCSTExp    = NaN(length(subjnrs),length(roisToPlot));
medianCSSExp    = NaN(length(subjnrs),length(roisToPlot));

% Resampled Data
resampledPRFSz        = NaN(length(subjnrs),length(roisToPlot),nboot);
resampledCSSExp       = NaN(length(subjnrs),length(roisToPlot),nboot);
resampledCSTExp       = NaN(length(subjnrs),length(roisToPlot),nboot);
resampledBetavalSust  = NaN(length(subjnrs),length(roisToPlot),nboot);
resampledBetavalTrans = NaN(length(subjnrs),length(roisToPlot),nboot);
resampledCVR2         = NaN(length(subjnrs),length(roisToPlot),3,nboot);
diff_resampledCVR2    = NaN(length(subjnrs),length(roisToPlot),3,nboot);

% Resampled Data Median
median_resampledPRFSz      = NaN(length(subjnrs),length(roisToPlot));
median_resampledCSSExp     = NaN(length(subjnrs),length(roisToPlot));
median_resampledCSTExp     = NaN(length(subjnrs),length(roisToPlot));
mean_resampledBetavalSust  = NaN(length(subjnrs),length(roisToPlot));
mean_resampledBetavalTrans = NaN(length(subjnrs),length(roisToPlot));
mean_diffCVR2              = NaN(length(subjnrs),length(roisToPlot),3);

allCVR2  = cell(length(subjnrs),length(roisToPlot),3);
diffCVR2 = allCVR2;

% Max noise ceiling (NC) -- split-half-correlation
maxNC = NaN(length(subjnrs),length(roisToPlot));


for idx = 1:length(roisToPlot)
    for sj = 1:length(subjnrs)
        szCSTToPlot = ds.pRFCSTsize(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));
        expCSTToPlot = ds.pRFCSTexp(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));
        
         expCSSToPlot = ds.pRFCSSexp(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));

        if ~isempty(szCSTToPlot)
            medianPRFSz(sj,idx)  = median(szCSTToPlot,'omitnan');
            medianCSTExp(sj,idx) = median(expCSTToPlot,'omitnan');
            medianCSSExp(sj,idx) = median(expCSSToPlot,'omitnan');
            
            resampledPRFSz(sj,idx,:)  = randsample(szCSTToPlot, nboot, true);
            resampledCSTExp(sj,idx,:) = randsample(expCSTToPlot, nboot, true);
            resampledCSSExp(sj,idx,:) = randsample(expCSSToPlot, nboot, true);

            median_resampledPRFSz(sj,idx,:) = median(squeeze(resampledPRFSz(sj,idx,:)),'omitnan');
            median_resampledCSTExp(sj,idx,:) = median(squeeze(resampledCSTExp(sj,idx,:)),'omitnan');
            median_resampledCSSExp(sj,idx,:) = median(squeeze(resampledCSSExp(sj,idx,:)),'omitnan');
        end 
        
        %% Cross-validated R2
         lss = 100.*ds.R2LSS(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));
        
        css = 100.*ds.R2CSS(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));
        
        cst = 100.*ds.R2CST(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));
        
        nc = 100.*ds.NC(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));
        
        if ~isempty(lss) && sum(isnan(lss))~=length(lss)
            allCVR2{sj,idx,1} = lss;
            allCVR2{sj,idx,2} = css;
            allCVR2{sj,idx,3} = cst;
            
%             subjR2mat{idx,1} = cat(1,subjR2mat{idx,1},lss);
%             subjR2mat{idx,2} = cat(1,subjR2mat{idx,2},css);
%             subjR2mat{idx,3} = cat(1,subjR2mat{idx,3},cst);
            
%             subjR2adj{sj,idx,1} = 100.*adjustedR2forNumParams(lss./100, 648,1);
%             subjR2adj{sj,idx,2} = 100.*adjustedR2forNumParams(css./100, 648,1);
%             subjR2adj{sj,idx,3} = 100.*adjustedR2forNumParams(cst./100, 648,2);
            
            diffCVR2{sj,idx,1} = (css-lss);
            diffCVR2{sj,idx,2} = (cst-lss);
            diffCVR2{sj,idx,3} = (cst-css);
            
            mean_diffCVR2(sj,idx,1) = mean(css-lss,'omitnan');
            mean_diffCVR2(sj,idx,2) = mean(cst-lss,'omitnan');
            mean_diffCVR2(sj,idx,3) = mean(cst-css,'omitnan');

            resampledCVR2(sj,idx,1,:) = randsample(lss, 1000, true);
            resampledCVR2(sj,idx,2,:) = randsample(css, 1000, true);
            resampledCVR2(sj,idx,3,:) = randsample(cst, 1000, true);
            maxNC(sj,idx) = max(nc,[], 'omitnan');
            
            diff_resampledCVR2(sj,idx,1,:) = randsample(css-lss, 1000, true);
            diff_resampledCVR2(sj,idx,2,:) = randsample(cst-lss, 1000, true);
            diff_resampledCVR2(sj,idx,3,:) = randsample(cst-css, 1000, true);
            
%             mean_diff_resampledCVR2(sj,idx,1) = mean(squeeze(diff_resampledCVR2(sj,idx,1,:)),'omitnan');
            
        else
            resampledCVR2(sj,idx,:,:) = NaN(1,1,3,1000);
            maxNC(sj,idx) = NaN;
        end
        
        %% SUSTAINED & TRANSIENT BETA WEIGHTS
        beta_s = ds.BetaCST_s(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));
        beta_t = ds.BetaCST_t(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));
        
        if ~isempty(beta_s) || length(beta_s) > 1
            resampledBetavalSust(sj, idx,:) = randsample(beta_s, nboot, true);
            resampledBetavalTrans(sj, idx,:) = randsample(beta_t, nboot, true);
            
            mean_resampledBetavalSust(sj,idx) = mean(squeeze(resampledBetavalSust(sj, idx,:)),'omitnan');
            mean_resampledBetavalTrans(sj,idx) = mean(squeeze(resampledBetavalTrans(sj, idx,:)),'omitnan');
        end

    end
end

