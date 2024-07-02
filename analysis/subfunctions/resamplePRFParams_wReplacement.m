function output = resamplePRFParams_wReplacement(ds, roisToPlot, useSTRetParams, temporalModel, spatialModel, subjnrs)
% Function to resample pRF parameters
% For example, effective size, CSS exponent, CST exponent, cv-R2, etc.
% 
% Title:   Rethinking simultaneous suppression in visual cortex via 
%          compressive spatiotemporal population receptive fields.
% Authors: Kupers, Kim, Grill-Spector (2024).
% Journal: Nature Communications
% DOI:     XXX
%
% Main OSF storage URL: https://osf.io/rpuhs/
% Supplemental OSF storage URL: https://osf.io/e83az/
%
% Code written by E.R. Kupers (2024) Stanford University
%
% INPUTS (required):
% - ds              : dataset
% - roisToPlot      : cell with ROI names
% - useSTRetParams  : boolean flag to indicate if we use supplementary
%                       spatiotemporal retinotopy PRF parameters
% - spatialModel    : spatial components of pRF models
% - temporalModel   : temporal components of pRF models
% - subjnrs         : subject nrs from the fMRI experiment
%
% OUTPUTS:
% - output          : struct with sample and resampled pRF parameters.
%
%% Define params
if ~exist('useSTRetParams','var') || isempty(useSTRetParams)
    useSTRetParams = false;
end

if ~exist('subjnrs','var') || isempty(subjnrs)
    subjnrs     = double(unique(ds.Subject)');
end

if useSTRetParams
    allSubjnrs = [1:3,7:10];
else
    allSubjnrs  = [1:3,7:13];
end
nboot       = 1000; % nr of bootstraps for resampling
bins = linspace(0,nboot,50);

% Allocate space for pRF median
if useSTRetParams
    % Sample median & resampled Data
    medianPRFSz     = NaN(length(subjnrs),length(roisToPlot));
    resampledPRFSz  = NaN(length(subjnrs),length(roisToPlot),nboot);
    resampledCVR2   = NaN(length(subjnrs),length(roisToPlot),nboot);
    
    if strcmp(temporalModel,'3ch-stLN')
        medianCST_exp         = NaN(length(subjnrs),length(roisToPlot));
        medianCST_tau         = NaN(length(subjnrs),length(roisToPlot));
        modeCST_tau           = NaN(length(subjnrs),length(roisToPlot));
        resampledCST_exp      = NaN(length(subjnrs),length(roisToPlot),nboot);
        resampledCST_tau      = NaN(length(subjnrs),length(roisToPlot),nboot);
        resampledBetavalSust  = NaN(length(subjnrs),length(roisToPlot),nboot);
        resampledBetavalTrans = NaN(length(subjnrs),length(roisToPlot),nboot);
        
        % Resampled Data Median
        median_resampledPRFSz      = NaN(length(subjnrs),length(roisToPlot));
        median_resampledCST_exp    = NaN(length(subjnrs),length(roisToPlot));
        median_resampledCST_tau    = NaN(length(subjnrs),length(roisToPlot));
        mode_resampledCST_tau      = NaN(length(subjnrs),length(roisToPlot));
        mean_resampledBetavalSust  = NaN(length(subjnrs),length(roisToPlot));
        mean_resampledBetavalTrans = NaN(length(subjnrs),length(roisToPlot));
        
    elseif strcmp(temporalModel,'1ch-dcts')
        medianDNST_n         = NaN(length(subjnrs),length(roisToPlot));
        medianDNST_tau1      = NaN(length(subjnrs),length(roisToPlot));
        medianDNST_tau2      = NaN(length(subjnrs),length(roisToPlot));
        medianDNST_semisat   = NaN(length(subjnrs),length(roisToPlot));
        
        resampledDNST_n         = NaN(length(subjnrs),length(roisToPlot),nboot);
        resampledDNST_tau1      = NaN(length(subjnrs),length(roisToPlot),nboot);
        resampledDNST_tau2      = NaN(length(subjnrs),length(roisToPlot),nboot);
        resampledDNST_semisat   = NaN(length(subjnrs),length(roisToPlot),nboot);
        resampledBetaval         = NaN(length(subjnrs),length(roisToPlot),nboot);
        
        median_resampledPRFSz           = NaN(length(subjnrs),length(roisToPlot));
        median_resampledDNST_n         = NaN(length(subjnrs),length(roisToPlot));
        median_resampledDNST_tau1      = NaN(length(subjnrs),length(roisToPlot));
        median_resampledDNST_tau2      = NaN(length(subjnrs),length(roisToPlot));
        median_resampledDNST_semisat   = NaN(length(subjnrs),length(roisToPlot));
        mean_resampledBetaval           = NaN(length(subjnrs),length(roisToPlot));
    end
    
    % Max noise ceiling (NC) -- split-half-correlation
    maxNC = NaN(length(subjnrs),length(roisToPlot));
    
    for idx = 1:length(roisToPlot)
        for sj = 1:length(subjnrs)
            if strcmp(temporalModel,'3ch-stLN')
                if any(ismember(ds.Properties.VarNames,'pRFCSTsize'))
                    szCST{sj,idx} = ds.pRFCSTsize(ds.Subject==sj & ...
                        ds.ROI==nominal(roisToPlot(idx)) & ...
                        ds.Condition==nominal(1));
                    expCSS{sj,idx} = ds.pRFCSSexp(ds.Subject==sj & ...
                        ds.ROI==nominal(roisToPlot(idx)) & ...
                        ds.Condition==nominal(1));
                    expCST{sj,idx} = ds.pRFCSTexp(ds.Subject==sj & ...
                        ds.ROI==nominal(roisToPlot(idx)) & ...
                        ds.Condition==nominal(1));
                    cst{sj,idx} = 100.*ds.R2CST(ds.Subject==sj & ...
                        ds.ROI==nominal(roisToPlot(idx)) & ...
                        ds.Condition==nominal(1));
                    tauCST{sj,idx} = [];
                    % SUSTAINED & TRANSIENT BETA WEIGHTS
                     beta_s{sj,idx} = ds.BetaCST_s(ds.Subject==sj & ...
                         ds.ROI==nominal(roisToPlot(idx)) & ...
                         ds.Condition==nominal(1));
                     beta_t{sj,idx} = ds.BetaCST_t(ds.Subject==sj & ...
                         ds.ROI==nominal(roisToPlot(idx)) & ...
                         ds.Condition==nominal(1));
                elseif any(ismember(ds.Properties.VarNames,'pRFCST_ST_size'))
                    szCST{sj,idx} = ds.pRFCST_ST_size(ds.Subject==sj & ...
                        ds.ROI==nominal(roisToPlot(idx)) & ...
                        ds.Condition==nominal(1));
                    expCST{sj,idx} = ds.pRFCST_ST_exp(ds.Subject==sj & ...
                        ds.ROI==nominal(roisToPlot(idx)) & ...
                        ds.Condition==nominal(1));
                    tauCST{sj,idx} = ds.pRFCST_ST_tau(ds.Subject==sj & ...
                        ds.ROI==nominal(roisToPlot(idx)) & ...
                        ds.Condition==nominal(1));
                    cst{sj,idx} = 100.*ds.R2CST_ST(ds.Subject==sj & ...
                        ds.ROI==nominal(roisToPlot(idx)) & ...
                        ds.Condition==nominal(1));
                     % SUSTAINED & TRANSIENT BETA WEIGHTS
                     beta_s{sj,idx} = ds.BetaCST_ST_s(ds.Subject==sj & ...
                         ds.ROI==nominal(roisToPlot(idx)) & ...
                         ds.Condition==nominal(1));
                     beta_t{sj,idx} = ds.BetaCST_ST_t(ds.Subject==sj & ...
                         ds.ROI==nominal(roisToPlot(idx)) & ...
                         ds.Condition==nominal(1));
                end
                    
                nc{sj,idx} = 100.*ds.NC(ds.Subject==sj & ...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(1));
                
                % Take sample median and mode, resampled median/mode
                if ~isempty(szCST{sj,idx})
                    medianPRFSz(sj,idx)   = median(szCST{sj,idx},'omitnan');
                    resampledPRFSz(sj,idx,:)   = randsample(szCST{sj,idx}, nboot, true);
                    median_resampledPRFSz(sj,idx,:)      = median(squeeze(resampledPRFSz(sj,idx,:)),'omitnan');
                    
                    medianCST_exp(sj,idx) = median(expCST{sj,idx},'omitnan');
                    resampledCST_exp(sj,idx,:) = randsample(expCST{sj,idx}, nboot, true);
                    median_resampledCST_exp(sj,idx,:) = median(squeeze(resampledCST_exp(sj,idx,:)),'omitnan');
                    
                    if ~isempty(tauCST{sj,idx})
                        medianCST_tau(sj,idx)       = median(tauCST{sj,idx},'omitnan');
                        modeCST_tau(sj,idx)         = mode(tauCST{sj,idx});
                        resampledCST_tau(sj,idx,:)  = randsample(tauCST{sj,idx}, nboot, true);
                        median_resampledCST_tau(sj,idx,:) = median(squeeze(resampledCST_tau(sj,idx,:)),'omitnan');
                        
                        [f,xi] = ksdensity(squeeze(resampledCST_tau(sj,idx,:)), bins);
                        [~,tauMode_idx] = max(f);
                        mode_resampledCST_tau(sj,idx,:) = xi(tauMode_idx);
                    end
                end

                if ~isempty(cst{sj,idx})
                    resampledCVR2(sj,idx,:) = randsample(cst{sj,idx}, 1000, true);
                    maxNC(sj,idx) = max(nc{sj,idx},[], 'omitnan');
                else
                    resampledCVR2(sj,idx,:,:) = NaN(1,1,1000);
                    maxNC(sj,idx) = NaN;
                end
                if ~isempty(beta_s{sj,idx}) || length(beta_s{sj,idx}) > 1
                    resampledBetavalSust(sj, idx,:) = randsample(beta_s{sj,idx}, nboot, true);
                    resampledBetavalTrans(sj, idx,:) = randsample(beta_t{sj,idx}, nboot, true);
                    
                    mean_resampledBetavalSust(sj,idx) = mean(squeeze(resampledBetavalSust(sj, idx,:)),'omitnan');
                    mean_resampledBetavalTrans(sj,idx) = mean(squeeze(resampledBetavalTrans(sj, idx,:)),'omitnan');
                end
                

            elseif strcmp(temporalModel,'1ch-dcts')
               
                szDNST{sj,idx} = ds.pRFDN_ST_size(ds.Subject==sj & ...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(1));
                nDNST{sj,idx} = ds.pRFDN_ST_n(ds.Subject==sj & ...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(1));
                tau1DNST{sj,idx} = ds.pRFDN_ST_tau1(ds.Subject==sj & ...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(1));
                tau2DNST{sj,idx} = ds.pRFDN_ST_tau2(ds.Subject==sj & ...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(1));
                semisatDNST{sj,idx} = ds.pRFDN_ST_semisat(ds.Subject==sj & ...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(1));
                
                if ~isempty(szDNST{sj,idx})
                    medianPRFSz(sj,idx)         = median(szDNST{sj,idx},'omitnan');
                    medianDNST_n(sj,idx)       = median(nDNST{sj,idx},'omitnan');
                    medianDNST_tau1(sj,idx)    = median(tau1DNST{sj,idx},'omitnan');
                    medianDNST_tau2(sj,idx)    = median(tau2DNST{sj,idx},'omitnan');
                    medianDNST_semisat(sj,idx) = median(semisatDNST{sj,idx},'omitnan');

                    resampledPRFSz(sj,idx,:)        = randsample(szDNST{sj,idx}, nboot, true);
                    resampledDNST_n(sj,idx,:)      = randsample(nDNST{sj,idx}, nboot, true);
                    resampledDNST_tau1(sj,idx,:)   = randsample(tau1DNST{sj,idx}, nboot, true);
                    resampledDNST_tau2(sj,idx,:)   = randsample(tau2DNST{sj,idx}, nboot, true);
                    resampledDNST_semisat(sj,idx,:)= randsample(semisatDNST{sj,idx}, nboot, true);

                    
                    median_resampledPRFSz(sj,idx,:)         = median(squeeze(resampledPRFSz(sj,idx,:)),'omitnan');
                    median_resampledDNST_n(sj,idx,:)       = median(squeeze(resampledDNST_n(sj,idx,:)),'omitnan');
                    median_resampledDNST_tau1(sj,idx,:)    = median(squeeze(resampledDNST_tau1(sj,idx,:)),'omitnan');
                    median_resampledDNST_tau2(sj,idx,:)    = median(squeeze(resampledDNST_tau2(sj,idx,:)),'omitnan');
                    median_resampledDNST_semisat(sj,idx,:) = median(squeeze(medianDNST_semisat(sj,idx,:)),'omitnan');
                end
                
                dnst{sj,idx} = 100.*ds.R2DN_ST(ds.Subject==sj & ...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(1));
                
                nc{sj,idx} = 100.*ds.NC(ds.Subject==sj & ...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(1));
                
                if ~isempty(dnst{sj,idx})
                    resampledCVR2(sj,idx,:) = randsample(dnst{sj,idx}, 1000, true);
                    maxNC(sj,idx) = max(nc{sj,idx},[], 'omitnan');
                else
                    resampledCVR2(sj,idx,:,:) = NaN(1,1,1000);
                    maxNC(sj,idx) = NaN;
                end
                
                
                %% BETA WEIGHT
                betaval{sj,idx} = ds.BetaDN_ST_s(ds.Subject==sj & ...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(1));
                
                if ~isempty(betaval{sj,idx}) || length(betaval{sj,idx}) > 1
                    resampledBetaval(sj, idx,:) = randsample(betaval{sj,idx}, nboot, true);
                    
                    mean_resampledBetaval(sj,idx) = mean(squeeze(resampledBetaval(sj, idx,:)),'omitnan');
                end
            
            end
        end
    end
    
    output = struct();
    output.resampledPRFSz = resampledPRFSz;
    output.median_resampledPRFSz = median_resampledPRFSz;
    output.resampledCVR2 = resampledCVR2;
    output.maxNC = maxNC;
    
    if strcmp(temporalModel,'3ch-stLN')
        output.szCST = szCST;
        output.expCST = expCST;
        output.tauCST = tauCST;
        output.R2CST = cst;
        output.NC = nc;
        output.beta_s = beta_s;
        output.beta_t = beta_t;
        output.median_resampledCST_exp = median_resampledCST_exp;
        output.median_resampledCST_tau = median_resampledCST_tau;
        output.mode_resampledCST_tau = mode_resampledCST_tau;
        output.resampledCST_exp = resampledCST_exp;
        output.resampledCST_tau = resampledCST_tau;
        output.mean_resampledBetavalSust = mean_resampledBetavalSust;
        output.mean_resampledBetavalTrans = mean_resampledBetavalTrans;
        output.resampledBetavalSust = resampledBetavalSust;
        output.resampledBetavalTrans = resampledBetavalTrans;
    
    elseif strcmp(temporalModel,'1ch-dcts')
        output.szDNST          = szDNST;
        output.expDNST         = nDNST;
        output.tau1DNST        = tau1DNST;
        output.tau2DNST        = tau2DNST;
        output.semisatDNST     = semisatDNST;
        output.R2DNST          = dnst;
        output.NC               = nc;
        output.beta             = betaval;
        output.mean_resampledBetaval        = mean_resampledBetaval;
        output.median_resampledDNST_n      = median_resampledDNST_n;
        output.median_resampledDNST_tau1   = median_resampledDNST_tau1;
        output.median_resampledDNST_tau2   = median_resampledDNST_tau2;
        output.median_resampledDNST_semisat = median_resampledDNST_semisat;
        output.resampledDNST_n         = resampledDNST_n;
        output.resampledDNST_tau1      = resampledDNST_tau1;
        output.resampledDNST_tau2      = resampledDNST_tau2;
        output.resampledDNST_semisat   = resampledDNST_semisat;
    end

elseif strcmp(spatialModel,'differenceOfGaussiansFit')
    medianPRFSz_center     = NaN(length(subjnrs),length(roisToPlot));
    resampledPRFSz_center  = NaN(length(subjnrs),length(roisToPlot),nboot);
    median_resampledPRFSz_center = NaN(length(subjnrs),length(roisToPlot));
    
    medianPRFSz_surround     = NaN(length(subjnrs),length(roisToPlot));
    resampledPRFSz_surround  = NaN(length(subjnrs),length(roisToPlot),nboot);
    median_resampledPRFSz_surround = NaN(length(subjnrs),length(roisToPlot));
    
    resampledBetaval      = NaN(length(subjnrs),length(roisToPlot),nboot);
    resampledCVR2         = NaN(length(subjnrs),length(roisToPlot),nboot);
    mean_resampledBetaval = NaN(length(subjnrs),length(roisToPlot));
    
    for idx = 1:length(roisToPlot)
        for sj = 1:length(subjnrs)
    
        szCenter{sj,idx} = ds.pRF_sizeCenter(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));
        szSurround{sj,idx} = ds.pRF_sizeSurround(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));
        
        if ~isempty(szCenter{sj,idx})
            medianPRFSz_center(sj,idx)        = median(szCenter{sj,idx},'omitnan');
            medianPRFSz_surround(sj,idx)      = median(szSurround{sj,idx},'omitnan');
            
            resampledPRFSz_center(sj,idx,:)   = randsample(szCenter{sj,idx}, nboot, true);
            resampledPRFSz_surround(sj,idx,:) = randsample(szSurround{sj,idx}, nboot, true);
            
            median_resampledPRFSz_center(sj,idx,:)   = median(squeeze(resampledPRFSz_center(sj,idx,:)),'omitnan');
            median_resampledPRFSz_surround(sj,idx,:) = median(squeeze(resampledPRFSz_surround(sj,idx,:)),'omitnan');
        end
        
        % Cross-validated R2
        dog{sj,idx} = 100.*ds.R2DoG(ds.Subject==find(subjnrs(sj)==allSubjnrs) & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));
        
        betaval{sj,idx} = ds.BetaDoG(ds.Subject==find(subjnrs(sj)==allSubjnrs) & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));

        resampledCVR2(sj,idx,:) = randsample(dog{sj,idx}, nboot, true);
        
        if ~isempty(betaval{sj,idx}) || length(betaval{sj,idx}) > 1
            resampledBetaval(sj, idx,:) = randsample(betaval{sj,idx}, nboot, true);            
            mean_resampledBetaval(sj,idx) = mean(squeeze(resampledBetaval(sj, idx,:)),'omitnan');
        end
        end
    end
    
    output = struct();
    output.szCenter = szCenter;
    output.szSurround = szSurround;
    output.R2DoG = dog;
    output.betaval = betaval;
    output.median_resampledPRFSz_center = median_resampledPRFSz_center;
    output.median_resampledPRFSz_surround = median_resampledPRFSz_surround;
    output.mean_resampledBetaval = mean_resampledBetaval;
    output.resampledCVR2 = resampledCVR2;
    
else
    
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
            
            
            szCST{sj,idx} = ds.pRFCSTsize(ds.Subject==sj & ...
                ds.ROI==nominal(roisToPlot(idx)) & ...
                ds.Condition==nominal(1));
            expCST{sj,idx} = ds.pRFCSTexp(ds.Subject==sj & ...
                ds.ROI==nominal(roisToPlot(idx)) & ...
                ds.Condition==nominal(1));
            expCSS{sj,idx} = ds.pRFCSSexp(ds.Subject==sj & ...
                ds.ROI==nominal(roisToPlot(idx)) & ...
                ds.Condition==nominal(1));
            
            if ~isempty(szCST{sj,idx})
                medianPRFSz(sj,idx)  = median(szCST{sj,idx},'omitnan');
                medianCSTExp(sj,idx) = median(expCST{sj,idx},'omitnan');
                medianCSSExp(sj,idx) = median(expCSS{sj,idx},'omitnan');
                
                resampledPRFSz(sj,idx,:)  = randsample(szCST{sj,idx}, nboot, true);
                resampledCSTExp(sj,idx,:) = randsample(expCST{sj,idx}, nboot, true);
                resampledCSSExp(sj,idx,:) = randsample(expCSS{sj,idx}, nboot, true);
                
                median_resampledPRFSz(sj,idx,:) = median(squeeze(resampledPRFSz(sj,idx,:)),'omitnan');
                median_resampledCSTExp(sj,idx,:) = median(squeeze(resampledCSTExp(sj,idx,:)),'omitnan');
                median_resampledCSSExp(sj,idx,:) = median(squeeze(resampledCSSExp(sj,idx,:)),'omitnan');
            end
            
            %% Cross-validated R2
            lss{sj,idx} = 100.*ds.R2LSS(ds.Subject==sj & ...
                ds.ROI==nominal(roisToPlot(idx)) & ...
                ds.Condition==nominal(1));
            
            css{sj,idx} = 100.*ds.R2CSS(ds.Subject==sj & ...
                ds.ROI==nominal(roisToPlot(idx)) & ...
                ds.Condition==nominal(1));
            
            cst{sj,idx} = 100.*ds.R2CST(ds.Subject==sj & ...
                ds.ROI==nominal(roisToPlot(idx)) & ...
                ds.Condition==nominal(1));
            
            nc{sj,idx} = 100.*ds.NC(ds.Subject==sj & ...
                ds.ROI==nominal(roisToPlot(idx)) & ...
                ds.Condition==nominal(1));
            
            if ~isempty(lss{sj,idx}) && sum(isnan(lss{sj,idx}))~=length(lss{sj,idx})
                allCVR2{sj,idx,1} = lss{sj,idx};
                allCVR2{sj,idx,2} = css{sj,idx};
                allCVR2{sj,idx,3} = cst{sj,idx};
                
                %             subjR2mat{idx,1} = cat(1,subjR2mat{idx,1},lss);
                %             subjR2mat{idx,2} = cat(1,subjR2mat{idx,2},css);
                %             subjR2mat{idx,3} = cat(1,subjR2mat{idx,3},cst);
                
                %             subjR2adj{sj,idx,1} = 100.*adjustedR2forNumParams(lss./100, 648,1);
                %             subjR2adj{sj,idx,2} = 100.*adjustedR2forNumParams(css./100, 648,1);
                %             subjR2adj{sj,idx,3} = 100.*adjustedR2forNumParams(cst./100, 648,2);
                
                diffCVR2{sj,idx,1} = (css{sj,idx}-lss{sj,idx});
                diffCVR2{sj,idx,2} = (cst{sj,idx}-lss{sj,idx});
                diffCVR2{sj,idx,3} = (cst{sj,idx}-css{sj,idx});
                
                mean_diffCVR2(sj,idx,1) = mean(css{sj,idx}-lss{sj,idx},'omitnan');
                mean_diffCVR2(sj,idx,2) = mean(cst{sj,idx}-lss{sj,idx},'omitnan');
                mean_diffCVR2(sj,idx,3) = mean(cst{sj,idx}-css{sj,idx},'omitnan');
                
                resampledCVR2(sj,idx,1,:) = randsample(lss{sj,idx}, 1000, true);
                resampledCVR2(sj,idx,2,:) = randsample(css{sj,idx}, 1000, true);
                resampledCVR2(sj,idx,3,:) = randsample(cst{sj,idx}, 1000, true);
                maxNC(sj,idx) = max(nc{sj,idx},[], 'omitnan');
                
                diff_resampledCVR2(sj,idx,1,:) = randsample(css{sj,idx}-lss{sj,idx}, 1000, true);
                diff_resampledCVR2(sj,idx,2,:) = randsample(cst{sj,idx}-lss{sj,idx}, 1000, true);
                diff_resampledCVR2(sj,idx,3,:) = randsample(cst{sj,idx}-css{sj,idx}, 1000, true);
                
                %             mean_diff_resampledCVR2(sj,idx,1) = mean(squeeze(diff_resampledCVR2(sj,idx,1,:)),'omitnan');
                
            else
                resampledCVR2(sj,idx,:,:) = NaN(1,1,3,1000);
                maxNC(sj,idx) = NaN;
            end
            
            %% SUSTAINED & TRANSIENT BETA WEIGHTS
            beta_s{sj,idx} = ds.BetaCST_s(ds.Subject==sj & ...
                ds.ROI==nominal(roisToPlot(idx)) & ...
                ds.Condition==nominal(1));
            beta_t{sj,idx} = ds.BetaCST_t(ds.Subject==sj & ...
                ds.ROI==nominal(roisToPlot(idx)) & ...
                ds.Condition==nominal(1));
            
            if ~isempty(beta_s{sj,idx}) || length(beta_s{sj,idx}) > 1
                resampledBetavalSust(sj, idx,:) = randsample(beta_s{sj,idx}, nboot, true);
                resampledBetavalTrans(sj, idx,:) = randsample(beta_t{sj,idx}, nboot, true);
                
                mean_resampledBetavalSust(sj,idx) = mean(squeeze(resampledBetavalSust(sj, idx,:)),'omitnan');
                mean_resampledBetavalTrans(sj,idx) = mean(squeeze(resampledBetavalTrans(sj, idx,:)),'omitnan');
            end
            
        end
    end
    
    output = struct();
    output.szCST        = szCST;
    output.expCST       = expCST;
    output.expCSS       = expCSS;
    output.R2CST        = cst;
    output.R2CSS        = css;
    output.R2LSS        = lss;
    output.NC           = nc;
    output.beta_s       = beta_s;
    output.beta_t       = beta_t;
    output.median_resampledPRFSz     = median_resampledPRFSz;
    output.median_resampledCSSExp    = median_resampledCSSExp;
    output.median_resampledCSTExp    = median_resampledCSTExp;
    output.mean_resampledBetavalSust = mean_resampledBetavalSust;
    output.mean_resampledBetavalTrans = mean_resampledBetavalTrans;
    output.resampledCVR2 = resampledCVR2;
    output.diff_resampledCVR2 = diff_resampledCVR2;
    output.mean_diffCVR2 = mean_diffCVR2;
    output.maxNC = maxNC;
    
end





