function allDS = createDataTableMatchedSTRet(projectDir, subjnrs)

%% Load data
roiOrder               = [1:6,8:11,7,12,13];
roiNames               = {'V1','V2','V3','hV4','VO1','VO2','V3AB','LO1','LO2','TO1','TO2','IPS0','IPS1'};
combROIs               = [5,8,10,12]; %VO1/2 (5,6), LO1/2(8,9), TO1/2 (10,11), IPS0/1 (12,13)
conditionNamesSimSeq   = {'Small-0.2s','Small-1s','Large-0.2s','Large-1s'};
nrConditions           = length(conditionNamesSimSeq);
conditionOrderSimSeq   = [1,2,3,4];

% Get paths to make figure dir
pths = getSubjectPaths(projectDir,1);

modelNames = {'CST_fix','CST_opt','DN_ST'};

% Get CST fix data
sNames           = sprintf('S%i_',subjnrs);
saveFolderNameDS = fullfile(projectDir,'experiments/simseq/results/group/');

% Get top 250 voxels based on split half reliability (i.e. noise ceiling)
nc_thresh      = 0.1; % threshold in percentage noise ceiling from simseq exp
ve_thresh      = 0.2; % threshold in percentage of pRF var expl. from toonotopy
tau_thresh     = 999; % ms
tauScaleFactor = 10; % go from centisecond to millisecond

% Load f data file
d = dir(fullfile(saveFolderNameDS,subDir, sprintf('*stRet_params.mat')));
dat = load(fullfile(d(end).folder,d(end).name));

%% Create data table
clear cstfix_coords cstfix_selectedVoxels
rNC = cell(length(subjnrs), length(roiOrder));
matchingVoxels_coords = rNC;
matchingVoxels_idx = rNC;

for sj = 1:length(subjnrs)
    for idx = 1:length(roiOrder)
        if ~isempty(dat.alignMeROI{sj,roiOrder(idx),1})
            clear ve rNC_tmp
            
            ve  = [dat.paramsROI{sj, roiOrder(idx),1}.varexpl.both];
            rNC_tmp = [dat.noiseceilingRun{sj, roiOrder(idx),1}];
            
            thresh = ((rNC_tmp>=nc_thresh) & (ve>=ve_thresh) & (dat.paramsROI{sj, roiOrder(idx),2}.temporal.tau_s.*tauScaleFactor<=tau_thresh));

            matchingVoxels_idx{sj,roiOrder(idx)} = thresh;
            if size(dat.alignMeROI{sj,roiOrder(idx),1}.commonCoords.lh,2)~=0
                thresh_lh = thresh(1:size(dat.alignMeROI{sj,roiOrder(idx),1}.commonCoords.lh,2));
                matchingVoxels_coords{sj,roiOrder(idx)} = dat.alignMeROI{sj,roiOrder(idx),1}.commonCoords.lh(:,thresh_lh);
            else
                thresh_lh = [];
                matchingVoxels_coords{sj,roiOrder(idx)} = [];
            end
            if size(dat.alignMeROI{sj,roiOrder(idx),1}.commonCoords.rh,2)~=0
                thresh_rh = thresh(length(thresh_lh)+[1:size(dat.alignMeROI{sj,roiOrder(idx),1}.commonCoords.rh,2)]);
                matchingVoxels_coords{sj,roiOrder(idx)} = ...
                    cat(2,matchingVoxels_coords{sj,roiOrder(idx)},...
                    dat.alignMeROI{sj,roiOrder(idx),1}.commonCoords.rh(:,thresh_rh));
            else
                thresh_lh = 0;
            end
            
        end
    end
end

%% Create data table for LMM
clear ds allDS

for mm = 1:length(modelNames)
    T = table();
    T_tmp = [];
    for sj = 1:length(subjnrs)
        for idx = 1:length(roiOrder)
            
            % Preallocate tables
            T_sim = [];
            T_seq = [];
            T_seqModel = [];
            T_simModel = [];
            T_params = [];
            
            if ismember(roiOrder(idx),combROIs) && ~isempty(matchingVoxels_idx{sj,roiOrder(idx)})
               coords = cat(2,matchingVoxels_coords{sj,roiOrder(idx)},matchingVoxels_coords{sj,roiOrder(idx+1)});
               sv_idx = [matchingVoxels_idx{sj,roiOrder(idx)}, matchingVoxels_idx{sj,roiOrder(idx+1)}];
               sv_idx = logical(sv_idx);
               ncToAdd = NaN(sum(sv_idx),1);
            elseif ~ismember(roiOrder(idx),combROIs+1) && ~isempty(matchingVoxels_idx{sj,roiOrder(idx)})
                coords = matchingVoxels_coords{sj,roiOrder(idx)};
                sv_idx = matchingVoxels_idx{sj,roiOrder(idx)};
                sv_idx = logical(sv_idx);
                ncToAdd = NaN(sum(sv_idx),1);
            else
                sv_idx = [];
                coords = {};
                ncToAdd = [];
            end

            if ~isempty(sv_idx)
                if ismember(roiOrder(idx),combROIs) && ~isempty(dat.paramsROI{sj,roiOrder(idx+1),mm})
 
                    if strcmp(modelNames{mm}, 'CST_opt')
                        R2 = [dat.r2Runs{sj,roiOrder(idx),mm},dat.r2Runs{sj,roiOrder(idx+1),mm}];
                        szPRF = [dat.paramsROI{sj,roiOrder(idx),mm}.effectiveSize.both,dat.paramsROI{sj,roiOrder(idx+1),mm}.effectiveSize.both];
                        vePRF = [dat.paramsROI{sj,roiOrder(idx),mm}.varexpl.both,dat.paramsROI{sj,roiOrder(idx+1),mm}.varexpl.both];
                        betaCST_ST_sus = [dat.betaRuns{sj,roiOrder(idx),mm}(1,:),dat.betaRuns{sj,roiOrder(idx+1),mm}(1,:)];
                        betaCST_ST_trans = [dat.betaRuns{sj,roiOrder(idx),mm}(2,:),dat.betaRuns{sj,roiOrder(idx+1),mm}(2,:)];
                        tauPRF = [dat.paramsROI{sj,roiOrder(idx),mm}.temporal.tau_s.*tauScaleFactor,dat.paramsROI{sj,roiOrder(idx+1),mm}.temporal.tau_s.*tauScaleFactor];
                        expPRF = [dat.paramsROI{sj,roiOrder(idx),mm}.temporal.exponent,dat.paramsROI{sj,roiOrder(idx+1),mm}.temporal.exponent];
                        nc = [dat.noiseceilingRun{sj, roiOrder(idx),mm},dat.noiseceilingRun{sj, roiOrder(idx+1),mm}];
                    elseif strcmp(modelNames{mm}, 'DN_ST')
                        R2 = [dat.r2Runs{sj,roiOrder(idx),mm},dat.r2Runs{sj,roiOrder(idx+1),mm}];
                        szPRF = [dat.paramsROI{sj,roiOrder(idx),mm}.effectiveSize.both,dat.paramsROI{sj,roiOrder(idx+1),mm}.effectiveSize.both];
                        vePRF = [dat.paramsROI{sj,roiOrder(idx),mm}.varexpl.both,dat.paramsROI{sj,roiOrder(idx+1),mm}.varexpl.both];
                        betaDN_ST = [dat.betaRuns{sj,roiOrder(idx),mm}(1,:),dat.betaRuns{sj,roiOrder(idx+1),mm}(1,:)];
                        tau1PRF = [dat.paramsROI{sj,roiOrder(idx),mm}.temporal.tau1,dat.paramsROI{sj,roiOrder(idx+1),mm}.temporal.tau1];
                        tau2PRF = [dat.paramsROI{sj,roiOrder(idx),mm}.temporal.tau2,dat.paramsROI{sj,roiOrder(idx+1),mm}.temporal.tau2];
                        nPRF = [dat.paramsROI{sj,roiOrder(idx),mm}.temporal.n,dat.paramsROI{sj,roiOrder(idx+1),mm}.temporal.n];
                        semisatpRF = [dat.paramsROI{sj,roiOrder(idx),mm}.temporal.sigma,dat.paramsROI{sj,roiOrder(idx+1),mm}.temporal.sigma];
                        nc = [dat.noiseceilingRun{sj, roiOrder(idx),mm},dat.noiseceilingRun{sj, roiOrder(idx+1),mm}];
                    else
                        expCSTPRF = [dat.paramsROI{sj,roiOrder(idx),mm}.exp_temporal.both',dat.paramsROI{sj,roiOrder(idx+1),mm}.exp_temporal.both'];
                        expCSSPRF = [dat.paramsROI{sj,roiOrder(idx),mm}.exp_spatial.both,dat.paramsROI{sj,roiOrder(idx+1),mm}.exp_spatial.both];
                        szPRF = [dat.paramsROI{sj,roiOrder(idx),mm}.effectiveSize.both,dat.paramsROI{sj,roiOrder(idx+1),mm}.effectiveSize.both];
                        vePRF = [dat.paramsROI{sj,roiOrder(idx),mm}.varexpl.both,dat.paramsROI{sj,roiOrder(idx+1),mm}.varexpl.both];
                        R2CST = [dat.r2Runs{sj,roiOrder(idx),mm},dat.r2Runs{sj,roiOrder(idx+1),mm}];
                        betaCST_sus = [dat.betaRuns{sj,roiOrder(idx),mm}(1,:),dat.betaRuns{sj,roiOrder(idx+1),mm}(1,:)];
                        betaCST_trans = [dat.betaRuns{sj,roiOrder(idx),mm}(2,:),dat.betaRuns{sj,roiOrder(idx+1),mm}(2,:)];
                        nc = [dat.noiseceilingRun{sj, roiOrder(idx),mm},dat.noiseceilingRun{sj, roiOrder(idx+1),mm}];
                        
                    end
                elseif ~isempty(dat.paramsROI{sj,roiOrder(idx),mm})
                    if strcmp(modelNames{mm}, 'CST_opt')
                        szPRF = dat.paramsROI{sj,roiOrder(idx),mm}.effectiveSize.both;
                        vePRF = dat.paramsROI{sj,roiOrder(idx),mm}.varexpl.both;
                        R2 = dat.r2Runs{sj,roiOrder(idx),mm};
                        betaCST_ST_sus = dat.betaRuns{sj,roiOrder(idx),mm}(1,:);
                        betaCST_ST_trans = dat.betaRuns{sj,roiOrder(idx),mm}(2,:);
                        tauPRF = dat.paramsROI{sj,roiOrder(idx),mm}.temporal.tau_s.*tauScaleFactor;
                        expPRF = dat.paramsROI{sj,roiOrder(idx),mm}.temporal.exponent;
                        nc = dat.noiseceilingRun{sj, roiOrder(idx),mm};
                    elseif strcmp(modelNames{mm}, 'DN_ST')
                        szPRF = dat.paramsROI{sj,roiOrder(idx),mm}.effectiveSize.both;
                        vePRF = dat.paramsROI{sj,roiOrder(idx),mm}.varexpl.both;
                        R2 = dat.r2Runs{sj,roiOrder(idx),mm};
                        betaDN_ST = dat.betaRuns{sj,roiOrder(idx),mm}(1,:);
                        tau1PRF = dat.paramsROI{sj,roiOrder(idx),mm}.temporal.tau1;
                        tau2PRF = dat.paramsROI{sj,roiOrder(idx),mm}.temporal.tau2;
                        nPRF = dat.paramsROI{sj,roiOrder(idx),mm}.temporal.n;
                        semisatpRF = dat.paramsROI{sj,roiOrder(idx),mm}.temporal.sigma;
                        nc = dat.noiseceilingRun{sj, roiOrder(idx),mm};
                    else
                        expCSTPRF = dat.paramsROI{sj,roiOrder(idx),mm}.exp_temporal.both';
                        expCSSPRF = dat.paramsROI{sj,roiOrder(idx),mm}.exp_spatial.both;
                        szPRF = dat.paramsROI{sj,roiOrder(idx),mm}.effectiveSize.both;
                        vePRF = dat.paramsROI{sj,roiOrder(idx),mm}.varexpl.both;
                        R2CST = dat.r2Runs{sj,roiOrder(idx),mm};
                        betaCST_sus = dat.betaRuns{sj,roiOrder(idx),mm}(1,:);
                        betaCST_trans = dat.betaRuns{sj,roiOrder(idx),mm}(2,:);
                        nc = dat.noiseceilingRun{sj, roiOrder(idx),mm};
                    end
                end
                
                
                if ~isempty(sv_idx)
                    if strcmp(modelNames{mm}, 'CST_opt')
                        T_params = array2table([...
                            double(szPRF(sv_idx))', ... size
                            double(vePRF(sv_idx))', ... ve
                            double(tauPRF(sv_idx))', ... tau
                            double(expPRF(sv_idx))',... exp
                            double(R2(sv_idx))', ... R2
                            double(betaCST_ST_sus(sv_idx))',... beta s
                            double(betaCST_ST_trans(sv_idx))',... beta t
                            double(nc(sv_idx))'], ... nc],...
                            'VariableNames',{'pRFCST_ST_size',...
                            'pRFCST_ST_ve','pRFCST_ST_tau', 'pRFCST_ST_exp'...
                            'R2CST_ST','BetaCST_ST_s','BetaCST_ST_t','NC'});
                        
                    elseif strcmp(modelNames{mm}, 'DN_ST')
                        T_params = array2table([...
                            double(szPRF(sv_idx))', ... size
                            double(vePRF(sv_idx))', ... ve
                            double(tau1PRF(sv_idx))', ... tau1
                            double(tau2PRF(sv_idx))',... tau2
                            double(nPRF(sv_idx))', ... n
                            double(semisatpRF(sv_idx))',... sigma temp
                            double(R2(sv_idx))', ... R2
                            double(betaDN_ST(sv_idx))', ... beta
                            double(nc(sv_idx))'],...
                            'VariableNames',{'pRFDN_ST_size','pRFDN_ST_ve',...
                            'pRFDN_ST_tau1','pRFDN_ST_tau2',...
                            'pRFDN_ST_n','pRFDN_ST_semisat', ...
                            'R2DN_ST','BetaDN_ST_s','NC'});
                    else
                        T_params = array2table([...
                            double(szPRF(sv_idx))', ...
                            double(expCSTPRF(sv_idx))', ...
                            double(expCSSPRF(sv_idx))', ...
                            double(vePRF(sv_idx))', ...
                            double(R2CST(sv_idx))', ...
                            double(betaCST_sus(sv_idx))',...
                            double(betaCST_trans(sv_idx))',...
                            double(nc(sv_idx))'], ...
                            'VariableNames',{'pRFCSTsize','pRFCSTexp','pRFCSSexp','pRFCSTve', ...
                            'R2CST','BetaCST_s','BetaCST_t','NC'});
                    end
                    T_coord = table(mat2cell(coords',ones(1,size(coords,2)),3),'VariableNames',{'Coords'});
                    
                    % Get SEQ and SIM Data
                    % Array shape: 10 subjects x 8 conditions (seq: 1-4,
                    % sim:5-8) x 13 ROIs (x 3 models)
                    
                    for c = conditionOrderSimSeq
                        simModel = []; seqModel = [];
                        
                        if ~isempty(squeeze(dat.amplTrial{sj, c, roiOrder(idx), mm}))

                            if ismember(roiOrder(idx),combROIs) && ~isempty(dat.amplTrial{sj,c, roiOrder(idx+1), mm})
                                seqData = [dat.amplTrial{sj, c, roiOrder(idx), mm},dat.amplTrial{sj, c, roiOrder(idx+1), mm}];
                                simData = [dat.amplTrial{sj, c+nrConditions, roiOrder(idx), mm},dat.amplTrial{sj, c+nrConditions, roiOrder(idx+1), mm}];

                                seqModel = [dat.modelTrial{sj,c,roiOrder(idx),mm}, dat.modelTrial{sj,c,roiOrder(idx+1),mm}];
                                simModel = [dat.modelTrial{sj,c+nrConditions,roiOrder(idx),mm}, dat.modelTrial{sj,c+nrConditions,roiOrder(idx+1),mm}];
                            else
                                seqData = dat.amplTrial{sj, c, roiOrder(idx),mm};
                                simData = dat.amplTrial{sj, c+nrConditions, roiOrder(idx),mm};

                                seqModel = dat.modelTrial{sj,c,roiOrder(idx),mm};
                                simModel = dat.modelTrial{sj,c+nrConditions,roiOrder(idx),mm};
                            end
                            
                            if ismember(roiOrder(idx),combROIs)
                                roiName = {[pths.allROIs{roiOrder(idx)} '/' pths.allROIs{roiOrder(idx+1)}]};
                            else
                                roiName = pths.allROIs(roiOrder(idx));
                            end
                            
                            
                            T_seq = array2table(double(seqData(sv_idx))','VariableNames',{'MeanSeqAmp'});
                            T_sim = array2table(double(simData(sv_idx))','VariableNames',{'MeanSimAmp'});
                            
                            if strcmp(modelNames{mm}, 'CST_opt')
                                T_seqModel = table(double(seqModel(:,sv_idx))','VariableNames',{'MeanSeqAmpModelCST_ST'});
                                T_simModel = table(double(simModel(:,sv_idx))','VariableNames',{'MeanSimAmpModelCST_ST'});
                            elseif strcmp(modelNames{mm}, 'DN_ST')
                                T_seqModel = table(double(seqModel(:,sv_idx))','VariableNames',{'MeanSeqAmpModelDN_ST'});
                                T_simModel = table(double(simModel(:,sv_idx))','VariableNames',{'MeanSimAmpModelDN_ST'});
                            else
                                T_seqModel = table(double(seqModel(:,sv_idx))','VariableNames',{'MeanSeqAmpModelCST'}); 
                                T_simModel = table(double(simModel(:,sv_idx))','VariableNames',{'MeanSimAmpModelCST'}); 
                            end
                            
                            sz = size(T_seq,1);
                            
                            labels_Tall = table( ...
                                ones(sz,1).*sj,...
                                repmat(roiName,sz,1), ...
                                ones(sz,1).*c, ...
                                'VariableNames',{'Subject','ROI','Condition'});
                            
                            T_tmp = [labels_Tall, T_seq, T_sim, T_params, T_seqModel, T_simModel, T_coord]; % T_params has (nc,R2 LSS/CSS/CST, pRF sz/exp/ve)
                            T = [T; T_tmp];
                        end % existing amplitudes
                    end % condition
                end % selected voxels index
            end % noise ceiling data
        end % roi idx
    end % subject idx
    
    ds = table2dataset(T);
    ds.Subject        = ds.Subject;
    ds.Condition      = nominal(ds.Condition);
    ds.ROI            = nominal(ds.ROI);
    allDS{mm}.ds = ds;
    allDS{mm}.modelName = modelNames{mm};
end % model idx

