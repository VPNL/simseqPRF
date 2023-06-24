function ds = createDataTable(projectDir)

%% Load data
subjnrs      = [1,2,3,7,8,9,10,11,12,13]; 
roiOrder     = [1:6,8:11,7,12,13]; 
roiNames     = {'V1','V2','V3','hV4','VO1','VO2','V3AB','LO1','LO2','TO1','TO2','IPS0','IPS1'};
combROIs     = [5,8,10,12]; %VO1/2 (5,6), LO1/2(8,9), TO1/2 (10,11), IPS0/1 (12,13)
conditionNamesSimSeq        = {'Small-0.2s','Small-1s','Large-0.2s','Large-1s'};
nrConditions                = length(conditionNamesSimSeq);
conditionOrderSimSeq        = [1,2,3,4];
modelNames = {'LSS','CSS','ST'};

% Define path to data
sNames         = sprintf('S%i_',subjnrs);
load(fullfile(projectDir,'data','simseq/group',...
    sprintf('%sBlockAmps_pRFparams_modelFit_results_cv_variableBlockOnset.mat',sNames))); 

% Get paths to make figure dir
pths = getSubjectPaths(projectDir,1,3);

% Get top 250 voxels based on split half reliability (i.e. noise ceiling)
nc_thresh      = 0.1; % threshold in percentage noise ceiling from simseq exp
ve_thresh      = 0.2; % threshold in percentage of pRF var expl. from toonotopy
[selectedVoxels, rNC, ve] = selectVoxels(noiseceilingRun, paramsROI, nc_thresh, ve_thresh, roiOrder, combROIs);

%% Create data table for LMM

T = table();

for sj = 1:length(subjnrs)
    for idx = 1:length(roiOrder)
        T_sim = [];
        T_seq = [];
        T_seqModel = [];
        T_simModel = [];
        T_params = [];
        
        if ismember(roiOrder(idx),setdiff(roiOrder,combROIs+1))
            % selectVoxels function already combined ROIs, so no need to do
            % again
            sv_idx = [selectedVoxels{sj,roiOrder(idx)}];

            ncToAdd = NaN(sum(sv_idx),1);
            expCSTPRFToAdd = ncToAdd;
            expCSSPRFToAdd = ncToAdd;
            szPRFToAdd = ncToAdd;
            vePRFToAdd = ncToAdd;
            R2LSSToAdd = ncToAdd;
            R2CSSToAdd = ncToAdd;
            R2CSTToAdd = ncToAdd;
            BetaLSSToAdd = ncToAdd;
            BetaCSSToAdd = ncToAdd;
            BetaCSTSustToAdd = ncToAdd;
            BetaCSTTransToAdd = ncToAdd;
     
            % Combine ROI pRF data, betas and R^2 data 
            if ismember(roiOrder(idx),combROIs) && ~isempty(paramsROI{sj,roiOrder(idx+1)})         
                expCSTPRF = [paramsROI{sj,roiOrder(idx),3}.exp_temporal.both,paramsROI{sj,roiOrder(idx+1),3}.exp_temporal.both];
                expCSSPRF = [paramsROI{sj,roiOrder(idx),2}.exp_spatial.both,paramsROI{sj,roiOrder(idx+1),2}.exp_spatial.both];
                szPRF = [paramsROI{sj,roiOrder(idx),3}.effectiveSize.both,paramsROI{sj,roiOrder(idx+1),3}.effectiveSize.both];
                vePRF = [paramsROI{sj,roiOrder(idx),3}.varexpl.both,paramsROI{sj,roiOrder(idx+1),3}.varexpl.both];
                
                R2LSS = [r2Runs{sj,roiOrder(idx),1},r2Runs{sj,roiOrder(idx+1),1}];
                R2CSS = [r2Runs{sj,roiOrder(idx),2},r2Runs{sj,roiOrder(idx+1),2}];
                R2CST = [r2Runs{sj,roiOrder(idx),3},r2Runs{sj,roiOrder(idx+1),3}];
                
                betaLSS = [betaRuns{sj,roiOrder(idx),1},betaRuns{sj,roiOrder(idx+1),1}];
                betaCSS = [betaRuns{sj,roiOrder(idx),2},betaRuns{sj,roiOrder(idx+1),2}];
                betaCST_sus = [betaRuns{sj,roiOrder(idx),3}(1,:),betaRuns{sj,roiOrder(idx+1),3}(1,:)];
                betaCST_trans = [betaRuns{sj,roiOrder(idx),3}(2,:),betaRuns{sj,roiOrder(idx+1),3}(2,:)];
                
                if ~isequal(paramsROI{sj,roiOrder(idx),1}.effectiveSize.both,paramsROI{sj,roiOrder(idx),3}.effectiveSize.both) || ...
                         ~isequal(paramsROI{sj,roiOrder(idx),2}.effectiveSize.both,paramsROI{sj,roiOrder(idx),3}.effectiveSize.both)
                    warning('[%s]: CST pRF size is not the same as CSS or LSS pRF model size')
                end
                
            elseif ~isempty(paramsROI{sj,roiOrder(idx)})             
                expCSTPRF = paramsROI{sj,roiOrder(idx),3}.exp_temporal.both;
                expCSSPRF = paramsROI{sj,roiOrder(idx),2}.exp_spatial.both;
                szPRF = paramsROI{sj,roiOrder(idx),3}.effectiveSize.both;
                vePRF = [paramsROI{sj,roiOrder(idx),3}.varexpl.both];
                
                R2LSS = [r2Runs{sj,roiOrder(idx),1}];
                R2CSS = [r2Runs{sj,roiOrder(idx),2}];
                R2CST = [r2Runs{sj,roiOrder(idx),3}]; 
                
                betaLSS = betaRuns{sj,roiOrder(idx),1};
                betaCSS = betaRuns{sj,roiOrder(idx),2};
                betaCST_sus = betaRuns{sj,roiOrder(idx),3}(1,:);
                betaCST_trans = betaRuns{sj,roiOrder(idx),3}(2,:);
            end
            nc = rNC{sj,roiOrder(idx)};

            % Put pRF params into table
            if ~isempty(selectedVoxels{sj,roiOrder(idx)})
                T_params = array2table([...
                    double(szPRF(sv_idx))', double(expCSTPRF(sv_idx))', ...
                    double(expCSSPRF(sv_idx))', double(vePRF(sv_idx))', ...
                    double(R2LSS(sv_idx))', double(R2CSS(sv_idx))', ...
                    double(R2CST(sv_idx))', ...
                    double(betaLSS(sv_idx))', ...
                    double(betaCSS(sv_idx))', ...
                    double(betaCST_sus(sv_idx))',double(betaCST_trans(sv_idx))',...
                    double(nc(sv_idx))'], ...
                     'VariableNames',{'pRFCSTsize','pRFCSTexp','pRFCSSexp','pRFCSTve', ...
                     'R2LSS','R2CSS','R2CST','BetaLSS','BetaCSS','BetaCST_s','BetaCST_t','NC'});
%                 

                % Get SEQ and SIM Data 
                % Array shape: 10 subjects x 8 conditions (seq: 1-4,
                % sim:5-8) x 13 ROIs (x 3 models)
                for c = conditionOrderSimSeq  
                    simModel = []; seqModel = []; 
                    if ~isempty(squeeze(amplTrial{sj, c, roiOrder(idx)}))
                        
                        if ismember(roiOrder(idx),combROIs) && ~isempty(amplTrial{sj,c, roiOrder(idx+1)})
                            seqData = [amplTrial{sj, c, roiOrder(idx)},amplTrial{sj, c, roiOrder(idx+1)}];
                            simData = [amplTrial{sj, c+nrConditions, roiOrder(idx)},amplTrial{sj, c+nrConditions, roiOrder(idx+1)}];
                            
                            for mm = 1:3
                                seqModel(mm,:) = [modelTrial{sj,c,roiOrder(idx),mm}, ...
                                                  modelTrial{sj,c,roiOrder(idx+1),mm}];
                                simModel(mm,:) = [modelTrial{sj,c+nrConditions,roiOrder(idx),mm}, ...
                                                  modelTrial{sj,c+nrConditions,roiOrder(idx+1),mm}];         
                            end
                        else
                            seqData = amplTrial{sj, c, roiOrder(idx)};
                            simData = amplTrial{sj, c+nrConditions, roiOrder(idx)};
                            
                            for mm = 1:3
                                seqModel(mm,:) = modelTrial{sj,c,roiOrder(idx),mm};
                                simModel(mm,:) = modelTrial{sj,c+nrConditions,roiOrder(idx),mm};        
                            end
                        end
                        if ismember(roiOrder(idx),combROIs)
                            roiName = {[pths.allROIs{roiOrder(idx)} '/' pths.allROIs{roiOrder(idx+1)}]};
                        else
                            roiName = pths.allROIs(roiOrder(idx));
                        end
                        
                    end
                    
                    % Create observed and modeled amplitudes
                    T_seq = array2table(double(seqData(sv_idx))','VariableNames',{'MeanSeqAmp'});
                    T_sim = array2table(double(simData(sv_idx))','VariableNames',{'MeanSimAmp'});
                    T_seqModel = array2table(double(seqModel(:,sv_idx))','VariableNames',{'MeanSeqAmpModelLSS','MeanSeqAmpModelCSS','MeanSeqAmpModelCST'});
                    T_simModel = array2table(double(simModel(:,sv_idx))','VariableNames',{'MeanSimAmpModelLSS','MeanSimAmpModelCSS','MeanSimAmpModelCST'});
                    
                    sz = size(T_seq,1);
                    
                    % Create labels
                    labels_Tall = table( ...
                        ones(sz,1).*sj,...
                        repmat(roiName,sz,1), ...
                        ones(sz,1).*c, ...
                        'VariableNames',{'Subject','ROI','Condition'});
                    
                    T_tmp = [labels_Tall, T_seq, T_sim, T_params, T_seqModel, T_simModel]; % T_params has (nc,R2 LSS/CSS/CST, pRF sz/exp/ve)
                    T = [T; T_tmp];
                end
            end
        end
    end
end

%% Convert to a dataset
ds = table2dataset(T);
ds.Subject        = ds.Subject;
ds.Condition      = nominal(ds.Condition);
ds.ROI            = nominal(ds.ROI);

