function [] = downloadDataTableFromOSF(dataType, varargin)
% Function to download data presented in the paper:
% "Rethinking simultaneous suppression in visual cortex via compressive
%   spatiotemporal population receptive fields."
% by Kupers, Kim, Grill-Spector (2024).
% Journal: bioRxiv
% DOI: https://doi.org/10.1101/2023.06.24.546388
%
% Main OSF storage URL: https://osf.io/rpuhs/
% Supplemental OSF storage URL: https://osf.io/e83az/
%
% Subjects 1 through 10 have the following labels:
% subj01, subj02, subj03, subj07, subj08 subj09, subj10, subj11, subj12, subj13.
%
% INPUTS:
% dataType : Referring to different processed versions of the data
% 1: Summary data table main paper
% 2: Summary data table and time series for suppl figure 7 (DoG)
% 3: Summary data table and time series for suppl figure 8-9 (CSTopt, DN-ST)
% 4: Subject's preproc simseq fMRI modelfits (main)
% 5: Subject's preproc simseq fMRI data (main)
% 6: Subject's stimulus runs
% 7: Subject's behavioral data
% 8: DoG suppl material: Subject's simulated fMRI data and model voxel time series
% 9: stRet suppl material: Subject's preprocessed fMRI data and model voxel time series
%
% OUTPUTS:
% None
%
% Examples:
%   downloadDataTableFromOSF(9, 'subjToLoad', 3)

% Parse inputs
p = inputParser;
p.addRequired('dataType', @isnumeric);
p.addParameter('subjToLoad', 3, @isnumeric); % Subject nrs: 1-10, corresponding to subj[1,2,3,7,8,9,10,11,12,13]
p.parse(dataType,varargin{:});
subjToLoad =  p.Results.subjToLoad;

% Set OSF URL
urlBase = 'https://osf.io';

% s1-s10 preproc data
subjnrs = [1,2,3,7,8,9,10,11,12,13];
for ii = 1:length(subjnrs); subjStr{ii} = sprintf('subj%02d',subjnrs(ii)); end

% Folder to download files
writeDir = fullfile(simseqRootPath,'data','simseq');
if ~exist(writeDir,'dir'), mkdir(writeDir); end

fprintf('[%s]: Start downloading data...\n',mfilename)
drawnow;

switch dataType
    case 1 % Main results: Summary data table
        fprintf('Main results: Summary data table')
        fName = 'SIMSEQ_dataTable_20230426.mat';
        writePth = fullfile(writeDir,'group',fName);
        dataStr = '7648g';
        readPth = fullfile(urlBase,dataStr,'?action=download&version=1');
        websave(writePth,readPth);
        
    case 2 % DoG Suppl material: Summary data table and time series for suppl fig 7
        fprintf('DoG Suppl material: Summary data table')
        fName = 'SIMSEQ_dataTable_differenceOfGaussiansFit_20231128.mat';
        writePth = fullfile(writeDir,'group',fName);
        dataStr = 't7qd3';
        readPth = fullfile(urlBase,dataStr,'?action=download&version=1');
        websave(writePth,readPth);
        
    case 3 % stRet Suppl material (CSTopt, DN-ST): Summary data table and time series for suppl fig 8-9
        fprintf('stRet Suppl material: Summary data table')
        fName = 'SIMSEQ_dataTable_stRetParams_matchVoxels_20240310.mat';
        writePth = fullfile(writeDir,'group',fName);
        dataStr = 'g73dz';
        readPth = fullfile(urlBase,dataStr,'?action=download&version=1');
        websave(writePth,readPth);
        
    case 4 % Main results: Subject's preproc simseq fMRI modelfits
        fprintf('Main results: Subject''s preproc simseq fMRI modelfits..')
        dataStr = {'wmxya','daws7','e8yhj','xq2gc','6fj5d',... subj01-08
            'urj83','3baqw','dr7uj','53gwm','4ewv2'}; % subj09-13
        if ~exist(writeDir,'dir'), mkdir(writeDir); end
        for s = 1:length(subjToLoad)
            fprintf('s%d.',subjToLoad(s)); drawnow;
            readPth = fullfile(urlBase,dataStr{subjToLoad(s)},'?action=download&version=1');
            fName = sprintf('%s_main_cvfits.zip',subjStr{subjToLoad(s)});
            writePth = fullfile(writeDir,fName);
            websave(writePth,readPth);
            unzip(writePth,fullfile(writeDir,subjStr{subjToLoad(s)}))
        end
        
    case 5 % Main results: Subject's preproc simseq fMRI data
        fprintf('Main results: Subject''s preproc simseq fMRI data..')
        dataStr = {'ga4w7','wrgqm','mec4j','cbsgn','7v3n2',... subj01-08
            '4b6sf','gjtxb','69ras','m6kz7','g75kf'}; % subj09-13
        for s = 1:length(subjToLoad)
            fprintf('s%d.',subjToLoad(s)); drawnow;
            readPth = fullfile(urlBase,dataStr{subjToLoad(s)},'?action=download&version=1');
            fName = sprintf('%s_main_preprocData.zip',subjStr{subjToLoad(s)});
            writePth = fullfile(writeDir,fName);
            websave(writePth,readPth);
            unzip(writePth,fullfile(writeDir,subjStr{subjToLoad(s)}))
        end
        
    case 6 % Subject's stimulus runs
        fprintf('Subject''s stimulus runs..')
        dataStr = {'mgkfh','c4kyd'; '2t6q5','tp3sf'; 'hnykp','df9rk'; ... subj01, subj02, subj03
            'aqmfc','49rq3';'tsw8v','34fgx'; 'g852a','kysuq'; ... subj07, subj07, subj09
            '7hgbf','rfng2';'54t7c','x7srj';'mdbv3','jqxvp'; ... subj10, subj11, subj12
            '96uxt','9p6yh'};
        % Loop over subjects
        for s = 1:length(subjToLoad)
            fprintf('s%d.',subjToLoad(s)); drawnow;
            subjnr = subjToLoad(s);
            pths = getSubjectPaths(fullfile(simseqRootPath), subjnr);
            writeDir = fullfile(simseqRootPath,'data','stimuli',subjStr{subjnr});
            if ~exist(writeDir,'dir'), mkdir(writeDir); end
            for r = [1,2]
                readPth = fullfile(urlBase,dataStr{subjnr,r},'?action=download&version=1');
                fName = sprintf('stim_simseq_run%d_v%d.mat',r,pths.expversionNr);
                writePth = fullfile(writeDir,fName);
                websave(writePth,readPth);
            end
        end
        
    case 7 % Subject's behavioral data
        fprintf('All subjects behavioral data..')
        dataStr = {'3yvwr';'n67za'}; % Toon & SimSeq
        projStr = {'toon','simseq'};
        % Loop over projects
        for pp = 1:length(dataStr)
            fName = sprintf('%sBehavior.zip',projStr{pp});
            writeDir = fullfile(simseqRootPath,'data',projStr{pp},'behavior');
            if ~exist(writeDir,'dir'), mkdir(writeDir); end
            writePth = fullfile(writeDir,fName);
            readPth = fullfile(urlBase,dataStr{pp},'?action=download&version=1');
            websave(writePth,readPth);
            unzip(writePth,fullfile(writeDir));
        end
        
    case 8 % DoG suppl material: Subject's simulated voxel time series and fMRI modelfits
        fprintf('DoG suppl material: Subject''s time series and modelfits..')
        dataStr = {'8dbfj', 'nrk8w','r35h9','8txs5','x2efk', ... s1-s5
            'w2yvx','u4tb3','c34xf','5hbsu','hbgjw'}; % s6-s10
        if ~exist(writeDir,'dir'), mkdir(writeDir); end
        for s = 1:length(subjToLoad)
            fprintf('s%d.',subjToLoad(s)); drawnow;
            readPth = fullfile(urlBase,dataStr{subjToLoad(s)},'?action=download&version=1');
            fName = sprintf('%s_DoGPRFSimulation.zip',subjStr{subjToLoad(s)});
            writePth = fullfile(writeDir,fName);
            websave(writePth,readPth);
            unzip(writePth,fullfile(writeDir,subjStr{subjToLoad(s)}))
        end
        
    case 9 % stRet suppl material: Subject's preproc and model time series
        fprintf('stRet suppl material: Subject''s time series and modelfits..')
        dataStr = {'tv572','at372','pg7j4','ek6tn','vuwn5',... s1-s5
            'n4th5','fa4jg'}; % s6, s7
        modelStrCSTfix = {'mxe8j', 'py64c','h58uw','kf25c','s6h2e',... s1-s5
            't9gru', '3t78s'}; % s6, s7
        modelStrCSTopt = {'werfc', '2j8r3','fvgue','8nmz6','tm4qx',... s1-s5
            'ub5j8', 'rtxmy'}; % s6, s7
        modelStrDNST = {'fprza', '84m6x','q6ybv','qpbez','n37yb',... s1-s5
            'xcuqt', '79mcu'}; % s6, s7
        for s = 1:length(subjToLoad)
            fprintf('s%d.',subjToLoad(s)); drawnow;
            for mm = 1:4
                switch mm
                    case 1
                        loadStr = dataStr{subjToLoad(s)};
                        fName = sprintf('%s_stRet_preprocData.zip',subjStr{subjToLoad(s)});
                        versionNr = 2;
                    case 2
                        loadStr = modelStrCSTfix{subjToLoad(s)};
                        fName = sprintf('%s_CSTfix_cvfits.zip',subjStr{subjToLoad(s)});
                        versionNr = 1;
                    case 3
                        loadStr = modelStrCSTopt{subjToLoad(s)};
                        fName = sprintf('%s_CSTopt_cvfits.zip',subjStr{subjToLoad(s)});
                        versionNr = 1;
                    case 4
                        loadStr = modelStrDNST{subjToLoad(s)};
                        fName = sprintf('%s_DNST_cvfits.zip',subjStr{subjToLoad(s)});
                        versionNr = 1;
                end
                readPth = fullfile(urlBase,loadStr,sprintf('?action=download&version=%d',versionNr));
                writePth = fullfile(writeDir,fName);
                websave(writePth,readPth);
                unzip(writePth,fullfile(writeDir,subjStr{subjToLoad(s)}))
            end
        end
        
        if exist(writePth,'file')
            fprintf('\n[%s]: Download complete\n',mfilename)
        else
            error('\n[%s]: Download unsuccessful, please check\n',mfilename)
        end
        
end