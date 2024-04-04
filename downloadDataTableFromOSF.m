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
% INPUTS:
% dataType : Referring to different processed versions of the data
% 1: Summary data table main paper
% 2: Summary data table and time series for suppl figure 7 (DoG)
% 3: Summary data table and time series for suppl figure 8-9 (CSTopt, DN-ST)
% 4: Subject's preproc simseq fMRI modelfits (main)
% 5: Subject's preproc simseq fMRI data (main)
% 6: Subject's stimulus runs
% 7: Subject's behavioral data
%
% OUTPUTS:
% None

if nargin == 1 && ~exist('subjToLoad','var')
    subjToLoad = 1;
end

urlBase = 'https://osf.io';

% s1-s10 preproc data
subjnrs = [1,2,3,7,8,9,10,11,12,13];
for ii = 1:length(subjnrs); subjStr{ii} = sprintf('subj%02d',subjnrs(ii)); end

% Folder to download files
writeDir = fullfile(simseqRootPath,'data','simseq');
if ~exist(writeDir,'dir'), mkdir(writeDir); end

fprintf('[%s]: Start downloading data %s...\n',mfilename)

switch dataType
    case 1 % Summary data table main paper
        fName = 'SIMSEQ_dataTable_20230426.mat';
        writePth = fullfile(writeDir,'group',fName);
        dataStr = '7648g';
        readPth = fullfile(urlBase,dataStr,'?action=download&version=1');
        websave(writePth,readPth);
        
    case 2 % Summary data table and time series for suppl figure 7 (DoG)
        fName = 'SIMSEQ_dataTable_differenceOfGaussiansFit_20231128.mat';
        writePth = fullfile(writeDir,'group',fName);
        dataStr = 't7qd3';
        readPth = fullfile(urlBase,dataStr,'?action=download&version=1');
        websave(writePth,readPth);
        
    case 3 % Summary data table and time series for suppl figure 8-9 (CSTopt, DN-ST)
       fName = 'SIMSEQ_dataTable_stRetParams_matchVoxels_20240310.mat';
        writePth = fullfile(writeDir,'group',fName);
        dataStr = 'g73dz';
        readPth = fullfile(urlBase,dataStr,'?action=download&version=1');
        websave(writePth,readPth);
        
    case 4 % Subject's preproc simseq fMRI modelfits (main)
        dataStr = {'wmxya','daws7','e8yhj','xq2gc','6fj5d',... subj01-08
            'urj83','3baqw','dr7uj','53gwm','4ewv2'}; % subj09-13
        if ~exist(writeDir,'dir'), mkdir(writeDir); end
        for s = 1:length(subjToLoad)
            readPth = fullfile(urlBase,dataStr{subjToLoad(s)},'?action=download&version=1');
            fName = sprintf('%s_cvfits.zip',subjStr{subjToLoad(s)});
            writePth = fullfile(writeDir,fName);
            websave(writePth,readPth);
            unzip(writePth,fullfile(writeDir,subjStr{subjToLoad(s)}))
        end
        
    case 5 % Subject's preproc simseq fMRI data (main)
        dataStr = {'ga4w7','wrgqm','mec4j','cbsgn','7v3n2',... subj01-08
                    '4b6sf','gjtxb','69ras','m6kz7','g75kf'}; % subj09-13
        for s = 1:length(subjToLoad)
            readPth = fullfile(urlBase,dataStr{subjToLoad(s)},'?action=download&version=1');
            fName = sprintf('%s_preprocData.zip',subjStr{subjToLoad(s)});
            writePth = fullfile(writeDir,fName);
            websave(writePth,readPth);
            unzip(writePth,fullfile(writeDir,subjStr{subjToLoad(s)}))
        end
        
    case 6 % Subject's stimulus runs
        dataStr = {'mgkfh','c4kyd'; '2t6q5','tp3sf'; 'hnykp','df9rk'; ... subj01, subj02, subj03
            'aqmfc','49rq3';'tsw8v','34fgx'; 'g852a','kysuq'; ... subj07, subj07, subj09
            '7hgbf','rfng2';'54t7c','x7srj';'mdbv3','jqxvp'; ... subj10, subj11, subj12
            '96uxt','9p6yh'};
        % Loop over subjects
        for s = 1:length(subjToLoad)
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
end

if exist(writePth,'file')
    fprintf('[%s]: Download complete\n',mfilename)
else
    error('[%s]: Download unsuccessful, please check\n',mfilename)
end

end