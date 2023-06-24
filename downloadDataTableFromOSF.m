function [] = downloadDataTableFromOSF(dataType, subjToLoad)
% Function to download data presented in the paper:
% "Rethinking simultaneous suppression in visual cortex via compressive
%   spatiotemporal population receptive fields."
% by Kupers, Kim, Grill-Spector (2023). Journal: XXX. DOI: XXX

% INPUTS:
% dataType : Referring to different processed versions of the data
%               - option 1: Summary data table
%               - option 2: Subject's preproc simseq fMRI modelfits
%               - option 3: Subject's preproc simseq fMRI data
%               - option 4: Subject's stimulus runs
%               - option 5: Subject's behavioral data

% OUTPUTS:
% None

if isempty(subjToLoad) || ~exist('subjToLoad','var')
    subjToLoad = 1;
end

urlBase = 'https://osf.io';
subjStr = {'subj01','subj02','subj03','subj07','subj08','subj09',...
    'subj10','subj11','subj12','subj13'}; % s1-s10 preproc data
subjnrs = [1,2,3,7,8,9,10,11,12,13];

switch dataType
    case 1
        fName = 'SIMSEQ_dataTable_20230426.mat';
        writeDir = fullfile(simseqRootPath,'data','simseq','group');
        if ~exist(writeDir,'dir'), mkdir(writeDir); end
        writePth = fullfile(writeDir,fName);
        dataStr = '7648g';
        readPth = fullfile(urlBase,dataStr,'?action=download&version=1');
        websave(writePth,readPth);
    case 2
        dataStr = {'wmxya','daws7','e8yhj','xq2gc','6fj5d',... subj01-08
            'urj83','3baqw','dr7uj','53gwm','4ewv2'}; % subj09-13
        writeDir = fullfile(simseqRootPath,'data','simseq');
        if ~exist(writeDir,'dir'), mkdir(writeDir); end
        for s = 1:length(subjToLoad)
            readPth = fullfile(urlBase,dataStr{subjToLoad(s)},'?action=download&version=1');
            fName = sprintf('%s_cvfits.zip',subjStr{subjToLoad(s)});
            writePth = fullfile(writeDir,fName);
            websave(writePth,readPth);
            unzip(writePth,fullfile(writeDir,subjStr{subjToLoad(s)}))
        end
    case 3
        dataStr = {'ga4w7','wrgqm','mec4j','cbsgn','7v3n2',... subj01-08
                    '4b6sf','gjtxb','69ras','m6kz7','g75kf'}; % subj09-13
        writeDir = fullfile(simseqRootPath,'data','simseq');
        if ~exist(writeDir,'dir'), mkdir(writeDir); end
        for s = 1:length(subjToLoad)
            readPth = fullfile(urlBase,dataStr{subjToLoad(s)},'?action=download&version=1');
            fName = sprintf('%s_preprocData.zip',subjStr{subjToLoad(s)});
            writePth = fullfile(writeDir,fName);
            websave(writePth,readPth);
            unzip(writePth,fullfile(writeDir,subjStr{subjToLoad(s)}))
        end
    case 4
        dataStr = {'mgkfh','c4kyd'; '2t6q5','tp3sf'; 'hnykp','df9rk'; ... subj01, subj02, subj03
            'aqmfc','49rq3';'tsw8v','34fgx'; 'g852a','kysuq'; ... subj07, subj07, subj09
            '7hgbf','rfng2';'54t7c','x7srj';'mdbv3','jqxvp'; ... subj10, subj11, subj12
            '96uxt','9p6yh'};
        
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
    case 5
        dataStr = {'3yvwr';'n67za'}; % Toon & SimSeq
        projStr = {'toon','simseq'};
        
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
    fprintf('[%s]: Download complete',mfilename)
else
    error('[%s]: Download unsuccessful, please check',mfilename)
end

end