function data = simseq_loadEventRelatedData(simSeqSessionDir, roiName, varargin)
% data = simseq_loadEventRelatedData(simSeqSessionDir, roiName, ...
%             [dataType], [scansToLoad], [recomputeFlag], [preserveCoordsFlag])

%% Parse inputs
p = inputParser;
p.addRequired('simSeqSessionDir', @ischar);
p.addRequired('roiName',@ischar);
p.addParameter('dataType','MotionComp_RefScan1', @ischar);
p.addParameter('scansToLoad',[1:8], @isnumeric);
p.addParameter('recomputeFlag', false, @islogical);
p.addParameter('preserveCoordsFlag', false, @islogical);
p.parse(simSeqSessionDir,roiName,varargin{:});

%% rename variables
fnToRename = fieldnames(p.Results);
for ff = 1:length(fnToRename)
    eval([sprintf('%s = p.Results.(fnToRename{ff});', fnToRename{ff})]);
end
clear fnToRename ff 

%% Remember current path, so we can change it back later
currPth = pwd;

% Go to mrVista session folder with simseq data
cd(simSeqSessionDir);

% Load single sim & seq TRIALs from GLM / tc* structs 
setVAnatomyPath('./3DAnatomy/t1.nii.gz');
hvol    = initHiddenGray;
hvol    = viewSet(hvol, 'curdt', dataType);
roiPth  = fullfile('simseqROIs', roiName);

load('mrSESSION.mat')
for ii = 1:length(dataTYPES)
    if strcmp(dataType,dataTYPES(ii).name)
        dtIdx = ii;
        break;
    end
end

for s = scansToLoad
    dataTYPES(dtIdx).eventAnalysisParams(s).detrend = 0;
    dataTYPES(dtIdx).eventAnalysisParams(s).inhomoCorrect = 1;
end
saveSession;

try
    [data, ~] = er_voxelData(hvol, roiPth, scansToLoad, recomputeFlag, preserveCoordsFlag);
catch ME
    warning(ME.message);
    data = [];
end
        
%% Reset path
cd(currPth);

end