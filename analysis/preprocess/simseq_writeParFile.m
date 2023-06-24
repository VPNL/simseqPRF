function [] = simseq_writeParFile(sessionDir)
% Function to create a parfile for mrVista GLM purposes
%
% simseq_writeParFile(subjnr, sessionDir)
%
% Example:
% sessionDir = pwd;
% simseq_writeParFile(subjnr, sessionDir)

cd(sessionDir);

% Get number of runs/files
d = dir(fullfile(sessionDir,'/Stimuli/simseqPRF*.mat'));

colors = {[1 1 1]; [0.545 0 0]; [0 0 0.545]}; % white, red, blue
condNames = {'blank', 'seq', 'sim'};

for run = 1:length(d)
    
    par = struct();
    par.onset = []; % onset in seconds
    par.code  = []; % condition code (blank=0, seq = 1, sim=2)
    
    load(fullfile(d(run).folder, d(run).name));
    stimTR = trial(2).onset-trial(1).onset;

    % Create vectors instead of structs
    for ii = 1:length(trial)
        onset(ii) = (ii-1)*stimTR;
        code(ii) = trial(ii).cond;
    end

    codeBinary = code;
    codeBinary(code>0 & code<5) = 1;
    codeBinary(code==5) = 1;

    % Define initial blank period
    par.onset(1) = onset(1);
    par.code(1)  = code(1);
    par.cond{1}  = condNames{1};
    par.color{1} = colors{1};
    
    % Define first condition onset
    transitions = find(diff(codeBinary>0));
    startStim = transitions(1);
    
    allTransitions = startStim:32:length(code);
    
    for jj = 1:length(allTransitions)-1
        idx = allTransitions(jj)+1;
        par.onset(jj+1) = onset(idx);
        if code(idx)==0
            par.code(jj+1)  = 0;
            par.cond{jj+1}  = condNames{1};
            par.color{jj+1} = colors{1};
        elseif code(idx)<5
            par.code(jj+1)  = 1;
            par.cond{jj+1}  = condNames{2};
            par.color{jj+1} = colors{2};
        elseif code(idx)==5
            par.code(jj+1)  = 2;
            par.cond{jj+1}  = condNames{3};
            par.color{jj+1} = colors{3};
        end
    end
    
    outFile = sprintf('simseq_run%d.par', run);
    saveDir = fullfile(sessionDir,'Stimuli', 'Parfiles');
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    
    fidout = fopen(fullfile(saveDir,outFile),'w');
    for pp=1:length(par.onset)
        fprintf(fidout,'%d \t %d \t', par.onset(pp), par.code(pp));
        fprintf(fidout,'%s \t', par.cond{pp});
        fprintf(fidout,'%1.3f %1.3f %1.3f \n', par.color{pp});
    end
    fclose(fidout);
    fclose('all');
    

%     save(fullfile(saveDir,outFile));
    fprintf('Wrote parfile %s successfully.\n', fullfile(saveDir,outFile));

end
