function trial = createTrialSet_varyArea(params, ims, scan,trialSetsDir)
% function to create specific trial structure for simseq PRF experiment
% outputs struct called 'trial' with params for each trial, which is
% each screen flip, which is 60 Hz (so 1/60 = 0.016667 s)  
  
trial = struct('cond',num2cell(length(params.seq)), ... % condition nr
    'onset',[],...      % Onset of trial (s), relative to start = 0
    'location',[], ...  % location in degrees of visual angle
    'rsvptarg', [], ... % if there is a 1-back repeat (i.e. a target)
    'rsvpID',[], ...    % RSVP letter identity
    'pos', [], ...      % position of stimulus square in grid?
    'block', [], ...    % block number, conditions SEQ=1-6,SIM=8,9,10,12
    'condName', [], ... % Name of the condition
    'dur', []);         % Duration of trial (s)

% Selects condition, block, onset, duration
imageNr = 1;
for n = 1:length(params.seq)
    
    trial(n).cond     = params.seq(:,n);
    trial(n).block    = params.blockNr(n);
    trial(n).onset    = params.seqTiming(n);
    trial(n).dur      = params.refreshRateHz;
    if trial(n).block==0
        trial(n).condName = 'blank';
    else
        trial(n).condName = params.condNames{trial(n).block};
    end
    
    %%%% Set squares
    switch (trial(n).cond(1))
        case {2,3,4,5}  % SEQ 4 squares each 4 deg2
            idx1 = (trial(n).cond(1)==[2:5]);
            idx2 = (trial(n).cond(2)==[2:5]);
            
            trial(n).pos        = cat(2, params.pos{1}(idx1),...   % upper right
                                         params.pos{2}(idx2));          % lower left
            trial(n).location   = cat(2,params.loc{1}(:,idx1,1), ... % upper right
                                        params.loc{1}(:,idx2,2));        % lower left
            
            if trial(n).cond == trial(n-1).cond
                trial(n).IDs = trial(n-1).IDs;
            else
                nrSquaresL     = length(trial(n).pos(1));
                nrSquaresR     = length(trial(n).pos(2));
                imageNrL       = imageNr:(imageNr+nrSquaresL-1);
                imageNr        = imageNr+nrSquaresL;
                
                imageNrR       = imageNr:(imageNr+nrSquaresR-1);
                trial(n).IDs   = cat(2,ims.order(imageNrL),ims.order(imageNrR))';
                imageNr        = imageNr+nrSquaresR;
            end
            
        case {8,9,10,11}
            idx1 = (trial(n).cond(1)==[8:11]);
            idx2 = (trial(n).cond(2)==[8:11]);
            
            trial(n).pos        = cat(2, params.pos{1}(idx1), ... % upper right
                                         params.pos{2}(idx2));    % lower left
            trial(n).location   = cat(2,params.loc{2}(:,idx1,1), ... % upper right
                                        params.loc{2}(:,idx2,2));    % lower left
            
            if trial(n).cond == trial(n-1).cond
                trial(n).IDs = trial(n-1).IDs;
            else
                nrSquaresL     = length(trial(n).pos(1));
                nrSquaresR     = length(trial(n).pos(2));
                imageNrL       = imageNr:(imageNr+nrSquaresL-1);
                imageNr        = imageNr+nrSquaresL;
                imageNrR       = imageNr:(imageNr+nrSquaresR-1);
                trial(n).IDs   = cat(2,ims.order(imageNrL),ims.order(imageNrR))';
                imageNr        = imageNr+nrSquaresR;
            end
            
        case 6  % SIM 4 squares
            trial(n).pos    = cat(2,params.pos{1},...  % upper right & % lower left
                                    params.pos{2});  % upper right & % lower left
            
            trial(n).location = cat(2,params.loc{1}(:,:,1),params.loc{1}(:,:,2));  % upper right &  % lower left
            if trial(n).cond == trial(n-1).cond
                trial(n).IDs = trial(n-1).IDs;
            else
                nrSquaresL     = length(params.pos{1});
                nrSquaresR     = length(params.pos{2});
                imageNrL       = imageNr:(imageNr+nrSquaresL-1);
                imageNr        = imageNr+nrSquaresL;
                imageNrR       = imageNr:(imageNr+nrSquaresR-1);
                trial(n).IDs   = cat(2,ims.order(imageNrL),ims.order(imageNrR))';
                imageNr        = imageNr+nrSquaresR;
            end
            
        case 12  % SIM 4 squares
            trial(n).pos    = cat(2, params.pos{1},...  % upper right
                                     params.pos{2});  %  lower left
            trial(n).location = cat(2,params.loc{2}(:,:,1),params.loc{2}(:,:,2));  % upper right &  % lower left
            if trial(n).cond == trial(n-1).cond
                trial(n).IDs = trial(n-1).IDs;
            else
                nrSquaresL     = length(params.pos{1});
                nrSquaresR     = length(params.pos{2});
                imageNrL       = imageNr:(imageNr+nrSquaresL-1);
                imageNr        = imageNr+nrSquaresL;
                imageNrR       = imageNr:(imageNr+nrSquaresR-1);
                trial(n).IDs   = cat(2,ims.order(imageNrL),ims.order(imageNrR))';
                imageNr        = imageNr+nrSquaresR;
            end
            
        case 1 % BLANKS
            trial(n).IDs = [0;0];
            trial(n).pos = [0;0];
    end
    
    %%%% Set RSVP target
    % Pre-set the first trial RSVP
    if n == 1
        trial(n).rsvptarg  = 0;
        trial(n).rsvpID    = datasample(1:length(params.str),1,'Replace',false);
        trial(n).rsvpColor = params.taskColor(1,:,:);
        
        % Update RSVP stream depending on task freq
    else
        if mod(n,params.taskFreq)==0
            % Set RSVP color (alternating black and white for visibility)
            [~,bi] = intersect(params.taskColor,trial(n-1).rsvpColor, 'rows');
            trial(n).rsvpColor = params.taskColor(setdiff(1:size(params.taskColor,1),bi),:);
            if rand<params.taskProb % if 1-back target repeat letter
                trial(n).rsvptarg = 1;               % Set target to 1
                trial(n).rsvpID = trial(n-1).rsvpID; % and set rsvp ID to previous RSVP letter ID
            else % no 1-back target: regular update
                trial(n).rsvptarg = 0;
                % get previous rsvp ID to avoid getting triple doubles
                prevRSVPID = trial(n-1).rsvpID;
                trial(n).rsvpID = datasample(setdiff(1:length(params.str),prevRSVPID),1,'Replace',false);
            end
        else % Repeat previous RSVP letter, make it a non-target
            trial(n).rsvptarg  = 0;
            trial(n).rsvpID    = trial(n-1).rsvpID;
            trial(n).rsvpColor = trial(n-1).rsvpColor;
        end
    end
end    
    % Save trialSet 
    save(fullfile(trialSetsDir, ...
        [datestr(now,30) '_trials' num2str(scan.runNum) '_v' num2str(scan.versionNr) '.mat']),'trial');
end