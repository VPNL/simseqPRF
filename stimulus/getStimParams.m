function params = getStimParams(trialType, varargin)
% Function to generate stimulus parameters and sequence for a variety of
% sim-seq experiments.
%       params = getStimParams(trialType, varargin)
%
% INPUTS
% trialType       : (str) what kind of trial sequence do you want? Choose:
%                   'Kastner98original':    4x squares, 216ms, each 2x2 deg
%                   'varyTransients':       4x squares, 216ms, 2x2 deg, 1-4x
%                                           more transients by dividing 216
%                                           by 1-4 for SEQ
%                   'varyIndivStimArea':    4x squares, 216ms, each 1x1, 2x2,
%                                           or 2.5x2.5 deg
%                   'varyStimDur':          4x squares, 54, 216, 864, 3456ms,
%                                           each 1x1 deg
%                   'varyNrStim':            1, 2, 3, or 4 squares
%                   'varyStimDur_x_Area':   vary total area (4 or 16 1x1 deg
%                                           squares) x 4 durations
%                                           (54-3456ms)
%                   'varyStimDur_x_Transients': 4 durations (54-3456ms) x
%                                            1-4 transients for 4x 1x1 deg
%                   'varyTransients_x_Area': 1-4 transients for 4, 9 or
%                                            16 squares of 1x1 deg
%                   'varyTransients_x_Area_x_StimDur': 1 or 4 transients,
%                                           4 or 16 squares, 216 or 3456 ms
% [nrGridPoints]  : (int) nr of points along single axes to construct grid (default=11)
% [gridSpaceDeg]  : (int) space between grid points (default=sqrt(2))
% [hemifield]     : (str) what hemifield to present stimuli,
%                       choose: 'llhf' for lower left,
%                               'urhf' for upper right (=default)
%                               'bothf' for both lower left and upper right
% [verbose]       : (bool) plot figures or not (default=true)

% OUTPUTS
% params          : (struct) parameters with stimulus size, location,
%                     duration
%
%
% Written by EK, 2021, Stanford U

%% Parse inputs
p = inputParser;
p.addRequired('trialType', @ischar);
p.addParameter('nrGridPoints',11, @isnumeric);
p.addParameter('gridSpaceDeg',sqrt(2), @isnumeric);
p.addParameter('hemifield','urhf', @(x) any(ismember(x,{'llhf','urhf','bothf'})));
p.addParameter('verbose', true, @islogical);
p.parse(trialType,varargin{:});

% Rename variables
trialType           = p.Results.trialType;
nrGridPoints        = p.Results.nrGridPoints;
gridSpaceDeg        = p.Results.gridSpaceDeg;
hemifield           = p.Results.hemifield;
verbose             = p.Results.verbose;

%% Create grid
gr = gridSpaceDeg*([1:nrGridPoints]-ceil(nrGridPoints/2));
[X,Y]=(meshgrid(gr));
X=X';Y=-Y';
XY = [X(:), Y(:)];

%% Set defaults given Kastner et al 98 paper

% Stimulus durations
params.durations.simSingleFrameDurS = 0.216; % Duration of stimulus frame with all squares (in seconds)
params.durations.seqSingleFrameDurS = 0.216; % Duration of stimulus frame with single square (in seconds)

% Show squares only once in sequential condition, same duration as
% individual squares in simultaneous condition
params.nrSeqRepetitions = 1;

% nr of trials
params.totalTrialsPerBlock = 1; % count

params.nrOfBlocks = 2; % one Seq, one Sim;

%number of locations of individual squares, [x,y] degrees from central fixation
params.locationIdx = [10, 11, 21, 22]; % for 4x 1 deg in upper right corner at ~8+ deg eccentricity

switch hemifield
    case 'urhf'
        % upper right hemi field, 4x 1 deg2 at ~8+ deg eccentricity:
        % use default
    case 'llhf'
        % lower left hemi field, 4x 1 deg2 at ~8+ deg eccentricity:
        params.locationIdx = params.locationIdx+90; %[100, 101, 111, 112];
    case 'bothf'
        % both upper right and lower left, 4x 1 deg2 at ~8+ deg eccen:
        params.locationIdx = [params.locationIdx, params.locationIdx+90];
end

params.locations = XY(params.locationIdx,:);

% Dimensions of each square, in degrees visual angle
params.stimWidthDeg  = 1; %deg
params.stimHeightDeg = 1; %deg
params.stimAreaDeg2  = params.stimWidthDeg*params.stimHeightDeg; % deg2

% Attach string to filename
params.postFix = '';

%% Update default params depending on trial type
switch trialType
    case 'Kastner98original'
        params.locationIdx = [33, 31, 11, 9]; % equal to eccen deg [[4.2426,4.2426]; [4.2426,7.0711]; [7.0711,4.2426]; [7.0711,7.0711]];
        if strcmp(hemifield, 'llhf'), params.locationIdx = params.locationIdx+80; end
        if strcmp(hemifield, 'bothf'), params.locationIdx = [params.locationIdx, params.locationIdx+80]; end
        params.locations = XY(params.locationIdx,:);
        params.stimWidthDeg  = 2; %deg
        params.stimHeightDeg = 2; %deg
        params.stimAreaDeg2  = params.stimWidthDeg*params.stimHeightDeg; % deg2
        
    case 'varyTransients'
        params.stimWidthDeg  = 2; %deg
        params.stimHeightDeg = 2; %deg
        params.stimAreaDeg2  = params.stimWidthDeg*params.stimHeightDeg; % deg2
        params.nrSeqRepetitions = [1,2,3,4]; % show single square 1x (=216ms), 2x (=108 ms), ... etc
        params.durations.simSingleFrameDurS = [0.216, 0.216, 0.216, 0.216]; % sec
        params.durations.seqSingleFrameDurS = params.durations.simSingleFrameDurS./params.nrSeqRepetitions;
        
    case 'varyIndivStimArea'
        params.locationIdx = [33, 31, 11, 9];      % for 3 deg squares choose [44, 41, 11, 8];
        if strcmp(hemifield, 'llhf'), params.locationIdx = params.locationIdx+80; end
        if strcmp(hemifield, 'bothf'), params.locationIdx = [params.locationIdx, params.locationIdx+80]; end
        params.locations   = XY(params.locationIdx,:); %deg
        params.stimWidthDeg  = [1, 2, 2.5]; %deg
        params.stimHeightDeg = [1, 2, 2.5]; %deg
        params.stimAreaDeg2 = params.stimWidthDeg.*params.stimHeightDeg; % deg2
        
    case 'varyStimDur'
        params.durations.simSingleFrameDurS = [0.054, 0.216, 0.864, 3.456]; %(in seconds)
        params.durations.seqSingleFrameDurS = params.durations.simSingleFrameDurS./params.nrSeqRepetitions;
        
    case 'varyStimDur_x_Area16'
        params.locationIdx = {[10, 11, 21, 22], ... 4x
            [41:44, 30:33, 19:22, 8:11]}; %16x 1 deg
        if strcmp(hemifield, 'llhf')
            params.locationIdx{1} = params.locationIdx{1}+90;
            tmp = reshape(params.locationIdx{2},4,4)' + [37;59;81;103];
            params.locationIdx{2} = tmp(:)';
        elseif strcmp(hemifield, 'bothf')
            params.locationIdx{1} = [params.locationIdx{1}, params.locationIdx{1}+90];
            tmp = reshape(params.locationIdx{2},4,4)' + [37;59;81;103];
            params.locationIdx{2} = [params.locationIdx{2}, tmp(:)'];
        end
        params.locations = {};
        params.locations{1}   = XY(params.locationIdx{1},:);
        params.locations{2}   = XY(params.locationIdx{2},:);
        params.durations.simSingleFrameDurS = [0.054, 0.216, 0.864, 3.456]; % seconds
        params.durations.seqSingleFrameDurS = params.durations.simSingleFrameDurS./params.nrSeqRepetitions; % seconds
        
    case 'varyStimDur_x_Transients'
        params.nrSeqRepetitions = [1,2,3,4]; % show single square 1x (=216 ms), 2x (=108 ms), 3x (=72 ms), 4x (=54 ms)
        params.durations.simSingleFrameDurS = [0.216, 0.864, 3.456]; % seconds
        params.durations.seqSingleFrameDurS = params.durations.simSingleFrameDurS./params.nrSeqRepetitions'; % seconds
        params.durations.simSingleFrameDurS = repmat(params.durations.simSingleFrameDurS, 4,1); % seconds
        
    case 'varyTransients_x_Area16'
        params.locationIdx = {[10, 11, 21, 22], ... 4x
            [31:33, 20:22, 9:11], ... % 9
            [41:44, 30:33, 19:22, 8:11]}; %16x 1 deg
        if strcmp(hemifield, 'llhf')
            params.locationIdx{1} = params.locationIdx{1}+90;
            tmp = reshape(params.locationIdx{2},3,3)' + [58;80;102];
            params.locationIdx{2} = tmp(:)';
            tmp = reshape(params.locationIdx{3},4,4)' + [37;59;81;103];
            params.locationIdx{3} = tmp(:)';
        elseif strcmp(hemifield, 'bothf')
            params.locationIdx{1} = [params.locationIdx{1}, params.locationIdx{1}+90];
            tmp = reshape(params.locationIdx{2},3,3)' + [58;80;102];
            params.locationIdx{2} = [params.locationIdx{2}, tmp(:)'];
            tmp = reshape(params.locationIdx{3},4,4)' + [37;59;81;103];
            params.locationIdx{3} = [params.locationIdx{3}, tmp(:)'];
        end
        params.locations = {};
        params.locations{1}   = XY(params.locationIdx{1},:);
        params.locations{2}   = XY(params.locationIdx{2},:);
        params.locations{3}   = XY(params.locationIdx{3},:);
        
        % Option 1: for 4 x 1 deg + 2 deg space
        %         locationIdx = [41, 44, 11, 8];
        %         postFix = '_w2degspace';
        
        % Option 2: for 4 x 1 deg + 1 deg space
        %         locationIdx = [33, 31, 11, 9];
        %         postFix = '_w1degspace';
        %         locations   = XY(locationIdx,:);
        
        % Show squares only once in sequential condition, same duration as
        % individual squares in simultaneous condition
        params.nrSeqRepetitions = [1,4];
        
        % Duration of stimulus frame with all squares (in seconds)
        params.durations.simSingleFrameDurS = [0.216, 0.216];
        
        % Duration of stimulus frame with single square (in seconds)
        params.durations.seqSingleFrameDurS = params.durations.simSingleFrameDurS./params.nrSeqRepetitions;
        
    case 'varyNrStim'
        
        %number of locations of individual squares
        % [x,y] degrees from central fixation
        % Option 1: for 1 x 1 deg + 2 deg space
        params.locationIdx = {params.locationIdx(1), ... % 1 stim
            params.locationIdx(1:2), ...     % 2 stim
            params.locationIdx(1:3), ...     % 3 stim
            params.locationIdx(1:4)};        % 4 stim
        if strcmp(hemifield,'bothf')
            params.locationIdx = {params.locationIdx(1,5), ... % 1 stim
                [params.locationIdx(1:2),params.locationIdx(5:6)], ...     % 2 stim
                [params.locationIdx(1:3),params.locationIdx(5:7)], ...     % 3 stim
                [params.locationIdx(1:4),params.locationIdx(5:8)]};        % 4 stim
        end
        params.locations = {};
        params.locations{1}   = XY(params.locationIdx{1},:);
        params.locations{2}   = XY(params.locationIdx{2},:);
        params.locations{3}   = XY(params.locationIdx{3},:);
        params.locations{4}   = XY(params.locationIdx{4},:);
        
        params.postFix = {'_1stim','_2stim','_3stim','_4stim'};
        params.stimWidthDeg  = 2; %deg
        params.stimHeightDeg = 2; %deg
        params.stimAreaDeg2 = params.stimWidthDeg*params.stimHeightDeg; % deg2
        
    case 'varyTransients_x_Area_x_StimDur'
        params.locationIdx      = {[31, 32, 20, 21], ... 4x 1 deg
            [41:44, 30:33, 19:22, 8:11]}; %16x 1 deg
        if strcmp(hemifield, 'llhf')
            params.locationIdx{1} = params.locationIdx{1}+70;
            tmp = reshape(params.locationIdx{2},4,4)' + [37;59;81;103];
            params.locationIdx{2} = tmp(:)';
        elseif strcmp(hemifield, 'bothf')
            params.locationIdx{1} = [params.locationIdx{1}, params.locationIdx{1}+70];
            tmp = reshape(params.locationIdx{2},4,4)' + [37;59;81;103];
            params.locationIdx{2} = [params.locationIdx{2}, tmp(:)'];
        end
        params.locations = {};
        params.totalTrialsPerBlock = {};
        params.nrOfBlocks = {};
        
        params.locations{1}     = XY(params.locationIdx{1},:);   % deg
        params.locations{2}     = XY(params.locationIdx{2},:);   % deg
        params.nrSeqRepetitions = [1,4,1]; % counts
        params.durations.simSingleFrameDurS = [0.2, 0.2,0.2*16]; % seconds
        params.durations.seqSingleFrameDurS = params.durations.simSingleFrameDurS./params.nrSeqRepetitions;  % seconds
        params.totalTrialsPerBlock{1} = [8, 8, 1];
        params.totalTrialsPerBlock{2} = [2, 2, 1];
        params.nrOfBlocks{1} = [7, 5, 4; 7, 5, 16]; % seq; sim -4 squares
        params.nrOfBlocks{2} = [7, 5, 1; 7, 5, 16]; % seq; sim -16 squares
        
    case 'varyStimDur_x_Area9'
        params.stimWidthDeg  = 1.5; %deg
        params.stimHeightDeg = 1.5; %deg
        params.stimAreaDeg2  = params.stimWidthDeg*params.stimHeightDeg; % deg2
        
        params.locationIdx      = {[71,75,11,15], ...%{[20, 21, 31, 32], ... 4x 2 deg
            [71:2:75, 41:2:45, 11:2:15]}'; %[41:43, 30:32, 19:21]}; % 9x 2 deg
        if strcmp(hemifield, 'llhf')
            params.locationIdx{1} = params.locationIdx{1}+140;
            tmp = reshape(params.locationIdx{2},3,3)' + 140; %[38;60;82];
            params.locationIdx{2} = tmp(:)';
        elseif strcmp(hemifield, 'bothf')
            params.locationIdx{1} = [params.locationIdx{1}, params.locationIdx{1}+140];
            tmp = reshape(params.locationIdx{2},3,3)' + 140;% [38;60;82];
            params.locationIdx{2} = [params.locationIdx{2}, tmp(:)'];
        end
        params.locations = {};
        params.totalTrialsPerBlock = {};
        params.nrOfBlocks = {};
        
        params.locations{1}     = XY(params.locationIdx{1},:);   % deg
        params.locations{2}     = XY(params.locationIdx{2},:);   % deg
        params.nrSeqRepetitions = [1,1,1]; % counts
        params.durations.simSingleFrameDurS = [0.1, 0.2, 0.4]; % seconds
        params.durations.seqSingleFrameDurS = params.durations.simSingleFrameDurS./params.nrSeqRepetitions;  % seconds
        params.totalTrialsPerBlock{1} = [20, 10, 6]; % 4 squares, 200 ms, 100ms amd 400ms presentation
        params.totalTrialsPerBlock{2} = [10, 5, 3]; % 9 squares, 200 ms, 100ms amd 400ms presentation
        params.nrOfBlocks{1} = [1, 1, 1; 1, 1, 1; 1 1 1]; % seq orig dur; sim orig dur -4 squares
        params.nrOfBlocks{2} = [1, 1, 1; 1, 1, 1; 1 1 1]; % seq orig dur; sim orig dur -9 squares
        
    case 'varyStimDur_x_Area9_Inner4'
        params.stimWidthDeg  = 1.5; %deg
        params.stimHeightDeg = 1.5; %deg
        params.stimAreaDeg2  = params.stimWidthDeg*params.stimHeightDeg; % deg2
        
        params.locationIdx      = {[57,59,27,29], ...%{[20, 21, 31, 32], ... 4x 2 deg
            [71:2:75, 41:2:45, 11:2:15]}'; %[41:43, 30:32, 19:21]}; % 9x 2 deg
        if strcmp(hemifield, 'llhf')
            params.locationIdx{1} = params.locationIdx{1}+140;
            tmp = reshape(params.locationIdx{2},3,3)' + 140; %[38;60;82];
            params.locationIdx{2} = tmp(:)';
        elseif strcmp(hemifield, 'bothf')
            params.locationIdx{1} = [params.locationIdx{1}, params.locationIdx{1}+140];
            tmp = reshape(params.locationIdx{2},3,3)' + 140;% [38;60;82];
            params.locationIdx{2} = [params.locationIdx{2}, tmp(:)'];
        end
        params.locations = {};
        params.totalTrialsPerBlock = {};
        params.nrOfBlocks = {};
        
        params.locations{1}     = XY(params.locationIdx{1},:);   % deg
        params.locations{2}     = XY(params.locationIdx{2},:);   % deg
        params.nrSeqRepetitions = [1,1,1]; % counts
        params.durations.simSingleFrameDurS = [0.1, 0.2, 0.4]; % seconds
        params.durations.seqSingleFrameDurS = params.durations.simSingleFrameDurS./params.nrSeqRepetitions;  % seconds
        params.totalTrialsPerBlock{1} = [20, 10, 6]; % 4 squares, 200 ms, 100ms amd 400ms presentation
        params.totalTrialsPerBlock{2} = [10, 5, 3]; % 9 squares, 200 ms, 100ms amd 400ms presentation
        params.nrOfBlocks{1} = [1, 1, 1; 1, 1, 1; 1 1 1]; % seq orig dur; sim orig dur -4 squares
        params.nrOfBlocks{2} = [1, 1, 1; 1, 1, 1; 1 1 1]; % seq orig dur; sim orig dur -9 squares
        
    case 'varyStimDur_x_Area4'
        
        % Get center from, we want the squares extend size laterally, i.e.
        % moving center position outward.
         params.locationIdx = {};
         params.locations = {};
        if strcmp(hemifield, 'llhf')
            mainCenter(1,:) = XY(:,57); % LH
            params.locationIdx{1} = [1:4]; % LH
        elseif strcmp(hemifield, 'bothf')
            mainCenter(1,:) = XY(43,:); % LH
            mainCenter(2,:) = XY(183,:); % RH
            params.locationIdx{1} = [1:4]; % LH
            params.locationIdx{2} = [5:8]; % RH
        end
        
        gapFromCenter = 0.41; % deg;
        
        params.stimWidthDeg  = [2, 4]; %deg
        params.stimHeightDeg = [2, 4]; %deg
        params.stimAreaDeg2  = params.stimWidthDeg.*params.stimHeightDeg; % deg2
        for ll = 1:size(mainCenter,1)
            
            for ii = 1:length(params.stimWidthDeg)
                
                squareCenterX = 0.5*params.stimWidthDeg(ii);
                squareCenterY = 0.5*params.stimHeightDeg(ii);
                
                % upper right
                ctrX1 = mainCenter(ll,1) + gapFromCenter + squareCenterX;
                ctrY1 = mainCenter(ll,2) + gapFromCenter + squareCenterY;
                % lower right
                ctrX2 = mainCenter(ll,1) + gapFromCenter + squareCenterX;
                ctrY2 = mainCenter(ll,2) - gapFromCenter - squareCenterY;
                % upper left
                ctrX3 = mainCenter(ll,1) - gapFromCenter - squareCenterX;
                ctrY3 = mainCenter(ll,2) + gapFromCenter + squareCenterY;
                % lower left
                ctrX4 = mainCenter(ll,1) - gapFromCenter - squareCenterX;
                ctrY4 = mainCenter(ll,2) - gapFromCenter - squareCenterY;
                
                % Select pixels for stimulus
                params.locations{ll}(ii,:,1) = [ctrX1,ctrY1]; % deg
                params.locations{ll}(ii,:,2) = [ctrX2,ctrY2]; % deg
                params.locations{ll}(ii,:,3) = [ctrX3,ctrY3]; % deg
                params.locations{ll}(ii,:,4) = [ctrX4,ctrY4]; % deg
            end
            
        end
        
        params.totalTrialsPerBlock = {};
        params.nrOfBlocks = {};
        
        params.nrSeqRepetitions = [1,1,1]; % counts
        params.durations.simSingleFrameDurS = [0.1, 0.2, 1.0]; % seconds
        params.durations.seqSingleFrameDurS = params.durations.simSingleFrameDurS./params.nrSeqRepetitions;  % seconds
        params.totalTrialsPerBlock{1} = [23, 13, 3]; % 4 squares 2 deg2, 100, 200, 400, 800 & 1000ms presentation
        params.totalTrialsPerBlock{2} = [23, 13, 3]; % 4 squares 5 deg2, 100, 200, 400, 800 & 1000ms presentation
        params.nrOfBlocks{1} = [3, 3, 3; 3, 3, 3]; % seq orig dur; sim orig dur -4x2deg2 squares
        params.nrOfBlocks{2} = [3, 3, 3; 3, 3, 3]; % seq orig dur; sim orig dur -4x5deg2 squares
        
    case 'varyStimDur_x_Area4_2x2design'
        % Get center from, we want the squares extend size laterally, i.e.
        % moving center position outward.
         params.locationIdx = {};
         params.locations = {};
        if strcmp(hemifield, 'llhf')
            mainCenter(1,:) = XY(:,57); % LH
            params.locationIdx{1} = [1:4]; % LH
        elseif strcmp(hemifield, 'bothf')
            mainCenter(1,:) = XY(43,:); % LH
            mainCenter(2,:) = XY(183,:); % RH
            params.locationIdx{1} = [1:4]; % LH
            params.locationIdx{2} = [5:8]; % RH
        end
        
        params.gapFromCenter = 0.41; % deg; same as for first pilot
        
        params.stimWidthDeg  = [2, 4]; % deg
        params.stimHeightDeg = [2, 4]; % deg
        params.stimAreaDeg2  = params.stimWidthDeg.*params.stimHeightDeg; % deg2
        for ll = 1:size(mainCenter,1)
            
            for ii = 1:length(params.stimWidthDeg)
                
                squareCenterX = 0.5*params.stimWidthDeg(ii);
                squareCenterY = 0.5*params.stimHeightDeg(ii);
                
                % upper right
                ctrX1 = mainCenter(ll,1) +  params.gapFromCenter + squareCenterX;
                ctrY1 = mainCenter(ll,2) +  params.gapFromCenter + squareCenterY;
                % lower right
                ctrX2 = mainCenter(ll,1) +  params.gapFromCenter + squareCenterX;
                ctrY2 = mainCenter(ll,2) -  params.gapFromCenter - squareCenterY;
                % upper left
                ctrX3 = mainCenter(ll,1) -  params.gapFromCenter - squareCenterX;
                ctrY3 = mainCenter(ll,2) +  params.gapFromCenter + squareCenterY;
                % lower left
                ctrX4 = mainCenter(ll,1) -  params.gapFromCenter - squareCenterX;
                ctrY4 = mainCenter(ll,2) -  params.gapFromCenter - squareCenterY;
                
                % Select pixels for stimulus
                params.locations{ll}(ii,:,1) = [ctrX1,ctrY1]; % deg
                params.locations{ll}(ii,:,2) = [ctrX2,ctrY2]; % deg
                params.locations{ll}(ii,:,3) = [ctrX3,ctrY3]; % deg
                params.locations{ll}(ii,:,4) = [ctrX4,ctrY4]; % deg
            end
            
        end
        
        params.totalTrialsPerBlock = {};
        params.nrOfBlocks = {};
        
        params.nrSeqRepetitions = [1,1]; % counts
        params.durations.simSingleFrameDurS = [0.2, 1.0]; % seconds
        params.durations.seqSingleFrameDurS = params.durations.simSingleFrameDurS./params.nrSeqRepetitions;  % seconds
        params.totalTrialsPerBlock{1} = [8, 2]; %[13, 3]; % 4 squares 2 deg2, 200 & 1000ms presentation
        params.totalTrialsPerBlock{2} = [8, 2]; %[13, 3]; % 4 squares 16 deg2, 200 & 1000ms presentation
        params.nrOfBlocks{1} = [4 4; 4 4]; % seq orig dur; sim orig dur -4x2deg2 squares
        params.nrOfBlocks{2} = [4 4; 4 4]; % seq orig dur; sim orig dur -4x5deg2 squares
end

% Plot square locations
if verbose
    if strcmp(trialType,'varyStimDur_x_Area4') || strcmp(trialType,'varyStimDur_x_Area4_2x2design')
        figure(10); clf; set(gcf, 'Position', [227,387,1200,575]); hold all;
        for ll = 1:length(params.locations)
            for ii = 1:length(params.stimWidthDeg)
                subplot(1,length(params.locations),ii); hold all;
                scatter(params.locations{ll}(ii,1,:), params.locations{ll}(ii,2,:),'ko', 'lineWidth',2); hold on; axis square;
                plot([0 0], [-10,10],'k','lineWidth',0.25); plot([-10,10],[0 0],'k', 'lineWidth',0.25);
                text(0.2+params.locations{ll}(ii,1,:),params.locations{ll}(ii,2,:), ...
                    cellstr(string(params.locationIdx{ll})), 'FontSize',12);
                ylabel('Y locations (deg)'); xlabel('X locations (deg)');
                for jj = 1:size(params.locations{ll},3)
                    rectangle('Position',[params.locations{ll}(ii,1,jj)-(params.stimWidthDeg(ii)/2), ...
                        params.locations{ll}(ii,2,jj)-(params.stimHeightDeg(ii)/2), ...
                        params.stimWidthDeg(ii),params.stimHeightDeg(ii)])
                end
            end
        end
    
        
    else
        figure(10); clf; set(gcf, 'Position', [227,387,1200,575]); hold all;
        for ii = 1:length(params.locations)
            subplot(1,length(params.locations),ii);
            scatter(XY(:,1), XY(:,2),'ko', 'lineWidth',2); hold on; axis square;
            plot([0 0], [-10,10],'k','lineWidth',0.25); plot([-10,10],[0 0],'k', 'lineWidth',0.25);
            text(0.2+XY(:,1),XY(:,2),cellstr(string([1:nrGridPoints*nrGridPoints])), 'FontSize',12);
            ylabel('Y Grid locations'); xlabel('X Grid locations');
            
            scatter(params.locations{ii}(:,1),params.locations{ii}(:,2),'rx', 'LineWidth',3)
            
            for jj = 1:length(params.stimWidthDeg)
                for ll = 1:size(params.locations{ii},1)
                    rectangle('Position',[params.locations{ii}(ll,1)-(params.stimWidthDeg(jj)/2), ...
                        params.locations{ii}(ll,2)-(params.stimHeightDeg(jj)/2), ...
                        params.stimWidthDeg(jj),params.stimHeightDeg(jj)])
                end
            end
        end
    end
end