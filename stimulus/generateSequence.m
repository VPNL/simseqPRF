function seqMaster = generateSequence(locations, dt, durations, ...
                        nrSeqRepetitions, conditionOrder, totalTrialsPerBlock, ...
                        hemifield)
% Stimulus order:
% Simultaneous trial order is one frame of 6 = all stim simultaneous, 
%   rest of frames are 1 = blank count, 1x4 squares 0.217s+33 ms gap/trial
% Sequential trial order: randomize single square presentation across 4 
% locations [2-5]: 4x 1 square of 0.217 s  + 33 ms gap per trial
                    
framesPerSIMStim = durations.simSingleFrameDurS/dt;
framesPerSEQStim = durations.seqSingleFrameDurS/dt;
framesPerBlank   = durations.blankDurS/dt;
framesPerGap     = durations.gap/dt;

% Number of positions
nLocs = size(locations);
% nLocs = size(locations,1);
% if strcmp(hemifield,'bothf')
%     tmp(:,1) = find(locations(:,1)>0)';
%     tmp(:,2) = setdiff(1:nLocs,tmp(:,1));
%     nLocs = size(tmp);
% end

% Randomize sequential order
seqOrder  = @(x,y) Shuffle(repmat((1:x)+1,1,y));

% Start after initial blanks
framesOfInitialBlank = durations.initialBlankDurS/dt;

% Run sequence generation twice if both hemifield
if strcmp(hemifield,'bothf'), rep = [1,2]; else rep = 1; end

for r = rep
    sequence = [];
    sequence(1:framesOfInitialBlank) = 1; % blank
    counter  = framesOfInitialBlank+1;
    
    for cc = conditionOrder
        % Simultaneous, show all squares at once
        if cc == 2
            if r ==1
                for tt = 1:totalTrialsPerBlock
                    simOrder = ones(1,nLocs(1));
                    if tt == 1
                        onsetSim = 1;
                    else
                        if (totalTrialsPerBlock>3)
                            onsetSim = randperm(nLocs(1),1);
                        else
                            onsetSim = 1;
                        end
                    end
                    simOrder(onsetSim) = nLocs(1)+2;
                    tmpSIM = repmat(simOrder', 1, framesPerSIMStim);
                    tmpSIM = [tmpSIM, ones(size(tmpSIM,1),nrSeqRepetitions*framesPerGap)];
                    tmpSIM = tmpSIM';
                    sequence(counter : (counter+length(tmpSIM(:))-1)) = tmpSIM(:);
                    counter = counter+length(tmpSIM(:));                
                end
            else
                sequence = seqMaster(1,:);
            end
               
            % Sequential, show single square per frame
        elseif cc == 1
            for tt = 1:totalTrialsPerBlock
                tmp = repmat([seqOrder(nLocs(1),nrSeqRepetitions)]', 1, framesPerSEQStim);
                tmp = [tmp, ones(size(tmp,1),framesPerGap)];
                tmp = tmp';
                sequence(counter : (counter+length(tmp(:))-1)) = tmp(:);
                counter = counter+length(tmp(:));
            end
            
        % Simultaneous, show all squares at once
        elseif cc == 3
            for tt = 1:totalTrialsPerBlock
                simCondNr = nLocs(1)+2;
                tmp = simCondNr*ones(nLocs(1), framesPerSIMStim);
                tmp = [tmp, ones(size(tmp,1),framesPerGap)];
                tmp = tmp';
                sequence(counter : (counter+length(tmp(:))-1)) = tmp(:);
                counter = counter+length(tmp(:));
            end
        end

        % add 8 s of blank after block presentation
        sequence(counter:(counter+framesPerBlank-1))  = ones(1,framesPerBlank);
        counter = counter+framesPerBlank;
    end 
    seqMaster(r,:) = sequence;
end

return