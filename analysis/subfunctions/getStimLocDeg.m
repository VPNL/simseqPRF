function loc = getStimLocDeg(projectDir, subjnr, sessionNr, useStimFlag, forcedGridLoc, hemi, relativeSize, offset)
% Define pRF XY centers that fall within stimulus squares

if ~exist('useStimFlag','var') || isempty(useStimFlag)
    useStimFlag = true;
    forcedGridLoc = [];
end

if ~exist('relativeSize','var') || isempty(relativeSize)
    relativeSize = 'small';
end

if ~exist('offset','var') || isempty(offset)
    offset = 0;
end


if strcmp(relativeSize,'small')
    sizeIdx = 1;
else
    sizeIdx = 2;
end

%% Set up directories
pths       = getSubjectPaths(projectDir,subjnr,sessionNr);
flipLR     = true;
if strcmp(hemi, 'lh'), hemiIdx = 1; else hemiIdx = 2; end


if useStimFlag
    % Load stimulus file
    load(pths.stimFiles{1})
    
    % Get square size
    halfSquareSize = params.squareSizeDeg(sizeIdx)/2;

    % Create stimulus center grid
    gr = params.gridSpaceDeg*([1:params.grid]-ceil(params.grid/2));
    [X,Y]=(meshgrid(gr));X=X';Y=-Y';
    X = X(:); Y = Y(:);
    % If we only select one center point per hemi, retrieve it here..
    if isfield(params,'locationsIdx') && ~isempty(params.locationsIdx)
        X = X(params.locationsIdx(hemiIdx));
        Y = Y(params.locationsIdx(hemiIdx));  
    end
    
    % If we stored location in deg per square, use that
%     if isfield(params, 'loc') && ~isempty(params.loc)
%         stimXYCenter = repmat([X; Y],1,size(params.loc{sizeIdx},1));
%         stimDiagCenter = [params.loc{sizeIdx}(:,[1,2],hemiIdx);
%                          params.loc{sizeIdx}(:,[3,4],hemiIdx)];
%         stimXbounds = params.loc{sizeIdx}(:,[1,3],hemiIdx) + offset + ([-1 1].*params.squareSizeDeg/2);
%         stimYbounds = params.loc{sizeIdx}(:,[2,4],hemiIdx) + offset + ([-1 1].*params.squareSizeDeg/2);
%     else
    pos = [];
    if size(flip.pos{end},1)>1
        if strcmp(hemi, 'lh'), hemiIdx = 1; else hemiIdx = 2; end
        for ii =1:length(flip.pos) 
            tmp = flip.pos{ii};
            if ~isempty(tmp)
                if length(tmp)==2
                    pos(ii) = tmp(hemiIdx);
                elseif length(tmp)>2
%                     start = (idx-1)*(length(tmp)/2);
%                     stop = idx*(length(tmp)/2);
%                     pos(ii,:) = tmp(start:stop);
                    pos(ii) = NaN;  
                elseif ~isempty(tmp) 
                    pos(ii) = tmp;
                end
            else
                pos(ii) = NaN;
            end
        end
        pos(isnan(pos))=[];
        tmp=unique(pos);
        if length(X(:))>10
            if strcmp(hemi, 'lh'), 
                gridloc=tmp(tmp>0 & tmp<60); else
                gridloc=tmp(tmp>60); end
        else
           gridloc = 1:length(tmp(tmp>0));
           X = repmat(X,1,length(gridloc));
           Y = repmat(Y,1,length(gridloc));
           
           relCenterX = ones(1,length(gridloc));
           relCenterY = ones(1,length(gridloc));
           relCenterX([1:2:end])= -1;
           relCenterY([1:(length(gridloc)/2)]) = -1;
           X = X + relCenterX.*(offset + halfSquareSize);
           Y = Y + relCenterY.*(offset + halfSquareSize);       
        end
    else
        gridloc = unique(cell2mat(flip.pos));
        gridloc = gridloc(gridloc>0);
    end
    
else % if forcing grid (pilot 2 only)
    gridSpaceDeg = sqrt(1.9);
    grid = 11;
    gr = gridSpaceDeg*([1:grid]-ceil(grid/2));
    [X,Y]=(meshgrid(gr));X=X';Y=-Y';
    params.squareSizeDeg = 1;
    halfSquareSize = params.squareDeg/2;
    gridloc = forcedGridLoc;
end

for ii = 1:length(gridloc)
    stimXYCenter(ii,:) = [X(gridloc(ii));Y(gridloc(ii))];
    stimXbounds(ii,:)  = stimXYCenter(ii,1) + ([-1 1].*halfSquareSize);
    stimYbounds(ii,:)  = stimXYCenter(ii,2) + ([-1 1].*halfSquareSize);
    stimDiagCenter(ii) = sqrt(sum([stimXYCenter(ii,1).^2,stimXYCenter(ii,2).^2]));
end


stimAngleCenter = rad2deg(atan(abs(stimXYCenter(:,2))./abs(stimXYCenter(:,1))))';
% stimAngleBound = rad2deg(atan(abs(5)./abs(stimYbounds)));

loc = struct();
loc.units = 'deg';
loc.grid = gr;
loc.gridloc = gridloc;
loc.xyCenter = stimXYCenter;
loc.xBounds = stimXbounds;
loc.yBounds = abs(stimYbounds);
loc.eccCenter  = stimDiagCenter;
loc.angCenter  = stimAngleCenter;

return
        