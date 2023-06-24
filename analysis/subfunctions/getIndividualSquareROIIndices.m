function squareIdx = getIndividualSquareROIIndices(pths, prfParams, roiType, hemi)

squareLoc   = [];
useStimFlag = true;
offset      = 0.41; % deg

tmp = strsplit(roiType,'_');
nrSquares = str2double(tmp{1}(regexp(tmp{1},'[0-9]')));
sqArea = tmp{2}(regexp(tmp{2},'[0-9]'));
if strcmp(sqArea,'2')
    relativeSize = 'small';
elseif strcmp(sqArea,'4')
    relativeSize = 'large';
end

loc = getStimLocDeg(pths.projectDir, pths.subjnr, str2double(pths.session(end)), useStimFlag, ...
    squareLoc, hemi, relativeSize, offset);

if useStimFlag && strcmp(hemi, 'rh')
    loc.yBounds = -1.*loc.yBounds;
end

minX = min(loc.xBounds(:));
maxX = max(loc.xBounds(:));
minY = min(loc.yBounds(:));
maxY = max(loc.yBounds(:));

if ~useStimFlag && strcmp(hemi, 'rh')
    minX = -minX;
    minY = -minY;
    maxX = -maxX;
    maxY = -maxY;
    
    loc.xBounds = -1.*loc.xBounds;
    loc.yBounds = -1.*loc.yBounds;
end


for sq = 1:size(loc.xBounds,1)
    squareIdx{sq} = (prfParams.(hemi).x0 >= loc.xBounds(sq,1)) & ...
        (prfParams.(hemi).x0 <= loc.xBounds(sq,2)) & ...
        (prfParams.(hemi).y0 >= loc.yBounds(sq,1)) & ...
        (prfParams.(hemi).y0 <= loc.yBounds(sq,2));
end


end