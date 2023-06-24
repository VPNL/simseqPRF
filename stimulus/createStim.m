function createStim()

sizeStimBlockWDeg = 4; % deg
sizeStimBlockHDeg = sizeStimBlockWDeg;
% screenWidth       = 38; %104; % cm (old CNI = 104cm at both 16ch and 32ch)
screenHeight      = 20; % cm 

scanCoil = 16;

if scanCoil == 32
    viewingDist = 272;   % cm; CNI = 270-273cm at 32ch, 265 at 16ch;
else
    viewingDist = 25; % cm (old measurement from Sonia's exp 265);
end

rng('default')
rng('shuffle')
maxCartoonImg = 63;

nCartoons = maxCartoonImg;
nSamples  = 500;
rmsThresh = 0.1;
blankIdx = zeros(nSamples,nCartoons);

normPix = @(x) double(x)./255;
rms = @(x) (sqrt(mean((x(:)-mean(x(:))).^2)));

for cr = 1:nCartoons
    
    cartoonNrToLoad = randi(maxCartoonImg);
    img = imread(fullfile(simseqRootPath, 'stimulus', 'stim_mri', 'cartoon',sprintf('pic%d.jpg',cartoonNrToLoad)));
    [imgH, imgW, ~] = size(img);
        
    pixPerDeg = pi*imgH / (atan(screenHeight/viewingDist/2)) / 360;    
    nrPixBlock = floor(sizeStimBlockHDeg*pixPerDeg);
    
    imgCut = NaN(nrPixBlock,nrPixBlock, size(img,3), nSamples);
    
    for ii = 1:nSamples
        pxStartH = randi(imgH - nrPixBlock);
        pxStartW = randi(imgW - nrPixBlock);
        
        pxCutH = (pxStartH:(pxStartH+nrPixBlock-1));
        pxCutW = (pxStartW:(pxStartW+nrPixBlock-1));
        
        thisCut = img( pxCutH, pxCutW, :);

        % Visualize squares:
        % figure(1); clf; 
        % imagesc(thisCut); 
        % title(sprintf('RMS: %1.3f',rms(normPix(rgb2gray(thisCut))))); 
        % axis image
        
        if rms(normPix(rgb2gray(thisCut)))>rmsThresh
            imgCut(:,:,:,ii) = thisCut;
        else
            blankIdx(ii,cr) = true;
        end
    end
    
    imgCut = uint8(imgCut);
    allImg(:,:,:,:, cr) = imgCut;
    
end

allImgRshp = reshape(allImg, size(allImg,1),size(allImg,2),size(allImg,3), []);

numel(find(~blankIdx(:)))

images = allImgRshp(:,:,:,~blankIdx(:));
images = images(:,:,:,randperm(size(images,4)));

figure(1); clf; set(gcf, 'color', 'w')
for jj = 1:100
    subplot(10,10,jj);
    imshow(squeeze(images(:,:,:,jj)), []); drawnow;
end

% Visualize images
% figure(1); clf; set(gcf, 'color', 'w')
% for jj = 1:size(images,4)
%     clf;
% 
%     imagesc(squeeze(images(:,:,:,jj))); title(jj); pause(0.001); drawnow;
% end

save(fullfile(simseqRootPath, 'stimulus', 'stim_mri', 'cartoon', sprintf('cartoonChopped_%ddeg_FOV38by29cm_rms%1.1f.mat', sizeStimBlockWDeg, rmsThresh)), 'images', '-v7.3')

return


    

