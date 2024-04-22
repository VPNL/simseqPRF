function ax = plotpRFlocOnImage(ax,prfData, stim, mask, hemi)

if strcmp(hemi, 'lh')
    xylim = [-1,12];
    xyticks = [0 5 10];
elseif strcmp(hemi, 'rh')
    xylim = [-12,1];
    xyticks = [-10 -5 0];
elseif strcmp(hemi, 'both')
    xylim = [-12,12];
    xyticks = [-10 -5 0 5 10];
end

% Plot stimulus squares
X = linspace(-12,12,size(stim,1)); 
Y = linspace(12,-12,size(stim,2)); 
% imagesc(X,Y,double(~stim));
imagesc(X,Y,stim);
set(gca, 'yDir','normal', 'CLim', [0 max(stim(:))])
hold all; colormap gray;

hold all;
gxPos = get(ax,'Position'); set(ax,'Units','normalized')
set(ax,'Position',[gxPos(1) gxPos(2) gxPos(3) gxPos(4)])
plot(xylim, [0,0],'k', [0,0],xylim, 'k'); hold all;

axis equal;  xlim(xylim); ylim(xylim);

set(ax, 'FontSize', 15, 'TickDir', 'out', ...
    'XTick',xyticks, 'XTickLabel',xyticks, ...
    'YTick',xyticks, 'YTickLabel',xyticks);
xlabel({'x position (deg)',''});
ylabel('y position (deg)');
% title(sprintf('%s: All pRF locations', roiName))


ax.Units = 'Points';
axSizePix = ax.Position(3:4);
ax.Units = 'normalized';
pixPerUnit = ceil(max(axSizePix ./ [range(xlim(ax)),range(ylim(ax))]));
radiusToPlot = 2.5; % radius (in deg)
if strcmp(hemi, 'both')
    if isfield(prfData, 'lh') || isfield(prfData, 'rh')
        voxelsToPlot = 1:(size(prfData.lh.x0,2)+size(prfData.rh.x0,2));
        x0 = cat(2,prfData.lh.x0,prfData.rh.x0);
        y0 = cat(2,prfData.lh.y0,prfData.rh.y0);
        sz = cat(2,prfData.lh.effectiveSize,prfData.rh.effectiveSize);
        if isfield(prfData.lh,'sigmaSurround') || isfield(prfData.rh,'sigmaSurround')
            sz_c = cat(2,prfData.lh.sigmaMajor,prfData.rh.sigmaMajor);
            sz_s = cat(2,prfData.lh.sigmaSurround,prfData.rh.sigmaSurround);
        end
    else
        x0 = prfData.x0;
        y0 = prfData.y0;
        sz = prfData.effectiveSize;
        voxelsToPlot = 1:length(sz);
        if isfield(prfData,'sigmaSurround')
            sz_c = prfData.sigmaMajor;
            sz_s = prfData.sigmaSurround;
        end
    end
elseif strcmp(hemi, 'lh')
    voxelsToPlot = 1:size(prfData.lh.x0,2);
    x0 = prfData.lh.x0;
    y0 = prfData.lh.y0;
    sz = prfData.lh.effectiveSize;
    if isfield(prfData.lh,'sigmaSurround')
        sz_c = prfData.lh.sigmaMajor;
        sz_s = prfData.lh.sigmaSurround;
    end
elseif strcmp(hemi,'rh')
    voxelsToPlot = 1:size(prfData.rh.x0,2);
    x0 = prfData.rh.x0;
    y0 = prfData.rh.y0;
    sz = prfData.rh.effectiveSize;
    if isfield(prfData.rh,'sigmaSurround')
        sz_c = prfData.rh.sigmaMajor;
        sz_s = prfData.rh.sigmaSurround;
    end
end

if ~isempty(mask)
    voxelsToPlot = voxelsToPlot(mask);
else
    voxelsToPlot= voxelsToPlot;
end

if length(voxelsToPlot)==1
    alpha = 0.5;
else
    alpha = 0.05;
end
for vv = voxelsToPlot
    if sz(vv) > 10.5
        xlim([-1,20]); ylim([-1,20]);
    end
    if ~isnan(x0(vv))
    % Plot pRF (scatter size = area in pixels, given that area = pi*radius^2
%     scatter(x0(vv),y0(vv),pi*(radiusToPlot*sz(vv)*pixPerUnit).^2, ...
%         'MarkerEdgeColor',[0 0 0],...
%         'MarkerFaceColor',[1 1 .5],...
%         'MarkerFaceAlpha',alpha, ...
%         'LineWidth',0.5)
        if exist('sz_s','var')
            rectangle('Position',[x0(vv)-((radiusToPlot*sz_s(vv))/2),...
                y0(vv)-((radiusToPlot*sz_s(vv))/2),...
                radiusToPlot*sz_s(vv),radiusToPlot*sz_s(vv)], 'Curvature',[1 1], ...
                'EdgeColor',[0 0 0], 'FaceColor',[1 0 0 alpha]); hold on;
            rectangle('Position',[x0(vv)-((radiusToPlot*sz_c(vv))/2),...
                y0(vv)-((radiusToPlot*sz_c(vv))/2),...
                radiusToPlot*sz_c(vv),radiusToPlot*sz_c(vv)], 'Curvature',[1 1], ...
                'EdgeColor',[0 0 0], 'FaceColor',[1 1 .5 alpha]);
        else
            rectangle('Position',[x0(vv)-((radiusToPlot*sz(vv))/2),...
                y0(vv)-((radiusToPlot*sz(vv))/2),...
                radiusToPlot*sz(vv),radiusToPlot*sz(vv)], 'Curvature',[1 1], ...
                'EdgeColor',[0 0 0], 'FaceColor',[1 1 .5 alpha]);
        end
    end
    xlim(xylim); ylim(xylim);
end

ax = gca;
end