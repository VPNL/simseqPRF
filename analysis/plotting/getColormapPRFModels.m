function cmap = getColormapPRFModels(type)

switch type
    case 0 % Dark gray LSS
        darkgray = [0.5 0.5 0.5];     % LSS
        orange   = [255 140 0]./255;  % CSS
        blue     = [26 115 225]./255; % CST
        cmap     = [darkgray;orange;blue];
    case 1 % Light gray LSS
        darkgray = [0.7 0.7 0.7];     % LSS
        orange   = [255 140 0]./255;  % CSS
        blue     = [26 115 225]./255; % CST
        cmap     = [darkgray;orange;blue];
end