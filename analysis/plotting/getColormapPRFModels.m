function cmap = getColormapPRFModels(type)

switch type
    case 0 % Dark gray LSS
        gray     = [0.5 0.5 0.5];     % LSS
        orange   = [255 140 0]./255;  % CSS
        blue     = [26 115 225]./255; % CST
        cmap     = [gray;orange;blue];
    
    case 1 % Light gray LSS
        gray     = [0.7 0.7 0.7];     % LSS
        orange   = [255 140 0]./255;  % CSS
        blue     = [26 115 225]./255; % CST
        cmap     = [gray;orange;blue];
        
    case 2 % DoG
        purple = [148 0 211]./255;
        cmap = purple;
        
    case 3 % Supplementary spatiotemporal pRF models
        darkblue = [26 115 225]./255; % CST_fix
        cyan     = [0 174 239]./255; % CST_ST
        pink     = [236 0 150]./255; % DN_ST
        cmap     = [darkblue;cyan;pink];
end