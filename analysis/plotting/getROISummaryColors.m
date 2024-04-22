function [colors1, colors2] = getROISummaryColors(type)

switch type
    case 0
        colors1 = [...
            38,60,142; ... % V1: Dark blue
            73,115,184; ... % V2: Medium dark blue
            142, 189, 228; ... % V3: Light blue
            15,104,53; ... % hV4: Dark green
            132,197, 81; ... % VO1/2: Light green
            120, 91, 167; ... % V3AB: Dark purlple
            207, 162, 203; % IPS0/1: light purlple
            213, 77, 40; ... % LO1/2: Red/orange
            213, 183, 24]./256; ... % TO1/2: dark yellow
        colors2 = [];
    case 1
        satValues = 1-linspace(0,0.9,4);
        turbomap = cmapturbo(5);
        cmap = varysat(turbomap,satValues);

        colors1 = NaN(20,1,3);
        colors2 = [];
        
        colors1(1,:,:)    = cmap(1,1,:); % V1
        colors1(2,:,:)    = cmap(1,2,:); % V2
        colors1(3,:,:)    = cmap(1,3,:); % V3

        colors1(4,:,:)    = cmap(2,1,:); % Ventral areas (hV4)
        colors1(5,:,:)    = cmap(2,2,:); % Ventral areas (VO1)
        colors1(6,:,:)    = cmap(2,3,:); % Ventral areas (VO2)

        colors1(7,:,:)    = cmap(3,1,:); % Lateral areas (V3AB)
  
        colors1(8,:,:)   = cmap(4,1,:); % Dorsal areas (LO1, LO2)
        colors1(9,:,:)   = cmap(4,2,:); % Dorsal areas (LO1, LO2)

        colors1(10,:,:) = cmap(5,1,:); % Dorsal areas (TO1, TO2)
        colors1(11,:,:) = cmap(5,2,:); % Dorsal areas (TO1, TO2)
        
        colors1(12,:,:)   = cmap(3,2,:); % Lateral areas (IPS0)
        colors1(13,:,:)   = cmap(3,3,:); % Lateral areas (IPS1)
        colors1(14,:,:)   = cmap(3,4,:); % Lateral areas (IPS2)
        colors1(15,:,:)   = cmap(3,4,:); % Lateral areas (IPS3)
        colors1(16,:,:)   = cmap(3,4,:); % Lateral areas (IPS4)
        colors1(17,:,:)   = cmap(3,4,:); % Lateral areas (IPS5)

        colors1(18,:,:) = cmap(2,1,:); % VO
        colors1(19,:,:) = cmap(4,1,:); % LO
        colors1(20,:,:) = cmap(5,1,:); % TO
    case 2
        satValues = 1-linspace(0,0.9,3);
        turbomap = cmapturbo(5);
        cmap = varysat(turbomap,satValues);

        colors1 = NaN(20,1,3); % full color (unsaturated)
        colors2 = NaN(20,1,3); % saturated color
        
        colors1(1,:,:)    = cmap(1,1,:); % V1
        colors1(2,:,:)    = cmap(1,1,:); % V2
        colors1(3,:,:)    = cmap(1,1,:); % V3
        colors1(4,:,:)    = cmap(2,1,:); % Ventral areas (hV4)
        colors1(5,:,:)    = cmap(2,1,:); % Ventral areas (VO1)
        colors1(6,:,:)    = cmap(2,1,:); % Ventral areas (VO2)
        colors1(7,:,:)    = cmap(3,1,:); % Lateral areas (V3AB)
        
        colors1(8,:,:)   = cmap(4,1,:); % Dorsal areas (LO1, LO2)
        colors1(9,:,:)   = cmap(4,1,:); % Dorsal areas (LO1, LO2)
        colors1(10,:,:)  = cmap(5,1,:); % Dorsal areas (TO1, TO2)
        colors1(11,:,:)  = cmap(5,1,:); % Dorsal areas (TO1, TO2)
        
        colors1(12,:,:)   = cmap(3,1,:); % Lateral areas (IPS0)
        colors1(13,:,:)   = cmap(3,1,:); % Lateral areas (IPS1)
        colors1(14,:,:)   = cmap(3,1,:); % Lateral areas (IPS2)
        colors1(15,:,:)   = cmap(3,1,:); % Lateral areas (IPS3)
        colors1(16,:,:)   = cmap(3,1,:); % Lateral areas (IPS4)
        colors1(17,:,:)   = cmap(3,1,:); % Lateral areas (IPS5)

        colors1(18,:,:) = cmap(2,1,:); % VO
        colors1(19,:,:) = cmap(4,1,:); % LO
        colors1(20,:,:) = cmap(5,1,:); % TO
        
        colors2(1,:,:)    = cmap(1,2,:); % V1
        colors2(2,:,:)    = cmap(1,2,:); % V2
        colors2(3,:,:)    = cmap(1,2,:); % V3
        colors2(4,:,:)    = cmap(2,2,:); % Ventral areas (hV4)
        colors2(5,:,:)    = cmap(2,2,:); % Ventral areas (VO1)
        colors2(6,:,:)    = cmap(2,2,:); % Ventral areas (VO2)
        colors2(7,:,:)    = cmap(3,2,:); % Lateral areas (V3AB)
           
        colors2(8,:,:)    = cmap(4,2,:); % Dorsal areas (LO1)
        colors2(9,:,:)    = cmap(4,2,:); % Dorsal areas (LO2)
        colors2(10,:,:)   = cmap(5,2,:); % Dorsal areas (TO1)
        colors2(11,:,:)   = cmap(5,2,:); % Dorsal areas (TO2)
        
        colors2(12,:,:)   = cmap(3,2,:); % Lateral areas (IPS0)
        colors2(13,:,:)   = cmap(3,2,:); % Lateral areas (IPS1)
        colors2(14,:,:)   = cmap(3,2,:); % Lateral areas (IPS2) 
        colors2(15,:,:)   = cmap(3,2,:); % Lateral areas (IPS3) 
        colors2(16,:,:)   = cmap(3,2,:); % Lateral areas (IPS4)
        colors2(17,:,:)   = cmap(3,2,:); % Lateral areas (IPS5) 

        colors2(18,:,:)   = cmap(2,2,:); % VO
        colors2(19,:,:)   = cmap(4,2,:); % LO
        colors2(20,:,:)   = cmap(5,2,:); % TO
        
end
        
return