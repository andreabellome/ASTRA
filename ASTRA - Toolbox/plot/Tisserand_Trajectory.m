function [r_p_pl , r_a_pl ] = Tisserand_Trajectory ( r_p_1 , r_a_1 , tj , Fb_planet , hFB , Plot_Points )
% ------------------------------------------------------------------------------------------------
%
% This functions calculated the positions for a maximum FB (considered in
% Tisserand analysis of CASTAway) dividing the fly by angle in small steps
% to allow plotting of the results in Tisserand Map ( see
% Main_CASTAway_2FB.m)

% INPUTS:
%
%          r_p_1 : initial orbit radius of perigee (AU)
%     
%          r_a_1 : initial orbit radius of apogee (AU)
%     
%          tj : type of fly by
%     
%                 0 = in front the planet fly by to decrease orbital energy
% 
%                 1 = behind the planet fly by to increase orbital energy
%     
%          Fb_planet : fly by planet
%
%          Plot_Points : number of points to divide the FB in
%
% OUTPUTS :
%
%          r_p_pl : vector containing the progress of radius of perigee over the FB (AU)
%     
%          r_a_pl : vector containing the progress of radius of apogee over the FB (AU)
%
% Author: Ramiro Gallego Fernández
%
% Funcions used: ATD_Planets_Constants.m 
%
% ------------------------------------------------------------------------------------------------

MuSun = 132724487690;

    [a_planet, R_Planet, MuPlanet] = astroConstantsj2000(Fb_planet);


    a_pl             = (r_p_1 + r_a_1)/2;
    ecc_pl           = 1 - r_p_1/a_pl;
    p_pl             = a_pl*(1-ecc_pl^2);
    V_arrival_pl     = sqrt(MuSun*(2/a_planet-1/a_pl));
    Fl_angle_pl      = acos(sqrt(p_pl*MuSun)/(a_planet*V_arrival_pl));
    V_planet         = sqrt(MuSun/a_planet);
    V_inf_pl         = sqrt(V_arrival_pl^2+V_planet^2-2*V_planet*V_arrival_pl*cos(Fl_angle_pl));
    Beta_pl          = acos((V_arrival_pl^2 - V_inf_pl^2 - V_planet^2)/(-2*V_inf_pl*V_planet))*180/pi;
    Alpha_initial_pl = 180 - Beta_pl;

    Perigee_FlyBy_Altitude   = hFB;
    DeltaMax_FB              = 2*asin(MuPlanet/(MuPlanet + V_inf_pl^2*(R_Planet + Perigee_FlyBy_Altitude)))*180/pi;
   
    Delta_plot        = linspace(0,DeltaMax_FB,Plot_Points);
    
    r_p_pl = zeros(1,Plot_Points);  r_a_pl = zeros(1,Plot_Points);

 for q = 1:Plot_Points
   
       if tj == 0 % in front fly by --> energy is decreased
        
        Alpha_end_pl = Alpha_initial_pl + Delta_plot(1,q);
        
        if Alpha_end_pl > 180   
        Alpha_end_pl = 180;
        end
        
           elseif tj == 1 % behind fly by --> energy is increased
               
         Alpha_end_pl = Alpha_initial_pl - Delta_plot(1,q);
         
        if Alpha_end_pl < 0   
        Alpha_end_pl = 0;
        end
        
           end
           
   v_Heliocentric_Out_pl    = sqrt(V_planet^2+V_inf_pl^2-2*V_inf_pl*V_planet*cos((180-Alpha_end_pl)*pi/180)); % New heliocentric velocity
   a_new_pl                 = 1/(2/a_planet-v_Heliocentric_Out_pl^2/MuSun); % New semi major axis
    
    if a_new_pl < 0
        r_a_pl(1,q) = NaN; r_p_pl(1,q) = NaN; %Eliminates hyperbolic transfers
    else
    
        if Alpha_end_pl == 180
            New_Flight_path_pl = 0 ; % Corrects rounding errors
        else
            New_Flight_path_pl = acos((v_Heliocentric_Out_pl^2+V_planet^2-V_inf_pl^2)/(2*V_planet*v_Heliocentric_Out_pl));
            
        end
    p_new_pl              = (a_planet*v_Heliocentric_Out_pl*cos(New_Flight_path_pl))^2/MuSun;
    ecc_new_pl            = sqrt(1-p_new_pl/a_new_pl);
            
    r_a_pl(1,q)           = a_new_pl*(1+ecc_new_pl);
    r_p_pl(1,q)           = (2*a_new_pl-r_a_pl(1,q));
    
    end
 end
   
end
