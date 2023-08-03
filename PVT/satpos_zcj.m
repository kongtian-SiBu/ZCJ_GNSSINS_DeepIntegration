function [satPositions, satVelocities, satClkCorr, satClkDrift] = satpos_zcj(transmitTimeatSat, ...
                                        prnList, eph)
                                    
numOfSatellites = size(prnList, 2);
gpsPi          = 3.1415926535898;   
Omegae_dot     = 7.2921151467e-5;  % Earth rotation rate, [rad/s]
GM             = 3.986005e14;      % Universal gravitational constant times
                                   % the mass of the Earth, [m^3/s^2]
F              = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]

satClkCorr     = zeros(1, numOfSatellites);    % [s]
satPositions   = zeros(3, numOfSatellites);    
satVelocities  = zeros(3, numOfSatellites);    
satClkDrift    = zeros(1, numOfSatellites);    % [s/s]

for satNr = 1 : numOfSatellites    
    prn = prnList(satNr);
    
    dt = check_t(transmitTimeatSat(satNr) - eph(prn).t_oc);   % revise here1

    %--- Calculate clock correction ---------------------------------------
    satClkCorr(satNr) = (eph(prn).a_f2 * dt + eph(prn).a_f1) * dt + ...
                         eph(prn).a_f0 - ...
                         eph(prn).T_GD;
                     
    time = transmitTimeatSat(satNr) - satClkCorr(satNr);      % revise here2
    
    %% Find satellite's position ----------------------------------------------

    %Restore semi-major axis
    a   = eph(prn).sqrtA * eph(prn).sqrtA;

    %Time correction
    tk  = check_t(time - eph(prn).t_oe);

    %Initial mean motion
    n0  = sqrt(GM / a^3);
    %Mean motion
    n   = n0 + eph(prn).deltan;                           % XieGang 3.52

    %Mean anomaly
    M   = eph(prn).M_0 + n * tk;                          % XieGang 3.53
    %Reduce mean anomaly to between 0 and 360 deg
    M   = rem(M + 2*gpsPi, 2*gpsPi);
    
    M_dot = n;                                            % XieGang 3.76
    
    %Initial guess of eccentric anomaly
    E   = M;

    %--- Iteratively compute eccentric anomaly ----------------------------
    for ii = 1:10                                         % XieGang 3.39
        E_old   = E;
        E       = M + eph(prn).e * sin(E);
        dE      = rem(E - E_old, 2*gpsPi);

        if abs(dE) < 1.e-12
            % Necessary precision is reached, exit from the loop
            break;
        end
    end

    %Reduce eccentric anomaly to between 0 and 360 deg
    E   = rem(E + 2*gpsPi, 2*gpsPi);

    E_dot = M_dot / (1 - eph(prn).e * cos(E));            % XieGang 3.75
    
    %Compute relativistic correction term
    dtr = F * eph(prn).e * eph(prn).sqrtA * sin(E);       % XieGang 4.23

    %Calculate the true anomaly
    nu   = atan2(sqrt(1 - eph(prn).e^2) * sin(E), cos(E)-eph(prn).e);  % XieGang 3.46
    
    nu_dot = sqrt(1 - eph(prn).e ^ 2) * E_dot / (1 - eph(prn).e * cos(E));  % XieGang 3.74
    
    %Compute angle phi
    phi = nu + eph(prn).omega;                            % XieGang 3.54
    %Reduce phi to between 0 and 360 deg
    phi = rem(phi, 2*gpsPi);   
    
    phi_dot = nu_dot;                                     % XieGang 3.73
    
    %Correct argument of latitude                         % XieGang 3.58 3.59 3.60 
    u = phi + ...
        eph(prn).C_uc * cos(2*phi) + ...
        eph(prn).C_us * sin(2*phi);
    %Correct radius
    r = a * (1 - eph(prn).e*cos(E)) + ...
        eph(prn).C_rc * cos(2*phi) + ...
        eph(prn).C_rs * sin(2*phi);
    %Correct inclination
    i = eph(prn).i_0 + eph(prn).iDot * tk + ...
        eph(prn).C_ic * cos(2*phi) + ...
        eph(prn).C_is * sin(2*phi);
    
    % XieGang 3.66 3.67 3.68 && 3.70 3.71 3.72
    u_dot = phi_dot + 2 * phi_dot * ( eph(prn).C_us * cos(2*phi) - ...
            eph(prn).C_uc * sin(2*phi));
    r_dot = a * eph(prn).e * E_dot * sin(E) + 2 * phi_dot * ( eph(prn).C_rs * cos(2*phi) - ...
            eph(prn).C_rc * sin(2*phi));    
    i_dot = eph(prn).iDot + 2 * phi_dot * ( eph(prn).C_is * cos(2*phi) - ...
            eph(prn).C_ic * sin(2*phi));
    
    %Compute the angle between the ascending node and the Greenwich meridian
    Omega = eph(prn).omega_0 + (eph(prn).omegaDot - Omegae_dot)*tk - ...
            Omegae_dot * eph(prn).t_oe;              % XieGang 3.62
    %Reduce to between 0 and 360 deg 
    Omega = rem(Omega + 2*gpsPi, 2*gpsPi);
    
    Omega_dot = eph(prn).omegaDot - Omegae_dot;      % XieGang 3.69

    %--- Compute satellite coordinates ------------------------------------
    satPositions(1, satNr) = cos(u)*r * cos(Omega) - sin(u)*r * cos(i)*sin(Omega);
    satPositions(2, satNr) = cos(u)*r * sin(Omega) + sin(u)*r * cos(i)*cos(Omega);
    satPositions(3, satNr) = sin(u)*r * sin(i);

    %--- Compute satellite velocities ------------------------------------
    x1dot = r_dot * cos(u) - r * u_dot * sin(u);  % 谢刚《GPS原理与接收机设计》P64公式3.65A
    y1dot = r_dot * sin(u) + r * u_dot * cos(u);  % 谢刚《GPS原理与接收机设计》P64公式3.65B
    
    % 谢刚《GPS原理与接收机设计》P63公式3.64ABC
    satVelocities(1, satNr) = -satPositions(2,satNr)*Omega_dot - (y1dot*cos(i) - ...
                              satPositions(3,satNr)*i_dot)*sin(Omega) + x1dot*cos(Omega);
    satVelocities(2, satNr) = satPositions(1,satNr)*Omega_dot + (y1dot*cos(i) - ...
                              satPositions(3,satNr)*i_dot)*cos(Omega) + x1dot*sin(Omega);
    satVelocities(3, satNr) = y1dot * sin(i) + sin(u)*r * i_dot * cos(i);
    
%% Include relativistic correction in clock correction --------------------
    satClkCorr(satNr) = (eph(prn).a_f2 * dt + eph(prn).a_f1) * dt + ...
                         eph(prn).a_f0 - ...
                         eph(prn).T_GD + dtr;             % XieGang 4.22

    dtr_dot = F * eph(prn).e * eph(prn).sqrtA * E_dot * cos(E);    % XieGang 4.27            
    satClkDrift(satNr)= eph(prn).a_f1 + 2 * eph(prn).a_f2 * dt + dtr_dot;   % XieGang 4.26
                     
end

end

