function [pos, vel, el, az, dop] = leastSquarePos_zcj(satpos, satvel, obs, obs_dot, settings)
%Function calculates the Least Square Solution.
%
%[pos, el, az, dop] = leastSquarePos(satpos, obs, settings);
%
%   Inputs:
%       satpos      - Satellites positions (in ECEF system: [X; Y; Z;] -
%                   one column per satellite)
%       obs         - Observations - the pseudorange measurements to each
%                   satellite:
%                   (e.g. [20000000 21000000 .... .... .... .... ....])
%       settings    - receiver settings
%
%   Outputs:
%       pos         - receiver position and receiver clock error 
%                   (in ECEF system: [X, Y, Z, dt]) 
%       el          - Satellites elevation angles (degrees)
%       az          - Satellites azimuth angles (degrees)
%       dop         - Dilutions Of Precision ([GDOP PDOP HDOP VDOP TDOP])

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%--------------------------------------------------------------------------
%Based on Kai Borre
%Copyright (c) by Kai Borre
%Updated by Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
%
% CVS record:
% $Id: leastSquarePos.m,v 1.1.2.12 2006/08/22 13:45:59 dpl Exp $
%==========================================================================

%=== Initialization =======================================================
nmbOfIterations = 7;

dtr     = pi/180;
pos     = zeros(4, 1);
X       = satpos;
Rot_satpos = satpos;
nmbOfSatellites = size(satpos, 2);

A       = zeros(nmbOfSatellites, 4);
omc     = zeros(nmbOfSatellites, 1);     % XieGang 5.36 G*[delta_pos]=omc
az      = zeros(1, nmbOfSatellites);
el      = az;

%% === Iteratively find receiver position ===================================
for iter = 1:nmbOfIterations

    for i = 1:nmbOfSatellites
        if iter == 1
            %--- Initialize variables at the first iteration --------------
            Rot_X = X(:, i);
            trop = 2;
        else
            %--- Update equations -----------------------------------------
            rho2 = (X(1, i) - pos(1))^2 + (X(2, i) - pos(2))^2 + ...
                   (X(3, i) - pos(3))^2;
            traveltime = sqrt(rho2) / settings.c ;

            %--- Correct satellite position (do to earth rotation) --------
            Rot_X = e_r_corr(traveltime, X(:, i));  % Rot_X: 3*1

            %--- Find the elevation angel of the satellite ----------------
            [az(i), el(i), dist] = topocent(pos(1:3, :), Rot_X - pos(1:3, :));

            if (settings.useTropCorr == 1)
                %--- Calculate tropospheric correction --------------------
                trop = tropo(sin(el(i) * dtr), ...
                             0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
            else
                % Do not calculate or apply the tropospheric corrections
                trop = 0;
            end
        end % if iter == 1 ... ... else 

        %--- Apply the corrections ----------------------------------------
        omc(i) = (obs(i) - norm(Rot_X - pos(1:3), 'fro') - pos(4) - trop);

        %--- Construct the A matrix ---------------------------------------
        A(i, :) =  [ (-(Rot_X(1) - pos(1))) / obs(i) ...
                     (-(Rot_X(2) - pos(2))) / obs(i) ...
                     (-(Rot_X(3) - pos(3))) / obs(i) ...
                     1 ];
    
        if iter == nmbOfIterations    % record satellite true position at the last epoch
            Rot_satpos(:, i) = Rot_X; 
        end         
                 
    end % for i = 1:nmbOfSatellites
       
    % These lines allow the code to exit gracefully in case of any errors
    if rank(A) ~= 4
        pos     = zeros(1, 4);
        vel     = zeros(1, 4);
        return
    end

    %--- Find position update ---------------------------------------------
    x   = A \ omc;                          % 此处没使用最小二乘，不知道为什么；不过改成最小二乘后，似乎也没啥区别
    % x   = (A' * A)^(-1) * A' * omc;       % 可能因为我这组数据是低动态无噪声吧
    
    %--- Apply position update --------------------------------------------
    pos = pos + x;    
    
end % for iter = 1:nmbOfIterations

pos = pos';

%% === Iteratively find receiver velocity ===================================
omc_dot = zeros(nmbOfSatellites, 1);     % XieGang 5.77 G*[vel]=omc_dot
II = zeros(nmbOfSatellites, 4);          % 速度最小二乘法雅可比矩阵重新计算一次，直接用A误差很大

for ii = 1 : nmbOfSatellites
    distance = sqrt((Rot_satpos(1, ii) - pos(1))^2 + ...
                    (Rot_satpos(2, ii) - pos(2))^2 + ...
                    (Rot_satpos(3, ii) - pos(3))^2);
    II(ii, 1:3) = -[Rot_satpos(1, ii) - pos(1), Rot_satpos(2, ii) - pos(2), Rot_satpos(3, ii) - pos(3)] ./ distance; 
    II(ii, 4)   = 1;
    omc_dot(ii, 1) = obs_dot(ii) + II(ii, 1:3) * e_r_corr(distance/settings.c, satvel(:, ii));   % 扣除了卫星运动引起的doppler
end
vel = II \ omc_dot;    vel = vel';


%% === Calculate Dilution Of Precision ======================================
if nargout  >= 4
    %--- Initialize output ------------------------------------------------
    dop     = zeros(1, 5);
    
    %--- Calculate DOP ----------------------------------------------------
    Q       = inv(A'*A);
    
    dop(1)  = sqrt(trace(Q));                       % GDOP    
    dop(2)  = sqrt(Q(1,1) + Q(2,2) + Q(3,3));       % PDOP
    dop(3)  = sqrt(Q(1,1) + Q(2,2));                % HDOP
    dop(4)  = sqrt(Q(3,3));                         % VDOP
    dop(5)  = sqrt(Q(4,4));                         % TDOP
    
    dop     = dop';
end

end