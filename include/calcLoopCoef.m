function [tau1, tau2] = calcLoopCoef(LBW, zeta, k)
%Function finds loop coefficients. The coefficients are used then in PLL-s
%and DLL-s.
%
%[tau1, tau2] = calcLoopCoef(LBW, zeta, k)
%
%   Inputs:
%       LBW           - Loop noise bandwidth
%       zeta          - Damping ratio
%       k             - Loop gain
%
%   Outputs:
%       tau1, tau2   - Loop filter coefficients 
%--------------------------------------------------------------------------

% 与谢钢书上二阶环路滤波公式严格对应

% Solve natural frequency
Wn = LBW*8*zeta / (4*zeta.^2 + 1);

% solve for t1 & t2
tau1 = k / (Wn * Wn);
tau2 = 2.0 * zeta / Wn;

end