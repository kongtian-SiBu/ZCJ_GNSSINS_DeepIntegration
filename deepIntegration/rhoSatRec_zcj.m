function [rho, LOS, AzEl, v_r_s, BLH, Cen] = rhoSatRec_zcj(satPos, recPos, rho0, satVel, recVeln)
% v_r_s  : velocity between Satellite and Receiver (n*1)
% rho    : distance between Satellite and Receiver
% recPos : receiver position in ECEF
% recVeln: receiver velocity in n frame

% satPos: n*3
% recPos: 3*1
% rho0  : n*1
% satVel: n*3
% recVeln: 3*1

global ggps
    recPos = repmat(recPos(1:3)',size(satPos,1),1);
	if nargin<3, dpos=satPos-recPos; rho0=sqrt(dpos(:,1).^2+dpos(:,2).^2+dpos(:,3).^2); end
    wtau = ggps.wie*rho0/ggps.c;
    sw = sin(wtau);  cw = cos(wtau);
    satRot = [cw.*satPos(:,1)+sw.*satPos(:,2), -sw.*satPos(:,1)+cw.*satPos(:,2), satPos(:,3)];
    dpos = satRot - recPos;
    rho = sqrt(dpos(:,1).^2+dpos(:,2).^2+dpos(:,3).^2);
    LOS = dpos./[rho,rho,rho];
    
    % topocent
    [BLH, Cen] = xyz2blh(recPos(1,1:3));
    Ln = LOS*Cen;  % Ln: LOS in receiver n-frame
    AzEl = [atan2(Ln(:,1), Ln(:,2)), atan(Ln(:,3)./sqrt(Ln(:,1).^2+Ln(:,2).^2))];
    idx = find(AzEl(:,1)<0);
    AzEl(idx,1) = AzEl(idx,1) + 2*pi;
    
    % vel ZCJ Added
    recVele   = Cen * recVeln;   % recVel in ECEF
    recVele   = repmat(recVele(1:3)',size(satPos,1),1); 
    satVelRot = [cw.*satVel(:,1)+sw.*satVel(:,2), -sw.*satVel(:,1)+cw.*satVel(:,2), satVel(:,3)];
    vrs       = satVelRot - recVele;
    v_r_s     = zeros(size(satPos,1), 1);
    for ii = 1 : size(satPos,1)
        v_r_s(ii,1) = vrs(ii,:) * LOS(ii,:)'; 
    end
end


