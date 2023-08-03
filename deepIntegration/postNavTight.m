function navSolut = postNavTight(trackDeepIn, settings, eph, TOW)
%% GNSS紧组合数据准备（不对仰角进行限制）
%
% 输入参数:
%        - trackDeepIn: 跟踪结构体，通道之间进行了准同步，即距离给定时刻不足一个相干积分周期
%        - settings: 接收机相关参数
%        - eph: 星历
%        - TOW: 首个子帧帧头的发射时刻
%
% 输出参数:
%        - navSolut: GNSS观测值结构体，伪距和伪距率。不定位
%
%--------------------------------------------------------------------------
%% 
numActChnList = length(trackDeepIn);

if (isempty(trackDeepIn) || (size(trackDeepIn, 2) < 1))
    disp('Too few satellites with ephemeris data for postion calculations. Exiting!');
    navSolut = [];
    return
end

%% 计算信号发射时刻
launchTime = zeros(1, numActChnList);
rawP = zeros(1, numActChnList);
doppler = zeros(1, numActChnList);
remSampleNum = zeros(1, numActChnList);
totalSampleNum = zeros(1, numActChnList);

for ii = 1 : numActChnList
    remSampleNum(ii) = fix((settings.recvTime - trackDeepIn(ii).recvTime) * settings.samplingFreq);
    totalSampleNum(ii) = trackDeepIn(ii).SamplePos + remSampleNum(ii);  
end

maxTotal = max(totalSampleNum);
for ii = 1 : numActChnList
    remSampleNum(ii) = remSampleNum(ii) + (maxTotal - totalSampleNum(ii));
end

for ii = 1 : numActChnList   
    codePhaseStep = trackDeepIn(ii).codeFreq / settings.samplingFreq;
    num_Cyclic = trackDeepIn(ii).numOfCoInt;

    codePhaseTao = trackDeepIn(ii).remCodePhase + codePhaseStep * ( remSampleNum(ii) - 1 );  % -1?

    launchTime(ii) = TOW + num_Cyclic * 1e-3 + codePhaseTao / 1023 * 1e-3;   % [s]

    rawP(ii) = (settings.recvTime - launchTime(ii)) * settings.c; 

    doppler(ii) = trackDeepIn(ii).carrFreq - settings.IF;
end
rawP_dot = -doppler .* settings.c / 1575.42e6;

%% 
[satPositions, satVelocity, satClkCorr, satClkDrift] = satpos_zcj(launchTime, ...
                                [trackDeepIn.PRN], eph);

%% 
navSolut.rawP = rawP;  navSolut.rawP_dot = rawP_dot;  
navSolut.satClkCorr = satClkCorr;  navSolut.satClkDrift = satClkDrift;  navSolut.satPositions = satPositions;
navSolut.satVelocity = satVelocity;

end
