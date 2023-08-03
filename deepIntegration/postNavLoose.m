function [navSolut, satPositions, satVelocity, satClkCorr, rawP, rawP_dot] = postNavLoose(trackDeepIn, settings, eph, TOW)
%% 最小二乘定位解算，为松组合做数据准备（不对仰角进行限制）?
%
% 输入参数:
%        - trackDeepIn: 跟踪结构体，通道之间进行了准同步，即距离给定时刻不足一个相干积分周期?
%        - settings: 接收机相关参数?
%        - eph: 星历
%        - TOW: 首个子帧帧头的发射时刻?
%
% 杈撳嚭鍙傛暟:
%        - navSolut: GNSS观测值结构体
%
%--------------------------------------------------------------------------
%% 成功跟踪的通道数不足4个，则失败?
numActChnList = length(trackDeepIn);

if (isempty(trackDeepIn) || (size(trackDeepIn, 2) < 4))
    % Show error message and exit
    disp('Too few satellites with ephemeris data for postion calculations. Exiting!');
    navSolut = [];
    return
end

%% 计算给定时刻伪距和多普勒
launchTime = zeros(1, numActChnList);
rawP = zeros(1, numActChnList);
doppler = zeros(1, numActChnList);
remSampleNum = zeros(1, numActChnList);
totalSampleNum = zeros(1, numActChnList);

for ii = 1 : numActChnList
    % 该通道距离给定时刻还差多少个采样点
    remSampleNum(ii) = fix((settings.recvTime - trackDeepIn(ii).recvTime) * settings.samplingFreq);
    totalSampleNum(ii) = trackDeepIn(ii).SamplePos + remSampleNum(ii);  
end
% 利用采样点计算时间，此处可能存在舍入误差，误差一个点将造成3e8/15e6=20m的巨大误差
maxTotal = max(totalSampleNum);  % 向最多的点数看齐
for ii = 1 : numActChnList
    remSampleNum(ii) = remSampleNum(ii) + (maxTotal - totalSampleNum(ii));
end

for ii = 1 : numActChnList   
    % 该通道对应的码相位步长
    codePhaseStep = trackDeepIn(ii).codeFreq / settings.samplingFreq;
    
    % 该位置距离首次出现帧头的相干积分次数
    num_Cyclic = trackDeepIn(ii).numOfCoInt;
        
    % 不足一整周的伪码码相位
    codePhaseTao = trackDeepIn(ii).remCodePhase + codePhaseStep * ( remSampleNum(ii) - 1 );  % -1?
    
    % 计算给定时刻的卫星信号发射时间
    launchTime(ii) = TOW + num_Cyclic * 1e-3 + codePhaseTao / 1023 * 1e-3;   % [s]

    % 计算伪距
    rawP(ii) = (settings.recvTime - launchTime(ii)) * settings.c; 

    % 计算多普勒 [Hz]
    doppler(ii) = trackDeepIn(ii).carrFreq - settings.IF;
end
rawP_dot = -doppler .* settings.c / 1575.42e6;

%% 计算卫星位置等信息
[satPositions, satVelocity, satClkCorr, satClkDrift] = satpos_zcj(launchTime, ...
                                [trackDeepIn.PRN], eph);

%% 最小二乘定位?
[xyzdt, vxyzdf, el, az, DOP] = ...
            leastSquarePos_zcj(satPositions, satVelocity, ...                         % 卫星位置和速度
            rawP + satClkCorr * settings.c, rawP_dot + satClkDrift * settings.c, ...  % 补偿掉卫星引起的误差后的伪距和伪距率
            settings);

% 记录结果
navSolut.X = xyzdt(1,1);    navSolut.Y = xyzdt(1,2);    navSolut.Z = xyzdt(1,3);    navSolut.dt = xyzdt(1,4);
navSolut.VX = vxyzdf(1,1);  navSolut.VY = vxyzdf(1,2);  navSolut.VZ = vxyzdf(1,3);  navSolut.df = vxyzdf(1,4);
navSolut.DOP = DOP;    navSolut.el = el;    navSolut.az = az;
navSolut.recvTime = settings.recvTime - navSolut.dt / settings.c;

end
