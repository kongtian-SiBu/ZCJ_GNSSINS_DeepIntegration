function [pseudoranges, transmitTimeatSat, Doppler] = calculatePseudoranges_zcj(trackResults, ...
                                                msOfTheSignal, ...
                                                activeChnList, settings, TOW)

%% 将首次定位时刻设置为最晚到达的子帧到达接收机的时刻
SampleNum = zeros(1, length(activeChnList));
Doppler   = zeros(1, length(activeChnList));

for i = 1 : length(activeChnList)
    k = activeChnList(i);
    SampleNum(i) = trackResults(k).absoluteSample(msOfTheSignal(k));
end
maxSampleNum = max(SampleNum);  % 最晚到达的信号位置对应的采样点数量

num_Cyclic = zeros(1, length(activeChnList));   % 其他子帧帧头距离最晚到达的帧头之间的CA码整周数
remSapmle = zeros(1, length(activeChnList));    % 不足一整周的剩余的采样点数
codePhaseTao = zeros(1, length(activeChnList)); % 不足一整周的最后一个采样点对应的码相位
remTime = zeros(1, length(activeChnList));      % 对应的发射时间
j = 1;
for i = 1 : length(activeChnList)
    k = activeChnList(i);
    while trackResults(k).absoluteSample(msOfTheSignal(k) + j) < maxSampleNum
        j = j + 1;
    end
    num_Cyclic(i) = j - 1;
    j = 1;
    remSapmle(i)    = maxSampleNum - trackResults(k).absoluteSample(msOfTheSignal(k) + num_Cyclic(i));
    
    Doppler(i)      = trackResults(k).carrFreq(msOfTheSignal(k) + num_Cyclic(i)) - settings.IF;

    remCodePhase    = trackResults(k).remCodePhase(msOfTheSignal(k) + num_Cyclic(i));
    codeFreq        = trackResults(k).codeFreq(msOfTheSignal(k) + num_Cyclic(i));
    codePhaseStep   = codeFreq / settings.samplingFreq;
    codePhaseTao(i) = remCodePhase + codePhaseStep * ( remSapmle(i) - 1 );
    remTime(i)      = num_Cyclic(i) * 1e-3 + codePhaseTao(i) / 1023 * 1e-3;   % 单位: s
end                                            

%%
travelTime = zeros(1, length(activeChnList));               
channelList = reshape(activeChnList, 1, length(activeChnList));
for channelNr = 1:length(channelList)
    
    %--- 计算信号传输时间 -----------------------------------------    
    travelTime(channelNr) = settings.recvTime - (TOW + remTime(channelNr));   % p = c (tr - tx)
end

pseudoranges = travelTime .* settings.c;     % 伪距
transmitTimeatSat = TOW + remTime;

end

