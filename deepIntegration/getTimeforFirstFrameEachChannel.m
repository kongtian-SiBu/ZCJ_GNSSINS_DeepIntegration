function recvTimeforFirstFrameperChannel = getTimeforFirstFrameEachChannel(settings, firstFrameSamplePos)
%% 计算每个通道首次出现帧头的时候对应的本地接收机时间估计值
%
% 输入参数: 
%        - settings: 接收机相关参数
%        - firstFrameSamplePos: 各通道首次出现帧头时对应的采样点数
%
% 输出参数:?
%        - recvTimeforFirstFrameperChannel: 各通道首次出现帧头时接收机的本地时间估计值
%
% -------------------------------------------------------------------------
%%
numActChnList = length(firstFrameSamplePos);
recvTimeforFirstFrameperChannel = zeros(1, numActChnList);
maxSamplePos = max(firstFrameSamplePos);

for ii = 1 : numActChnList
    remTime = (maxSamplePos - firstFrameSamplePos(ii) ) / settings.samplingFreq;
    recvTimeforFirstFrameperChannel(ii) = settings.recvTime - remTime;
end

end