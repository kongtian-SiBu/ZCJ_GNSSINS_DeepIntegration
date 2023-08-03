%% GNSS单点定位，不同于GNSS_SDR,本文件下所有通道是交替执行，而不是一个通道运行结束后再运行下一个通道
settings = initSettings();
[fid, message] = fopen(settings.fileName, 'rb');

%% 跟踪通道初始化，需要依赖传统跟踪环和接收机的一些信息，本代码认为信号成功跟踪的通道都能帧同步成功
% 记录首次出现帧头的时候，各通道的跟踪情况
activeChnList = find([trackResults.status] ~= '-');   % 成功跟踪信号的通道
numActChnList = length(activeChnList);                % 成功跟踪信号的通道数量
BitSyncTime = subFrameStart - 1;
ii = 1;
for ch = activeChnList      % 去除跟踪失败的通道，trackDeepIn大小等于成功跟踪的通道数[][][][][][][][][]
    % 基本信息
    % 1) PRN号
    trackDeepIn(ii).PRN = trackResults(ch).PRN;     
    % 2) 跟踪状态
    trackDeepIn(ii).status = trackResults(ch).status;
    % 3) 首次出现帧头的位置
    trackDeepIn(ii).SamplePos = trackResults(ch).absoluteSample(BitSyncTime(ch));
    
    % 跟踪信息
    % 1. 码频率
    trackDeepIn(ii).codeFreq = trackResults(ch).codeFreq(BitSyncTime(ch));
    % 2. 码相位残差
    trackDeepIn(ii).remCodePhase = trackResults(ch).remCodePhase(BitSyncTime(ch));
    % 3. 载波频率
    trackDeepIn(ii).carrFreq = trackResults(ch).carrFreq(BitSyncTime(ch));
    % 4. 载波相位差
    trackDeepIn(ii).remCarrPhase = trackResults(ch).remCarrPhase(BitSyncTime(ch));

    % 5. 码NCO
    trackDeepIn(ii).codeError = trackResults(ch).dllDiscr(BitSyncTime(ch));
    trackDeepIn(ii).codeNco = trackResults(ch).dllDiscrFilt(BitSyncTime(ch));
    % 6. 载波NCO
    trackDeepIn(ii).carrError = trackResults(ch).pllDiscr(BitSyncTime(ch));
    trackDeepIn(ii).carrNco = trackResults(ch).pllDiscrFilt(BitSyncTime(ch));
    
    % 7. 载波基准频率，该数值来自捕获
    trackDeepIn(ii).carrFreqBasis = channel(ch).acquiredFreq;
    
    % 8. 从子帧头开始算起，又进行了几次相干积分（用于计算发射时间）
    trackDeepIn(ii).numOfCoInt = 0;
    
    % 9. 保存一些中间结果，用于查看跟踪环情况
    trackProcess(ii).codeErrorList = [];
    trackProcess(ii).carrErrorList = [];
    trackProcess(ii).codeFreqList = [];
    trackProcess(ii).carrFreqList = []; 

    ii = ii + 1;
end

% 单独记录下首次出现帧头时的跟踪情况，便于调试
SamplePosatFirstFrame = zeros(1, numActChnList);
for ch = 1 : numActChnList
    SamplePosatFirstFrame(ch) = trackDeepIn(ch).SamplePos;
end

%% 开始单点定位
roundTime = 10;    % 定位次数，不要超过GNSS Obs的个数
navResults = [];
navResults.X = zeros(1, roundTime); navResults.Y = zeros(1, roundTime); navResults.Z = zeros(1, roundTime); navResults.dt = zeros(1, roundTime);
navResults.VX = zeros(1, roundTime); navResults.VY = zeros(1, roundTime); navResults.VZ = zeros(1, roundTime); navResults.df = zeros(1, roundTime);


for currMeasNr = 1 : roundTime
    currMeasNr  
    
%     if currMeasNr > 16                     % 在运行过程中修改带宽
%     settings.pllNoiseBandwidth = 5;
%     settings.dllNoiseBandwidth = 0.5;
%     end
    
    if currMeasNr == 1
        positioningTime = TOW + settings.navSolPeriod / 1000;   % 规定一个首次实现定位的时刻
        
        % 最晚的一个通道首次出现子帧帧头时的接收机本地时间估计值(在实现首次定位之前无法得到准确的接收机本地时间)
        settings.recvTime = TOW + (settings.startOffset)/1000;  

        % 每个通道首次出现子帧帧头时接收机本地时间估计值?
        recvTimeforFirstFrameperChannel = getTimeforFirstFrameEachChannel(settings, SamplePosatFirstFrame);
        for ii = 1 : numActChnList
            trackDeepIn(ii).recvTime = recvTimeforFirstFrameperChannel(ii);  % 将该信息合并入trackDeepIn结构体
        end
    end
    
    settings.recvTime = positioningTime; % 设置接收机本地时钟为定位时刻，跟踪环先运行至该时刻
    
    % 每个通道进行跟踪直到到达定位时刻之前（所谓通道间的准同步）
    for ii = 1 : numActChnList
        trackans = trackDeepIn(ii);
        while trackans.recvTime < positioningTime    % 若该通道距离给定时刻的时间差超过一次相干积分，则进行一次相干积分
            trackDeepIn(ii) = trackans;
            trackans = perChannelTrackOnce(trackans, settings, fid);
            
            if trackans.recvTime < positioningTime
                trackProcess(ii).codeErrorList = [trackProcess(ii).codeErrorList, trackans.codeError];
                trackProcess(ii).carrErrorList = [trackProcess(ii).carrErrorList, trackans.carrError];
                trackProcess(ii).codeFreqList = [trackProcess(ii).codeFreqList, trackans.codeFreq];
                trackProcess(ii).carrFreqList = [trackProcess(ii).carrFreqList, trackans.carrFreq];
        
            end
        end
    end
    
    navSolut = postNavLoose(trackDeepIn, settings, eph, TOW);
    navResults.X(1,currMeasNr) = navSolut.X; navResults.Y(1,currMeasNr) = navSolut.Y; navResults.Z(1,currMeasNr) = navSolut.Z;
    navResults.VX(1,currMeasNr) = navSolut.VX; navResults.VY(1,currMeasNr) = navSolut.VY; navResults.VZ(1,currMeasNr) = navSolut.VZ;
    navResults.dt(1,currMeasNr) = navSolut.dt;  navResults.df(1,currMeasNr) = navSolut.df;
    
    % 钟差修正
    for ii = 1 : numActChnList
        trackDeepIn(ii).recvTime = trackDeepIn(ii).recvTime - navSolut.dt / settings.c;   
    end
    
    % 计算下一次定位时刻
    positioningTime = positioningTime + settings.navSolPeriod / 1000;
end


