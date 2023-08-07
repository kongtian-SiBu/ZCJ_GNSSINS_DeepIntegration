%% GNSS\INS松组合，不同于GNSS_SDR,本文件下所有通道是交替执行，而不是一个通道运行结束后再运行下一个通道
%% 仅仅是在GNSS_SPP基础上嵌入了惯导模块
% trackResults, channel, TOW, eph, subFrameStart 需要提前准备好
settings = initSettings();
[fid, ~] = fopen(settings.fileName, 'rb');

positioningTime = TOW + settings.navSolPeriod / 1000;    % 首次松组合时刻

%% INS相关信息初始化，需提前安装好PSINS
glvs
psinstypedef(153);
trj = trjfile('trj_sc522.mat');
% initial settings
[nn, ts, nts] = nnts(2, trj.ts);
imuerr = imuerrset(5,1000,0.05,30);
imu = imuadderr(trj.imu, imuerr);
davp0 = avperrset([10;10;60], 0.5, [2; 2; 6]); 

ins = insinit(avpadderr(trj.avp0,davp0), ts);
% KF filter
rk = poserrset([1;1;3]);
kf = kfinit(ins, davp0, imuerr, rk);
kf.Pmin = [avperrset(0.01,1e-4,0.1); gabias(1e-3, [1,10])].^2;  kf.pconstrain=1;

%% 惯导运行至首次松组合之前
k = 1;
k1 = 1;            
t = imu(k1, end);   % INS时间戳
while t < positioningTime - 518400      % 518400是本次仿真的起始GPS周内时刻
    k1 = k+nn-1;
    wvm = imu(k:k1,1:6);  t = imu(k1,end);
    ins = insupdate(ins, wvm);
    kf.Phikk_1 = kffk(ins);
    kf = kfupdate(kf);

    k = k + nn;      % 用于标识惯导数据用到第几行了
end

%% 跟踪通道初始化，需要依赖传统跟踪环和接收机的一些信息，本代码认为信号成功跟踪的通道都能帧同步成功
activeChnList = find([trackResults.status] ~= '-');  
numActChnList = length(activeChnList);               
BitSyncTime = subFrameStart - 1;
ii = 1;
for ch = activeChnList     
    % 基本信息
    % 1) PRN号
    trackDeepIn(ii).PRN = trackResults(ch).PRN;    
    % 2) 跟踪状态
    trackDeepIn(ii).status = trackResults(ch).status;
    % 3) 首次出现帧头的位置
    trackDeepIn(ii).SamplePos = trackResults(ch).absoluteSample(BitSyncTime(ch));
    

    trackDeepIn(ii).codeFreq = trackResults(ch).codeFreq(BitSyncTime(ch));

    trackDeepIn(ii).remCodePhase = trackResults(ch).remCodePhase(BitSyncTime(ch));

    trackDeepIn(ii).carrFreq = trackResults(ch).carrFreq(BitSyncTime(ch));

    trackDeepIn(ii).remCarrPhase = trackResults(ch).remCarrPhase(BitSyncTime(ch));


    trackDeepIn(ii).codeError = trackResults(ch).dllDiscr(BitSyncTime(ch));
    trackDeepIn(ii).codeNco = trackResults(ch).dllDiscrFilt(BitSyncTime(ch));

    trackDeepIn(ii).carrError = trackResults(ch).pllDiscr(BitSyncTime(ch));
    trackDeepIn(ii).carrNco = trackResults(ch).pllDiscrFilt(BitSyncTime(ch));
    

    trackDeepIn(ii).carrFreqBasis = channel(ch).acquiredFreq;
    
    trackDeepIn(ii).numOfCoInt = 0;

    ii = ii + 1;
end


SamplePosatFirstFrame = zeros(1, numActChnList);
for ch = 1 : numActChnList
    SamplePosatFirstFrame(ch) = trackDeepIn(ch).SamplePos;
end

settings.recvTime = TOW + (settings.startOffset)/1000;  

recvTimeforFirstFrameperChannel = getTimeforFirstFrameEachChannel(settings, SamplePosatFirstFrame);
for ii = 1 : numActChnList
    trackDeepIn(ii).recvTime = recvTimeforFirstFrameperChannel(ii);  
end

%% GNSS运行至首次定位时刻之前
for ii = 1 : numActChnList
    trackans = trackDeepIn(ii);
    while trackans.recvTime < positioningTime    % 若该通道距离给定时刻的时间差超过一次相干积分，则进行一次相干积分
        trackDeepIn(ii) = trackans;
        trackans = perChannelTrackOnce(trackans, settings, fid);
    end
end

%% GNSS INS 松组合
roundTime = 70;    % 松组合次数
navResults = [];
navResults.X = zeros(1, roundTime); navResults.Y = zeros(1, roundTime); navResults.Z = zeros(1, roundTime); navResults.dt = zeros(1, roundTime);
navResults.VX = zeros(1, roundTime); navResults.VY = zeros(1, roundTime); navResults.VZ = zeros(1, roundTime);

%settings.pllNoiseBandwidth = 3;
%settings.dllNoiseBandwidth = 2;   

for currMeasNr = 1 : roundTime
    currMeasNr
    settings.recvTime = positioningTime;

    %% 松组合过程
    % 1. GNSS观测值计算
    navSolut = postNavLoose(trackDeepIn, settings, eph, TOW);
    [phi, lambda, h] = cart2geo(navSolut.X, navSolut.Y, navSolut.Z, 5);
    posGPS = [phi * pi/180; lambda * pi/180; h];  % pos in BLH
    
    % 2. EKF 
    kf = kfupdate(kf, ins.pos-posGPS, 'M');   % 实测距离-理论距离 ，此处不能弄反了  ins.pos - posGPS，反了之后不收敛
    [kf, ins] = kffeedback(kf, ins, 1, 'avp');

    % 3. 记录松组合的结果，将结果再转换回ECEF系
    [posX, posY, posZ] = geo2cart(ins.avp(7,1), ins.avp(8,1), ins.avp(9,1), 5);
    Cenu2xyz = [-sin(ins.pos(2))                  cos(ins.pos(2))                  0
                -sin(ins.pos(1))*cos(ins.pos(2)) -sin(ins.pos(1))*sin(ins.pos(2))  cos(ins.pos(1))
                 cos(ins.pos(1))*cos(ins.pos(2))  cos(ins.pos(1))*sin(ins.pos(2))  sin(ins.pos(1))];
    vxyz = Cenu2xyz' * ins.vn;
    navResults.X(1, currMeasNr) = posX;
    navResults.Y(1, currMeasNr) = posY;
    navResults.Z(1, currMeasNr) = posZ;
    navResults.dt(1,currMeasNr) = navSolut.dt;
    navResults.VX(1, currMeasNr)= vxyz(1);
    navResults.VY(1, currMeasNr)= vxyz(2);
    navResults.VZ(1, currMeasNr)= vxyz(3);

    % 4. 钟差修正
    for ii = 1 : numActChnList
        trackDeepIn(ii).recvTime = trackDeepIn(ii).recvTime - navSolut.dt / settings.c;  
    end
    
    % 5. 下一次松组合时刻
    positioningTime = positioningTime + settings.navSolPeriod / 1000;

    % 6. INS运行至下一次组合时刻
    while t < positioningTime - 518400     
        k1 = k+nn-1;
        wvm = imu(k:k1,1:6);  t = imu(k1,end);
        ins = insupdate(ins, wvm);
        kf.Phikk_1 = kffk(ins);
        kf = kfupdate(kf);
    
        k = k + nn;      
    end

    % 7. GNSS运行至下一次组合时刻
    for ii = 1 : numActChnList
        trackans = trackDeepIn(ii);
        while trackans.recvTime < positioningTime     
            trackDeepIn(ii) = trackans;
            trackans = perChannelTrackOnce(trackans, settings, fid);
        end
    end
    
end


