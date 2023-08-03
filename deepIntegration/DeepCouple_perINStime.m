%% GNSS INS 深组合，在紧组合的基础上对跟踪环路进行辅助
% trackResults, channel, TOW, eph, subFrameStart 需要提前准备好
close all;

settings = initSettings();
[fid, ~] = fopen(settings.fileName, 'rb');

positioningTime = TOW + settings.navSolPeriod / 1000;     

%% INS模块初始化
glvs
ggpsvars
psinstypedef('test_SINS_GPS_tightly_def');
trj = trjfile('trj_sc522.mat');
[nn, ts, nts] = nnts(2, diff(trj.imu(1:2,end)));    
avp0 = trj.avp0;  

davp = avperrset([10;10;60], 0.5, [2; 2; 6]);     

% ins = insinit(avpadderr(trj.avp0, davp), ts);
ins = insinit(trj.avp0, ts);

imuerr = imuerrset(5,1000,0.05,30);  
trj.imu = imuadderr(trj.imu, imuerr);

kf = kfinit(ins, davp, imuerr);

%% INS运行至首次组合时刻
k = 1;
k1 = 1;                 
t = trj.imu(k1, end);   
while t < positioningTime - 518400     
    k1 = k+nn-1;
    wvm = trj.imu(k:k1,1:6);  t = trj.imu(k1,end);
    ins = insupdate(ins, wvm);
    kf.Phikk_1 = kffk(ins);
    kf = kfupdate(kf);

    k = k + nn;      
end


%% GNSS跟踪结构体初始化
activeChnList = find([trackResults.status] ~= '-');  
numActChnList = length(activeChnList);               
BitSyncTime = subFrameStart - 1;
ii = 1;
for ch = activeChnList      
    trackDeepIn(ii).PRN = trackResults(ch).PRN;  
    
    trackDeepIn(ii).status = trackResults(ch).status;

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

    trackDeepIn(ii).carrFreqBasis = trackResults(ch).carrFreq(BitSyncTime(ch));
    trackDeepIn(ii).carrNco = 0;
    
    trackDeepIn(ii).numOfCoInt = 0;

    trackProcess(ii).codeErrorList = [];
    trackProcess(ii).carrErrorList = [];
    trackProcess(ii).codeFreqList = [];
    trackProcess(ii).carrFreqList = [];   
    trackProcess(ii).PLI = [];
    
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

%% GNSS运行至于首次组合时刻
I_P_1_list = [];
Q_P_1_list = [];


for ii = 1 : numActChnList
    trackans = trackDeepIn(ii);
    while trackans.recvTime < positioningTime    
        trackDeepIn(ii) = trackans;
        [trackans, I_P, Q_P] = perChannelTrackOnce(trackans, settings, fid);
        
        if ii == 1
            I_P_1_list = [I_P_1_list, I_P];
            Q_P_1_list = [Q_P_1_list, Q_P];
        end
        
        trackProcess(ii).codeErrorList = [trackProcess(ii).codeErrorList, trackans.codeError];
        trackProcess(ii).carrErrorList = [trackProcess(ii).carrErrorList, trackans.carrError];
        trackProcess(ii).codeFreqList = [trackProcess(ii).codeFreqList, trackans.codeFreq];
        trackProcess(ii).carrFreqList = [trackProcess(ii).carrFreqList, trackans.carrFreq];
        trackProcess(ii).PLI = [trackProcess(ii).PLI, (I_P^2-Q_P^2)/(I_P^2+Q_P^2)];
    end
end


%% 深组合
roundTime = 40;   
navResults = [];
navResults.X = zeros(1, roundTime); navResults.Y = zeros(1, roundTime); navResults.Z = zeros(1, roundTime); navResults.dt = zeros(1, roundTime);
navResults.VX = zeros(1, roundTime); navResults.VY = zeros(1, roundTime); navResults.VZ = zeros(1, roundTime);

settings.pllNoiseBandwidth = 3;    % 对于test_522这组数据，PLL的噪声带宽为3时刚刚失锁
settings.dllNoiseBandwidth = 2;

oldAidFreq = zeros(numActChnList, 1);

for currMeasNr = 1 : roundTime
    currMeasNr   
    settings.recvTime = positioningTime;    
    
    % 先进行一次单点定位，缩小时间差
    if currMeasNr < 2
        navSolut_1 = postNavLoose(trackDeepIn, settings, eph, TOW);
        positioningTime = positioningTime + settings.navSolPeriod / 1000;
        
        % 钟差修正
        for ii = 1 : numActChnList
            trackDeepIn(ii).recvTime = trackDeepIn(ii).recvTime - navSolut_1.dt / settings.c;  
        end
        
        % 惯导运行至于下一次组合时刻
        while t < positioningTime - 518400
            k1 = k+nn-1;
            wvm = trj.imu(k:k1,1:6);  t = trj.imu(k1,end);
            ins = insupdate(ins, wvm);
            kf.Phikk_1 = kffk(ins);
            kf = kfupdate(kf);

            k = k + nn;       
        end
        
        % GNSS 运行至于下一次组合时刻
        for ii = 1 : numActChnList
            trackans = trackDeepIn(ii);                
            
            oldAidFreq(ii,1) = trackans.carrFreq - settings.IF;      
            
            while trackans.recvTime < positioningTime     
                trackDeepIn(ii) = trackans;
                [trackans, I_P, Q_P] = perChannelTrackOnce(trackans, settings, fid);             

                if trackans.recvTime < positioningTime
                    trackProcess(ii).codeErrorList = [trackProcess(ii).codeErrorList, trackans.codeError];
                    trackProcess(ii).carrErrorList = [trackProcess(ii).carrErrorList, trackans.carrError];
                    trackProcess(ii).codeFreqList = [trackProcess(ii).codeFreqList, trackans.codeFreq];
                    trackProcess(ii).carrFreqList = [trackProcess(ii).carrFreqList, trackans.carrFreq];
                    trackProcess(ii).PLI = [trackProcess(ii).PLI, (I_P^2-Q_P^2)/(I_P^2+Q_P^2)];
                    
                    if ii == 1
                        I_P_1_list = [I_P_1_list, I_P];
                        Q_P_1_list = [Q_P_1_list, Q_P];
                    end
        
                end
            end
        end
        
    else
             
    %% 开始紧组合
    % 1. GNSS观测值
    navSolut = postNavTight(trackDeepIn, settings, eph, TOW);
    % 只保留一颗卫星观测值
%     navSolut.rawP = navSolut.rawP(1); navSolut.satPositions = navSolut.satPositions(:,1); 
%     navSolut.satVelocity = navSolut.satVelocity(:,1);  navSolut.satClkCorr = navSolut.satClkCorr(1); 
    
    [posxyz, ~] = blh2xyz(ins.pos)    
    [rho, LOS, AzEl] = rhoSatRec(navSolut.satPositions', posxyz, navSolut.rawP');
    el = AzEl(:,2); el(el<15*pi/180) = 1*pi/180;  P = diag(sin(el.^2));
    delta_rawP = navSolut.rawP' + settings.c * navSolut.satClkCorr' - rho;   
    
    % 2. EKF
    kf.Hk = kfhk(ins, LOS);  
    kf.Rk = P^-1 * 10^2;
    kf = kfupdate(kf, delta_rawP);
    [kf, ins] = kffeedback(kf, ins, 1, 'avp');  
    
    % 3. 将结果转换为ECEF坐标系
    [posX, posY, posZ] = geo2cart(ins.avp(7,1), ins.avp(8,1), ins.avp(9,1), 5);
    pos_XYZ = [posX; posY; posZ];
    Cenu2xyz = [-sin(ins.pos(2))                  cos(ins.pos(2))   0
                -sin(ins.pos(1))*cos(ins.pos(2)) -sin(ins.pos(1))*sin(ins.pos(2))  cos(ins.pos(1))
                 cos(ins.pos(1))*cos(ins.pos(2))  cos(ins.pos(1))*sin(ins.pos(2))  sin(ins.pos(1))];
    vxyz = Cenu2xyz' * ins.vn;
    navResults.X(1, currMeasNr) = posX;  navResults.Y(1, currMeasNr) = posY;
    navResults.Z(1, currMeasNr) = posZ;  navResults.dt(1, currMeasNr) = kf.xk(end-1);
    navResults.VX(1, currMeasNr)= vxyz(1);
    navResults.VY(1, currMeasNr)= vxyz(2);
    navResults.VZ(1, currMeasNr)= vxyz(3);
    
    % 3.5 构造反馈量
    % correctedP = navSolut.rawP + navSolut.satClkCorr .* settings.c - kf.xk(end-1); 
    % [~, ~, ~, vrs] = rhoSatRec_zcj(navSolut.satPositions', pos_XYZ, correctedP', navSolut.satVelocity', ins.vn);
    % 上面那种计算方法和下面这行没啥区别
    [~, ~, ~, vrs] = rhoSatRec_zcj(navSolut.satPositions', pos_XYZ, navSolut.rawP', navSolut.satVelocity', ins.vn);  % 根据接收机信息反推多普勒数值
    dopplerFeedback = -vrs / settings.c * 1575.42e6;   % numOfSat * 1
    
    % 4. 钟差修正
    for ii = 1 : numActChnList
        trackDeepIn(ii).recvTime = trackDeepIn(ii).recvTime - kf.xk(end-1) / settings.c;  
    end
          
    % 5. 下一次定位时刻
    positioningTime = positioningTime + settings.navSolPeriod / 1000;   
       
    aidFreq = dopplerFeedback;  
    
    %% 深组合反馈模块
    while t < positioningTime - 518400    
        
        % 沈聪的操作，可防止carrNco过大，我很难领悟其中奥妙
        if currMeasNr == 2
            for iii = 1 : numActChnList
                trackDeepIn(ii).carrNco = 0;
            end
        end
               
        % 6.1 跳过当前相干积分
        for ii = 1 : numActChnList    
            % 该方程可参考一些论文，总之正确性有待验证
            [trackDeepIn(ii), I_P, Q_P] = perChannelTrackOnce_DeepIn(trackDeepIn(ii), settings, fid, oldAidFreq(ii), aidFreq(ii) - oldAidFreq(ii));
            
            if ii == 1
                I_P_1_list = [I_P_1_list, I_P];
                Q_P_1_list = [Q_P_1_list, Q_P];
            end
            
            trackProcess(ii).codeErrorList = [trackProcess(ii).codeErrorList, trackDeepIn(ii).codeError];
            trackProcess(ii).carrErrorList = [trackProcess(ii).carrErrorList, trackDeepIn(ii).carrError];
            trackProcess(ii).codeFreqList = [trackProcess(ii).codeFreqList, trackDeepIn(ii).codeFreq];
            trackProcess(ii).carrFreqList = [trackProcess(ii).carrFreqList, trackDeepIn(ii).carrFreq];   
            trackProcess(ii).PLI = [trackProcess(ii).PLI, (I_P^2-Q_P^2)/(I_P^2+Q_P^2)];
        end
        
        % 6.3 惯导运行至下一次更新时刻，即惯导更新一次
        k1 = k+nn-1;
        wvm = trj.imu(k:k1,1:6);  t = trj.imu(k1,end);
        ins = insupdate(ins, wvm);
        kf.Phikk_1 = kffk(ins);
        kf = kfupdate(kf);

        k = k + nn;       
        
        % 6.4 跟踪环运行至惯导更新时刻
        for ii = 1 : numActChnList
            trackans = trackDeepIn(ii);
            
            kkk = 1;
            
            while trackans.recvTime < t + 518400         
                trackDeepIn(ii) = trackans;
                % 可参考某些论文，我也只是懂个大概
                [trackans, I_P, Q_P] = perChannelTrackOnce_DeepIn(trackans, settings, fid, oldAidFreq(ii), kkk * (aidFreq(ii) - oldAidFreq(ii)));
                
                kkk = kkk + 1;

                if trackans.recvTime < t + 518400
                    trackProcess(ii).codeErrorList = [trackProcess(ii).codeErrorList, trackans.codeError];
                    trackProcess(ii).carrErrorList = [trackProcess(ii).carrErrorList, trackans.carrError];
                    trackProcess(ii).codeFreqList = [trackProcess(ii).codeFreqList, trackans.codeFreq];
                    trackProcess(ii).carrFreqList = [trackProcess(ii).carrFreqList, trackans.carrFreq];
                    trackProcess(ii).PLI = [trackProcess(ii).PLI, (I_P^2-Q_P^2)/(I_P^2+Q_P^2)];
                    
                    if ii == 1
                        I_P_1_list = [I_P_1_list, I_P];
                        Q_P_1_list = [Q_P_1_list, Q_P];
                    end
                end
            end
        end
              
        
        % 6.5 重新计算惯导的反馈量
        settings.recvTime = settings.recvTime + nts;
        navSolut = postNavTight(trackDeepIn, settings, eph, TOW);  % GNSS观测值
        [posxyz, ~] = blh2xyz(ins.pos);
        % correctedP = navSolut.rawP + navSolut.satClkCorr .* settings.c - kf.xk(end-1);  
        correctedP = navSolut.rawP;  % 和上一行区别不大
        [~, ~, ~, vrs] = rhoSatRec_zcj(navSolut.satPositions', posxyz, correctedP', navSolut.satVelocity', ins.vn);
        dopplerFeedback = -vrs / settings.c * 1575.42e6;  
        
        oldAidFreq = aidFreq;
        aidFreq = dopplerFeedback;
    end
       
    
    end  % if currMeasNr == 1
end










