function [trackResults, channel]= trackpll2nd(fid, channel, settings)
% Performs code and carrier tracking for all channels.
%
%[trackResults, channel] = tracking(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute starting 
%                       positions of spreading codes, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------


%% 初始化信号跟踪结构体 ============================================

% 信号跟踪状态
trackResults.status         = '-';      % '-'表示该通道无信号

% 每一次相干积分结束后，数据文件指针的位置
trackResults.absoluteSample = zeros(1, settings.msToProcess);

% CA码频率
trackResults.codeFreq       = inf(1, settings.msToProcess);

% 载波频率
trackResults.carrFreq       = inf(1, settings.msToProcess);

% E,L,P支路I路
trackResults.I_P            = zeros(1, settings.msToProcess);
trackResults.I_E            = zeros(1, settings.msToProcess);
trackResults.I_L            = zeros(1, settings.msToProcess);

% E,L,P支路Q路
trackResults.Q_E            = zeros(1, settings.msToProcess);
trackResults.Q_P            = zeros(1, settings.msToProcess);
trackResults.Q_L            = zeros(1, settings.msToProcess);

% 鉴相器输出
trackResults.dllDiscr       = inf(1, settings.msToProcess);
trackResults.pllDiscr       = inf(1, settings.msToProcess);

% NCO输出
trackResults.dllDiscrFilt   = inf(1, settings.msToProcess);
trackResults.pllDiscrFilt   = inf(1, settings.msToProcess);

%------------------------ 为每一个通道初始化 -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% 初始化信号跟踪的相关参数 ==========================================

codePeriods = settings.msToProcess;     % For GPS one C/A code is one ms

%--- DLL 参数 --------------------------------------------------------
% 相关器间隔 (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

% 码环更新周期，在该函数中等于相干积分时间
PDIcode = 0.001;

% 2阶DLL环路滤波器系数计算
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
                                    settings.dllDampingRatio, ...
                                    1.0);

%--- PLL 参数 --------------------------------------------------------
% 载波环更新周期，在该函数中等于相干积分时间
PDIcarr = 0.001;

% 2阶PLL环路滤波器系数计算
a2 = 1.414;
wn = settings.pllNoiseBandwidth / 0.53;
kca1 = wn ^ 2 * PDIcarr / settings.pllLoopGain;
kca2 = a2 * wn / settings.pllLoopGain;

hwb = waitbar(0,'Tracking...');

%% 开始信号跟踪 ==============================================
for channelNr = 1:settings.numberOfChannels
    
    % 如果成功捕获，则该通道PRN号不为0
    if (channel(channelNr).PRN ~= 0)

        trackResults(channelNr).PRN     = channel(channelNr).PRN;
        
        % 移动到捕获到的码相位处，如果码相位计算准确，移动后码相位应该为0
        fseek(fid, ...
              settings.fileType * settings.dataFormat * (settings.skipNumberOfSamples + channel(channelNr).codePhase-1), ...
              'bof');

        % 生成CA码
        caCode = generateCAcode(channel(channelNr).PRN);
        % 前后各添加一个码片
        caCode = [caCode(1023) caCode caCode(1)];

        %--- 初始化一些变量用于捕获 ------------------------------
      
        codeFreq      = settings.codeFreqBasis;     % 码频率初始化为1.023e6即可  
        remCodePhase  = 0.0;                        % 下一轮的码起始相位
    
        carrFreq      = channel(channelNr).acquiredFreq;  % 载波频率初始化为捕获到的频率
        carrFreqBasis = channel(channelNr).acquiredFreq;  % 载波环的基准频率初始化为捕获到的频率
        remCarrPhase  = 0.0;                              % 下一轮的载波起始相位

        % 码环环路滤波器的中间量初始化
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;

        % 载波环环路滤波器的中间量初始化
        zkminus1     = 0.0;             % 对应于PLL信号流图第一个延时器之后的位置
        oldCarrError = 0.0;
        
        % 通道之间串行跟踪，即一个通道全部跟踪结束之后再跟踪下一通道
        for loopCnt =  1:codePeriods
            
%% 跟踪进度条 -------------------------------------------------------------
            % 进度条每50次相干积分更新一次
            if (rem(loopCnt, 50) == 0)
                try
                    waitbar(loopCnt/codePeriods, ...
                            hwb, ...
                            ['Tracking: Ch ', int2str(channelNr), ...
                            ' of ', int2str(settings.numberOfChannels), ...
                            '; PRN#', int2str(channel(channelNr).PRN), ...
                            '; Completed ',int2str(loopCnt), ...
                            ' of ', int2str(codePeriods), ' msec']);                       
                catch
                    % 进度条如果被关闭则自动退出跟踪
                    disp('Progress bar closed, exiting...');
                    return
                end
            end

%% 处理一个相干积分数据块 ------------------------------------------------                        
            
            % 对于当前码环计算得到的码频率来说，一个采样点对应的码相位步进量是多少 
            codePhaseStep = codeFreq / settings.samplingFreq;   
            
            % 对于当前码环计算得到的码频率、当前周期的起始码相位、每个采样点对应的相位步进量来说，
            % 生成一个周期的码对应多少采样点（该数值的大小在该函数中约等于1ms对应的采样点个数）
            blksize = ceil((settings.codeLength - remCodePhase) / codePhaseStep);  
            
            % 读取对应数量的采样点
            [rawSignal, samplesRead] = fread(fid, ...
                                             settings.fileType * blksize, settings.dataType);
            rawSignal = transpose(rawSignal);  % 转置，注意这里不要用'，'号表示共轭转置
            if settings.fileType == 2
                rawSignalI = rawSignal(1:2:end);
                rawSignalQ = rawSignal(2:2:end);
                rawSignal  = rawSignalI + 1j * rawSignalQ; Flag = 0; 
                % rawSignal  = rawSignalI + 0 * rawSignalQ;  Flag = 1;  % I,Q not combined               
            end
                                                               
            % 如果数据不够了，直接退出
            if (samplesRead ~= settings.fileType * blksize)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
 %               fclose(fid);
                return
            end

%% 生成E、L、P支路对应的码 ------------------------------------------
            % E支路
            tcode       = (remCodePhase-earlyLateSpc) : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            earlyCode   = caCode(tcode2);
            
            % L支路
            tcode       = (remCodePhase+earlyLateSpc) : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            lateCode    = caCode(tcode2);
            
            % P支路
            tcode       = remCodePhase : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase);
            tcode2      = ceil(tcode) + 1;
            promptCode  = caCode(tcode2);
            
            remCodePhase = (tcode(blksize) + codePhaseStep) - 1023.0;  % 计算下一轮CA码的起始相位，保证码相位连续

%% 生成载波 -------------------------------------------------------
            time    = (0:blksize) ./ settings.samplingFreq;
            
            % 生成载波相位
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;  % remCarrPhase作为本轮的起始相位，这样能保持每一轮的相位连续        
            remCarrPhase = rem(trigarg(blksize+1), (2 * pi));          % 计算下一轮载波的起始相位
    
            % 实数据和复数据分开处理，是否有必要？以后再说
            if settings.fileType == 1 || Flag == 1
                carrCos = cos(trigarg(1:blksize));
                carrSin = sin(trigarg(1:blksize));
                qBasebandSignal = carrCos .* rawSignal;
                iBasebandSignal = carrSin .* rawSignal;
                             
            else
                carr = exp(-1j* trigarg(1:blksize));
                qBasebandSignal = imag(carr .* rawSignal);
                iBasebandSignal = real(carr .* rawSignal);
            end
                        
            % 相干积分
            I_E = sum(earlyCode  .* iBasebandSignal);
            Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode   .* iBasebandSignal);
            Q_L = sum(lateCode   .* qBasebandSignal);            
            
%% PLL环路更新 -----------------------------------------------------
            % 鉴相器，atan(Q_P / I_P)的单位为弧度，除以2pi后单位变成 '周'，也等价于Hz
            carrError = atan(Q_P / I_P) / (2.0 * pi);
            
            % 环路滤波器的时域公式，是从Z域转换得到
            zk = zkminus1 + (carrError + oldCarrError) * kca1;  
            carrNco = zk / 2 + carrError * kca2;
            zkminus1 = zk;
            oldCarrError = carrError;
           
            % 更新载波频率
            carrFreq = carrFreqBasis + carrNco;
            trackResults(channelNr).carrFreq(loopCnt) = carrFreq; % 记录载波频率

%% DLL环路更新 -----------------------------------------------------
            % 码鉴相器，和书上相比少了个1/2，不过无所谓不影响
            codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
                (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));
            
            % 环路滤波器的时域公式，是从Z域转换得到，和载波环毫无区别
            codeNco = oldCodeNco + (tau2code/tau1code) * ...
                (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;
            
            % 注意此处应是减号。码相位完美同步时应该有E和L支路信号能量相等，都
            % 小于P支路信号能量。不妨设此轮E支路能量大于L支路，此时由鉴相器公式
            % 可得codeError > 0。此时环路认为本地码相位超前了输入信号的码相位，因此
            % 应降低本地码频率，稍微“等一等”输入信号。
            codeFreq = settings.codeFreqBasis - codeNco;
            % codeFreq = settings.codeFreqBasis - codeNco + (carrFreq - settings.IF) / 1540;   % 该公式为载波环辅助码环
            
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;

%% 记录信号跟踪的结果 ----------------------
            % 记录数据文件指针的位置，该数值的单位为采样点个数，该数值是为了从
            % 跟踪环提取观测值用的。而ftell的单位是字节，此处要注意单位匹配
            trackResults(channelNr).absoluteSample(loopCnt) = ftell(fid) / ...
                                    settings.dataFormat / settings.fileType;

            trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
            trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;
            
            % 用于画莉萨如图，观察E、P、L三支路信号能量变化
            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;
            
            trackResults(channelNr).remCodePhase(loopCnt)   = remCodePhase;
            trackResults(channelNr).remCarrPhase(loopCnt)   = remCarrPhase;
            
        end % for loopCnt

        % 简单的认为成功捕获到卫星信号的通道在跟踪环节都成功跟踪了信号
        % 这里其实应该进行信号跟踪检测，实时判断信号是否失锁
        % 但我比较懒
        trackResults(channelNr).status  = channel(channelNr).status;        
        
    end % if a PRN is assigned
end % for channelNr 

% 关闭进度条
close(hwb)
