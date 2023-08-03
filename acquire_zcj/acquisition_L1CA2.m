function acqResults = acquisition_L1CA2(longSignal, settings)
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 20 ms of raw signal from the front-end 
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the 
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number. 

% 本代码无法应对比特跳变，且不适用精确捕获

%% 初始化 =========================================================
tempFreq = settings.samplingFreq;
rawSignal = longSignal;     % 降采样之前的频率
downratio = 1;              % 降采样倍数初始化，若没有降采样则该数值为1
% 是否使用降采样
if settings.acqDownSample == 1
    downratio = settings.samplingFreq / settings.acqDownFreq;  % 降采样倍数
    longSignal = downsample(longSignal, downratio);
    settings.samplingFreq = settings.acqDownFreq;
end
signal0DC = rawSignal - mean(rawSignal);      % 用于精确捕获

% 一个CA码周期对应的采样点个数
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

% 创建两个连续的时间序列，由于限定了总相干积分时长不超过10ms，因此必然有一个时间序列没有比特翻转
signal1 = longSignal(1 : settings.acqnonCoIntime * settings.acqCoIntime * samplesPerCode);
signal2 = longSignal(settings.acqnonCoIntime * settings.acqCoIntime * samplesPerCode+1 : ...
                     2 * settings.acqnonCoIntime * settings.acqCoIntime * samplesPerCode);
 
% 采样率
ts = 1 / settings.samplingFreq;

% 载波的相位
phasePoints = (0 : (settings.acqnonCoIntime * settings.acqCoIntime * samplesPerCode - 1)) * 2 * pi * ts;

% 频率搜索单元个数
numberOfFrqBins = round(settings.acqSearchBand / settings.acqSearchStep) + 1;

% 根据采样率和相干积分时间生成CA码
caCodesTable = makeCaTable(settings);   % 其长度 = settings.acqCoIntime * samplesPerCode

%------------- 初始化数组用于存储捕获结果 ----------------------------------
% 弄到这里了，待会儿继续
% 存储搜索到的频率和码相位 (一颗卫星)
results = zeros(numberOfFrqBins, samplesPerCode); % 码相位数不会受积分时间的影响

% 载波频率网格
frqBins = zeros(1, numberOfFrqBins);

%--- acqResults结构体初始化 ------------------------------------------------
acqResults.carrFreq     = zeros(1, 32);     % 捕获到的频率
acqResults.codePhase    = zeros(1, 32);     % 捕获到的码相位

acqResults.peakMetric   = zeros(1, 32);     % Correlation peak ratios of the detected signals

fprintf('(');

% 开始捕获过程，使用基于FFT的并行码相位捕获算法
for PRN = settings.acqSatelliteList

%% 对信号进行相关操作 ======================================================   
    %---------- 对CA码进行FFT操作 ------------------------------------------
    caCodeFreqDom = conj(fft(caCodesTable(PRN, :)));
    
    for frqBinIndex = 1 : numberOfFrqBins       % 频率搜索
        
        % 对于实信号和复信号分开处理或许没有必要，以后再探究吧     
        acqRes1 = 0; acqRes2 = 0;       % 积分变量初始化，MATLAB自带广播机制，会将0变为0向量
        
        % 生成频率单元格
        frqBins(frqBinIndex) = settings.IF - (settings.acqSearchBand/2) + settings.acqSearchStep * (frqBinIndex - 1);
        
        for nonCoIntIndex = 1 : settings.acqnonCoIntime   % 非相干积分次数
            if settings.fileType == 1  % 处理实信号
                
                % 生成载波
                sinCarr = sin(frqBins(frqBinIndex) * phasePoints((nonCoIntIndex-1)*settings.acqCoIntime*samplesPerCode+1: ...
                                                                 nonCoIntIndex*settings.acqCoIntime*samplesPerCode));
                cosCarr = cos(frqBins(frqBinIndex) * phasePoints((nonCoIntIndex-1)*settings.acqCoIntime*samplesPerCode+1: ...
                                                                 nonCoIntIndex*settings.acqCoIntime*samplesPerCode));

                % 相干积分
                I1      = sinCarr .* signal1((nonCoIntIndex-1)*settings.acqCoIntime*samplesPerCode+1: ...
                                                                 nonCoIntIndex*settings.acqCoIntime*samplesPerCode);
                Q1      = cosCarr .* signal1((nonCoIntIndex-1)*settings.acqCoIntime*samplesPerCode+1: ...
                                                                 nonCoIntIndex*settings.acqCoIntime*samplesPerCode);
                I2      = sinCarr .* signal2((nonCoIntIndex-1)*settings.acqCoIntime*samplesPerCode+1: ...
                                                                 nonCoIntIndex*settings.acqCoIntime*samplesPerCode);
                Q2      = cosCarr .* signal2((nonCoIntIndex-1)*settings.acqCoIntime*samplesPerCode+1: ...
                                                                 nonCoIntIndex*settings.acqCoIntime*samplesPerCode);

                % 将时域信号转换到频域 
                IQfreqDom1 = fft(I1 + 1j*Q1) ./ length(I1);  % 防止数值过大
                IQfreqDom2 = fft(I2 + 1j*Q2) ./ length(I2);
                    
                % 频率乘积等价于时域卷积
                convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
                convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;
                
            else  % 处理复信号
                carr = exp(-1j * frqBins(frqBinIndex) * phasePoints((nonCoIntIndex-1)*settings.acqCoIntime*samplesPerCode+1: ...
                                                                 nonCoIntIndex*settings.acqCoIntime*samplesPerCode));
                IQ1 = carr .* signal1((nonCoIntIndex-1)*settings.acqCoIntime*samplesPerCode+1: ...
                                                                 nonCoIntIndex*settings.acqCoIntime*samplesPerCode);
                IQ2 = carr .* signal2((nonCoIntIndex-1)*settings.acqCoIntime*samplesPerCode+1: ...
                                                                 nonCoIntIndex*settings.acqCoIntime*samplesPerCode);
                                                             
                % 将时域信号转换到频域               
                IQfreqDom1 = fft(IQ1 ./ length(IQ1));  % 防止数值过大
                IQfreqDom2 = fft(IQ2 ./ length(IQ2));
                    
                % 频率乘积等价于时域卷积
                convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
                convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;                
            end
            
            acqRes1 = acqRes1 + abs(ifft(convCodeIQ1)) .^ 2;   % 在此处进行累加才是真正的非相干积分！abs()抹除相位，典型的非相干
            acqRes2 = acqRes2 + abs(ifft(convCodeIQ2)) .^ 2;          
        end  
      
        % 令积分数值较大的一个为捕获结果，因为值较大的一个更可能不涉及比特翻转
        if (max(acqRes1) > max(acqRes2))
            results(frqBinIndex, :) = acqRes1(1 : samplesPerCode);  % 捕获结果仅取第一个CA码周期即可
        else
            results(frqBinIndex, :) = acqRes2(1 : samplesPerCode);
        end
    end

    
%% 在acqRes中寻找自相关峰值 =========================================
    % 如果信号捕获成功，则第一峰值应明显大于第二峰值
    
    % 寻找捕获到的频率
    [~, frequencyBinIndex] = max(max(results, [], 2));
    
    % 寻找捕获到的码相位
    [peakSize, codePhase] = max(max(results));
    
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
    excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    excludeRangeIndex2 = codePhase + samplesPerCodeChip;

    % 去除第一峰值附近的采样点，继而寻找第二峰值
    if excludeRangeIndex1 < 2
        codePhaseRange = excludeRangeIndex2 : ...
                         (samplesPerCode + excludeRangeIndex1);
                         
    elseif excludeRangeIndex2 >= samplesPerCode
        codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                         excludeRangeIndex1;
    else
        codePhaseRange = [1:excludeRangeIndex1, ...
                          excludeRangeIndex2 : samplesPerCode];
    end

    % 寻找第二峰值
    secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));

    % 计算第一峰值与第二峰值的比值
    acqResults.peakMetric(PRN) = peakSize / secondPeakSize;
    
    % 记录捕获到的码相位
    if settings.acqDownSample == 1
        codePhase = (codePhase - 1) * downratio;  % 复原到降采样之前的码相位
    end
    acqResults.codePhase(PRN) = codePhase;
    
    % 如比值大于阈值，则说明捕获成功，继而进行精捕获过程
    if (peakSize/secondPeakSize) > settings.acqThreshold      
        fprintf('%02d ', PRN);   % 显示捕获成功的PRN号

%% 频率精确捕获。请注意此时的采样率必须是中频数据的两倍以上
        samplesPerCode_findacq = round(tempFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));
        ts_findacq = 1/tempFreq;
                    
        caCode = generateCAcode(PRN);
        codeValueIndex = floor((ts_findacq * (1:10*samplesPerCode_findacq)) / (1/settings.codeFreqBasis)); % 使用10ms长度的数据
        longCaCode = caCode((rem(codeValueIndex, 1023) + 1));     % 生成长度为10ms的CA码
        xCarrier = signal0DC(codePhase:(codePhase + 10*samplesPerCode_findacq-1)) .* longCaCode;   % 去除CA码，只保留载波
        fftNumPts = 8*(2^(nextpow2(length(xCarrier))));
        
        if settings.IF <= 0
            fftxc = abs(fft(xCarrier, fftNumPts)) ./ fftNumPts;         
            fftxc = fftshift(fftxc);
            ffs = settings.samplingFreq;
            fff = -ffs/2 : ffs/fftNumPts : ffs/2 - ffs/fftNumPts;
            [~, fftMaxIndex] = max(fftxc); 
            acqResults.carrFreq(PRN)  = fff(fftMaxIndex);
        else
            fftxc = abs(fft(xCarrier, fftNumPts)) ./ fftNumPts;         
            uniqFftPts = ceil((fftNumPts + 1) / 2); % focus on freq between(0,pi), exclude (pi,2pi)
            [~, fftMaxIndex] = max(fftxc);       
            fftFreqBins = (0 : uniqFftPts-1) * tempFreq / fftNumPts; 
            acqResults.carrFreq(PRN)  = fftFreqBins(fftMaxIndex);
        end
        
                           
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
end    % for PRN = satelliteList

%=========== 捕获过程结束 ==================================================
fprintf(')\n');
settings.samplingFreq = tempFreq;    % 恢复采样率
end
