% ZCJ重构GNSS_SDR项目，并修正其中的若干BUG

clear; close all; 
%%
format ('compact');
format ('long', 'g');

%-------- 添加一些文件夹到工作路径中 ---------------------------------------
addpath include             % 将一些子函数放在这里
addpath geoFunctions        % 一些和地球相关的子函数
addpath acquire_zcj         % 信号捕获
addpath track_zcj           % 信号跟踪函数
addpath PVT                 % 定位
addpath deepIntegration     % 各种组合导航实现

%% 初始化 ==========================================================
disp ('Starting processing...');

settings = initSettings();
[fid, message] = fopen(settings.fileName, 'rb');

% 数据文件成功加载
if (fid > 0)
    
    % 将文件指针移动到起始位置，注意第二个参数的单位是字节
    fseek(fid, settings.skipNumberOfSamples * settings.dataFormat * settings.fileType, 'bof');

%% Acquisition ============================================================

    % 开始捕获过程，如果捕获过程没被强制跳过的话
    if ((settings.skipAcquisition == 0) || ~exist('acqResults', 'var'))
        
        % 一个CA码周期对应的采样点个数
        samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));
        
        % 读取数据用于捕获. 使用长度为20ms的数据用于精确捕获(精确捕获这一过程
        % 不是必须的，应根据实际应用进行考虑)
        data = fread(fid, 20 * samplesPerCode * settings.fileType, settings.dataType)';
        if settings.fileType == 2
            dataI = data(1:2:end);
            dataQ = data(2:2:end);
            % data  = 1 * dataI + 0 * dataQ;   % using I only 
            data  = 1 * dataI + 1j * dataQ;  % using IQ 
        end

        %-------- 开始捕获过程 ---------------------------------------------
        disp ('   Acquiring satellites...');
        
        notUsingFineFreqAcq = 1;                     % 不使用频率精确捕获
        if notUsingFineFreqAcq == 1
            acqResults = acquisition_L1CA1(data, settings);   % 该算法无法抵抗比特翻转。不使用精捕获方法。不利于弱信号捕获
        else
            acqResults = acquisition_L1CA2(data, settings);   % 该算法无法抵抗比特翻转。使用精捕获方法。不利于弱信号捕获
        end
        
        plotAcquisition(acqResults, settings);
        clear data dataI dataQ
    end

%% 将捕获得到的结果设置成通道的初始状态 ===============================
    if (any(acqResults.carrFreq))   % 如果有通道捕获成功则进行通道初始化
        channel = preRun(acqResults, settings);    % 利用捕获到的数据初始化通道结构体
        showChannelStatus(channel, settings);
    else
        % 没检测到任何卫星信号
        disp('No GNSS signals detected, signal processing finished.');
        trackResults = [];
        return;
    end

%% 信号跟踪 =========================================================
    startTime = now;
    disp (['   Tracking started at ', datestr(startTime)]);
    
    % 对于给定长度的数据(由settings.msToProcess决定)进行信号跟踪
    % 跟踪方式为逐个通道跟踪，即通道之间不是并行的
    % 尝试写出多种跟踪函数
    
    % GNSS_SDR的跟踪函数，仅使用1ms相干积分，因此可以不考虑比特同步的问题
    % [trackResults, channel] = tracking(fid, channel, settings);
    
    % 二阶码环，二阶载波环。PLL环路滤波采用书本上信号流图的写法，仅使用1ms相干积分
    % [trackResults, channel] = trackpll2nd(fid, channel, settings);  
    
    % 二阶码环，三阶载波环。PLL环路滤波采用书本上信号流图的写法，仅使用1ms相干积分
    % [trackResults, channel] = trackpll3rd(fid, channel, settings);
    
    % 二阶码环，三阶载波环。PLL环路滤波采用另外一种写法，仅使用1ms相干积分
    % [trackResults, channel] = trackpll3rd2(fid, channel, settings);
    
    % 二阶码环，一阶FLL辅助二阶PLL。仅使用1ms相干积分，仿照gnss_sdrlib这份C语言代码
    % 该函数的正确性有待验证,虽然很像是对的(doge),仅使用1ms相干积分
    % isFLL = 1, 则断开PLL变为纯PLL；否则为FLL辅助PLL
    isFLL = 0; [trackResults, channel] = trackfll1stpll2nd(fid, channel, settings, isFLL);
    
    % 二阶码环，一阶FLL辅助二阶PLL。仅使用1ms相干积分，仿照gnss_sdrlib这份C语言代码
    % 该函数的正确性有待验证,虽然很像是对的(doge),仅使用1ms相干积分
    % 完全按照信号流图的写法
    % isFLL = 1, 则断开PLL变为纯PLL；否则为FLL辅助PLL
    % [trackResults, channel] = trackfll1stpll2nd2(fid, channel, settings, isFLL);
    
    % 二阶码环，二阶FLL辅助三阶PLL。仅使用1ms相干积分，仿照gnss_sdrlib这份C语言代码
    % 该函数的正确性有待验证,虽然很像是对的(doge),仅使用1ms相干积分
    % 完全按照信号流图的写法
    % isFLL = 1, 则断开PLL变为纯PLL；否则为FLL辅助PLL
    % [trackResults, channel] = trackfll2ndpll3rd(fid, channel, settings, isFLL);
    
    % 关闭数据文件
    fclose(fid);
    
    disp(['   Tracking is over (elapsed time ', ...
                                        datestr(now - startTime, 13), ')'])     

    % 跟踪会耗费很长时间，跟踪结束后将跟踪结果存储到本地
    disp('   Saving Acq & Tracking results to file "trackingResults.mat"')
    save('trackingResults', ...
                      'trackResults', 'settings', 'acqResults', 'channel');                  

%% 定位 ============================================================
    disp('   Calculating navigation solutions...');
    [navSolutions, eph, subFrameStart, TOW] = postNavigation_zcj(trackResults, settings);
    disp('   Processing is complete for this data block');
    
%% 展示结果 ===================================================
    disp ('   Ploting results...');
    if settings.plotTracking
        plotTracking(1:settings.numberOfChannels, trackResults, settings);
    end

%    plotNavigation(navSolutions, settings);

    disp('Post processing of the signal is over.');

else
    % 数据文件加载失败
    error('Unable to read file %s: %s.', settings.fileName, message);
end 
