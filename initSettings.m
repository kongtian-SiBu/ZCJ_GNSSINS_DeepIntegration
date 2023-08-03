function settings = initSettings()
% Functions initializes and saves settings. 
%
% All settings are described inside function code.
%
% settings = initSettings()
%
%   Inputs: none
%
%   Outputs:
%       settings     - Receiver settings (a structure). 

%% Processing settings ====================================================
settings.msToProcess          = 17000;        % 总的信号处理时间长度 [ms]
settings.numberOfChannels     = 6;            % 接收机通道数 (导频和数据各一半)
settings.skipNumberOfSamples  = 0;            % 信号处理的起始点 [sample number]

%% Raw signal file name and other parameter ===============================
settings.fileName           = ...
   'E:\zcj_masterDegree_code_assemble\test_522_B1C_4MSpan.bin';
settings.dataType           = 'int16';
settings.IF                 = 0e6;      % [Hz]
settings.samplingFreq       = 5e6;     % [Hz]
settings.codeFreqBasis      = 1.023e6;  % [Hz]
settings.fileType           = 2;        % 2 for IQ; 1 for I only
settings.dataFormat         = 2;        % 2 for int16; 1 for int8
settings.codeLength         = 1023;


%% Acquisition settings ===================================================
settings.skipAcquisition    = 0;        % 0, 捕获; 1, 跳过捕获过程
settings.acqSatelliteList   = 1:32;     % 卫星列表 [PRN numbers]
settings.acqSearchBand      = 20e3;     % 频率搜索范围 [Hz]，以中频为中心，±多少Hz
settings.acqSearchStep      = 500;      % 频率搜索步长 [Hz]

settings.acqThreshold       = 2.5;      % 捕获门限，与捕获方法有关，且需要根据实际情况谨慎选择

settings.acqDownSample      = 0;        % 0, 捕获时不降采样; 1, 捕获时降采样
settings.acqDownFreq        = 3e6;      % 对于捕获过程来说，可适当降低采样率以提高捕获速度。但不能低于码率的2倍。
                                        % 为了方便，将降采样后的频率设置为原采样率的整数分之1
                                        % 精确捕获使用的是FFT方法，但是采样率必须高于中频的2倍
                                        % 降采样不影响多普勒频率的计算，但影响绝对频率的计算

% 对于CA码来说，由于未实现位同步，总相干积分时间应小于10ms
settings.acqCoIntime        = 1;        % 捕获的相干积分时间 [ms]
settings.acqnonCoIntime     = 1;        % 捕获的非相干积分次数 [times]   
if settings.acqCoIntime * settings.acqnonCoIntime >= 10 || settings.acqCoIntime <= 0
    error('Too long Integration Time or Wrong acqCoIntime Time ! ');
end  


%% Tracking loops settings ================================================
% 本工程中PLL和DLL的环路更新周期一致

% 仅考虑二阶码环 ----------------------------------------------------------------- 
settings.dllDampingRatio         = 0.7;
settings.dllNoiseBandwidth       = 2;       % [Hz]

settings.dllCorrelatorSpacing    = 0.5;     % [chips]

% 载波环, 考虑二阶和三阶的两种情况 -------------------------------------------------         
settings.pllNoiseBandwidth       = 25;       % 对于三阶及以上PLL来说，带宽必须谨慎选择以防出现不稳定现象 [Hz]
settings.pllLoopGain             = 1;        % 环路增益，在工程中是由电路设计计算得到，不能随便调整
settings.pllDampingRatio         = 0.7;

% 对于CA码来说，相干积分时间 settings.trkCoIntime 应小于20ms
settings.trkCoIntime             = 2;       % 跟踪的相干积分时间 [ms]
settings.trknonCoIntime          = 2;       % 跟踪的非相干积分次数 [times]

% FLL
settings.fllNoiseBandwidth       = 17;
settings.fllLoopGain             = 1;

%% Navigation solution settings ===========================================
settings.navSolPeriod       = 500;          % 定位周期 [ms]

settings.elevationMask      = 10;           % 截止仰角[degrees 0 - 90]

% 启用/禁用对流层校正
settings.useTropCorr        = 0;            % 0 - Off
                                            % 1 - On

%% Plot settings ==========================================================
% 是否显示信号跟踪的效果
settings.plotTracking       = 1;            % 0 - Off
                                            % 1 - On

%% Constants ==============================================================
settings.c                  = 299792458;    % 光速 [m/s]
settings.startOffset        = 80.000;       % [ms] Initial sign. travel time
