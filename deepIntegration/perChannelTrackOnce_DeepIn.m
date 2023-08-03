function  [trackans, I_P, Q_P] = perChannelTrackOnce_DeepIn(trackans, settings, fid, dppl, delta_dppl)
%% 某通道进行一次相干积分，和原本的tracking函数无区别
%
% 输入参数:
%        - trackans: 一个通道的跟踪结构体
%        - settings: 接收机相关参数
%        - fid: 中频数据文件
%        - dppl: 上一时刻惯导反馈的频率值
%        - delta_dppl: 当前惯导反馈量减去上一时刻惯导的反馈量[Hz]
%
% 杈撳嚭鍙傛暟:
%        - trackans: 跟踪结构体
%
%--------------------------------------------------------------------------
%% 鍒濆鍖栦竴浜涘彉閲?
%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

% Summation interval
PDIcode = 0.001;

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
                                    settings.dllDampingRatio, ...
                                    settings.dllGain);

%--- PLL variables --------------------------------------------------------
% Summation interval
PDIcarr = 0.001;

% Calculate filter coefficient values
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
                                    settings.pllDampingRatio, ...
                                    0.25);

% Move the starting point of processing. skipNumberOfBytes 宸茬粡绠楀湪SamplePos涓簡锛屼笉蹇呭啀娆¤绠?
fseek(fid, settings.fileType * settings.dataFormat * (trackans.SamplePos), 'bof');

%--------------------------------------------------------------------------
% Get a vector with the C/A code sampled 1x/chip
caCode = generateCAcode(trackans.PRN);
% Then make it possible to do early and late versions
caCode = [caCode(1023) caCode caCode(1)];

% define initial code frequency basis of NCO
codeFreq = trackans.codeFreq;
% define residual code phase (in chips)
remCodePhase = trackans.remCodePhase;
% define carrier frequency which is used over whole tracking period
carrFreq = trackans.carrFreq;
% carrFreqBasis = trackans.carrFreqBasis;
% define residual carrier phase
remCarrPhase = trackans.remCarrPhase;

% code tracking loop parameters
oldCodeNco   = trackans.codeNco;
oldCodeError = trackans.codeError;

% carrier/Costas loop parameters
oldCarrNco   = trackans.carrNco;
oldCarrError = trackans.carrError;

%% 寮?濮嬭窡韪?
% Find the size of a "block" or code period in whole samples
codePhaseStep = codeFreq / settings.samplingFreq;            
blksize = ceil((settings.codeLength - remCodePhase) / codePhaseStep);

trackans.recvTime = trackans.recvTime + blksize / settings.samplingFreq; % [s]

% Read in the appropriate number of samples to process this interation
[rawSignal, samplesRead] = fread(fid, settings.fileType * blksize, settings.dataType);

rawSignal = transpose(rawSignal);  % transpose vector
if settings.fileType == 2
    rawSignalI = rawSignal(1:2:end);
    rawSignalQ = rawSignal(2:2:end);
    rawSignal  = rawSignalI + 1j * rawSignalQ; 
end

% If did not read in enough samples, then could be out of data - better exit 
if (samplesRead ~= settings.fileType * blksize)
    disp('Not able to read the specified number of samples  for tracking, exiting!')
    return
end

%--------------------------------------------------------------------------
% Define index into early code vector
tcode       = (remCodePhase-earlyLateSpc) : ...
              codePhaseStep : ...
              ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
tcode2      = ceil(tcode) + 1;
earlyCode   = caCode(tcode2);

% Define index into late code vector
tcode       = (remCodePhase+earlyLateSpc) : ...
              codePhaseStep : ...
              ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
tcode2      = ceil(tcode) + 1;
lateCode    = caCode(tcode2);

% Define index into prompt code vector
tcode       = remCodePhase : ...
              codePhaseStep : ...
              ((blksize-1)*codePhaseStep+remCodePhase);
tcode2      = ceil(tcode) + 1;
promptCode  = caCode(tcode2);

remCodePhase = (tcode(blksize) + codePhaseStep) - 1023.0;

% Generate the carrier frequency to mix the signal to baseband -----------
time        = (0:blksize) ./ settings.samplingFreq;

% Get the argument to sin/cos functions
trigarg     = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
remCarrPhase = rem(trigarg(blksize+1), (2 * pi));

carr = exp(-1j* trigarg(1:blksize));
qBasebandSignal = imag(carr .* rawSignal);
iBasebandSignal = real(carr .* rawSignal);

% Now get early, late, and prompt values for each
I_E = sum(earlyCode  .* iBasebandSignal);
Q_E = sum(earlyCode  .* qBasebandSignal);
I_P = sum(promptCode .* iBasebandSignal);
Q_P = sum(promptCode .* qBasebandSignal);
I_L = sum(lateCode   .* iBasebandSignal);
Q_L = sum(lateCode   .* qBasebandSignal);

% Find PLL error and update code NCO --------------------------------------
% Implement carrier loop discriminator (phase detector)
carrError = atan(Q_P / I_P) / (2.0 * pi);

% Implement carrier loop filter and generate NCO command
carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
    (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
% oldCarrNco   = carrNco;
% oldCarrError = carrError;

% Modify carrier freq based on NCO command; something has been changed here
% carrFreq = carrFreqBasis + carrNco;
carrFreq = carrNco + settings.IF + dppl + delta_dppl / 20;  % 仅在此处进行了改动

trackans.carrFreq = carrFreq;   

% Find DLL error and update code NCO --------------------------------------
codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
                (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));
            
% Implement code loop filter and generate NCO command
codeNco = oldCodeNco + (tau2code/tau1code) * ...
    (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
% oldCodeNco   = codeNco;
% oldCodeError = codeError;

% Modify code freq based on NCO command
codeFreq = settings.codeFreqBasis - codeNco + (carrFreq - settings.IF) / 1540;   %% PLL Aided DLL correct answer!
trackans.codeFreq = codeFreq;

trackans.SamplePos = ftell(fid) / settings.dataFormat / settings.fileType;

trackans.codeError          = codeError;
trackans.codeNco            = codeNco;
trackans.carrError          = carrError;
trackans.carrNco            = carrNco;

trackans.remCodePhase       = remCodePhase;
trackans.remCarrPhase       = remCarrPhase;

trackans.numOfCoInt         = trackans.numOfCoInt + 1; 
end
