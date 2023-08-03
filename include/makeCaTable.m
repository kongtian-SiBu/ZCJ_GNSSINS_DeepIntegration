function caCodesTable = makeCaTable(settings)
%Function generates CA codes for all 32 satellites based on the settings
%provided in the structure "settings". The codes are digitized at the
%sampling frequency specified in the settings structure.
%One row in the "caCodesTable" is one C/A code. The row number is the PRN
%number of the C/A code.
%
%caCodesTable = makeCaTable(settings)
%
%   Inputs:
%       settings        - receiver settings
%   Outputs:
%       caCodesTable    - an array of arrays (matrix) containing C/A codes
%                       for all satellite PRN-s


%---------- 一个CA码周期对应的采样点个数 -----------------------------------
samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));

%---------- 存储上采样后的CA码 ---------------------------------------------
caCodesTable = zeros(32, settings.acqCoIntime * samplesPerCode);
 
ts = 1 / settings.samplingFreq;   % 采样周期 [s]
tc = 1 / settings.codeFreqBasis;  % CA码码片周期, 约977.5ns
 
%=== 开始生成CA码 ...
for PRN = 1:32
    %------------ 生成CA码 ------------------------------------------------
    caCode = generateCAcode(PRN);  % 双极性
    for ii = 1 : settings.acqCoIntime - 1
        caCode = [caCode, caCode(1:1023)];
    end
  
    %======= 上采样 =======================================================
    
    %--- Make index array to read C/A code values -------------------------
    % The length of the index array depends on the sampling frequency -
    % number of samples per millisecond (because one C/A code period is one
    % millisecond).
    codeValueIndex = ceil((ts * (1 : settings.acqCoIntime * samplesPerCode)) / tc);
    
    %--- Correct the last index (due to number rounding issues) -----------
    codeValueIndex(end) = settings.acqCoIntime * 1023;
    
    %--- Make the digitized version of the C/A code -----------------------
    % The "upsampled" code is made by selecting values form the CA code
    % chip array (caCode) for the time instances of each sample.
    caCodesTable(PRN, :) = caCode(codeValueIndex);
    
end 

end