function [firstSubFrame, activeChnList, firstSubFrameSampleNum] = findPreambles(trackResults, settings)
% findPreambles finds the first preamble occurrence in the bit stream of
% each channel. The preamble is verified by check of the spacing between
% preambles (6sec) and parity checking of the first two words in a
% subframe. At the same time function returns list of channels, that are in
% tracking state and with valid preambles in the nav data stream.
%
%[firstSubFrame, activeChnList] = findPreambles(trackResults, settings)
%
%   Inputs:
%       trackResults    - output from the tracking function
%       settings        - Receiver settings.
%
%   Outputs:
%       firstSubframe   - the array contains positions of the first
%                       preamble in each channel. The position is ms count 
%                       since start of tracking. Corresponding value will
%                       be set to 0 if no valid preambles were detected in
%                       the channel.
%       activeChnList   - list of channels containing valid preambles

%--------------------------------------------------------------------------

% 帧同步不一定要从一开始就进行。由于跟踪环收敛较慢，刚开始的几个数据比特误码率
% 较高。因此可以延后一段时间再进行帧同步。不过此处忽略该问题。
searchStartOffset = 0;

% 用于存储帧同步结果
firstSubFrame = zeros(1, settings.numberOfChannels);

% 帧同步码
preamble_bits = [1 -1 -1 -1 1 -1 1 1];

% 将帧同步码上采样。由于本工程仅采用1ms相干积分，因此一个数据比特包含20周CA码
preamble_ms = kron(preamble_bits, ones(1, 20));

% 去除未能实现信号跟踪的通道
activeChnList = find([trackResults.status] ~= '-');

% 我自己定义的变量，不一定能用上
firstSubFrameSampleNum = zeros(1, settings.numberOfChannels);

% 开始帧同步过程
for channelNr = activeChnList

%% 将数据比特和帧同步码进行互相关 ================================
    % 读取数据比特
    bits = trackResults(channelNr).I_P(1 + searchStartOffset : end);   % 从跟踪环的I_P处获得

    % 二值化 
    bits(bits > 0)  =  1;
    bits(bits <= 0) = -1;

    % 互相关
    tlmXcorrResult = xcorr(bits, preamble_ms);   % xcorr between data bits stream and preamble bits

%% 寻找可能的帧头位置 ===============================================
    clear index
    clear index2

    xcorrLength = (length(tlmXcorrResult) +  1) /2;

    %--- 寻找可能的帧头位置 ------------------------
    index = find(...
        abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) > 153)' + ...
        searchStartOffset;   % 如果存在帧同步码，则互相关值应为160。考虑到可能存在误码，因此将该数值设置为153

%% 分析每个出现帧同步码的位置 ========================================
    % 寻找首次出现帧头的位置
    for i = 1 : size(index) 

        %--- Find distances in time between this occurrence and the rest of
        %preambles like patterns. If the distance is 6000 milliseconds (one
        %subframe), the do further verifications by validating the parities
        %of two GPS words      
        index2 = index - index(i);  % each matched pattern has 6000ms interval

        if (~isempty(find(index2 == 6000, 1))) % 帧头每6000ms出现一次

            %=== Re-read bit vales for preamble verification ==============
            % Preamble occurrence is verified by checking the parity of
            % the first two words in the subframe. Now it is assumed that
            % bit boundaries a known. Therefore the bit values over 20ms are
            % combined to increase receiver performance for noisy signals.
            % in Total 62 bits mast be read :
            % 2 bits from previous subframe are needed for parity checking;
            % 60 bits for the first two 30bit words (TLM and HOW words).
            % The index is pointing at the start of TLM word.
            bits = trackResults(channelNr).I_P(index(i)-40 : ...
                                               index(i) + 20 * 60 -1)'; % 数值的设置与信号体制和相干积分时间有关

            %--- Combine the 20 values of each bit ------------------------
            bits = reshape(bits, 20, (size(bits, 1) / 20));
            bits = sum(bits);

            % Now threshold and make it -1 and +1 
            bits(bits > 0)  = 1;
            bits(bits <= 0) = -1;

            %--- Check the parity of the TLM and HOW words ----------------
            % 能否通过奇偶校验
            if (navPartyChk(bits(1:32)) ~= 0) && ...
               (navPartyChk(bits(31:62)) ~= 0)
                % Parity was OK. Record the preamble start position. Skip
                % the rest of preamble pattern checking for this channel
                % and process next channel. 
                
                firstSubFrame(channelNr) = index(i);
                break;    
            end % if parity is OK ...
            
        end % if (~isempty(find(index2 == 6000)))
    end % for i = 1:size(index)

    % Exclude channel from the active channel list if no valid preamble was
    % detected
    if firstSubFrame(channelNr) == 0
        
        % Exclude channel from further processing. It does not contain any
        % valid preamble and therefore nothing more can be done for it.
        activeChnList = setdiff(activeChnList, channelNr);

        disp(['Could not find valid preambles in channel ', ...
                                                  num2str(channelNr),'!']);
    else    
        firstSubFrameSampleNum(channelNr) = trackResults(channelNr).absoluteSample(firstSubFrame(channelNr));
    end
    
end % for channelNr = activeChnList

end