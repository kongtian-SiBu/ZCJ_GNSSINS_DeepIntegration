function showChannelStatus(channel, settings)
%Prints the status of all channels in a table.
%
%showChannelStatus(channel, settings)
%
%   Inputs:
%       channel     - data for each channel. It is used to initialize and
%                   at the processing of the signal (tracking part).
%       settings    - receiver settings


fprintf('\n*=========*=====*===============*===========*=============*========*\n');
fprintf(  '| Channel | PRN |   Frequency   |  Doppler  | Code Offset | Status |\n');
fprintf(  '*=========*=====*===============*===========*=============*========*\n');

for channelNr = 1 : settings.numberOfChannels
    if (channel(channelNr).status ~= '-')
        fprintf('|      %2d | %3d |  %2.5e |   %5.0f   |    %6d   |     %1s  |\n', ...
                channelNr, ...
                channel(channelNr).PRN, ...
                channel(channelNr).acquiredFreq, ...
                channel(channelNr).acquiredFreq - settings.IF, ...
                channel(channelNr).codePhase, ...
                channel(channelNr).status);
    else
        fprintf('|      %2d | --- |  ------------ |   -----   |    ------   |   Off  |\n', ...
                channelNr);
    end
end

fprintf('*=========*=====*===============*===========*=============*========*\n\n');
