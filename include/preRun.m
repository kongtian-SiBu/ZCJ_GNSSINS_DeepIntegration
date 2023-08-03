function [channel] = preRun(acqResults, settings)
%Function initializes tracking channels from acquisition data. The acquired
%signals are sorted according to the signal strength. This function can be
%modified to use other satellite selection algorithms or to introduce
%acquired signal properties offsets for testing purposes.
%
%[channel] = preRun(acqResults, settings)
%
%   Inputs:
%       acqResults  - results from acquisition.
%       settings    - receiver settings
%
%   Outputs:
%       channel     - structure contains information for each channel (like
%                   properties of the tracked signal, channel status etc.). 

%% Initialize all channels ================================================
channel                 = [];   % Clear, create the structure

channel.PRN             = 0;    % PRN number of the tracked satellite
channel.acquiredFreq    = 0;    % Used as the center frequency of the NCO
channel.codePhase       = 0;    % Position of the C/A  start

channel.status          = '-';  % Mode/status of the tracking channel
                                % "-" - "off" - no signal to track
                                % "T" - Tracking state

%--- Copy initial data to all channels ------------------------------------
channel = repmat(channel, 1, settings.numberOfChannels);

%% Copy acquisition results ===============================================

%--- Sort peaks to find strongest signals, keep the peak index information
[junk, PRNindexes]          = sort(acqResults.peakMetric, 2, 'descend');

%--- Load information about each satellite --------------------------------
% Maximum number of initialized channels is number of detected signals, but
% not more as the number of channels specified in the settings.
for ii = 1:min([settings.numberOfChannels, sum(acqResults.carrFreq > 0)])
    channel(ii).PRN          = PRNindexes(ii);
    channel(ii).acquiredFreq = acqResults.carrFreq(PRNindexes(ii));
    channel(ii).codePhase    = acqResults.codePhase(PRNindexes(ii));
    
    % Set tracking into mode (there can be more modes if needed e.g. pull-in)
    channel(ii).status       = 'T';
end
