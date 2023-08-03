function plotAcquisition(acqResults, settings)
%Functions plots bar plot of acquisition results (acquisition metrics). No
%bars are shown for the satellites not included in the acquisition list (in
%structure SETTINGS). 
%
%plotAcquisition(acqResults)
%
%   Inputs:
%       acqResults    - Acquisition results from function acquisition.


%% Plot all results =======================================================
figure(101);

hAxes = newplot();

bar(hAxes, acqResults.peakMetric);

title (hAxes, 'Acquisition results');
xlabel(hAxes, 'PRN number (no bar - SV is not in the acquisition list)');
ylabel(hAxes, 'Acquisition Metric');

oldAxis = axis(hAxes);
axis  (hAxes, [0, 33, 0, oldAxis(4)]);
set   (hAxes, 'XMinorTick', 'on');
set   (hAxes, 'YGrid', 'on');

%% Mark acquired signals ==================================================
if settings.IF > 0
    acquiredSignals = acqResults.peakMetric .* (acqResults.carrFreq > 0);  
else
    acquiredSignals = acqResults.peakMetric .* (acqResults.carrFreq ~= 0);     % 针对0中频的情况单独画图
end

hold(hAxes, 'on');
bar (hAxes, acquiredSignals, 'FaceColor', [0 0.8 0]);
hold(hAxes, 'off');

legend(hAxes, 'Not acquired signals', 'Acquired signals');
