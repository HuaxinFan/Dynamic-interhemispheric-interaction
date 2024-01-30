function LCI_features = DynLagIdx(data,roi,win,step)
% Calculate the time-varing hemodynamic lag between two rois and the corresponding features.
% Input:
%         data: m * n array, m = number of time points, n = number of rois.
%         roi: the index of two rois with the second one as the reference. 
%         win: window length.
%         step: window step.
% Output:
%         Features of the time-varing hemodynamic lag analysis:
%         - time spent of delay state
%         - average delay
%         - time spent of lead state
%         - average lead
%         - time spent of synchronized state
%         - time spent of uncorrelated state

range = win/2;
ntime = size(data,1);
nwin = fix((ntime - win + step) / step); 
LCI = nan(nwin,1);
for t = 1:nwin
    tWIN = (t-1)*step + 1:(t-1)*step + win; 
    [LCI_sub,lags] = xcorr(data(tWIN,roi(1)), data(tWIN,roi(2)), range, 'normalized');
    [~,idx] = max(LCI_sub);
    LCI(t) = lags(idx);
end
LCInan = LCI; 
LCInan(abs(LCInan)==range) = nan; 
LCI_features = {sum(LCInan>0,1)/nwin, mean(LCInan.*(LCInan>0),1,'omitnan'), sum(LCInan<0,1)/nwin, mean(LCInan.*(LCInan<0),1,'omitnan'), sum(LCInan==0,1)/nwin, sum(isnan(LCInan),1)/nwin};
end
