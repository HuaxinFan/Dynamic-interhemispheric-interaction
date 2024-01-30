function dLI = DynLatIdx(data,roi,net,win,step)
% Calculate the dynamic laterality indexes
% Input:
%         data: m * n array, m = number of time points, n = number of rois.
%         roi: index of the two rois.
%         net: 2 * 8 array, indicating the index of the 8 brain networks, 
%              the first line is the network in the left hemisphere, 
%              the second line is the network in the right hemisphere.
%         win: window length.
%         step: window step.
% Output:
%         Dynamic laterality indexes

ntime = size(data,1);
nwin = fix((ntime - win + step) / step);
dLI = nan(2,8,nwin);
for t = 1:nwin
    tWIN = (t-1)*step + 1:(t-1)*step + win;
    dFC = corr(data(tWIN,:));
    dFC_net = nan(2,16); 
    m = 0;
    for jj = roi    
        m = m+1;
        n = 0;
        for j = net
            n = n+1;
            idx = find(Yeo7HMATHCPMMPASEG==j); % Yeo7HMATHCPMMPASEG is a column vector, indicating the network index of an roi in the whole brain
            idx(idx==jj) = [];
            dFC_net(m,n) = mean(dFC(jj,idx),2);
        end       
    end
    for n = 1:8
        dLI(:,n,t) = dFC_net(:,n) - dFC_net(:,n+8);
    end
end
end
