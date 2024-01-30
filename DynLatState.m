% Cluster the dynamic laterality indexes intro groups

% dLI_subjects is a 113 * 1 cell, 113 is the number of subjects, each cell is a 2 * 8 * nwin matrix indicating the dynamic laterality indexes
load('dLI_subjects.mat') 

nwin = size(dLI,3);
data = dLI_subjects;
for i = 1: 113
    % change the original data dimension from 2*8*nwin to 16*nwin
    data{i,1} = reshape(data{i,1},2*8,nwin)';   
end
data = cell2mat(data);


% Set the parameters for Kmeans clustering
maxk=5;
mink=1;
rangeK=mink:maxk;
Kmeans_results=cell(size(rangeK));
distance = 'sqEuclidean';
for k=1:length(rangeK)      
    disp(['- ' num2str(rangeK(k)) ' clusters'])
    [IDX, C]=kmeans(data,rangeK(k),'Replicates',10,'MaxIter',1000,'Distance',distance,'Display','off','Options',statset('UseParallel',1));
    Kmeans_results{k}.IDX=IDX;   % Cluster time course - numeric collumn vectos
    Kmeans_results{k}.C=C;       % Cluster centroids (lateralization patterns)
end


% Identify the optimal number of clusters
% Average Silhouette Coefficient
disp('Computing average Silhouette coefficient:')
avg_sil = zeros(length(rangeK),1);
for i = 1:length(rangeK)
    eva_sil = evalclusters(data,Kmeans_results{i}.IDX,'Silhouette','Distance',distance);
    avg_sil(i) = eva_sil.CriterionValues;
    clear eva_sil;
    disp(['- K = ' num2str(rangeK(i))])
end
[~, ind_maxsil] = max(avg_sil);
disp(['- Best clustering solution according to average Silhouette coefficient: ' num2str(rangeK(ind_maxsil)) ' clusters']);

% CH index
disp('Computing CH index:')
CH = zeros(length(rangeK),1);
for i = 1:length(rangeK)
    eva_CH = evalclusters(data,Kmeans_results{i}.IDX,'CalinskiHarabasz');
    CH(i) = eva_CH.CriterionValues;
    clear eva_CH;
    disp(['- K = ' num2str(rangeK(i))])
end
[~, ind_maxCH] = max(CH);
disp(['- Best clustering solution according to CH index: ' num2str(rangeK(ind_maxCH)) ' clusters']);

% the optimal number of clusters is 2 as indicated by the above two methods
Kmeans_IDX = reshape(Kmeans_results{1,2}.IDX,[nwin,113])';


% calculate the cluster centroids (lateralization patterns) of each state
Kmeans_centrl_1 = reshape(Kmeans_results{1,2}.C(1,:),2,8);
Kmeans_centrl_2 = reshape(Kmeans_results{1,2}.C(2,:),2,8);


% calculate the time spent of each states
for i = 1:2
    Kmeans_TS(:,i) = sum(Kmeans_IDX==i,2);
end
