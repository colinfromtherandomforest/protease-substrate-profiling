%% xCFH019 Filtering Analysis
%  VAMP1
%  Colin Hemez
%  2023 01 09

clear; close all;

data = readtable("aa_count_vamp1.csv");
dim = size(data);

% Define the referece substrate sequence
WT = 'LERDQKLSELDDRADAL';
WTd = double(WT);

% Define substitutions permitted in the library
pos = (1:length(WTd))';

%% Extract sets of sequences with defined Hamming values from WT
E0 = hammingExtractor(data, WT, 0);
E1 = hammingExtractor(data, WT, 1);

hammingFiltered = [E0 ; E1];
writetable(hammingFiltered, 'aa_count_xCFH019_vamp1_hammingFiltered.xlsx')

%% Load filtered table, normalize, and sort
hammingFiltered = readtable("aa_count_xCFH019_vamp1_hammingFiltered.xlsx");

dim = size(hammingFiltered);
normVals = zeros(dim(1),dim(2)-1);

% Normalization
for i = 2:dim(2)
    % normVals(:,i-1) = table2array(hammingFiltered(:,i)) ./ sum(table2array(hammingFiltered(:,i)));
    % Pseudocount
    normVals(:,i-1) = (table2array(hammingFiltered(:,i))+1) ./ sum((table2array(hammingFiltered(:,i))+1));
end

hammingFiltered = addvars(hammingFiltered, normVals);

% AA library
AA = {'*';'A';'C';'D';'E';'F';'G';'H';'I';'K';'L';'M';'N';'P';'Q';'R';'S';'T';'V';'W';'Y'};

hammingOrdered = hammingFiltered;
hammingOrdered(1:end, :) = [];

% Reorder table rows
for i = 1:length(WT)
    seqi = WT;
    for j = 1:length(AA)
        seqi(i) = AA{j};
        entry = find(strcmp(seqi,hammingFiltered.Sequence));
        Ri = hammingFiltered(entry,:);
        hammingOrdered = [hammingOrdered ; Ri];
    end
end

%% Calculate fold enrichment values
normVals = hammingOrdered.normVals;
CFH442fold  = log2(normVals(:,2) ./ normVals(:,1));
CFH443fold  = log2(normVals(:,4) ./ normVals(:,3));
CFH498fold = log2(normVals(:,6) ./ normVals(:,5));
CFH499fold = log2(normVals(:,8) ./ normVals(:,7));

hammingOrdered = addvars(hammingOrdered, CFH442fold, CFH443fold, CFH498fold, CFH499fold);

CFH442mat  = reshape(CFH442fold,[21,length(WT)]);
CFH443mat  = reshape(CFH443fold,[21,length(WT)]);
CFH498mat = reshape(CFH498fold,[21,length(WT)]);
CFH499mat = reshape(CFH499fold,[21,length(WT)]);

writematrix(CFH442mat, 'xCFH019_vamp1_CFH442.xlsx')
writematrix(CFH443mat, 'xCFH019_vamp1_CFH443.xlsx')
writematrix(CFH498mat, 'xCFH019_vamp1_CFH498.xlsx')
writematrix(CFH499mat, 'xCFH019_vamp1_CFH499.xlsx')