%% xCFH018 Filtering Analysis
%  Colin Hemez
%  2023 05 25

clear; close all;

data = readtable("aa_count_xCFH018.xlsx");
dim = size(data);

% Define the referece substrate sequence
WT = 'TSNRRLQQTQAQVEEVVDIIRVNVDKVLERDQKLSELDDRADALQAGASQFESSAAKLKR';
WTd = double(WT);

% Define substitutions permitted in the library
pos = (1:length(WTd))';

%% Extract sets of sequences with defined Hamming values from WT
E0 = hammingExtractor(data, WT, 0);
E1 = hammingExtractor(data, WT, 1);

hammingFiltered = [E0 ; E1];
writetable(hammingFiltered, 'aa_count_xCFH018_hammingFiltered.xlsx')

%% Load filtered table, normalize, and sort
hammingFiltered = readtable("aa_count_xCFH018_hammingFiltered.xlsx");

dim = size(hammingFiltered);
normVals = zeros(dim(1),dim(2)-1);

% Normalization
for i = 2:dim(2)
    normVals(:,i-1) = table2array(hammingFiltered(:,i)) ./ sum(table2array(hammingFiltered(:,i)));
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
Xfold  = log2(normVals(:,1) ./ normVals(:,2));
Ffold  = log2(normVals(:,3) ./ normVals(:,4));
dXfold = log2(normVals(:,5) ./ normVals(:,6));
dFfold = log2(normVals(:,7) ./ normVals(:,8));

hammingOrdered = addvars(hammingOrdered, Xfold, Ffold, dXfold, dFfold);

Xmat  = reshape(Xfold,[21,60]);
Fmat  = reshape(Ffold,[21,60]);
dXmat = reshape(dXfold,[21,60]);
dFmat = reshape(dFfold,[21,60]);

writematrix(Xmat, 'xCFH018_Xmat.xlsx')
writematrix(Fmat, 'xCFH018_Fmat.xlsx')
writematrix(dXmat, 'xCFH018_dXmat.xlsx')
writematrix(dFmat, 'xCFH018_dFmat.xlsx')