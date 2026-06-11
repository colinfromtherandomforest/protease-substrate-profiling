%% xCFH023 Filtering Analysis
%  Colin Hemez

clear; close all;

data = readtable("aa_count_annotated.csv");
dim = size(data);

% Define the referece substrate sequence
WT = 'ATGHKRSTSEGAWPQLPSGLSMMRCLHNFLTDGVPAEGAFT';
WTd = double(WT);

% Define substitutions permitted in the library
pos = (1:length(WTd))';

%% Extract sets of sequences with defined Hamming values from WT
E0 = hammingExtractor(data, WT, 0);
E1 = hammingExtractor(data, WT, 1);

hammingFiltered = [E0 ; E1];
writetable(hammingFiltered, 'aa_count_xCFH024_hammingFiltered.xlsx')

%% Load filtered table, normalize, and sort
hammingFiltered = readtable("aa_count_xCFH024_hammingFiltered.xlsx");

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
fold1  = log2(normVals(:,1) ./ normVals(:,2));
fold2  = log2(normVals(:,3) ./ normVals(:,4));
fold3 = log2(normVals(:,5) ./ normVals(:,6));
fold4 = log2(normVals(:,7) ./ normVals(:,8));
fold5 = log2(normVals(:,9) ./ normVals(:,10));
fold6 = log2(normVals(:,11) ./ normVals(:,12));

hammingOrdered = addvars(hammingOrdered, fold1, fold2, fold3, fold4, fold5, fold6);

fold1mat  = reshape(fold1,[21,length(WT)]);
fold2mat  = reshape(fold2,[21,length(WT)]);
fold3mat = reshape(fold3,[21,length(WT)]);
fold4mat = reshape(fold4,[21,length(WT)]);
fold5mat = reshape(fold5,[21,length(WT)]);
fold6mat = reshape(fold6,[21,length(WT)]);

writematrix(fold1mat, 'xCFH024_GSDMD_x3100_1_4.xlsx')
writematrix(fold2mat, 'xCFH024_GSDMD_dx3100_1_4.xlsx')
writematrix(fold3mat, 'xCFH024_GSDMD_x5139_14.xlsx')
writematrix(fold4mat, 'xCFH024_GSDMD_dx5139_14.xlsx')
writematrix(fold5mat, 'xCFH024_GSDMD_x5139_14c.xlsx')
writematrix(fold6mat, 'xCFH024_GSDMD_dx5139_14c.xlsx')