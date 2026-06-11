% sCFH030 Filtering Analysis
% Colin Hemez
% 2023 03 11

clear; close all;

data = readtable("aa_count_annotated.csv");

WT = 'NLSLPTTEEFEDDAIKK';
WTd = double(WT);

% Create empty data table
filtered = data;
filtered([1:end], :) = [];
discard = filtered;
dim = size(data);
hammingvals = zeros(dim(1), 1);

for i = 1:length(data.Sequence)
    Ti = data(i,:);

    substrate = data.Sequence{i};
    subd = double(substrate);

    hamming = 0;

    for j = 1:length(WTd)
        if subd(j) ~= WTd(j)
            hamming = hamming + 1;
        end
    end

    hammingvals(i) = hamming;

    if hamming < 2
        filtered = [filtered; Ti];
    else
        discard = [discard; Ti];
    end
end

% Tabulate the fraction of discarded reads
for colindex = 2:7
    counted = sum(filtered{:,colindex});
    discarded = sum(discard{:,colindex});
    fraction = discarded / (discarded + counted);
    disp(['Column ', num2str(colindex), ': ', num2str(discarded), ' reads with Hamming distance > 1 (', num2str(fraction*100), '%)'])
end
