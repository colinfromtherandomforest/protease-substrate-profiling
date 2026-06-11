function T = hammingExtractor(table, reference, H)
ref = double(reference);
T = table;
T(1:end, :) = [];

dim = size(table);

for i = 1:dim(1)
    Ti = table(i,:);
    Hi = 0;
    seq = double(table.Sequence{i});
    for j = 1:length(seq)
        if seq(j) ~= ref(j)
            Hi = Hi + 1;
        end
    end

    if Hi == H
        T = [T ; Ti];
    end
end
end