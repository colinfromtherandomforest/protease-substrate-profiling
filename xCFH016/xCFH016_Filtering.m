%% xCFH016 Filtering Analysis
%  Colin Hemez
%  2023 05 25

clear; close all;

data = readtable("aa_count.csv");
dim = size(data);

% Sort column names alphabetically
sortedNames = sort(data.Properties.VariableNames(2:end));
data = [data(:,1) data(:,sortedNames)];

% Define the referece substrate sequence
WT = 'NLSLPTTEEFEDDAIKK';
WTd = double(WT);

% Define substitutions permitted in the library
pos = (1:17)';
oligo =  {'N';
          'L';
          'S';
          'L';
          'P';
          'T';
          'T';
          'E';
          'E';
          'CFLMWY*';
          'E';
          'ACDEGSY*';
          'ACDEGFLMSVY*';
          'A';
          'ACDEGFHIKLMNPQRSTVWY*';
          'K';
          'K'};
active = {'N';
          'L';
          'S';
          'L';
          'P';
          'T';
          'T';
          'E';
          'E';
          'FLMWY';
          'E';
          'EACDGS';
          'ACDEGMSV';
          'A';
          'CFHILMPQTVWY';
          'K';
          'K'};

%% Filter for sequences covered within oligo pool
disp('Filtering sequences covered in oligo pool...')
oligoFiltered = data;
oligoFiltered(1:end, :) = [];

for i = 1:length(data.Sequence)
    Ri = data(i,:);
    substrate = double(data.Sequence{i});
    score = 1;
    for j = 1:length(substrate)
        if sum(double(oligo{j}) == substrate(j)) ~= 1
            score = score + 1;
        end
    end

    if score == 1
        oligoFiltered = [oligoFiltered; Ri];
    end

    if mod(i,10000) == 0
        disp(['Entries processed: ',num2str(i)])
    end

end
disp(['Sequences found: ',num2str(length(oligoFiltered.Sequence))])
disp(' ')

%% Calculate fold enrichment scores
% Assumes that +/- Abx conditions are adjacent columns (+Abx first)
disp('Calculating fold enrichment scores...')
disp(' ')

Enrichment = [];

num_conditions = (dim(2) - 1)/2;
for i = 1:num_conditions
    plus_abx = oligoFiltered(:,i*2);
    minus_abx = oligoFiltered(:,i*2+1);
    % logfold = log2(table2array(plus_abx) ./ table2array(minus_abx));
    logfold = log2((table2array(plus_abx)./sum(table2array(plus_abx))) ./ (table2array(minus_abx)./sum(table2array(minus_abx))));
    Enrichment = [Enrichment logfold];
end

oligoFiltered = addvars(oligoFiltered, Enrichment);

%% Filter for sequences covered within activity library
disp('Filtering sequences covered in activity library...')
activeFiltered = oligoFiltered;
activeFiltered(1:end, :) = [];
oligoNotActive = activeFiltered;

for i = 1:length(oligoFiltered.Sequence)
    Ri = oligoFiltered(i,:);
    substrate = double(oligoFiltered.Sequence{i});
    score = 1;
    for j = 1:length(substrate)
        if sum(double(active{j}) == substrate(j)) ~= 1
            score = score + 1;
        end
    end

    if score == 1
        activeFiltered = [activeFiltered; Ri];
    else
        oligoNotActive = [oligoNotActive; Ri];
    end
end
disp(['Sequences found: ',num2str(length(activeFiltered.Sequence))])
disp(' ')

%% Display fold enrichment for replicates
active_1 = activeFiltered.Enrichment(:,1);
active_2 = activeFiltered.Enrichment(:,3);
oligo_1 = oligoNotActive.Enrichment(:,1);
oligo_2 = oligoNotActive.Enrichment(:,3);

figure
plot(oligo_1, oligo_2, '.', active_1, active_2, '.')
legend('Oligo library','Active library','location','northwest')
grid on
xlabel('Replicate 1')
ylabel('Replicate 2')
title('log_2(fold enrichment) by replicate, ACTIVE PROTEASE')

dactive_1 = activeFiltered.Enrichment(:,2);
dactive_2 = activeFiltered.Enrichment(:,4);
doligo_1 = oligoNotActive.Enrichment(:,2);
doligo_2 = oligoNotActive.Enrichment(:,4);

figure
plot(doligo_1, doligo_2, '.', dactive_1, dactive_2, '.')
legend('Oligo library','Active library','location','northwest')
grid on
xlabel('Replicate 1')
ylabel('Replicate 2')
title('log_2(fold enrichment) by replicate, INACTIVE PROTEASE')

%% Filter sequences based on methionine content
activeMet = activeFiltered;
activeMet(1:end,:) = [];
activeNoMet = activeMet;
notActiveNoMet = activeMet;
notActiveMet = activeMet;

for i = 1:length(activeFiltered.Sequence)
    Ri = activeFiltered(i,:);
    substrate = double(activeFiltered.Sequence{i});
    if sum(substrate == double('M')) == 0
        activeNoMet = [activeNoMet ; Ri];
    else
        activeMet = [activeMet ; Ri];
    end
end

for i = 1:length(oligoNotActive.Sequence)
    Ri = oligoNotActive(i,:);
    substrate = double(oligoNotActive.Sequence{i});
    if sum(substrate == double('M')) == 0
        notActiveNoMet = [notActiveNoMet ; Ri];
    else
        notActiveMet = [notActiveMet ; Ri];
    end
end

activeWithMetTable = activeMet.Enrichment;
activeNoMetTable = activeNoMet.Enrichment;

notActiveWithMetTable = notActiveMet.Enrichment;
notActiveNoMetTable = notActiveNoMet.Enrichment;


figure
plot(activeNoMetTable(:,2),activeNoMetTable(:,1),'.',...
     activeWithMetTable(:,2),activeWithMetTable(:,1),'.');
grid on
legend('Without Methionine','With Methionine','location','northwest')
xlabel('dX(3015)8')
ylabel('X(3015)8')

figure
plot(notActiveNoMetTable(:,2),notActiveNoMetTable(:,1),'.',...
     notActiveWithMetTable(:,2),notActiveWithMetTable(:,1),'.');
grid on
legend('Without Methionine','With Methionine','location','northwest')
xlabel('dX(3015)8')
ylabel('X(3015)8')

writematrix(activeNoMetTable(:,1:2),'activeMet.csv')
writematrix(activeWithMetTable(:,1:2),'activeWithMet.csv')
writematrix(notActiveNoMetTable(:,1:2),'notActiveNoMet.csv')
writematrix(notActiveWithMetTable(:,1:2),'notActiveMet.csv')

%% Filter libraries for positive hits
activePositive = activeFiltered;
activePositive(1:end,:) = [];
activeNegative = activePositive;

notActivePositive = oligoNotActive;
notActivePositive(1:end,:) = [];
notActiveNegative = notActivePositive;

activeThreshold = log2(270/131) - 2; % 2 logs below IL1B enrichment
inactiveThreshold = -3;

for i = 1:length(activeFiltered.Sequence)
    Ri = activeFiltered(i,:);
    if Ri.Enrichment(1,1) >= activeThreshold
        if Ri.Enrichment(1,2) <= inactiveThreshold
            activePositive = [activePositive ; Ri];
        else
            activeNegative = [activeNegative ; Ri];
        end
    else
        activeNegative = [activeNegative ; Ri];
    end
end

for i = 1:length(oligoNotActive.Sequence)
    Ri = oligoNotActive(i,:);
    if Ri.Enrichment(1,1) >= activeThreshold
        if Ri.Enrichment(1,2) <= inactiveThreshold
            notActivePositive = [notActivePositive ; Ri];
        else
            notActiveNegative = [notActiveNegative ; Ri];
        end
    else
        notActiveNegative = [notActiveNegative ; Ri];
    end
end

writematrix(activePositive.Enrichment(:,1:2), 'sCFH035_activePositive.csv')
writematrix(activeNegative.Enrichment(:,1:2), 'sCFH035_activeNegative.csv')
writematrix(notActivePositive.Enrichment(:,1:2), 'sCFH035_notActivePositive.csv')
writematrix(notActiveNegative.Enrichment(:,1:2), 'sCFH035_notActiveNegative.csv')

writecell(activePositive.Sequence, 'sCFH035_activePositive.txt')
writecell(notActivePositive.Sequence, 'sCFH035_notActivePositive.txt')

% Create table of positive hits to pass into proteomeSearch
A = table(activePositive.Sequence,...
          activePositive.Enrichment(:,1),...
          activePositive.Enrichment(:,2),...
          activePositive.xCFH016_1,...
          activePositive.xCFH016_2,...
          activePositive.xCFH016_3,...
          activePositive.xCFH016_4,...
          'VariableNames',...
          ["Sequence","activeEnrichment","inactiveEnrichment",...
           "activeReadsAbx","activeReadsNoAbx","inactiveReadsAbx","inactiveReadsNoAbx"]);

% Save positive hits as csv file and matlab matrix file
writetable(A, 'xCFH016_sCFH035_activePositive.csv')
save('xCFH016_sCFH035_activePositive.mat','A');

