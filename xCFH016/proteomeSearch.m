% Read proteome FASTA file into memory
disp('Loading proteome file...')
disp(' ')

proteomefile = 'UP000005640_SwissProt.fasta';
proteome = fastaread(proteomefile);
dim = size(proteome);
n_protein = dim(1);

% Load enrichment data file
disp('Loading enrichment data file...')
disp(' ')

datafile = 'xCFH016_sCFH035_activePositive.mat';
load(datafile)

% Remove new line characters from protein sequences
disp('Removing new line characters from proteome sequences...')
disp(' ')

remove = '\';
for i = 1:n_protein
    seq = proteome(i).Sequence;
    seq = seq(~ismember(seq, remove));
    proteome(i).Sequence = seq;
end

% Search for sequences and add them into a list of hits
disp('Searching for occurrences of enriched sequences within proteome...')
disp(' ')

AA = 'ACDEFGHIKLMNPQRSTVWY';
motif_start = 10; % Start position of variable residues within substrate
motif_end = 15; % End position of variable residues within subatrate
pos_fixed = [2 5]; % Positions of fixed residues within variable substrate
cellind = ones(1,length(pos_fixed));
cellind = cellind.*length(AA);

match_count = 0;

% Data storage variables
proteinID = {};
querySequence = {};
queryStartPosition = [];
libraryHitSequence = {};
libraryHitSequenceReadsAbx = [];
libraryHitSequenceReadsNoAbx = [];
libraryHitSequenceReadsInactiveAbx = [];
libraryHitSequenceReadsInactiveNoAbx = [];

for hit = 1:length(A.Sequence)
    % Generate a set of searchstrings
    Si = A.Sequence{hit};
    motif = Si(motif_start:motif_end);
    motif_ij = motif;
    motif_search = cell(cellind);

    for i = 1:length(AA)
        fixi = AA(i);
        for j = 1:length(AA)
            fixj = AA(j);
            AAij = {fixi,fixj};
            for ij = 1:length(pos_fixed)
                motif_ij(pos_fixed(ij)) = AAij{ij};
            end
            motif_search{i,j} = motif_ij;
        end
    end

    motif_search = motif_search(:);

    % Search for each motif within each proteome sequence
    for prot = 1:n_protein
        prot_seq = proteome(prot).Sequence;
        for query = 1:length(motif_search)
            query_motif = motif_search{query};

            motif_ind = strfind(prot_seq,query_motif);

            if ~isempty(motif_ind)
                match_count = match_count + length(motif_ind);
                for dataindex = 1:length(motif_ind)
                    % Record data
                    proteinID = [proteinID ; proteome(prot).Header];
                    querySequence = [querySequence ; query_motif];
                    queryStartPosition = [queryStartPosition ; motif_ind(dataindex)];
                    libraryHitSequence = [libraryHitSequence ; Si];
                    libraryHitSequenceReadsAbx = [libraryHitSequenceReadsAbx ; A.activeReadsAbx(hit)];
                    libraryHitSequenceReadsNoAbx = [libraryHitSequenceReadsNoAbx ; A.activeReadsNoAbx(hit)];
                    libraryHitSequenceReadsInactiveAbx = [libraryHitSequenceReadsInactiveAbx ; A.inactiveReadsAbx(hit)];
                    libraryHitSequenceReadsInactiveNoAbx = [libraryHitSequenceReadsInactiveNoAbx ; A.inactiveReadsNoAbx(hit)];
                end
            end
        end
    end

    if mod(hit, 10) == 0
        disp(['Searched ',num2str(hit),' of ',num2str(length(A.Sequence)),' hits; Found ',num2str(match_count),' occurrences of motifs'])
    end
end
disp(['Searched ',num2str(hit),' of ',num2str(length(A.Sequence)),' hits'])

libraryHitActiveEnrichment = log2(libraryHitSequenceReadsAbx./libraryHitSequenceReadsNoAbx);
libraryHitInctiveEnrichment = log2(libraryHitSequenceReadsInactiveAbx./libraryHitSequenceReadsInactiveNoAbx);

M = table(proteinID,...
          querySequence,...
          queryStartPosition,...
          libraryHitSequence,...
          libraryHitSequenceReadsAbx,...
          libraryHitSequenceReadsNoAbx,...
          libraryHitSequenceReadsInactiveAbx,...
          libraryHitSequenceReadsInactiveNoAbx,...
          libraryHitActiveEnrichment,...
          libraryHitInctiveEnrichment);

n_uniqueProteins = length(unique(proteinID));

% Save table
writetable(M,'proteomeSearch_xCFH016_sCFH035.csv')
