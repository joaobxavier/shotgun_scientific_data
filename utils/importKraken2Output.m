function tblKraken2 = importKraken2Output(fn)

% Kraken 2's standard sample report format is tab-delimited with one line per taxon. The fields of the output, from left-to-right, are as follows:
% 
% Percentage of fragments covered by the clade rooted at this taxon
% Number of fragments covered by the clade rooted at this taxon
% Number of fragments assigned directly to this taxon
% A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. Taxa that are not at any of these 10 ranks have a rank code that is formed by using the rank code of the closest ancestor rank with a number indicating the distance from that rank. E.g., "G2" is a rank code indicating a taxon is between genus and species and the grandparent taxon is at the genus rank.
% NCBI taxonomic ID number
% Indented scientific name

tblKraken2 = readtable(fn, 'NumHeaderLines', 0, 'Delimiter', '\t');
tblKraken2.Properties.VariableNames = {'percentage' 'nFragsBelow' 'nFragsThis' 'rankCode' 'ncbiTaxonomy' 'name'};

%% The arbitrary denoomination D,K just complicates everything. 
% easy fic: replace all D with K
tblKraken2.rankCode = strrep(tblKraken2.rankCode, 'D', 'K');

%% Kraken2 has these intemediate ranks defined with a letter and a number
% that makes it harder to parse. 
biologicalRankOrder = 'URKPCOFGS';
tblRanks = table(unique(tblKraken2.rankCode, 'stable'));
tblRanks.Properties.VariableNames{1} = 'rankCode';
% get the order of ranks: 
tblRanks.oneLetterRank = cellfun(@(x) {x(1)}, tblRanks.rankCode); 
% special
tblRanks.oneLetterRankLevel = cellfun(@(x) find(biologicalRankOrder==x), tblRanks.oneLetterRank);

% check if this is uncertain rank
tblRanks.isUncertainRank = cellfun(@(x) length(x), tblRanks.rankCode) > 1;
tblRanks.uncertainLevel = zeros(height(tblRanks), 1);
tblRanks.uncertainLevel(tblRanks.isUncertainRank) =...
    cellfun(@(x) str2double(x(2)), tblRanks.rankCode(tblRanks.isUncertainRank));
% get the equivalent rank
tblRanks.rankLevel = tblRanks.oneLetterRankLevel + tblRanks.uncertainLevel;
tblRanks.Properties.RowNames = tblRanks.rankCode;
tblRanks = sortrows(tblRanks, {'rankLevel', 'uncertainLevel'});

%% add a column of the ranklevel
tblKraken2.rankLevel = tblRanks.rankLevel(tblKraken2.rankCode);

%% add the taxonomy columns to the kraken2 table
for i = 1:height(tblRanks)
    tblKraken2.(tblRanks.rankCode{i}) = repmat({''}, [height(tblKraken2), 1]);
    idx = strcmp(tblKraken2.rankCode, tblRanks.rankCode{i});
    tblKraken2(idx, tblRanks.rankCode(i)) = tblKraken2.name(idx);
end

%%
h = waitbar(0, ['Parsing Kraken2 file ' fn]);
h.Children.Title.Interpreter = 'none';
for i = 2:height(tblKraken2)
    waitbar(i/height(tblKraken2),h)
    rl = tblKraken2.rankLevel(i);
    % get the row above that comes from a lower taxonomic level
    idxParent = find(tblKraken2.rankLevel(1:i-1) < rl, 1, 'last');
    % copy the taxonomy of the parent, all the way to this rankLevel
    rc = tblRanks.rankCode(tblRanks.rankLevel< rl);
    tblKraken2(i, rc) = tblKraken2(idxParent, rc);
end
close(h)

end