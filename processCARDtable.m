function newTbl = processCARDtable(outputTbl)
% get the percentage of each gene in the table
readlength=75;
total = sum(outputTbl.Depth .* outputTbl.Template_length .* outputTbl.Template_Coverage);
outputTbl.RelavantPercentInCARD = outputTbl.Depth .* outputTbl.Template_length .* outputTbl.Template_Coverage ./ total;
ReadsOfCARD = sum((outputTbl.Template_length .* outputTbl.Depth)/readlength); % reads identified by VFDB
PercentageOfCARD = (ReadsOfCARD ./ outputTbl.shotgunReadcount);
outputTbl.PercentageInShotgun =  outputTbl.RelavantPercentInCARD .* PercentageOfCARD;


outputTbl.mutation = repmat({''},height(outputTbl),1);
%%% find genes with mutation
p = cellfun(@(X) contains(X, 'mutation'), outputTbl.x_Template);
p2 = cellfun(@(X) contains(X, 'mutant'), outputTbl.x_Template);
p=p+p2;
outputTbl.mutation(p>0) = {'mutation'};
outputTbl.mutation(p==0) = {'noMutation'};
clear p p2
% get bacteria names and gene names
b = outputTbl.x_Template(:);
b1 = cellfun(@(X) strrep(X, 'CARD|', ''), b, 'UniformOutput', false);
b2 = cellfun(@(X) strsplit(X, '['), b1, 'UniformOutput', false);
b3 = cellfun(@(X) X{1}, b2, 'UniformOutput', false);
b3 = cellfun(@(X) strsplit(X, ' '), b3,  'UniformOutput', false);
b3 = cellfun(@(X) X{1}, b3,  'UniformOutput', false);
outputTbl.Template = b3;
for i=1:length(b3)
   bidx = strfind(b3{i}, '_');
   if length(bidx) >= 2
       b3i = b3{i};
       b3{i} = b3i(1:(bidx(end-1)-1));
   end
end
outputTbl.Accession = b3;
b4 = cellfun(@(X) X{2}, b2, 'UniformOutput', false);
outputTbl.Genome = cellfun(@(X) strrep(X, ']', ''), b4, 'UniformOutput', false);
x = regexp(outputTbl.Genome, '^\S*\s\S*\s', 'match');
% ge = regexp(outputTbl.Genome, '^\S*', 'match');
% outputTbl.genus = cellfun(@(X) X{1}, ge, 'UniformOutput', false);
for i=1:length(x)
    if isempty(x{i})
        x{i} = outputTbl.Genome{i};
    else
        x(i) = x{i};
    end
end
outputTbl.Species = cellfun(@strtrim,x, 'UniformOutput', false);
clear x x1

%%
cardDB = readtable('aro_index_5.1.1.tsv', 'FileType', 'text');
cardDB=cardDB(:, {'AROAccession' 'AROName' 'ProteinAccession' 'DNAAccession' 'DrugClass' 'ResistanceMechanism'});

%% get gene name and drug class
outputTbl.drugClass = repmat({'unknown'}, height(outputTbl),1);
outputTbl.resistGene = repmat({'unknown'}, height(outputTbl),1);
outputTbl.resistMechanism = repmat({'unknown'}, height(outputTbl),1);
[~, ia1, ib1]=intersect(outputTbl.Accession, cardDB.DNAAccession);

[~, ia2, ib2]=intersect(outputTbl.Accession, cardDB.ProteinAccession);
outputTbl.drugClass(ia1) = cardDB.DrugClass(ib1);
outputTbl.resistGene(ia1) = cardDB.AROName(ib1);
outputTbl.resistMechanism(ia1) = cardDB.ResistanceMechanism(ib1);
outputTbl.drugClass(ia2) = cardDB.DrugClass(ib2);
outputTbl.resistGene(ia2) = cardDB.AROName(ib2);
outputTbl.resistMechanism(ia2) = cardDB.ResistanceMechanism(ib2);
%%
D = unique(cardDB.DrugClass, 'stable');
DrugList = {};
for i=1:length(D)
% for i=1:2
    d=D{i};
    d=strsplit(D{i}, ';');
    DrugList = [DrugList d];  
end
DrugList = unique(DrugList);

M = zeros(height(outputTbl), length(DrugList));
for i=1:length(DrugList)
    drugIdx=contains(outputTbl.drugClass, DrugList{i});
    M(:, i)= drugIdx;
end


%%
Mtbl = array2table(M, 'VariableNames', DrugList);
Mtbl.unknown = zeros(height(Mtbl),1);
Mtbl.unknown(ismember(outputTbl.drugClass, 'unknown'))=1;

outputTbl = [outputTbl Mtbl];
outputTbl.drugClass = [];

% newTbl = outputTbl;
newTbl=sortrows(outputTbl, 'RelavantPercentInCARD', 'descend');
idx = find(ismember(newTbl.Properties.VariableNames, 'Template'));
newTbl = [newTbl(:, idx:end) newTbl(:, 1:(idx-1))];
newTbl.x_Template =[];