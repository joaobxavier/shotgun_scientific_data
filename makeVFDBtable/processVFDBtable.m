function newTbl = processVFDBtable(outputTbl)
% get the percentage of each gene in the table
readlength=75;
total = sum(outputTbl.Depth .* outputTbl.Template_length .* outputTbl.Template_Coverage);
outputTbl.RelavantPercentInVF = outputTbl.Depth .* outputTbl.Template_length .* outputTbl.Template_Coverage ./ total;
ReadsOfVFDB = sum((outputTbl.Template_length .* outputTbl.Depth)/readlength); % reads identified by VFDB
PercentageOfVFDB = (ReadsOfVFDB ./ outputTbl.shotgunReadcount);
outputTbl.PercentageInShotgun =  outputTbl.RelavantPercentInVF .* PercentageOfVFDB;

b = outputTbl.x_Template(:);
b1 = cellfun(@(X) strrep(X, 'VFDB|', ''), b, 'UniformOutput', false);
b12 = regexp(b1,'\s','split','once');
outputTbl.Template = cellfun(@(X) X{1}, b12, 'UniformOutput', false);
b12 = cellfun(@(X) X{2}, b12, 'UniformOutput', false);
b13 =cellfun(@(X) strsplit(X, ' ['), b12, 'UniformOutput', false);
outputTbl.Function = cellfun(@(X) X{1}, b13, 'UniformOutput', false);
% b2 = cellfun(@(X) strsplit(X, '['), b1, 'UniformOutput', false);
b3 = cellfun(@(X) X{2}, b13, 'UniformOutput', false);

outputTbl.Genome = cellfun(@(X) strrep(X, ']', ''), b3,  'UniformOutput', false);
% x = regexp(outputTbl.Genome, '^\S*\s\S*\s', 'match');

newTbl=sortrows(outputTbl, 'RelavantPercentInVF', 'descend');
idx = find(ismember(newTbl.Properties.VariableNames, 'Template'));
newTbl = [newTbl(:, idx:end) newTbl(:, 1:(idx-1))];
newTbl.x_Template =[];
