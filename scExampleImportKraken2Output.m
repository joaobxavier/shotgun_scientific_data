fn = 'tempFiles/report.txt';
tblKraken2 = importKraken2Output(fn);

%% add all the fragments mapped to a genus and show the result
% ordered by the most adundant in sample
sumG = grpstats(tblKraken2, {'K' 'G' 'S'}, 'sum', 'DataVars', 'nFragsThis');
sumG.relativeAbundance = sumG.sum_nFragsThis ./ sum(sumG.sum_nFragsThis);
sumG = sortrows(sumG, 'sum_nFragsThis', 'descend');
head(sumG)
