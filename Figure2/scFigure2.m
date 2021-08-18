%% 
close all; 
clear;
% read the sample table
tblASVsamples = readtable('../deidentified_data_tables/samples/tblASVsamples.csv', ...
                        'Format', '%s%s%d%s%s%s%d%s');
% pick the patient with most shotgun sequenced samples
idxHasShotun = cellfun(@(x) ~isempty(x), tblASVsamples.AccessionShotgun);
tblASVsamples = tblASVsamples(idxHasShotun, :);
readcountTbl = readtable('../metagenome_data/tblShotgunReadcounts.csv','Format', '%s%f');
tblASVsamples=join(tblASVsamples, readcountTbl, 'Keys', 'SampleID');
rewriteShotgunAbundances = 0; % will rewrite everything if set to 1

%% plot the timeline for this patient
addpath('../utils');
data_path = '../deidentified_data_tables/'; % path to data
% get the 16S sample data
opts = detectImportOptions(strcat(data_path, 'samples/tblASVsamples.csv'));
opts = setvartype(opts,{'PatientID'},'categorical');
tblsamples = readtable(strcat(data_path, 'samples/tblASVsamples.csv'),opts);

tblsamples = tblsamples(ismember(tblsamples.SampleID, tblASVsamples.SampleID), :);
%tblsamples = tblsamples((tblsamples.PatientID==PatientID2Plot & ismember(tblsamples.SampleID, tblASVsamples.SampleID(idxHasShotun))), :);
tblsamples = sortrows(tblsamples, 'Timepoint'); % sort rows by time point of samples

% load counts table of 16S sequencing
tblcounts_16S = readtable(strcat(data_path, 'counts/tblcounts_asv_melt.csv'));
tblcounts_16S = tblcounts_16S(contains(tblcounts_16S.SampleID, tblsamples.SampleID), :);

% unstack counts table and normalize ASV counts to relative abundance
tblcounts_16S = unstack(tblcounts_16S, 'Count', 'ASV');
counts_matrix_16S = tblcounts_16S{:, 2:end}; % the first column is "SampleID"
counts_matrix_16S(isnan(counts_matrix_16S)) = 0; % missing count value is filled with 0
counts_matrix_16S = counts_matrix_16S ./ sum(counts_matrix_16S, 2); % convert to relative abundance
tblcounts_16S{:, 2:end} = counts_matrix_16S;
tblcounts_16S = innerjoin(tblsamples(:, {'SampleID', 'Timepoint', 'DayRelativeToNearestHCT'}), tblcounts_16S);
tblcounts_16S = sortrows(tblcounts_16S, 'Timepoint'); % sort rows by time point of samples

% load taxonomy table
tbltaxonomy = readtable(strcat(data_path,'taxonomy/tblASVtaxonomy_silva132_v4v5_filter.csv'));
tbltaxonomy = tbltaxonomy(ismember(tbltaxonomy.ASV,tblcounts_16S.Properties.VariableNames(4:end)), :);


% %#########################################################################
 
% plot samples

% the first 3 columns of tblcounts are SampleID, Timepoint and DayRelativeToNearestHCT
abundance_matrix_16S = tblcounts_16S{:, 4:end};

% calculate the cumulative sum of taxa with same color_order
% unique_color_order should be automatically sorted
[unique_color_order_16S,ia,ic] = unique(tbltaxonomy.ColorOrder);
uni_color_hex_16S = tbltaxonomy.HexColor(ia);

color_grouped_abundance_16S = zeros(size(abundance_matrix_16S,1), length(unique_color_order_16S));
for k = 1:length(unique_color_order_16S)
    currsum = sum(abundance_matrix_16S(:,ic==k),2);
    color_grouped_abundance_16S(:,k) = currsum;
end

ASV_cmap_16S = hex2rgb(uni_color_hex_16S);

recordkraken2file = 'samplesWithoutKraken2.txt';
if isfile(recordkraken2file)
    delete(recordkraken2file)
end

%% make the sample plot to shotgun
addpath('../');
if rewriteShotgunAbundances ==1
    tblShotgunAbundances = [];
    S=zeros(height(tblsamples),1);
    S2=zeros(height(tblsamples),1);
    for i = 1:height(tblsamples)
    % for i=1:5
        s = tblsamples.SampleID{i};
        fn  = sprintf('../PATRIC_output/kraken2/.%s_kraken2/report.txt', s);
        srr = tblASVsamples.AccessionShotgun{strcmp(tblASVsamples.SampleID, s)};
        % check if Kraken2 file exist at all
        if isfile(fn)
            gz  = sprintf('../PATRIC_output/kraken2/.%s_kraken2/output.txt.gz', s);
            if isfile(gz)
                system(['rm ' gz]); % remove the gz file because it takes too much space and we do not need it
            end
            % check if matlab table has been precomputed
            if isfile(sprintf('../PATRIC_output/kraken2/.%s_kraken2/tblKraken2.mat', s))
                fprintf('Sample %s has been processed by matlab. Loading...\n', s)
                load(sprintf('../PATRIC_output/kraken2/.%s_kraken2/tblKraken2.mat', s));
            else
                fprintf('Sample %s has not been processed by matlab yet.\n', s)
    %             fprintf('Doing it now...\n')
                srr = tblASVsamples.AccessionShotgun{strcmp(tblASVsamples.SampleID, s)};
                log = sprintf('../PATRIC_output/kraken2/%s_kraken2', s);
                %%%% check if the srr number can be found in the html file
                [status,result] = system(['cat ' log ' | grep ' srr]);
                if status == 1 || isempty(result)
    %                 error('Sample %s not confirmed by %s', s, log)
                    fprintf('Sample %s not confirmed by %s', s, log)
                    fileID = fopen(recordkraken2file,'a');
                    fprintf(fileID,'Sample %s not confirmed by %s, double check SRR number\n', s, log);
                    fclose(fileID);
                else
                    fprintf('Sample %s confirmed by %s\n', s, log)
                    % load the table for this file
                    fprintf('Doing it now...\n')
                    tblKraken2 = importKraken2Output(fn);
                    save(sprintf('../PATRIC_output/kraken2/.%s_kraken2/tblKraken2.mat',...
                        s), 'tblKraken2');
                end

            end
            % compute relative abundances
            % add all the fragments mapped to a genus and show the result
            % ordered by the most adundant in sample
            taxa = {'K' 'P' 'C' 'O' 'F' 'G'};
            taxaFull = {'Kingdom' 'Phylum' 'Class' 'Order' 'Family' 'Genus'};
            tblKrakenBugs = tblKraken2(ismember(tblKraken2.K, {'Bacteria' 'Archaea'}), :);
            sumG = grpstats(tblKrakenBugs, taxa, 'sum', 'DataVars', 'nFragsThis');
            sumG{:, {['s' s]}} = sumG.sum_nFragsThis ./ sum(sumG.sum_nFragsThis);
            % match shotgun taxa to best 16S
            for j = 2:length(taxa)
                t = sumG{:, taxa{j}};
                [~ , loc] = ismember(t, tbltaxonomy{:, taxaFull{j}});
                idx = loc>0;
                sumG.ColorOrder(idx) = tbltaxonomy.ColorOrder(loc(idx));
                sumG.HexColor(idx) = tbltaxonomy.HexColor(loc(idx));
%                 sumG.hit16S(idx) = tbltaxonomy{loc(idx), taxaFull{j}};
            end
            % if there's no match
            sumG.HexColor(sumG.ColorOrder == 0) = {'#000000'};
            % join the tables
            if isempty(tblShotgunAbundances)
                tblShotgunAbundances = sumG(:, [taxa {'ColorOrder'} {'HexColor'}  {['s' s]}]);
            else
                tblShotgunAbundances = outerjoin(tblShotgunAbundances, sumG(:, [taxa {'ColorOrder'} {'HexColor'} {['s' s]}]),...
                                        'MergeKeys',true);
            end
            S(i)=size(tblShotgunAbundances,1);
            S2(i)=size(sumG,1);
        else
            warning('Sample %s has no Kraken2 output yet...\n', s)
            fprintf('>>>>>> Go to PATRIC and request processing %s\n', s);
            fileID = fopen(recordkraken2file,'a');
            fprintf(fileID,'Sample %s has no Kraken2 output yet...>>>>>> Go PATRIC and request processing %s\n', s, s);
            fclose(fileID);
        end
        clear sumG
    end
    % the first 3 columns of tblcounts are SampleID, Timepoint and DayRelativeToNearestHCT
    %%%%%% only keep bacteria and Archea in abundance table
    tblShotgunAbundances=tblShotgunAbundances(ismember(tblShotgunAbundances.K, {'Archaea' 'Bacteria'}), :);
    abundance_matrix_shotgun = tblShotgunAbundances{:, 10:end};
    abundance_matrix_shotgun(isnan(abundance_matrix_shotgun))=0;
    tblShotgunAbundances{:, 10:end}=abundance_matrix_shotgun ;
    tblShotgunAbundances.Properties.VariableNames(1:6)={'Kindom', 'Phylum', 'Class', 'Order', 'Family', 'Genus'};
    writetable(tblShotgunAbundances, '../metagenome_data/tblShotgunRelAbundances.csv');
    abundance_matrix_shotgun =abundance_matrix_shotgun';
    figure(2)
    plot(S,'ro-')
    hold on
    plot(S2, 'bx-')
else
    tblShotgunAbundances=readtable('../metagenome_data/tblShotgunRelAbundances.csv');
    %%%%%% only keep bacteria and Archea in abundance table
    tblShotgunAbundances=tblShotgunAbundances(ismember(tblShotgunAbundances.Kindom, {'Archaea' 'Bacteria'}), :);
    abundance_matrix_shotgun = tblShotgunAbundances{:, 10:end};
    abundance_matrix_shotgun =abundance_matrix_shotgun';
end

%%
% figure(1)
% subplot(2, 1, 2)
% % calculate the cumulative sum of taxa with same color_order
% % unique_color_order should be automatically sorted
[unique_color_order_shotgun,ia,ic] = unique(tblShotgunAbundances.ColorOrder);
uni_color_hex_shotgun = tblShotgunAbundances.HexColor(ia);

color_grouped_abundance_shotgun = zeros(size(abundance_matrix_shotgun,1), length(unique_color_order_shotgun));
for k = 1:length(unique_color_order_shotgun)
    currsum = sum(abundance_matrix_shotgun(:,ic==k),2);
    color_grouped_abundance_shotgun(:,k) = currsum;
end

%% plot stacked bar of a given patient
sampleTblwShotgun = tblASVsamples(cellfun(@(X) ~isempty(X), tblASVsamples.AccessionShotgun),:);
sampleTblwShotgunCounts = groupcounts(sampleTblwShotgun, 'PatientID')
%%%%%% choose the patient with the most samples
patientID = sampleTblwShotgunCounts.PatientID( ...
                sampleTblwShotgunCounts.GroupCount==max(sampleTblwShotgunCounts.GroupCount) ...
                )
patientTbl = tblASVsamples(strcmp(tblASVsamples.PatientID, patientID) & ...
                                        cellfun(@(X) ~isempty(X), tblASVsamples.AccessionShotgun), :);
patientTbl = sortrows(patientTbl, 'DayRelativeToNearestHCT');
patientSamples = patientTbl.SampleID;            
%
[~, sampleIdx, ~] = intersect(tblcounts_16S.SampleID, patientSamples, 'stable');
% patient_tblcounts = tblcounts(sampleIdx, :);
patient_color_grouped_abundance_16S = color_grouped_abundance_16S(sampleIdx,:);

%%% plot the taxonomic composition of a chosen patient
figure()
subplot(2,1,1)
h=bar(patient_color_grouped_abundance_16S, 'stacked',  'BarWidth', height(tblcounts_16S)/(height(tblcounts_16S)+9));
% set barplot color
% ASV_cmap = hex2rgb(uni_color_hex);
for cc = 1:size(patient_color_grouped_abundance_16S,2)
    h(cc).FaceColor = 'flat';
    h(cc).CData = ASV_cmap_16S(cc,:);
end
colormap(gca, ASV_cmap_16S);
ylim([0 1]);
h = gca;
h.XTick = 1:size(patient_color_grouped_abundance_16S, 1);
% h.XTickLabel = tblcounts_16S.SampleID(sampleIdx);
L = cellstr(num2str(patientTbl.DayRelativeToNearestHCT));
L = cellfun(@(X) strrep(X, ' ', ''), L, 'UniformOutput', false);
h.XTickLabel = L;
h.XTickLabelRotation = 90;
ylabel('relative abundance' )
% title(sprintf('PatientID=%s', patientID{:}), 'fontSize', 14)
title('16S rRNA gene sequencing', 'fontsize', 14)
subplot(2,1,2)

V=tblShotgunAbundances.Properties.VariableNames(10:end);
V=cellfun(@(X) X(2:end), V, 'UniformOutput', false);
V=cellfun(@(X) strrep(X, '_', '.'), V, 'UniformOutput', false);
[~, sampleIdx2, ~] = intersect(V, patientSamples, 'stable');

patient_color_grouped_abundance_shotgun = color_grouped_abundance_shotgun(sampleIdx2,:);

% plot stacked bars
h=bar(patient_color_grouped_abundance_shotgun, 'stacked',  'BarWidth', height(tblcounts_16S)/(height(tblcounts_16S)+9));
% set barplot color
ASV_cmap_shotgun = hex2rgb(uni_color_hex_shotgun);
for cc = 1:size(patient_color_grouped_abundance_shotgun,2)
    h(cc).FaceColor = 'flat';
    h(cc).CData = ASV_cmap_shotgun(cc,:);
end
colormap(gca, ASV_cmap_shotgun);
ylim([0 1]);
h = gca;
h.XTick = 1:size(patient_color_grouped_abundance_shotgun,1);
h.XTickLabel = L;
h.XTickLabelRotation = 90;
ylabel('relative abundance')
title('Shotgun metagenomic sequencing', 'fontsize', 14)

%% correlation analysis to compare taxanomic composition
%%%%% 1. correlation forr each sample
tblcountsM=tblcounts_16S{:, 4:end};
a = strcat('s', tblcounts_16S.SampleID); % 's' was added to the sample ID to be compatible with the formating of table; the same for '.' replaced by '_'
a = cellfun(@(X) strrep(X, '.', '_'), a, 'UniformOutput', false);
tblcounts2_16s=array2table(tblcountsM', 'VariableNames', a);
tblcounts2_16s.ASV=tblcounts_16S.Properties.VariableNames(4:end)';
Tbl_16S = join(tblcounts2_16s, tbltaxonomy(:, [1, 3:8]), 'Keys', 'ASV');
clear a
%%
sSamples = tblShotgunAbundances.Properties.VariableNames(9:end);
Taxa = {'Phylum' 'Class' 'Order' 'Family' 'Genus'};
c = {'Cor', 'Rsq', 'Pvalue'};
[Tx,Cx] = ndgrid(1:numel(Taxa),1:numel(c));
colnames = strcat(Taxa(Tx(:)),'_',c(Cx(:)));
corrT = zeros(length(sSamples),length(colnames));
rownames = cellfun(@(X) X(2:end), sSamples, 'UniformOutput', false);
rownames = cellfun(@(X) strrep(X, '_', '.'), rownames, 'UniformOutput', false);
corrT = array2table(corrT, 'VariableNames', colnames, 'RowNames', rownames);
corrS= [];
abdcutoff = 0;
sampleWithSingleDominace = {};

shannon_Index_shotgun = zeros(size(abundance_matrix_shotgun,1), length(Taxa));
shannon_Index_16S = zeros(size(abundance_matrix_shotgun,1), length(Taxa));
for j=1:length(Taxa)
    corrS(j).name = Taxa(j);
%     X = [];
    T_compare =[];
for i=1:length(sSamples)
    s1 = sSamples{i};
%     % correct sample name
    s2 = s1(2:end);
    s2=strrep(s2, '_', '.');
    st_shotgun = tblShotgunAbundances(:, {Taxa{j}, s1}); % shotgun
    st_shotgun_g = grpstats(st_shotgun, Taxa{j},'sum', 'DataVars', s1);
    st_shotgun_g = st_shotgun_g(:, [1, 3]);
    p1 = st_shotgun_g{:,2};
    p1=p1(p1>0);
    shannon_Index_shotgun(i, j) = -sum(p1 .* log(p1));
    st_16S = Tbl_16S(:, {Taxa{j}, s1});  % 16S
    st_16S_g= grpstats(st_16S, Taxa{j}, 'sum', 'DataVars', s1);
    st_16S_g=st_16S_g(:, [1 3]);
    p2 = st_16S_g{:,2};
    p2 = p2(p2>0);
    shannon_Index_16S(i, j) = -sum(p2 .* log(p2));
    clear p1 p2
    T = outerjoin(st_shotgun_g, st_16S_g, 'Keys', Taxa{j});
    
    %%% get all the phylum name
    T(cellfun(@(X) isempty(X), T{:,1}), [Taxa{j} '_st_shotgun_g'])= T(cellfun(@(X) isempty(X), T{:,1}), [Taxa{j} '_st_16S_g']);
    T(:, [Taxa{j} '_st_16S_g']) =[];
    T.Properties.VariableNames = {Taxa{j}, ['shotgun-' s2], ['s16S-' s2]};    
    t = T{:, 2:3};
    t(isnan(t))=-10;  % nan as -10  
    T{:, 2:3}=t;

    if i>1
        T_compare = outerjoin(T_compare, T, 'Keys', Taxa{j}, 'MergeKeys', true);
    else
        T_compare = T;
    end
    % calculate correlation, pvalue and R^2
    if height(T)>1
        t = T{:, 2:3};
        t(t==-10)=0;  % nan as -10  
        T{:, 2:3}=t;
        [R, p] = corrcoef( T{:,3},T{:,2});
        corrT{i, [Taxa{j} '_Pvalue']} = p(1,2);
        corrT{i, [Taxa{j} '_Rsq']} = R(1,2)^2;
%         x = [ones(length(T{:,'shotgun'}),1) T{:,'shotgun'}];
%         b = x\T{:,'16S'};
%         corrT{i, [Taxa{j} '_Cor']} = b(2);
        corrT{i, [Taxa{j} '_Cor']}=corr( T{:,3},T{:,2}, 'Type','Pearson');

    else
        sampleWithSingleDominace{end+1} = [s2 ' ' Taxa{j}];    
    end
%     X =[X; T];
    clear T
end
corrS(j).TaxaCompare = T_compare;
clear X
end

shannon_Index_shotgun = array2table(shannon_Index_shotgun, 'VariableNames', Taxa, 'RowNames', sSamples);
shannon_Index_16S = array2table(shannon_Index_16S, 'VariableNames', Taxa, 'RowNames', sSamples);
shannon_Index_shotgun.SampleID = rownames';
shannon_Index_16S.SampleID = rownames';
corrT.SampleID = rownames';

%% get taxas that are always missing in shotgun or 16S
TaxaNotIn16S =[];
TaxaNotInShotgun =[];

for i=1:length(corrS)
    T = corrS(i).TaxaCompare;
    shotgun = T{:, contains(T.Properties.VariableNames, 'shotgun')};
    shotgun_s = sum(shotgun,2);
    s16s = T{:, contains(T.Properties.VariableNames, '16S')};
    s16s_s = sum(s16s,2);
    %%%% find not in shotgun
    idx1 = find(shotgun_s == -10*size(shotgun,2) & s16s_s > -10*size(s16s,2));   
    notInShotgun = T( idx1,1);    
    notInShotgun.Taxa = repmat(notInShotgun.Properties.VariableNames, height(notInShotgun),1);
    notInShotgun.Properties.VariableNames{1} ='TaxaName';
    for k=1:length(idx1)
        notInShotgun.frequency(k) = sum((s16s(idx1(k), :) >0)) / size(s16s,2);  
        
        notInShotgun.median(k) = median(s16s(idx1(k), s16s(idx1(k), :) >0));
    end
    TaxaNotInShotgun = [TaxaNotInShotgun; notInShotgun];
    clear id
   
    %%%% find taxa not in 16S
    idx2 = find(shotgun_s > -10*size(shotgun,2) & s16s_s == -10*size(s16s,2));
    notIn16S = T( idx2,1);
    notIn16S.Taxa = repmat(corrS(i).name, height(notIn16S),1);
    notIn16S.Properties.VariableNames{1} ='TaxaName';
    for k=1:length(idx2)
        notIn16S.frequency(k) = sum((shotgun(idx2(k), :) >0)) / size(shotgun,2);  
        notIn16S.median(k) = median(shotgun (idx2(k), shotgun(idx2(k), :) >0));
    end
    TaxaNotIn16S =[TaxaNotIn16S; notIn16S];
    clear notIn16S notInShotgun idx2
end

TaxaNotInShotgun(ismember(TaxaNotInShotgun.TaxaName, '<not present>'), :) = [];
TaxaNotIn16S(ismember(TaxaNotIn16S.TaxaName, '<not present>'), :) = [];
TaxaNotInShotgun = sortrows(TaxaNotInShotgun, 'median', 'descend');
TaxaNotIn16S = sortrows(TaxaNotIn16S, 'frequency', 'descend');
writetable(TaxaNotInShotgun, 'TaxaNotInShotgun.csv');
writetable(TaxaNotIn16S, 'TaxaNotIn16S.csv');
%%

sID = corrT.SampleID;
[~, ia, ic] = intersect(sID, tblASVsamples.SampleID);
sT = tblASVsamples(ic, {'SampleID', 'readcount'});
corrT.SampleID = sID;
corrT = join(corrT, sT, 'Keys', 'SampleID');
save('../savedMat/corrT', 'corrT');

%% scatter plot of stool form vs correlation
T1 = tblASVsamples(:, {'SampleID', 'Consistency'});
T1_corr = outerjoin(T1, corrT, 'Keys', 'SampleID', 'mergeKeys', true);

%% statistical test to see if consistency impact on shotgun sequence read counts
%%%% Mann-Whitney test
X = log10(T1_corr.readcount);
Xcat = T1_corr.Consistency;
[p12,h12,stats12] = ranksum(X(ismember(Xcat, {'formed'})), ...
                        X(ismember(Xcat, {'semi-formed'})))
[p13,h13,stats13] = ranksum(X(ismember(Xcat, {'formed'})), ...
                        X(ismember(Xcat, {'liquid'})))
[p23,h23,stats23] = ranksum(X(ismember(Xcat, {'semi-formed'})), ...
                        X(ismember(Xcat, {'liquid'})))
%% anova
[p,tbl,stats] = anova1(X, Xcat)
[c,~,~,gnames] = multcompare(stats)
% the 6th column of c is the pvalue
[gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))]

%% show diversity among stool consistency
T1 = tblASVsamples(:, {'SampleID', 'Consistency'});
T1_Shannon = outerjoin(T1, shannon_Index_shotgun, 'Keys', 'SampleID', 'mergeKeys', true);

Shannon_stool_consistency = zeros(3, 3,length(Taxa));
for i=1:length(Taxa)    
    A = T1_Shannon{:, Taxa{i}};
    A_cutoff = prctile(A, [33 67]);
    Shannon_stool_consistency(1, 1, i) = length(find(ismember(T1_Shannon.Consistency, 'formed') & A <= A_cutoff(1)));
    Shannon_stool_consistency(1, 2, i) = length(find(ismember(T1_Shannon.Consistency, 'formed') & ...
                                                                A > A_cutoff(1) & A <= A_cutoff(2)));
    Shannon_stool_consistency(1, 3, i) = length(find(ismember(T1_Shannon.Consistency, 'formed') & ...
                                                                    A > A_cutoff(2)));
    Shannon_stool_consistency(2, 1, i) = length(find(ismember(T1_Shannon.Consistency, 'semi-formed') & A <= A_cutoff(1)));
    Shannon_stool_consistency(2, 2, i) = length(find(ismember(T1_Shannon.Consistency, 'semi-formed') & ...
                                                                A > A_cutoff(1) & A <= A_cutoff(2)));
    Shannon_stool_consistency(2, 3, i) = length(find(ismember(T1_Shannon.Consistency, 'semi-formed') & ...
                                                                    A > A_cutoff(2)));   
    Shannon_stool_consistency(3, 1, i) = length(find(ismember(T1_Shannon.Consistency, 'liquid') & A <= A_cutoff(1)));
    Shannon_stool_consistency(3, 2, i) = length(find(ismember(T1_Shannon.Consistency, 'liquid') & ...
                                                                A > A_cutoff(1) & A <= A_cutoff(2)));
    Shannon_stool_consistency(3, 3, i) = length(find(ismember(T1_Shannon.Consistency, 'liquid') & ...
                                                                    A > A_cutoff(2)));
end


%%
Tcorr_shannon = outerjoin(T1_corr, shannon_Index_shotgun, 'Keys', 'SampleID', 'mergeKeys', true);
save('../savedMat/Tcorr_shannon','Tcorr_shannon');



%% find genera that are the most different between 16S and shotgun

P=[];
Ranksum=[];
DiffGenera = {};
%%
grporder = {'shotgun', '16S'};
dotsz = 40;
abundance_cutoff = 0.1;
figure
n=1;

cmap = [.4 .4 .4;
        .4 .4 .4];
for k=1:length(Taxa)
    gT=corrS(k).TaxaCompare;
    gT_m = gT{:, 2:end};
    gT_m(gT_m==-10)=0;
    gT{:, 2:end}=gT_m;
    taxaUnit = unique(gT.(Taxa{k}));
    for i=1:length(taxaUnit)
    % for i=503
        idx = ismember(gT.(Taxa{k}), taxaUnit{i});
        gT_taxa = gT(idx, :);
        a1 = gT_taxa{:, contains(gT.Properties.VariableNames, 'shotgun-')};
        a2 = gT_taxa{:, contains(gT.Properties.VariableNames, 's16S')};
        [p, h, stats] = ranksum(a1, a2);
        if p<=0.05
            P(end+1)=p;
            Ranksum(end+1)=stats.ranksum;
            DiffGenera{end+1} =  taxaUnit{i};
            if median(a1)>=abundance_cutoff & median(a2)>=abundance_cutoff
                subplot(4,4,n)
                n=n+1;
                vs = violinplot([a1; a2], ...
                    [repmat({'shotgun'}, length(a1),1); repmat({'16S'}, length(a2),1)], ...
                    'Bandwidth', 0.15, ...
                    'ViolinAlpha', 1, ...
                    'ShowNotches', false, ...
                    'ShowMean', false, ...
                    'GroupOrder', grporder, ...
                    'BoxColor', [.5 0.1 0.1]);
                vs(1).MedianPlot.MarkerFaceColor = [.5 0.1 0.1];
                vs(2).MedianPlot.MarkerFaceColor = [.5 0.1 0.1];
                vs(1).MeanPlot.LineWidth=2;
                vs(2).MeanPlot.LineWidth=2;
                vs(1).ViolinColor = cmap(1,:);
                vs(2).ViolinColor = cmap(2,:);
                vs(1).ScatterPlot.SizeData=dotsz;
                vs(2).ScatterPlot.SizeData=dotsz;
                set(gca, 'xticklabels', {'shotgun', '16S'}, 'fontsize', 12);
                ylabel('relative abundances', 'fontsize', 12);
                text([0.8;1.8], [1.1; 1.1], ...
                    {sprintf('median=%.3f', median(a1)), sprintf('median=%.3f', median(a2))}, ...
                    'fontsize', 12)
                title([taxaUnit{i} ' (' Taxa{k} ')'], 'fontsize', 14)
                set(gca, 'ylim', [-.1 1.3], 'ytick', 0:0.2:1)

            end

        end

    end
end
