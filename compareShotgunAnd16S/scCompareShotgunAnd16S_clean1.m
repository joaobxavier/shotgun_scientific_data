close all; clear;
% read the sample table
% tblASVsamples = readtable('../metagenome_data/tblASVsamplesUpdatedWithShotgun.csv', 'Format', '%s%s%d%s%s%s%d%s');
tblASVsamples = readtable('../metagenome_data/tblASVsamplesUpdatedWithShotgunWithReadcounts_final.xlsx');
% pick the patient with most shotgun sequenced samples
idxHasShotun = cellfun(@(x) ~isempty(x), tblASVsamples.AccessionShotgun);
pt = groupcounts(tblASVsamples(idxHasShotun, :), 'PatientID');
pt = sortrows(pt, 'GroupCount', 'descend');

rewriteShotgunAbundances = 2;
%% plot the timeline for this patient
% chensCodeBaseDir = '../../MSKCC_Microbiome_SD2021_Scripts/';
%PatientID2Plot = pt.PatientID{1};
addpath('../utils');
% data_path = [chensCodeBaseDir 'deidentified_data_tables/']; % path to data
data_path = '../deidentified_data_tables/'; % path to data
% get the 16S sample data
opts = detectImportOptions(strcat(data_path, 'samples/tblASVsamples.csv'));
opts = setvartype(opts,{'PatientID'},'categorical');
tblsamples = readtable(strcat(data_path, 'samples/tblASVsamples.csv'),opts);

tblsamples = tblsamples(ismember(tblsamples.SampleID, tblASVsamples.SampleID(idxHasShotun)), :);
%tblsamples = tblsamples((tblsamples.PatientID==PatientID2Plot & ismember(tblsamples.SampleID, tblASVsamples.SampleID(idxHasShotun))), :);
tblsamples = sortrows(tblsamples, 'Timepoint'); % sort rows by time point of samples

% load counts table
tblcounts = readtable(strcat(data_path, 'counts/tblcounts_asv_melt.csv'));
tblcounts = tblcounts(contains(tblcounts.SampleID, tblsamples.SampleID), :);

% unstack counts table and normalize ASV counts to relative abundance
tblcounts = unstack(tblcounts, 'Count', 'ASV');
counts_matrix = tblcounts{:, 2:end}; % the first column is "SampleID"
counts_matrix(isnan(counts_matrix)) = 0; % missing count value is filled with 0
counts_matrix = counts_matrix ./ sum(counts_matrix, 2); % convert to relative abundance
tblcounts{:, 2:end} = counts_matrix;
tblcounts = innerjoin(tblsamples(:, {'SampleID', 'Timepoint', 'DayRelativeToNearestHCT'}), tblcounts);
tblcounts = sortrows(tblcounts, 'Timepoint'); % sort rows by time point of samples

% load taxonomy table
tbltaxonomy = readtable(strcat(data_path,'taxonomy/tblASVtaxonomy_silva132_v4v5_filter.csv'));
tbltaxonomy = tbltaxonomy(ismember(tbltaxonomy.ASV,tblcounts.Properties.VariableNames(4:end)), :);

% %############# deal with <not present> in tbltaxonomy ##################
columns = {'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus'};
for i=2:length(columns)
    notpresent =  find(ismember(tbltaxonomy.(columns{i}), '<not present>'));
    for j=1:length(notpresent)
        if contains(tbltaxonomy.(columns{i-1})(notpresent(j)), columns{i-1})
            tbltaxonomy.(columns{i})(notpresent(j))=tbltaxonomy.(columns{i-1})(notpresent(j));
        else
            tbltaxonomy.(columns{i}){notpresent(j)} = ...
                [tbltaxonomy.(columns{i-1}){notpresent(j)} '[' columns{i-1} ']'];         
        end        
%         tbltaxonomy.(columns{i})(ismember(tbltaxonomy.(columns{i}), '<not present>')) = ...
%         tbltaxonomy.(columns{i-1})(ismember(tbltaxonomy.(columns{i}), '<not present>')); 
    end
end
% %#########################################################################
 
% plot samples

% the first 3 columns of tblcounts are SampleID, Timepoint and DayRelativeToNearestHCT
abundance_matrix_16S = tblcounts{:, 4:end};

% calculate the cumulative sum of taxa with same color_order
% unique_color_order should be automatically sorted
[unique_color_order_16S,ia,ic] = unique(tbltaxonomy.ColorOrder);
uni_color_hex_16S = tbltaxonomy.HexColor(ia);

color_grouped_abundance_16S = zeros(size(abundance_matrix_16S,1), length(unique_color_order_16S));
for k = 1:length(unique_color_order_16S)
    currsum = sum(abundance_matrix_16S(:,ic==k),2);
    color_grouped_abundance_16S(:,k) = currsum;
end

% plot stacked bars
figure(1)
subplot(2, 1, 1)
h=bar(color_grouped_abundance_16S, 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
% set barplot color
ASV_cmap_16S = hex2rgb(uni_color_hex_16S);
for cc = 1:size(color_grouped_abundance_16S,2)
    h(cc).FaceColor = 'flat';
    h(cc).CData = ASV_cmap_16S(cc,:);
end
colormap(gca, ASV_cmap_16S);
ylim([0 1]);
h = gca;
h.XTick = 1:height(tblsamples);
h.XTickLabel = tblsamples.SampleID;
h.XTickLabelRotation = 90;
ylabel('relative abundance from 16S')

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
    %             % load the table for this file
    %             fprintf('Doing it now...\n')
    %             tblKraken2 = importKraken2Output(fn);
    %             save(sprintf('../PATRIC_output/kraken2/.%s_kraken2/tblKraken2.mat',...
    %                 s), 'tblKraken2');
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
                sumG.hit16S(idx) = tbltaxonomy{loc(idx), taxaFull{j}};
            end
            % if there's no match
            sumG.HexColor(sumG.ColorOrder == 0) = {'#000000'};
            sumG.hit16S(sumG.ColorOrder == 0) = {'none'};
            % join the tables
            if isempty(tblShotgunAbundances)
                tblShotgunAbundances = sumG(:, [taxa {'ColorOrder'} {'HexColor'} {'hit16S'} {['s' s]}]);
            else
                tblShotgunAbundances = outerjoin(tblShotgunAbundances, sumG(:, [taxa {'ColorOrder'} {'HexColor'} {'hit16S'} {['s' s]}]),...
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
    writetable(tblShotgunAbundances, 'tblShotgunAbundances.xlsx');
    abundance_matrix_shotgun =abundance_matrix_shotgun';
    figure(2)
    plot(S,'ro-')
    hold on
    plot(S2, 'bx-')
else
    tblShotgunAbundances=readtable('tblShotgunAbundances.xlsx');
    %%%%%% only keep bacteria and Archea in abundance table
    tblShotgunAbundances=tblShotgunAbundances(ismember(tblShotgunAbundances.K, {'Archaea' 'Bacteria'}), :);
    abundance_matrix_shotgun = tblShotgunAbundances{:, 10:end};
    abundance_matrix_shotgun =abundance_matrix_shotgun';
end

%%
figure(1)
subplot(2, 1, 2)
% calculate the cumulative sum of taxa with same color_order
% unique_color_order should be automatically sorted
[unique_color_order_shotgun,ia,ic] = unique(tblShotgunAbundances.ColorOrder);
uni_color_hex_shotgun = tblShotgunAbundances.HexColor(ia);

color_grouped_abundance_shotgun = zeros(size(abundance_matrix_shotgun,1), length(unique_color_order_shotgun));
for k = 1:length(unique_color_order_shotgun)
    currsum = sum(abundance_matrix_shotgun(:,ic==k),2);
    color_grouped_abundance_shotgun(:,k) = currsum;
end

% plot stacked bars
h=bar(color_grouped_abundance_shotgun, 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
% set barplot color
ASV_cmap_shotgun = hex2rgb(uni_color_hex_shotgun);
for cc = 1:size(color_grouped_abundance_shotgun,2)
    h(cc).FaceColor = 'flat';
    h(cc).CData = ASV_cmap_shotgun(cc,:);
end
colormap(gca, ASV_cmap_shotgun);
ylim([0 1]);
h = gca;
h.XTick = 1:height(tblsamples);
h.XTickLabel = tblsamples.SampleID;
h.XTickLabelRotation = 90;
ylabel('relative abundance from shotgun')
return
%%
%% use regression to find which taxa can seperate shotgun and 16S



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
patientSamples = patientTbl.SampleID                     
%
[~, sampleIdx, ~] = intersect(tblcounts.SampleID, patientSamples, 'stable');
% patient_tblcounts = tblcounts(sampleIdx, :);
patient_color_grouped_abundance_16S = color_grouped_abundance_16S(sampleIdx,:);

%%% plot the taxonomic composition of a chosen patient
figure(3)
subplot(2,1,1)
h=bar(patient_color_grouped_abundance_16S, 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
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
h.XTickLabel = tblcounts.SampleID(sampleIdx);
h.XTickLabelRotation = 90;
ylabel('relative abundance from 16S' )
title(sprintf('PatientID=%s', patientID{:}), 'fontSize', 14)

subplot(2,1,2)

V=tblShotgunAbundances.Properties.VariableNames(10:end);
V=cellfun(@(X) X(2:end), V, 'UniformOutput', false);
V=cellfun(@(X) strrep(X, '_', '.'), V, 'UniformOutput', false);
[~, sampleIdx2, ~] = intersect(V, patientSamples, 'stable');

patient_color_grouped_abundance_shotgun = color_grouped_abundance_shotgun(sampleIdx2,:);

% plot stacked bars
h=bar(patient_color_grouped_abundance_shotgun, 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
% set barplot color
% ASV_cmap_shotgun = hex2rgb(uni_color_hex_shotgun);
for cc = 1:size(patient_color_grouped_abundance_shotgun,2)
    h(cc).FaceColor = 'flat';
    h(cc).CData = ASV_cmap_shotgun(cc,:);
end
colormap(gca, ASV_cmap_shotgun);
ylim([0 1]);
h = gca;
h.XTick = 1:size(patient_color_grouped_abundance_shotgun,1);
% h.XTickLabel = tblsamples.SampleID(sampleIdx);
h.XTickLabel = V(sampleIdx2);
h.XTickLabelRotation = 90;
ylabel('relative abundance from shotgun')


%% plot stacked bar of a given sample
lowreadsample ={'1042U' '1044R'};                 
%
[~, sampleIdx, ~] = intersect(tblcounts.SampleID, lowreadsample, 'stable');
% patient_tblcounts = tblcounts(sampleIdx, :);
lowReads_color_grouped_abundance_16S = color_grouped_abundance_16S(sampleIdx,:);

%%% plot the taxonomic composition of a chosen patient
figure(4)
subplot(2,1,1)
if size(lowReads_color_grouped_abundance_16S,1)==1
    bp=[lowReads_color_grouped_abundance_16S; zeros(1,size(lowReads_color_grouped_abundance_16S,2))];
else
    bp=lowReads_color_grouped_abundance_16S;
end

h=bar(bp, ...
        'stacked',  'BarWidth', 0.5);
% set barplot color
% ASV_cmap = hex2rgb(uni_color_hex);
for cc = 1:size(lowReads_color_grouped_abundance_16S ,2)
    h(cc).FaceColor = 'flat';
    h(cc).CData = ASV_cmap_16S(cc,:);
end
colormap(gca, ASV_cmap_16S);
ylim([0 1]);
h = gca;
h.XTick = 1:size(lowReads_color_grouped_abundance_16S , 1);
h.XTickLabel = tblcounts.SampleID(sampleIdx);
h.XTickLabelRotation = 90;
ylabel('relative abundance from 16S' )
title('low read shotgun sample=%s', 'fontSize', 14)

subplot(2,1,2)
V=tblShotgunAbundances.Properties.VariableNames(10:end);
[~, sampleIdx2, ~] = intersect(V, strcat('s' ,lowreadsample), 'stable');

lowread_color_grouped_abundance_shotgun = color_grouped_abundance_shotgun(sampleIdx2,:);

% plot stacked bars
if size(lowread_color_grouped_abundance_shotgun,1)==1
    bp=[lowread_color_grouped_abundance_shotgun; zeros(1, size(lowread_color_grouped_abundance_shotgun,2))];
else
    bp=lowread_color_grouped_abundance_shotgun;
end
h=bar(bp, ...
        'stacked',  'BarWidth', 0.5);
% set barplot color
% ASV_cmap_shotgun = hex2rgb(uni_color_hex_shotgun);
for cc = 1:size(lowread_color_grouped_abundance_shotgun,2)
    h(cc).FaceColor = 'flat';
    h(cc).CData = ASV_cmap_shotgun(cc,:);
end
colormap(gca, ASV_cmap_shotgun);
ylim([0 1]);
h = gca;
h.XTick = 1:size(lowread_color_grouped_abundance_shotgun,1);
% h.XTickLabel = tblsamples.SampleID(sampleIdx);
h.XTickLabel = V(sampleIdx2);
h.XTickLabelRotation = 90;
ylabel('relative abundance from shotgun')


%% correlation analysis to compare taxanomic composition
%%%%% 1. correlation forr each sample
tblShotgunAbundances.Properties.VariableNames(1:6)={'Kingdom' 'Phylum' 'Class' 'Order' 'Family' 'Genus'};
% Tbl2=join(tblShotgunAbundances, tbltaxonomy(:, [1, 3:end]), 'Keys', {'Kingdom' 'Phylum' 'Class' 'Order' 'Family' 'Genus'});
tblcountsM=tblcounts{:, 4:end};
tblcounts2=array2table(tblcountsM', 'VariableNames', tblcounts.SampleID);
tblcounts2.ASV=tblcounts.Properties.VariableNames(4:end)';
Tbl_16S = join(tblcounts2, tbltaxonomy(:, [1, 3:8]), 'Keys', 'ASV');
%%
sSamples = tblShotgunAbundances.Properties.VariableNames(10:end);
% for i=1:length(sSamples)
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
for j=1:length(Taxa)
    corrS(j).name = Taxa(j);
    X = [];
for i=1:length(sSamples)
    s1 = sSamples{i};
%     % correct sample name
    s2 = s1(2:end);
    s2=strrep(s2, '_', '.');

    st1 = tblShotgunAbundances(:, {Taxa{j}, s1}); % shotgun
    st1g = grpstats(st1, Taxa{j},'sum', 'DataVars', s1);
    st1g = st1g(:, [1, 3]);
    st2 = Tbl_16S(:, {Taxa{j}, s2});  % 16S
    st2g= grpstats(st2, Taxa{j}, 'sum', 'DataVars', s2);
    st2g=st2g(:, [1 3]);
    T = outerjoin(st1g, st2g, 'Keys', Taxa{j});
    
    %%% get all the phylum name
    T(cellfun(@(X) isempty(X), T{:,1}), [Taxa{j} '_st1g'])= T(cellfun(@(X) isempty(X), T{:,1}), [Taxa{j} '_st2g']);
    T(:, [Taxa{j} '_st2g']) =[];
    T.Properties.VariableNames = {Taxa{j}, 'shotgun', '16S'};    
    t = T{:, 2:3};
    t(isnan(t))=0;    
    T{:, 2:3}=t;
    %%%%% add cutoff %%%
    T{:,'shotgun'}(T{:,'shotgun'}<=abdcutoff) = 0;
    T{:,'16S'}(T{:,'16S'}<=abdcutoff)=0;
    %%%%%%%%%%%
    T(T{:,'shotgun'} ==0 & T{:,'16S'}==0,:)=[];
    T.NotInShotgun=zeros(height(T),1);
    T.NotInShotgun(T.shotgun==0 & T{:,'16S'}>0)=1;
    T.NotIn16S=zeros(height(T),1);
    T.NotIn16S(T.shotgun>0 & T{:,'16S'}==0)=1;
    
    T.SampleID = repmat({s2}, height(T),1);
    % calculate correlation, pvalue and R^2
    if height(T)>1
        [R, p] = corrcoef( T{:,'16S'},T{:,'shotgun'});
        corrT{i, [Taxa{j} '_Pvalue']} = p(1,2);
        corrT{i, [Taxa{j} '_Rsq']} = R(1,2)^2;
%         x = [ones(length(T{:,'shotgun'}),1) T{:,'shotgun'}];
%         b = x\T{:,'16S'};
%         corrT{i, [Taxa{j} '_Cor']} = b(2);
        corrT{i, [Taxa{j} '_Cor']}=corr( T{:,'16S'},T{:,'shotgun'}, 'Type','Pearson');

    else
        sampleWithSingleDominace{end+1} = [s2 ' ' Taxa{j}];    
    end
    X =[X; T];
    clear T
end
corrS(j).TaxaCompare = X;
clear X
end

corrT.SampleID = sSamples';
%% barplot of samples that are not correlated
p_cutoff = 0.05;
% Idx = find(contains(corrT.Properties.VariableNames, '_Pvalue'));
V2=cellfun(@(X) X(2:end), V, 'UniformOutput', false);
V2=cellfun(@(X) strrep(X, '_', '.'), V2,'UniformOutput', false);
tblASVsamples.NotCorWith16S = repmat({''},height(tblASVsamples),1);
for i=1:length(Taxa)
% for i=1
    nonSigSamples=corrT.SampleID(corrT{:,[Taxa{i} '_Pvalue']} > p_cutoff);
    nonSigSamples = cellfun(@(X) strrep(X, '_', '.'), nonSigSamples, 'UniformOutput', false);
    nonSigSamples=regexprep(nonSigSamples,'^s','');
    [~, nI, ~] = intersect(tblASVsamples.SampleID, nonSigSamples, 'stable');
    tblASVsamples.NotCorWith16S(nI) = cellfun(@(A, B) [A '_' B], ...
                                    tblASVsamples.NotCorWith16S(nI),repmat({Taxa{i}}, length(nI),1), ...
                                    'UniformOutput', false);
    clear nI
    
    [~, nonSigIdx, ~] = intersect(tblcounts.SampleID, nonSigSamples, 'stable');
    NonSig_16S=color_grouped_abundance_16S(nonSigIdx,:);    
    figure(i+100)    
    subplot(2,3,1)
    if length(nonSigIdx) == 1
        h=bar([NonSig_16S; zeros(1, length(NonSig_16S))], 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
    else
    	h=bar(NonSig_16S, 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
    end
    % set barplot color
    for cc = 1:size(NonSig_16S,2)
        h(cc).FaceColor = 'flat';
        h(cc).CData = ASV_cmap_16S(cc,:);
    end
    colormap(gca, ASV_cmap_16S);
    ylim([0 1]);
    h = gca;
    h.XTick = 1:size(NonSig_16S, 1);
    h.XTickLabel = tblcounts.SampleID(nonSigIdx);
    h.XTickLabelRotation = 90;
    ylabel('relative abundance from 16S' )
    title(sprintf('Non-significant samples in %s', Taxa{i}), 'fontSize', 14)

    subplot(2,3,2) % plot only at the tested taxa level
    % set barplot color
    NonSig_16S2= abundance_matrix_16S(nonSigIdx, :);
    [uf, ufi, ufic]=unique(tbltaxonomy.(Taxa{i}), 'stable');
    taxa_order = zeros(length(uf),1);
    for k=1:length(taxa_order)
        [taxa_order(k), I] = max(tbltaxonomy.ColorOrder(ufic==k));
    end
    [taxa_order2, ufi2]=sort(taxa_order);
%     reIdx = ufi(ufi2);
    newUf=uf(ufi2);
    color_uf = zeros(size(NonSig_16S2,1), length(uf));
    for k = 1:length(uf)
        idx = ismember(tbltaxonomy.(Taxa{i}), newUf{k});
        currsum = sum(NonSig_16S2(:,idx),2);
        color_uf(:,k) = currsum;
    end
    if length(nonSigIdx) == 1
        h=bar([color_uf; zeros(1, length(color_uf))], 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
    else
    	h=bar(color_uf, 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
    end
    ASV_cmap_t = hex2rgb(tbltaxonomy.HexColor(ufi));
    ASV_cmap_t =ASV_cmap_t(ufi2,:);
    for cc = 1:size(color_uf,2)
        h(cc).FaceColor = 'flat';
        h(cc).CData = ASV_cmap_t(cc,:);
    end    
    colormap(gca, ASV_cmap_t)
    ylim([0 1]);
    h = gca;
    h.XTick = 1:size(NonSig_16S2, 1);
    h.XTickLabel = tblcounts.SampleID(nonSigIdx);
    h.XTickLabelRotation = 90;
    ylabel('relative abundance from 16S' )
    title(sprintf('Color by %s', Taxa{i}), 'fontSize', 14)
    
    subplot(2,3,4)
    [~, nonSigIdx2, ~] = intersect(V2, nonSigSamples, 'stable');
    nonSig_shotgun = color_grouped_abundance_shotgun(nonSigIdx2,:);
    % plot stacked bars
    if length(nonSigIdx2) == 1
        h=bar([nonSig_shotgun; zeros(1, length(nonSig_shotgun))], 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
    else
        h=bar(nonSig_shotgun, 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
    end
    % set barplot color
    for cc = 1:size(nonSig_shotgun,2)
        h(cc).FaceColor = 'flat';
        h(cc).CData = ASV_cmap_shotgun(cc,:);
    end
    colormap(gca, ASV_cmap_shotgun);
    ylim([0 1]);
    h = gca;
    h.XTick = 1:size(nonSig_shotgun,1);
    h.XTickLabel = V2(nonSigIdx2);
    h.XTickLabelRotation = 90;
    ylabel('relative abundance from shotgun')    
        
    subplot(2,3,5)   
    NonSig_shotgun2= abundance_matrix_shotgun(nonSigIdx2, :);
    [uf, ufi, ufic]=unique(tblShotgunAbundances.(Taxa{i}), 'stable');
    taxa_order_s = zeros(length(uf),1);
    for k=1:length(taxa_order_s)
        [taxa_order_s(k), I] = max(tblShotgunAbundances.ColorOrder(ufic==k));
    end
    [taxa_order_s2, ufi2]=sort(taxa_order_s);
    newUf_s=uf(ufi2);
    color_uf_s = zeros(size(NonSig_shotgun2,1), length(uf));
    for k = 1:length(uf)
        idx = ismember(tblShotgunAbundances.(Taxa{i}), newUf_s{k});
        currsum = sum(NonSig_shotgun2(:,idx),2);
        color_uf_s(:,k) = currsum;
    end     
%     [~, ufi2]=sort(tblShotgunAbundances.ColorOrder(ufi));
%     reIdx2 = ufi(ufi2);
%     taxa2=tblShotgunAbundances.(Taxa{i})(reIdx2);
%     color_uf_s = zeros(size(NonSig_shotgun2,1), length(uf));
%     for k = 1:length(uf)
%         currsum = sum(NonSig_shotgun2(:,ufic==k),2);
%         color_uf_s(:,k) = currsum;
%     end
    if length(nonSigIdx) == 1
        h=bar([color_uf_s; zeros(1, length(color_uf_s))], 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
    else
    	h=bar(color_uf_s, 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
    end
    ASV_cmap_t2 = hex2rgb(tblShotgunAbundances.HexColor(ufi));
    ASV_cmap_t2 = ASV_cmap_t2(ufi2,:);
    for cc = 1:size(color_uf_s,2)
        h(cc).FaceColor = 'flat';
        h(cc).CData = ASV_cmap_t2(cc,:);
    end       
    colormap(gca, ASV_cmap_t2)
    ylim([0 1]);
    h = gca;
    h.XTick = 1:size(NonSig_shotgun2,1);
    h.XTickLabel = V2(nonSigIdx2);
    h.XTickLabelRotation = 90;
    ylabel('relative abundance from shotgun') 

    subplot(2,3, [3 6])
    legendIdx1=find(sum([color_uf; zeros(1, length(color_uf))])>=0.1);
    ASV_cmap_t_plot = ASV_cmap_t(legendIdx1,:);
    taxa_plot = newUf(legendIdx1);
    for k=1:size(ASV_cmap_t_plot,1)
        rectangle('Position',[1,k,5,1],'FaceColor',ASV_cmap_t_plot(k,:),'EdgeColor','k','LineWidth',1)     
        text(1+8, k+0.3, taxa_plot{k})        
    end
    legendIdx2=find(sum([color_uf_s; zeros(1, length(color_uf_s))])>=0.1);
    ASV_cmap_t2_plot = ASV_cmap_t2(legendIdx2,:);
    taxa2_plot = newUf_s(legendIdx2);
    for k=1:size(ASV_cmap_t2_plot,1)
        rectangle('Position',[70,k,5,1],'FaceColor',ASV_cmap_t2_plot(k,:),'EdgeColor','k','LineWidth',1)     
        text(70+8, k+0.7, taxa2_plot{k})        
    end
    set(gca, 'xlim', [0 150])
    box off

end

%%
%%
%%
%% barplot of the two samples that have the lowest reads
lowreadsample ={'1042U' '1044R'};
barwidth = 0.5;
for i=1:length(Taxa)
% for i=1    
    [~, lowreadIdx, ~] = intersect(tblcounts.SampleID, lowreadsample, 'stable');
    lowread_16S=color_grouped_abundance_16S(lowreadIdx,:);    
    figure(i+50)    
    subplot(2,3,1)
    if length(lowreadIdx) == 1
        h=bar([lowread_16S; zeros(1, length(lowread_16S))], 'stacked',  'BarWidth', barwidth);
    else
    	h=bar(lowread_16S, 'stacked',  'BarWidth', barwidth);
    end
    % set barplot color
    for cc = 1:size(lowread_16S,2)
        h(cc).FaceColor = 'flat';
        h(cc).CData = ASV_cmap_16S(cc,:);
    end
    colormap(gca, ASV_cmap_16S);
    ylim([0 1]);
    h = gca;
    h.XTick = 1:size(lowread_16S, 1);
    h.XTickLabel = tblcounts.SampleID(lowreadIdx);
    h.XTickLabelRotation = 90;
    ylabel('relative abundance from 16S' )
    title(sprintf('Non-significant samples in %s', Taxa{i}), 'fontSize', 14)

    subplot(2,3,2) % plot only at the tested taxa level
    % set barplot color
    lowread_16S2= abundance_matrix_16S(lowreadIdx, :);
    [uf, ufi, ufic]=unique(tbltaxonomy.(Taxa{i}), 'stable');
    taxa_order = zeros(length(uf),1);
    for k=1:length(taxa_order)
        [taxa_order(k), I] = max(tbltaxonomy.ColorOrder(ufic==k));
    end
    [taxa_order2, ufi2]=sort(taxa_order);
%     reIdx = ufi(ufi2);
    newUf=uf(ufi2);
    color_uf = zeros(size(lowread_16S2,1), length(uf));
    for k = 1:length(uf)
        idx = ismember(tbltaxonomy.(Taxa{i}), newUf{k});
        currsum = sum(lowread_16S2(:,idx),2);
        color_uf(:,k) = currsum;
    end
    if length(lowreadIdx) == 1
        h=bar([color_uf; zeros(1, length(color_uf))], 'stacked',  'BarWidth', barwidth);
    else
    	h=bar(color_uf, 'stacked',  'BarWidth', barwidth);
    end
    ASV_cmap_t = hex2rgb(tbltaxonomy.HexColor(ufi));
    ASV_cmap_t =ASV_cmap_t(ufi2,:);
    for cc = 1:size(color_uf,2)
        h(cc).FaceColor = 'flat';
        h(cc).CData = ASV_cmap_t(cc,:);
    end    
    colormap(gca, ASV_cmap_t)
    ylim([0 1]);
    h = gca;
    h.XTick = 1:size(lowread_16S2, 1);
    h.XTickLabel = tblcounts.SampleID(lowreadIdx);
    h.XTickLabelRotation = 90;
    ylabel('relative abundance from 16S' )
    title(sprintf('Color by %s', Taxa{i}), 'fontSize', 14)
    
    subplot(2,3,4)
    [~, lowreadIdx2, ~] = intersect(V2, lowreadsample, 'stable');
    lowread_shotgun = color_grouped_abundance_shotgun(lowreadIdx2,:);
    % plot stacked bars
    if length(lowreadIdx2) == 1
        h=bar([lowread_shotgun; zeros(1, length(lowread_shotgun))], 'stacked',  'BarWidth', barwidth);
    else
        h=bar(lowread_shotgun, 'stacked',  'BarWidth', barwidth);
    end
    % set barplot color
    for cc = 1:size(lowread_shotgun,2)
        h(cc).FaceColor = 'flat';
        h(cc).CData = ASV_cmap_shotgun(cc,:);
    end
    colormap(gca, ASV_cmap_shotgun);
    ylim([0 1]);
    h = gca;
    h.XTick = 1:size(lowread_shotgun,1);
    h.XTickLabel = V2(lowreadIdx2);
    h.XTickLabelRotation = 90;
    ylabel('relative abundance from shotgun')    
        
    subplot(2,3,5)   
    lowread_shotgun2= abundance_matrix_shotgun(lowreadIdx2, :);
    [uf, ufi, ufic]=unique(tblShotgunAbundances.(Taxa{i}), 'stable');
    taxa_order_s = zeros(length(uf),1);
    for k=1:length(taxa_order_s)
        [taxa_order_s(k), I] = max(tblShotgunAbundances.ColorOrder(ufic==k));
    end
    [taxa_order_s2, ufi2]=sort(taxa_order_s);
    newUf_s=uf(ufi2);
    color_uf_s = zeros(size(lowread_shotgun2,1), length(uf));
    for k = 1:length(uf)
        idx = ismember(tblShotgunAbundances.(Taxa{i}), newUf_s{k});
        currsum = sum(lowread_shotgun2(:,idx),2);
        color_uf_s(:,k) = currsum;
    end     
%     [~, ufi2]=sort(tblShotgunAbundances.ColorOrder(ufi));
%     reIdx2 = ufi(ufi2);
%     taxa2=tblShotgunAbundances.(Taxa{i})(reIdx2);
%     color_uf_s = zeros(size(NonSig_shotgun2,1), length(uf));
%     for k = 1:length(uf)
%         currsum = sum(NonSig_shotgun2(:,ufic==k),2);
%         color_uf_s(:,k) = currsum;
%     end
    if length(lowreadIdx2) == 1
        h=bar([color_uf_s; zeros(1, length(color_uf_s))], 'stacked',  'BarWidth', barwidth);
    else
    	h=bar(color_uf_s, 'stacked',  'BarWidth', barwidth);
    end
    ASV_cmap_t2 = hex2rgb(tblShotgunAbundances.HexColor(ufi));
    ASV_cmap_t2 = ASV_cmap_t2(ufi2,:);
    for cc = 1:size(color_uf_s,2)
        h(cc).FaceColor = 'flat';
        h(cc).CData = ASV_cmap_t2(cc,:);
    end       
    colormap(gca, ASV_cmap_t2)
    ylim([0 1]);
    h = gca;
    h.XTick = 1:size(lowread_shotgun2,1);
    h.XTickLabel = V2(lowreadIdx2);
    h.XTickLabelRotation = 90;
    ylabel('relative abundance from shotgun') 

    subplot(2,3, [3 6])
    color_cutoff =0.03;
    legendIdx1=find(sum([color_uf; zeros(1, length(color_uf))])>=color_cutoff);
    ASV_cmap_t_plot = ASV_cmap_t(legendIdx1,:);
    taxa_plot = newUf(legendIdx1);
    for k=1:size(ASV_cmap_t_plot,1)
        rectangle('Position',[1,k,5,1],'FaceColor',ASV_cmap_t_plot(k,:),'EdgeColor','k','LineWidth',1)     
        text(1+8, k+0.3, taxa_plot{k})        
    end
    legendIdx2=find(sum([color_uf_s; zeros(1, length(color_uf_s))])>=color_cutoff);
    ASV_cmap_t2_plot = ASV_cmap_t2(legendIdx2,:);
    taxa2_plot = newUf_s(legendIdx2);
    for k=1:size(ASV_cmap_t2_plot,1)
        rectangle('Position',[70,k,5,1],'FaceColor',ASV_cmap_t2_plot(k,:),'EdgeColor','k','LineWidth',1)     
        text(70+8, k+0.7, taxa2_plot{k})        
    end
    set(gca, 'xlim', [0 150])
    box off

end
%% get taxas that are always missing in shotgun or 16S
TaxaNotIn16S =[];
TaxaNotInShotgun =[];
for i=1:length(corrS)
    missingTaxa = grpstats(corrS(i).TaxaCompare,  corrS(i).name{:}, 'sum','DataVar', {'NotIn16S', 'NotInShotgun'});
    missIn16S = missingTaxa(missingTaxa.sum_NotIn16S == missingTaxa.GroupCount & missingTaxa.sum_NotInShotgun==0 , ...
                            corrS(i).name{:});
    missIn16S.Properties.VariableNames = {'missingTaxa'};
    missIn16S.Taxa = repmat(corrS(i).name, height(missIn16S), 1);
    missIn16S.Properties.RowNames={};
    TaxaNotIn16S = [TaxaNotIn16S; missIn16S];
    clear missIn16S
    missInShotgun = missingTaxa(missingTaxa.sum_NotInShotgun == missingTaxa.GroupCount & missingTaxa.sum_NotIn16S==0 , ...
                            corrS(i).name{:});
    missInShotgun.Properties.VariableNames = {'missingTaxa'};
    missInShotgun.Taxa = repmat(corrS(i).name, height(missInShotgun), 1);
    missInShotgun.Properties.RowNames ={};
    TaxaNotInShotgun = [TaxaNotInShotgun; missInShotgun];
    clear missInShotgun
end
TaxaNotInShotgun(ismember(TaxaNotInShotgun.Taxa, '<not present>'), :) = [];
writetable(TaxaNotInShotgun, 'TaxaNotInShotgun.xlsx');
writetable(TaxaNotIn16S, 'TaxaNotIn16S.xlsx');
%%
figure
color1=[.2 .2 .8];
color2 = [.8, .2 .2];

for i=1:length(Taxa)
    subplot(5,2, i*2-1)
    histogram(corrT{:,[Taxa{i} '_Cor']}, 0:0.1:1.2, ...
                'EdgeColor', color1, 'FaceColor', color1)
    set(gca, 'xlim', [-0.2 1.2], 'xtick', 0:0.5:1, ...
                'ylim', [0 350], 'ytick', 0:100:300)
    if i==length(Taxa)
        ylabel('counts of samples')
        xlabel('Correlation (shotgun/16S)')
    end
%     title(Taxa{i})
    subplot(5,2, i*2)
    histogram(corrT{:,[Taxa{i} '_Rsq']},0:0.1:1.2, ...
                'EdgeColor', color2, 'FaceColor', color2)
    set(gca, 'xlim', [-0.2 1.2], 'xtick', 0:0.5:1, ...
                'ylim', [0 300], 'ytick', 0:100:300)
    if i==length(Taxa)
        ylabel('counts of samples')
        xlabel('R^2')
    end
    u=sum(corrT{:,[Taxa{i} '_Pvalue']}>0.05);
    corrT(corrT{:,[Taxa{i} '_Pvalue']}>0.05, {'SampleID' [Taxa{i} '_Cor'] [Taxa{i} '_Rsq'],[Taxa{i} '_Pvalue']})
    title(sprintf('%s \n N-nonsig = %i', Taxa{i}, u))

end
sID = corrT.SampleID;
sID = cellfun(@(X) X(2:end), sID, 'UniformOutput', false);
sID = cellfun(@(X) strrep(X, '_', '.'), sID,'UniformOutput', false);
[~, ia, ic] = intersect(sID, tblASVsamples.SampleID);
sT = tblASVsamples(ic, {'SampleID', 'readcount'});
corrT.SampleID = sID;
corrT = join(corrT, sT, 'Keys', 'SampleID');
%% scatter plot of read counts vs correlation
figure(10)
cmap = [.2 .2 .8; ...
        .4 .8 .8; ...
        .8, .2 .2;...
        .8 .8 .2];
dotS = 25;
for i=1:length(Taxa)
    subplot(5, 2, i*2-1)
    scatter(corrT{corrT.([Taxa{i} '_Pvalue'])<=0.05,[Taxa{i} '_Cor']}, ...
            log10(corrT.readcount(corrT.([Taxa{i} '_Pvalue'])<=0.05)), ...
            dotS, 'MarkerFaceColor', cmap(1, :),'MarkerEdgeColor', cmap(1,:)*0.6)
    n1=length(corrT{corrT.([Taxa{i} '_Pvalue'])<=0.05,[Taxa{i} '_Cor']});
    hold on
    scatter(corrT{corrT.([Taxa{i} '_Pvalue'])>0.05,[Taxa{i} '_Cor']}, ...
            log10(corrT.readcount(corrT.([Taxa{i} '_Pvalue'])>0.05)), ...
            dotS, 'MarkerFaceColor', cmap(2, :), 'MarkerEdgeColor', cmap(2, :)*0.6)
    n2=length(corrT{corrT.([Taxa{i} '_Pvalue'])>0.05,[Taxa{i} '_Cor']});
    legend({sprintf('significant (n=%i)', n1), ...
                sprintf('nonsignificant (n=%i)', n2)}, ...
                'Location','WestOutside')
    title(Taxa{i}, 'FontSize', 14)    
    ylabel('log10(readscount)')
    if i==length(Taxa)
        xlabel('Correlation', 'FontSize', 14)
    end
    set(gca, 'xlim', [-0.2 1])
    subplot(5, 2, i*2)
    scatter(corrT{corrT.([Taxa{i} '_Pvalue'])<=0.05,[Taxa{i} '_Rsq']}, ...
            log10(corrT.readcount(corrT.([Taxa{i} '_Pvalue'])<=0.05)), ...
            dotS, 'MarkerFaceColor', cmap(3, :),'MarkerEdgeColor', cmap(3,:)*0.6)
    hold on
    scatter(corrT{corrT.([Taxa{i} '_Pvalue'])>0.05,[Taxa{i} '_Rsq']}, ...
            log10(corrT.readcount(corrT.([Taxa{i} '_Pvalue'])>0.05)), ...
            dotS, 'MarkerFaceColor', cmap(4, :), 'MarkerEdgeColor', cmap(4,:)*0.6)
    title(Taxa{i}, 'FontSize', 14)
    legend({sprintf('significant (n=%i)', n1), ...
                sprintf('nonsignificant (n=%i)', n2)}, ...
                'Location','EastOutside')
    
    ylabel('log10(readscount)')
    if i==length(Taxa)
        xlabel('Rsquare', 'FontSize', 14)
    end
    set(gca, 'xlim', [0 1])
end
%% scatter plot of stool form vs correlation
T1 = tblASVsamples(:, {'SampleID', 'Consistency'});
T1 = outerjoin(T1, corrT, 'Keys', 'SampleID', 'mergeKeys', true);
%
figure(11)
cmap = [.2 .2 .8; ...
        .4 .8 .8; ...
        .8, .2 .2;...
        .8 .8 .2];
dotS = 25;
for i=1:length(Taxa)
    subplot(5, 2, i*2-1)
    scatter(T1{ismember(T1.Consistency, 'formed'),[Taxa{i} '_Cor']}, ...
            log10(T1.readcount(ismember(T1.Consistency, 'formed'))), ...
            dotS, 'MarkerFaceColor', cmap(1, :),'MarkerEdgeColor', cmap(1,:)*0.6)
    n1=length(T1{ismember(T1.Consistency, 'formed'),[Taxa{i} '_Cor']});
    hold on
    scatter(T1{ismember(T1.Consistency, 'semi-formed'),[Taxa{i} '_Cor']}, ...
            log10(T1.readcount(ismember(T1.Consistency, 'semi-formed'))), ...
            dotS, 'MarkerFaceColor', cmap(2, :), 'MarkerEdgeColor', cmap(2, :)*0.6)
    n2=length(corrT{ismember(T1.Consistency, 'semi-formed'),[Taxa{i} '_Cor']});
        hold on
    scatter(T1{ismember(T1.Consistency, 'liquid'),[Taxa{i} '_Cor']}, ...
            log10(T1.readcount(ismember(T1.Consistency, 'liquid'))), ...
            dotS, 'MarkerFaceColor', cmap(3, :), 'MarkerEdgeColor', cmap(3, :)*0.6)
    n3=length(corrT{ismember(T1.Consistency, 'liquid'),[Taxa{i} '_Cor']});
        
    legend({sprintf('formed (n=%i)', n1), ...
                sprintf('semi-formed (n=%i)', n2), ...
                sprintf('liquid (n=%i)', n3)}, ...
                'Location','WestOutside')
    title(Taxa{i}, 'FontSize', 14)    
    ylabel('log10(readscount)')
    if i==length(Taxa)
        xlabel('Correlation', 'FontSize', 14)
    end
    set(gca, 'xlim', [-0.1 1])
    subplot(5, 2, i*2)
    scatter(T1{ismember(T1.Consistency, 'formed'),[Taxa{i} '_Rsq']}, ...
            log10(T1.readcount(ismember(T1.Consistency, 'formed'))), ...
            dotS, 'MarkerFaceColor', cmap(1, :),'MarkerEdgeColor', cmap(1,:)*0.6)
    hold on
    scatter(T1{ismember(T1.Consistency, 'semi-formed'),[Taxa{i} '_Rsq']}, ...
            log10(T1.readcount(ismember(T1.Consistency, 'semi-formed'))), ...
            dotS, 'MarkerFaceColor', cmap(2, :), 'MarkerEdgeColor', cmap(2,:)*0.6)
        hold on
    scatter(T1{ismember(T1.Consistency, 'liquid'),[Taxa{i} '_Rsq']}, ...
            log10(T1.readcount(ismember(T1.Consistency, 'liquid'))), ...
            dotS, 'MarkerFaceColor', cmap(3, :), 'MarkerEdgeColor', cmap(3,:)*0.6)
    title(Taxa{i}, 'FontSize', 14)
    legend({sprintf('formed (n=%i)', n1), ...
                sprintf('semi-formed (n=%i)', n2), ...
                sprintf('liquid (n=%i)', n3)}, ...
                'Location','WestOutside')
    
    ylabel('log10(readscount)')
    if i==length(Taxa)
        xlabel('Rsquare', 'FontSize', 14)
    end
    set(gca, 'xlim', [0 1])
end

%%
addpath('../util');
grouporder = {'formed', 'semi-formed', 'liquid'};
stoolC=cellstr(T1.Consistency);
figure
vs = violinplot(log10(T1.readcount), stoolC, ...
                'Bandwidth', 0.2, ...
                'ViolinAlpha', 0.35, ...
                'ShowNotches', false, ...
                'ShowMean', true, ...
                'GroupOrder', grouporder);
vs(1).MedianPlot.MarkerFaceColor = [0 0 0];
vs(2).MedianPlot.MarkerFaceColor = [0 0 0];
vs(3).MedianPlot.MarkerFaceColor = [0 0 0];
vs(1).MeanPlot.LineWidth=3;
vs(2).MeanPlot.LineWidth=3;
vs(3).MeanPlot.LineWidth=3;
% vs(1).WhiskerPlot.MarkerFaceColor = [ 0 0 0];
% vs(1).WhiskerPlot.MarkerFaceColor = [ 0 0 0];
vs(1).ViolinColor = [.2 0.8 .3];
vs(2).ViolinColor = [.5 0.5 1];
vs(3).ViolinColor = [.8 0.3 0.2];
vs(1).WhiskerPlot.Color = [.2 .2 .2];
vs(1).ScatterPlot.SizeData=50;
vs(2).ScatterPlot.SizeData=50;
vs(3).ScatterPlot.SizeData=50;
xlabel('Stool Consistancy');
ylabel('log10(counts in shotgun samples)');
set(gca, 'xticklabel', ...
    {sprintf('formed (n=%i)', sum(ismember(tblASVsamples.Consistency, {'formed'}))), ...
            sprintf('semi-formed (n=%i)', sum(ismember(tblASVsamples.Consistency, {'semi-formed'}))), ...
            sprintf('liquid (n=%i)', sum(ismember(tblASVsamples.Consistency, {'liquid'})))}, ...
    'Fontsize',14);
xlim([0.2, 3.8]);
title('impact of stool consistency on shotgun sequencing reads')
%% statistical test to see if consistency impact on shotgun sequence read counts
%%%% Mann-Whitney test
X = log10(T1.readcount);
Xcat = T1.Consistency;
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

%%
%% check if consistency impacts on correlation
cmap2 = [.2 0.8 .3;
         .5 0.5 1;
         .8 0.3 0.2];
dotsz = 30;
stoolC=cellstr(T1.Consistency);
figure
for i=1:length(Taxa)
    subplot(5, 2, i*2-1)
vs = violinplot(T1{:, [Taxa{i} '_Cor']}, stoolC, ...
                'Bandwidth', 0.2, ...
                'ViolinAlpha', 0.3, ...
                'ShowNotches', false, ...
                'ShowMean', true, ...
                'GroupOrder', grouporder);
vs(1).MedianPlot.MarkerFaceColor = [0 0 0];
vs(2).MedianPlot.MarkerFaceColor = [0 0 0];
vs(3).MedianPlot.MarkerFaceColor = [0 0 0];
vs(1).MeanPlot.LineWidth=3;
vs(2).MeanPlot.LineWidth=3;
vs(3).MeanPlot.LineWidth=3;
vs(1).MeanPlot.Color = [.1 .1 .1];
% vs(1).MeanPlot.Line.Color=[0 0 0]
vs(1).ViolinColor = cmap2(1,:);
vs(2).ViolinColor = cmap2(2,:);
vs(3).ViolinColor = cmap2(3,:);
vs(1).WhiskerPlot.Color = [.2 .2 .2];
vs(1).ScatterPlot.SizeData=dotsz;
vs(2).ScatterPlot.SizeData=dotsz;
vs(3).ScatterPlot.SizeData=dotsz;
vs(1).MeanPlot.Color = cmap2(1,:)*0.6;
vs(2).MeanPlot.Color = cmap2(2,:)*0.6;
vs(3).MeanPlot.Color = cmap2(3,:)*0.6;
xlabel('Stool Consistancy');
ylabel('correlation');
n1n = sum(ismember(T1.Consistency, {'formed'}) & T1.([Taxa{i} '_Pvalue'])>0.05);
n2n = sum(ismember(T1.Consistency, {'sem-formed'}) & T1.([Taxa{i} '_Pvalue'])>0.05);
n3n = sum(ismember(T1.Consistency, {'liquid'}) & T1.([Taxa{i} '_Pvalue'])>0.05);
text([1:3]'-0.1, repmat(1.3, 3,1), {sprintf('%i', n1n), ...
            sprintf('%i', n2n),  ...
            sprintf('%i non-sig', n3n)}, ...
            'fontsize', 12)
set(gca, 'xticklabel', ...
    {sprintf('formed (n=%i)', sum(ismember(T1.Consistency, {'formed'}))), ...
            sprintf('semi-formed (n=%i)', sum(ismember(T1.Consistency, {'semi-formed'}))), ...
            sprintf('liquid (n=%i)', sum(ismember(T1.Consistency, {'liquid'})))}, ...
    'Fontsize',14);
set(gca, 'Ytick', 0:0.5:1)
xlim([0.2, 3.8]);
ylim([-.1 1.5])
title(Taxa{i})

subplot(5, 2, i*2)
vs = violinplot(T1{:, [Taxa{i} '_Rsq']}, stoolC, ...
                'Bandwidth', 0.2, ...
                'ViolinAlpha', 0.3, ...
                'ShowNotches', false, ...
                'ShowMean', true, ...
                'GroupOrder', grouporder);
vs(1).MedianPlot.MarkerFaceColor = [0 0 0];
vs(2).MedianPlot.MarkerFaceColor = [0 0 0];
vs(3).MedianPlot.MarkerFaceColor = [0 0 0];
vs(1).MeanPlot.LineWidth=3;
vs(2).MeanPlot.LineWidth=3;
vs(3).MeanPlot.LineWidth=3;
vs(1).ViolinColor = cmap2(1,:);
vs(2).ViolinColor = cmap2(2,:);
vs(3).ViolinColor = cmap2(3,:);
vs(1).WhiskerPlot.Color = [.2 .2 .2];
vs(1).ScatterPlot.SizeData=dotsz;
vs(2).ScatterPlot.SizeData=dotsz;
vs(3).ScatterPlot.SizeData=dotsz;
vs(1).MeanPlot.Color = cmap2(1,:)*0.5;
vs(2).MeanPlot.Color = cmap2(2,:)*0.5;
vs(3).MeanPlot.Color = cmap2(3,:)*0.5;
xlabel('Stool Consistancy');
ylabel('Rsquare');
text([1:3]'-0.1, repmat(1.3, 3,1), {sprintf('%i', n1n), ...
            sprintf('%i', n2n),  ...
            sprintf('%i non-sig', n3n)}, ...
      'fontsize', 12)
set(gca, 'xticklabel', ...
    {sprintf('formed (n=%i)', sum(ismember(T1.Consistency, {'formed'}))); ...
            sprintf('semi-formed (n=%i)', sum(ismember(T1.Consistency, {'semi-formed'}))); ...
            sprintf('liquid (n=%i)', sum(ismember(T1.Consistency, {'liquid'})))}, ...
    'Fontsize',14);
set(gca, 'Ytick', 0:0.5:1)
xlim([0.2, 3.8]);
ylim([-0.1 1.5])
title(Taxa{i})
end

%% find genera that are the most different between 16S and shotgun

P=[];
Ranksum=[];
DiffGenera = {};
%%
grporder = {'shotgun', '16S'};
dotsz = 40;
figure
n=1;
for k=1:length(Taxa)
    gT=corrS(k).TaxaCompare;
    taxaUnit = unique(gT.(Taxa{k}));
for i=1:length(taxaUnit)
% for i=503
    idx = ismember(gT.(Taxa{k}), taxaUnit{i});
    a1 = gT.shotgun(idx);
    a2 = gT.('16S')(idx);
    [p, h, stats] = ranksum(a1, a2);
    if p<=0.05
        P(end+1)=p;
        Ranksum(end+1)=stats.ranksum;
        DiffGenera{end+1} =  taxaUnit{i};
        if median(a1)>0.01 | median(a2)>0.01
            subplot(4,4,n)
            n=n+1;
            vs = violinplot([a1; a2], ...
                [repmat({'shotgun'}, length(a1),1); repmat({'16S'}, length(a2),1)], ...
                'Bandwidth', 0.3, ...
                'ViolinAlpha', 0.3, ...
                'ShowNotches', false, ...
                'ShowMean', true, ...
                'GroupOrder', grporder);
            vs(1).MedianPlot.MarkerFaceColor = [0 0 0];
            vs(2).MedianPlot.MarkerFaceColor = [0 0 0];
            vs(1).MeanPlot.LineWidth=2;
            vs(2).MeanPlot.LineWidth=2;
            vs(1).ViolinColor = [.9 0.7 .3];
            vs(2).ViolinColor = [.5 0.1 1];
            vs(1).ScatterPlot.SizeData=dotsz;
            vs(2).ScatterPlot.SizeData=dotsz;
            set(gca, 'xticklabels', {'shotgun', '16S'}, 'fontsize', 12);
            ylabel('relative abundances', 'fontsize', 12);
            text([0.8;1.8], [1.1; 1.1], ...
                {sprintf('median=%.3f', median(a1)), sprintf('median=%.3f', median(a2))}, ...
                'fontsize', 12)
            title([taxaUnit{i} ' (' Taxa{k} ')'], 'fontsize', 14)
            set(gca, 'ylim', [-.1 1.3], 'ytick', 0:0.2:1)
            
%             figure
%             counta1=0;
%             counta2=0;
%             for j=1:length(a1)
%                 p = plot([1 2], [a1(j) a2(j)], 'o', ...
%                     'MarkerSize', 8, 'MarkerEdgeColor', [.2 .2 .2], ...
%                     'MarkerFaceColor', [.7 .7 .7]);
%                 p.Color(4) = 0.2;
%                 hold on
%                 if a1(j) > a2(j)
%                     clr = [.1 .4 1];
%                     counta1=counta1+1;
%                 elseif a1(j) < a2(j)
%                     clr = [1 .4 .1];
%                     counta2=counta2+1;
%                 elseif a1(j) == a2(j)
%                     clr = [.7 .7 .7];
%                 end
%                 line([1 2], [a1(j) a2(j)], 'Color', clr)
%                 hold on                
%             end
%             hold off
%             ylabel('relative abundances', 'fontsize', 14);
%             set(gca, 'xlim', [0.5 2.5], 'ylim', [-.1 1.1], 'ytick', 0:0.2:1);
%             set(gca, 'xtick', 1:1:2,'xticklabels', {'shotgun', '16S'}, 'fontsize', 14);
%             title([taxaUnit{i} '(' Taxa{k} ')'], 'fontsize', 16)
        end
        
    end
    
end
end
%% scatter plot of taxonomic composition for each unit
% 
% cutoff2 = 0.3;
% cmap = [.3 .3 .3;
%         .1 .6 .9;
%         .9 .6 .1];
% dotS = 45;
% lineColor = [1 1 1];
% Tall = [];
% Counts = sprintf('Names      NotIn16S     NotInShotgun\n\n');
% disp(Counts)
% fig6=figure(6);
% for i=1:length(corrS)
% % for i=1
%     
%     X = corrS(i).TaxaCompare;
%     X(strcmp(X{:, corrS(i).name}, '<not present>'),:)=[];
%     X = join(X, sT);
%     X.logreads = log10(X.readcount);
% %     size(X)
%     subplot(1, length(corrS), i)
%     scatter(X.('16S')(X.NotIn16S==0 & X.NotInShotgun==0) , ...
%                 X.shotgun(X.NotIn16S==0 & X.NotInShotgun==0), dotS, ...
%                 'MarkerEdgeColor', lineColor, ...
%                 'MarkerFaceColor', cmap(1,:));
%     l0=length(X.('16S')(X.NotIn16S==0 & X.NotInShotgun==0));
%     hold on
%     scatter(X.('16S')(X.NotIn16S==1 & X.NotInShotgun==0) , ...
%         X.shotgun(X.NotIn16S==1 & X.NotInShotgun==0), dotS, ...
%         'MarkerEdgeColor', lineColor, ...
%         'MarkerFaceColor', cmap(2,:));
%     l1=length(X.('16S')(X.NotIn16S==1 & X.NotInShotgun==0));
%     hold on
%     scatter(X.('16S')(X.NotIn16S==0 & X.NotInShotgun==1) , ...
%         X.shotgun(X.NotIn16S==0 & X.NotInShotgun==1), dotS, ...
%         'MarkerEdgeColor', lineColor, ...
%         'MarkerFaceColor', cmap(3,:));
%     l2=length(X.('16S')(X.NotIn16S==0 & X.NotInShotgun==1));
%     R = corrcoef(X.('16S'), X.shotgun);
%     R2 = R(1,2)^2;
%     xlabel('16S')
%     ylabel('shotgun')
%     title(sprintf('%s\nR^2=%.2f', corrS(i).name{:}, R2))
% %     legend({sprintf('In both (n=%d)', l0), ...
% %             sprintf('Not in 16S (n=%d)', l1), ...
% %             sprintf('Not in shotgun (n=%d)', l2)}, ...
% %             'Location','SouthOutside')
%     if i==length(corrS)
%         legend({'In common', ...
%                 'Not in 16S', ...
%                 'Not in shotgun'}, ...
%                 'Location','SouthOutside')
%     end
%     TnotIn16S = X(X.shotgun>cutoff2 & X.('16S')==0 , :);
%     TnotInShotgun = X(X.shotgun==0 & X.('16S')>cutoff2 , :);
%     T = [TnotIn16S; TnotInShotgun];
%     T.taxa = repmat(corrS(i).name, height(T),1);
%     T.Properties.VariableNames(1) = {'Names'};
%     Tall = [Tall; T];
% 
%     Counts = sprintf('%s\t%d\t%d\n', ...
%                      corrS(i).name{:}, height(TnotIn16S), height(TnotInShotgun) );
%     disp(Counts)
%      
%     clear TnotIn16S TnotInShotgun 
% end
% 
% fig6.Renderer='Painters';
% writetable(Tall, 'uniqueTaxa.xlsx');

%%
% %% scatter plot of taxonomic composition for each unit using reads as color
% % 
% 
% colormap(jet(10));
% dotS = 25;
% 
% fig7=figure(7);
% for i=1:length(corrS)
% % for i=1
%     
%     X = corrS(i).TaxaCompare;
%     X(strcmp(X{:, corrS(i).name}, '<not present>'),:)=[];
%     X = join(X, sT);
%     X.logreads = log10(X.readcount);
%     C = X.logreads;
%     C(isinf(C))=min(C)-1;
%     %     size(X)
%     subplot(1, length(corrS), i)
%     hScat=scatter(X.('16S') , ...
%                 X.shotgun, dotS, ...
%                 C, 'filled');
%     xlabel('16S')
%     ylabel('shotgun')
%     title(corrS(i).name{:})
%     hchar=colorbar
%     caxis([5 8.5])
%             
% %     l0=length(X.('16S')(X.NotIn16S==0 & X.NotInShotgun==0));
% %     hold on
% %     scatter(X.('16S')(X.NotIn16S==1 & X.NotInShotgun==0) , ...
% %         X.shotgun(X.NotIn16S==1 & X.NotInShotgun==0), dotS, ...
% %         'MarkerEdgeColor', lineColor, ...
% %         'MarkerFaceColor', cmap(2,:));
% %     l1=length(X.('16S')(X.NotIn16S==1 & X.NotInShotgun==0));
% %     hold on
% %     scatter(X.('16S')(X.NotIn16S==0 & X.NotInShotgun==1) , ...
% %         X.shotgun(X.NotIn16S==0 & X.NotInShotgun==1), dotS, ...
% %         'MarkerEdgeColor', lineColor, ...
% %         'MarkerFaceColor', cmap(3,:));
% %     l2=length(X.('16S')(X.NotIn16S==0 & X.NotInShotgun==1));
% %     R = corrcoef(X.('16S'), X.shotgun);
% %     R2 = R(1,2)^2;
% 
% %     title(sprintf('%s\nR^2=%.2f', corrS(i).name{:}, R2))
% % %     legend({sprintf('In both (n=%d)', l0), ...
% % %             sprintf('Not in 16S (n=%d)', l1), ...
% % %             sprintf('Not in shotgun (n=%d)', l2)}, ...
% % %             'Location','SouthOutside')
% %     if i==length(corrS)
% %         legend({'In common', ...
% %                 'Not in 16S', ...
% %                 'Not in shotgun'}, ...
% %                 'Location','SouthOutside')
% %     end
% %     TnotIn16S = X(X.shotgun>cutoff2 & X.('16S')==0 , :);
% %     TnotInShotgun = X(X.shotgun==0 & X.('16S')>cutoff2 , :);
% %     T = [TnotIn16S; TnotInShotgun];
% %     T.taxa = repmat(corrS(i).name, height(T),1);
% %     T.Properties.VariableNames(1) = {'Names'};
% %     Tall = [Tall; T];
% % 
% %     Counts = sprintf('%s\t%d\t%d\n', ...
% %                      corrS(i).name{:}, height(TnotIn16S), height(TnotInShotgun) );
% %     disp(Counts)
% %      
% %     clear TnotIn16S TnotInShotgun 
% end
% 
% fig7.Renderer='Painters';
