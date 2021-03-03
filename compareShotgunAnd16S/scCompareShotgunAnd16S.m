% read the sample table
tblASVsamples = readtable('../tblASVsamplesUpdatedWithShotgun.csv', 'Format', '%s%s%d%s%s%s%d%s');

% pick the patient with most shotgun sequenced samples
idxHasShotun = cellfun(@(x) ~isempty(x), tblASVsamples.AccessionShotgun);
pt = groupcounts(tblASVsamples(idxHasShotun, :), 'PatientID');
pt = sortrows(pt, 'GroupCount', 'descend');


%% plot the timeline for this patient
chensCodeBaseDir = '../../MSKCC_Microbiome_SD2021_Scripts/';
%PatientID2Plot = pt.PatientID{1};
addpath([chensCodeBaseDir 'utils']);
data_path = [chensCodeBaseDir 'deidentified_data_tables/']; % path to data

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

% plot samples

% the first 3 columns of tblcounts are SampleID, Timepoint and DayRelativeToNearestHCT
abundance_matrix = tblcounts{:, 4:end};

% calculate the cumulative sum of taxa with same color_order
% unique_color_order should be automatically sorted
[unique_color_order,ia,ic] = unique(tbltaxonomy.ColorOrder);
uni_color_hex = tbltaxonomy.HexColor(ia);

color_grouped_abundance = zeros(size(abundance_matrix,1), length(unique_color_order));
for k = 1:length(unique_color_order)
    currsum = sum(abundance_matrix(:,ic==k),2);
    color_grouped_abundance(:,k) = currsum;
end

% plot stacked bars
figure(1)
subplot(3, 1, 1)
h=bar(color_grouped_abundance, 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
% set barplot color
ASV_cmap = hex2rgb(uni_color_hex);
for cc = 1:size(color_grouped_abundance,2)
    h(cc).FaceColor = 'flat';
    h(cc).CData = ASV_cmap(cc,:);
end
colormap(gca, ASV_cmap);
ylim([0 1]);
h = gca;
h.XTick = 1:height(tblsamples);
h.XTickLabel = tblsamples.SampleID;
h.XTickLabelRotation = 90;
ylabel('relative abundance from 16S')



%% make the sample plot to shotgun
addpath('../');
tblShotgunAbundances = [];

for i = 1:height(tblsamples)
    s = tblsamples.SampleID{i};
    fn  = sprintf('../PATRIC_output/kraken2/.%s_kraken2/report.txt', s);
    srr = tblASVsamples.AccessionShotgun{strcmp(tblASVsamples.SampleID, s)};
    % check if matlab table has been precomputed
    if isfile(sprintf('../PATRIC_output/kraken2/.%s_kraken2/tblKraken2.mat', s))
        fprintf('Sample %s has been processed by matlab. Loading...\n', s)
        load(sprintf('../PATRIC_output/kraken2/.%s_kraken2/tblKraken2.mat', s));
        % compute relative abundances
        % add all the fragments mapped to a genus and show the result
        % ordered by the most adundant in sample
        taxa = {'K' 'P' 'C' 'O' 'F' 'G'};
        taxaFull = {'Kingdom' 'Phylum' 'Class' 'Order' 'Family' 'Genus'};
        tblKrakenBugs = tblKraken2(ismember(tblKraken2.K, {'Bacteria' 'Archaea'}), :);
        sumG = grpstats(tblKrakenBugs, taxa, 'sum', 'DataVars', 'nFragsThis');
        sumG{:, {['s' s]}} = sumG.sum_nFragsThis ./ sum(sumG.sum_nFragsThis);
        % join the tables
        if isempty(tblShotgunAbundances)
            tblShotgunAbundances = sumG(:, [taxa {['s' s]}]);
        else
            tblShotgunAbundances = outerjoin(tblShotgunAbundances,...
                sumG(:, [taxa {['s' s]}]), 'MergeKeys', true);
        end
    else
        str =...
            sprintf('!!!!Sample %s has no Kraken2 output yet. Go to PATRIC and request processing %s\n',...
            s, srr);
        warning(str);
        fid = fopen('samplesWithoutKraken2.txt', 'a+');
        fprintf(fid, '%s\n', str);
        fclose(fid);
    end
    
end
%% replace NaN with 0s
m = tblShotgunAbundances{:, 7:end}; 
m(isnan(m)) = 0;
tblShotgunAbundances{:, 7:end} = m;
%% find the matches with the 16S classification
% match shotgun taxa to best 16S
for i = 2:length(taxa)
    t = tblShotgunAbundances{:, taxa{i}};
    [~ , loc] = ismember(t, tbltaxonomy{:, taxaFull{i}});
    idx = loc>0;
    tblShotgunAbundances.ColorOrder(idx) = tbltaxonomy.ColorOrder(loc(idx));
    tblShotgunAbundances.HexColor(idx) = tbltaxonomy.HexColor(loc(idx));
    tblShotgunAbundances.hit16S(idx) = tbltaxonomy{loc(idx), taxaFull{i}};
end
% if there's no match
tblShotgunAbundances.HexColor(tblShotgunAbundances.ColorOrder == 0) = {'#000000'};
tblShotgunAbundances.hit16S(tblShotgunAbundances.ColorOrder == 0) = {'none'};


%%
figure(1)
subplot(3, 1, 2)
% the first 3 columns of tblcounts are SampleID, Timepoint and DayRelativeToNearestHCT
abundance_matrix = tblShotgunAbundances{:, 10:end-3}';

% calculate the cumulative sum of taxa with same color_order
% unique_color_order should be automatically sorted
[unique_color_order,ia,ic] = unique(tblShotgunAbundances.ColorOrder);
uni_color_hex = tblShotgunAbundances.HexColor(ia);

color_grouped_abundance = zeros(size(abundance_matrix,1), length(unique_color_order));
for k = 1:length(unique_color_order)
    currsum = sum(abundance_matrix(:,ic==k),2);
    color_grouped_abundance(:,k) = currsum;
end

% plot stacked bars
h=bar(color_grouped_abundance, 'stacked',  'BarWidth', height(tblcounts)/(height(tblcounts)+9));
% set barplot color
ASV_cmap = hex2rgb(uni_color_hex);
for cc = 1:size(color_grouped_abundance,2)
    h(cc).FaceColor = 'flat';
    h(cc).CData = ASV_cmap(cc,:);
end
colormap(gca, ASV_cmap);
ylim([0 1]);
h = gca;
h.XTick = 1:height(tblsamples);
h.XTickLabel = tblsamples.SampleID;
h.XTickLabelRotation = 90;
ylabel('relative abundance from shotgun')

%% make a tbl16SAbundances with a similar format to tblShotgunAbundances
tbl16SAbundances = tblcounts;
tbl16SAbundances.Timepoint = [];
tbl16SAbundances.DayRelativeToNearestHCT = [];
tbl16SAbundances.SampleID = cellfun(@(x) ['s' x], tbl16SAbundances.SampleID, 'UniformOutput', false);
SampleID = tbl16SAbundances.SampleID;
ASV = tbl16SAbundances.Properties.VariableNames(2:end);
m = tbl16SAbundances{:, 2:end}';
tbl16SAbundances = array2table(m);
tbl16SAbundances.Properties.VariableNames = SampleID;
tbl16SAbundances.ASV = ASV';
tbl16SAbundances = innerjoin(tbl16SAbundances, tbltaxonomy(:, [{'ASV'} taxaFull]));

%% correlaiton plots
% phylum level

for i = 2:length(taxaFull)
    sum16S = grpstats(tbl16SAbundances, taxaFull{i}, 'sum', 'DataVars', SampleID);
    sumShotgun = grpstats(tblShotgunAbundances, taxa{i}, 'sum', 'DataVars', SampleID);
    taxaInBoth = intersect(sumShotgun{:, taxa{i}},...
        sum16S{:, taxaFull{i}});
    x16S = sum16S{taxaInBoth, 3:end};
    yShotgun = sumShotgun{taxaInBoth, 3:end};

    figure(1)
    subplot(3, 5, 10+(i-1))
    plot(x16S, yShotgun, 'ko', 'MarkerSize', 2)
    %set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('16S')
    ylabel('shotgun')
    title(sprintf('%s\nR^2 = %0.2f', taxaFull{i}, corr(x16S(:), yShotgun(:))^2));
    axis square equal tight
end

%%
set(gcf, 'Position', [44         119        1397         679]);
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, 'comparison16sShotgun.eps', '-depsc')