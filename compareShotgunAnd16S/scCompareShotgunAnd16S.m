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
subplot(2, 1, 1)
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
    % check if Kraken2 file exist at all
    if isfile(fn)
        % check if matlab table has been precomputed
        if isfile(sprintf('../PATRIC_output/kraken2/.%s_kraken2/tblKraken2.mat', s))
            fprintf('Sample %s has been processed by matlab. Loading...\n', s)
            load(sprintf('../PATRIC_output/kraken2/.%s_kraken2/tblKraken2.mat', s));
        else
            fprintf('Sample %s has not been processed by matlab yet.\n', s)
            fprintf('Doing it now...\n')
            srr = tblASVsamples.AccessionShotgun{strcmp(tblASVsamples.SampleID, s)};
            log = sprintf('../PATRIC_output/kraken2/%s_kraken2', s);
            [status,result] = system(['cat ' log ' | grep ' srr]);
            if status == 1 || isempty(result)
                error('Sample %s not confirmed by %s', s, log)
            else
                fprintf('Sample %s confirmed by %s\n', s, log)
            end
            % load the table for this file
            tblKraken2 = importKraken2Output(fn);
            save(sprintf('../PATRIC_output/kraken2/.%s_kraken2/tblKraken2.mat',...
                s), 'tblKraken2');
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
        for i = 2:length(taxa)
            t = sumG{:, taxa{i}};
            [~ , loc] = ismember(t, tbltaxonomy{:, taxaFull{i}});
            idx = loc>0;
            sumG.ColorOrder(idx) = tbltaxonomy.ColorOrder(loc(idx));
            sumG.HexColor(idx) = tbltaxonomy.HexColor(loc(idx));
            sumG.hit16S(idx) = tbltaxonomy{loc(idx), taxaFull{i}};
        end
        % if there's no match
        sumG.HexColor(sumG.ColorOrder == 0) = {'#000000'};
        sumG.hit16S(sumG.ColorOrder == 0) = {'none'};
        % join the tables
        if isempty(tblShotgunAbundances)
            tblShotgunAbundances = sumG(:, [taxa {'ColorOrder'} {'HexColor'} {'hit16S'} {['s' s]}]);
        else
            tblShotgunAbundances = innerjoin(tblShotgunAbundances, sumG(:, [taxa {'ColorOrder'} {'HexColor'} {'hit16S'} {['s' s]}]));
        end
    else
        warning('Sample %s has no Kraken2 output yet...\n', s)
        fprintf('>>>>>> Go PATRIC and request processing %s\n', s);
    end
end

%%
figure(1)
subplot(2, 1, 2)
% the first 3 columns of tblcounts are SampleID, Timepoint and DayRelativeToNearestHCT
abundance_matrix = tblShotgunAbundances{:, 10:end}';

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