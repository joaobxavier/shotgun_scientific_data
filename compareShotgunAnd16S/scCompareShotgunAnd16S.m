% read the sample table
tblASVsamples = readtable('../tblASVsamplesUpdatedWithShotgun.csv', 'Format', '%s%s%d%s%s%s%d%s');

% pick the patient with most shotgun sequenced samples
idxHasShotun = cellfun(@(x) ~isempty(x), tblASVsamples.AccessionShotgun);
pt = groupcounts(tblASVsamples(idxHasShotun, :), 'PatientID');
pt = sortrows(pt, 'GroupCount', 'descend');


%% plot the timeline for this patient
chensCodeBaseDir = '../../MSKCC_Microbiome_SD2021_Scripts/';
PatientID2Plot = pt.PatientID{1};
addpath([chensCodeBaseDir 'utils']);
data_path = [chensCodeBaseDir 'deidentified_data_tables/']; % path to data

% get the 16S sample data
opts = detectImportOptions(strcat(data_path, 'samples/tblASVsamples.csv'));
opts = setvartype(opts,{'PatientID'},'categorical');
tblsamples = readtable(strcat(data_path, 'samples/tblASVsamples.csv'),opts);
tblsamples = tblsamples((tblsamples.PatientID==PatientID2Plot & ismember(tblsamples.SampleID, tblASVsamples.SampleID(idxHasShotun))), :);
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


