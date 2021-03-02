% download the samples table from figshare
url = 'https://ndownloader.figshare.com/files/26070716';
filename = 'tblASVsamples.csv';
outfilename = websave(['tempFiles/' filename] , url);
%% load the table
tblASVsamples = readtable(['tempFiles/' filename], 'Format', '%s%s%d%s%s%s%d');

%% download the latest metadata on all runs submitted to the SRA
% beleonging to MSKCC allo-HCT projects
% 1. nagigate to https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA394877%2CPRJNA491657%2CPRJNA479462%2CPRJNA545312%2CPRJNA548153%2CPRJNA576440%2CPRJNA607574%2CPRJNA606262%2CPRJNA701319%20&o=acc_s%3Aa
% 2. click on 'Metadata' to donwload a file called 'SraRunTable.txt'
% 3. save the file 'SraRunTable.txt' to 'tempFiles'
SraRunTable = readtable('tempFiles/metadata-8924168-processed-ok.tsv',...
    'FileType', 'text', 'Delimiter', '\t');

% keep only the shotgun sequences submitted in 2021 or late 2020
% bacause those the high quality ones
%idx = strcmp(SraRunTable.library_strategy, 'WGS') & contains(SraRunTable.ReleaseDate, {'2021-' '2020-12'});
%SraRunTable(~idx, :) = [];

%% keep only the Libray Name and the Run
SraRunTable = SraRunTable(:, {'library_ID' 'accession'});
% extract the SampleID from the LibraryName (remove the trailing '_shotgun')
SraRunTable.SampleID = strrep(SraRunTable.library_ID, '_shotgun', '');
% remove the columsn library_ID and rename Run to accession
SraRunTable.library_ID = [];
SraRunTable.Properties.VariableNames{1} = 'AccessionShotgun';

%% join the tables
tblASVsamplesUpdatedWithShotgun = outerjoin(tblASVsamples, SraRunTable, 'Type', 'left', 'MergeKeys', true);

%% write the updated table
writetable(tblASVsamplesUpdatedWithShotgun, 'tblASVsamplesUpdatedWithShotgun.csv');

