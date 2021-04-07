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
SraRunTable = readtable('tempFiles/SraRunTableAllMskProjects04062021.txt',...
    'FileType', 'text', 'Delimiter', ',');

% keep only the shotgun sequences submitted in 2021 or late 2020
% bacause those anre high quality ones
%idx = strcmp(SraRunTable.library_strategy, 'WGS') & contains(SraRunTable.ReleaseDate, {'2021-' '2020-12'});
idx =  strcmp(SraRunTable.AssayType, 'WGS') &...
     contains(SraRunTable.ReleaseDate, {'2021-' '2020-12'});
SraRunTable(~idx, :) = [];

%% keep only the Libray Name and the Run
SraRunTable = SraRunTable(:, {'LibraryName' 'Run'})

%% make a table with shotgun sequeced samples and save
shotgunSamples = SraRunTable(contains(SraRunTable.LibraryName, '_shotgun'),....
    {'LibraryName' 'Run' 'ReleaseDate'});
shotgunSamples.SampleID = strrep(shotgunSamples.LibraryName, '_shotgun', '');
% remove the columsn library_ID and rename Run to accession
shotgunSamples.LibraryName = [];
shotgunSamples.Properties.VariableNames{1} = 'AccessionShotgun';

% add an extra columns, kneadbatch, to help Jinyuan find out which 
% files were already processed in PATRIC
shotgunSamples.kneadbatch = repmat({'before_JY'}, [height(shotgunSamples), 1]);
shotgunSamples.kneadbatch(strcmp(shotgunSamples.ReleaseDate, '2021-02-12T00:00:00Z')) = {'kneadbatch1'};
shotgunSamples.kneadbatch(strcmp(shotgunSamples.ReleaseDate, '2021-03-29T00:00:00Z')) = {'kneadbatch2'};
shotgunSamples.kneadbatch(strcmp(shotgunSamples.ReleaseDate, '2021-04-05T00:00:00Z')) = {'kneadbatch3'};

% join the tables
tblASVsamplesUpdatedWithShotgun = innerjoin(tblASVsamples, shotgunSamples);
% write the updated table with shotgun sequences
writetable(tblASVsamplesUpdatedWithShotgun, 'tempFiles/tblASVsamplesWithShotgun04072021.csv');

%% make a table with isolate sequences and save
isolateSamples = SraRunTable(contains(SraRunTable.LibraryName, 'isolate'),....
    {'LibraryName' 'Run' 'ReleaseDate'});
s = cellfun(@(x) strsplit(x, '_isolate_'), isolateSamples.LibraryName, 'UniformOutput', false);
isolateSamples.SampleID = cellfun(@(x) x(1), s);
isolateSamples.isolate = cellfun(@(x) x(2), s);
isolateSamples.LibraryName = [];
isolateSamples.Properties.VariableNames{1} = 'AccessionIsolate';
% join the tables
isolateSamples = innerjoin(tblASVsamples, isolateSamples);
% write the updated table with WGS of isolates from those stool samples
writetable(isolateSamples, 'tempFiles/tblIsolateSamples04072021.csv');
