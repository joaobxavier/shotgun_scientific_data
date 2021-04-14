%% download the latest metadata on all runs submitted to the SRA
% beleonging to MSKCC allo-HCT projects
% 1. nagigate to https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA394877%2CPRJNA491657%2CPRJNA479462%2CPRJNA545312%2CPRJNA548153%2CPRJNA576440%2CPRJNA607574%2CPRJNA606262%2CPRJNA701319%20&o=acc_s%3Aa
% 2. click on 'Metadata' to donwload a file called 'SraRunTable.txt'
% 3. save the file 'SraRunTable.txt' to 'tempFiles'
SraRunTable = readtable('tempFiles/SraRunTableAllMskProjects04062021.txt',...
    'FileType', 'text', 'Delimiter', ',');

%%
gc = groupcounts(SraRunTable, {'BioProject' 'AssayType'});
gc = unstack(gc, 'GroupCount', 'AssayType');

%%
idx = contains(SraRunTable.LibraryName, '_VRE');
SraRunTable.LibraryName(idx)
unique(SraRunTable.BioSample(idx))
%% extract the patientID
libs = SraRunTable.LibraryName(idx);
samps = cellfun(@(x) strsplit(x, '_'), libs, 'UniformOutput', false);
pts = cellfun(@(x) x{1}(1:end-1), samps, 'UniformOutput', false);
unique(pts)