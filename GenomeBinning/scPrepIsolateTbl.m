close all; clear

sampletbl =readtable('metadata-9426271-processed-ok.tsv', ...
                    'FileType','text', ...
                        'delimiter', '\t', 'readVariableNames', true);

sampletbl = sampletbl(:,{'accession' 'library_ID' 'title'});
sampletbl.Properties.VariableNames = {'Isolate_accession' 'isolates' 'SampleID'};
%%
[sampleWIsolates, ia, ic] = unique(sampletbl.SampleID);
Isolates = repmat({''}, length(sampleWIsolates),1);
for i=1:length(sampleWIsolates)
    a = sampletbl.Isolate_accession(ic==i);
    b = sampletbl.isolates(ic==i);
    c = cellfun(@(X, Y) strcat(X, ' (', Y, ')'), a, b, 'UniformOutput', false);
    C = strjoin(c, ';');
    Isolates{i} = C;
end

T2 = cell2table([sampleWIsolates Isolates], 'VariableNames', {'SampleID' 'SRA_for_Isolates'});

%%
shotguntbl = readtable('../metagenome_data/tblASVsamplesWithShotgun04072021.csv', ...
                        'readVariableNames', true);

shotguntbl = shotguntbl(:, {'SampleID' 'AccessionShotgun'});

%%
isolates = innerjoin(shotguntbl, T2,'Keys', 'SampleID')
writetable(isolates, 'VREisolates_with_shotgun_2021May1.xlsx');