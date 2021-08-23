close all; clear

% read the sample table
tblASVsamples = readtable('../deidentified_data_tables/samples/tblASVsamples.csv', ...
                        'Format', '%s%s%d%s%s%s%d%s');

% pick the patient with most shotgun sequenced samples
idxHasShotgun = cellfun(@(x) ~isempty(x), tblASVsamples.AccessionShotgun);
tblShotgun = tblASVsamples(idxHasShotgun,:);
readcounts = readtable('../metagenome_data/tblShotgunReadcounts.csv', 'Format', '%s%f');
tblShotgun = join(tblShotgun, readcounts, 'Keys', 'SampleID');

processPatricOutput = 0;  % if set as 1, then re-process CARD output
%% read card output table and extract accession number
if processPatricOutput == 1
    recordCARDfile = 'samplesWithoutCARD.txt';
    if isfile(recordCARDfile)
        delete(recordCARDfile)
    end
    carddir = '../PATRIC_output/CARD';
    dinfo = dir(carddir);
    dinfo = dinfo(contains( {dinfo.name}, {'CARD'}));
    folders = {dinfo.name};
    sizefactor = ones(length(folders),1);
    cardTbl = [];
    for i=1:height(tblShotgun)    
        s = tblShotgun.SampleID{i};
        cardfile =sprintf('../PATRIC_output/CARD/.%s_CARD/kma.txt', s);
        cardfile2 =sprintf('../PATRIC_output/CARD/.%s_CARD/kma.res', s);
        if isfile(cardfile2)
            system(['mv ' cardfile2 ' ' cardfile])
        end

        if isfile(cardfile)
            gz  = sprintf('../PATRIC_output/CARD/.%s_CARD/kma.frag.gz', s);
            if isfile(gz)
                system(['rm ' gz]); % remove the gz file because it takes too much space and we do not need it
            end

            srr = tblShotgun.AccessionShotgun{i};
            log = sprintf('../PATRIC_output/CARD/%s_CARD', s);
            [status,result] = system(['cat ' log ' | grep ' srr]);
            if status == 1 || isempty(result)
                fprintf('Sample %s not matched with SRR in %s\n', s, log)
                fileID = fopen(recordCARDfile,'a');
                fprintf(fileID,'Sample %s matched with SRR in %s\n', s, log);
                fclose(fileID);
            else
                fprintf('Sample %s confirmed by %s\n', s, log)

                props = dir(cardfile);

                if props.bytes<130
                    sprintf('sample %s has empty output', s);
                    fileID = fopen(recordCARDfile,'a');
                    fprintf(fileID,'sample %s has empty output\n', s);
                    fclose(fileID);
                    continue;
                else
                    tbl = readtable(cardfile, 'delimiter', '\t');
                    tbl.shotgunReadcount = repmat(tblShotgun.readcount(i), height(tbl),1);
                    tbl2 = processCARDtable(tbl);
                    tbl2.SampleID = repmat({s}, height(tbl2), 1);
                    cardTbl = [cardTbl; tbl2];    
                    clear tbl tbl2
                end

            end


        else
            warning('Sample %s has no CARD output yet...\n', s)
            fprintf('>>>>>> Go PATRIC and request processing %s\n', s);
            fileID = fopen(recordCARDfile,'a');
            fprintf(fileID,'Sample %s has no CARD output yet...>>>>>> Go PATRIC and request processing %s\n', s, s);
            fclose(fileID);
        end
    end
    cardTbl.("N/A") =[];
    writetable(cardTbl, '../metagenome_data/cardTbl.csv');
else
    mycardTbl = '../metagenome_data/cardTbl.csv';
    cardTbl = readtable(mycardTbl);
end




