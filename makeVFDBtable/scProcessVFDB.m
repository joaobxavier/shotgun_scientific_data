close all; clear

% read the sample table
tblASVsamples = readtable('../deidentified_data_tables/samples/tblASVsamples.csv', ...
                        'Format', '%s%s%d%s%s%s%d%s');

% pick the patient with most shotgun sequenced samples
idxHasShotgun = cellfun(@(x) ~isempty(x), tblASVsamples.AccessionShotgun);
tblShotgun = tblASVsamples(idxHasShotgun,:);
readcounts = readtable('../metagenome_data/tblShotgunReadcounts.csv', 'Format', '%s%f');
tblShotgun = join(tblShotgun, readcounts, 'Keys', 'SampleID');
processPatricOutput = 0;
%% read VFDB output table and extract accession number

if processPatricOutput == 1
    recordVFDBfile = 'samplesWithoutVFDB.txt';
    if isfile(recordVFDBfile)
        delete(recordVFDBfile)
    end
    VFDBdir = '../PATRIC_output/VFDB';
    dinfo = dir(VFDBdir);
    dinfo = dinfo(contains( {dinfo.name}, {'VFDB'}));
    folders = {dinfo.name};

    VFDBtbl = [];

    for i=1:height(tblShotgun)   
% for i=1
        s = tblShotgun.SampleID{i};
        VFDBfile =sprintf('../PATRIC_output/VFDB/.%s_VFDB/kma.txt', s);
        VFDBfile2 =sprintf('../PATRIC_output/VFDB/.%s_VFDB/kma.res', s);
        if isfile(VFDBfile2)
            system(['mv ' VFDBfile2 ' ' VFDBfile])
        end

        if isfile(VFDBfile)
            %%%% This block of code for deleting files is only for saving
            %%%% storage space and files such as kma.aln can be kept for
            %%%% other analysis
            gz  = sprintf('../PATRIC_output/VFDB/.%s_VFDB/kma.frag.gz', s);
            aln = sprintf('../PATRIC_output/VFDB/.%s_VFDB/kma.aln', s);
            html = sprintf('../PATRIC_output/VFDB/.%s_VFDB/MetagenomicReadMappingReport.html', s);
            fsa = sprintf('../PATRIC_output/VFDB/.%s_VFDB/kma.fsa', s);
            if isfile(gz)
                system(['rm ' gz]); % remove the gz file because it takes too much space and we do not need it
            end
            if isfile(aln)
                system(['rm ' aln]); % remove the gz file because it takes too much space and we do not need it
            end
            if isfile(html)
                system(['rm ' html]); % remove the gz file because it takes too much space and we do not need it
            end
            if isfile(fsa)
                system(['rm ' fsa]); % remove the gz file because it takes too much space and we do not need it
            end
    %         % check if matlab table has been precomputed
            srr = tblShotgun.AccessionShotgun{i};
            log = sprintf('../PATRIC_output/VFDB/%s_VFDB', s);
            [status,result] = system(['cat ' log ' | grep ' srr]);
            if status == 1 || isempty(result)
    %                 error('Sample %s not confirmed by %s', s, log)
                fprintf('Sample %s not matched with SRR in %s\n', s, log)
                fileID = fopen(recordVFDBfile,'a');
                fprintf(fileID,'Sample %s matched with SRR in %s\n', s, log);
                fclose(fileID);
            else
                fprintf('Sample %s confirmed by %s\n', s, log)

                props = dir(VFDBfile);

                if props.bytes<130
                    sprintf('sample %s has empty output', s);
                    fileID = fopen(recordVFDBfile,'a');
                    fprintf(fileID,'sample %s has empty output\n', s);
                    fclose(fileID);
                    continue;
                else
                    tbl = readtable(VFDBfile, 'delimiter', '\t');
%                     tbl = tbl(:, {'x_Template', 'Template_length', 'Depth'});
                    tbl.shotgunReadcount = repmat(tblShotgun.readcount(i), height(tbl),1);
                    tbl2 = processVFDBtable(tbl);
                    tbl2.SampleID = repmat({s}, height(tbl2), 1);
                    VFDBtbl = [VFDBtbl; tbl2];    
                    clear tbl tbl2
                end

            end


        else
            warning('Sample %s has no VFDB output yet...\n', s)
            fprintf('>>>>>> Go PATRIC and request processing %s\n', s);
            fileID = fopen(recordVFDBfile,'a');
            fprintf(fileID,'Sample %s has no VFDB output yet...>>>>>> Go PATRIC and request processing %s\n', s, s);
            fclose(fileID);
        end
    end

    writetable(VFDBtbl, '../metagenome_data/vfdbTbl_2021.csv');
else
    VFDBtblname = '../metagenome_data/vfdbTbl_2021.csv';
    VFDBtbl = readtable(VFDBtblname);
end

