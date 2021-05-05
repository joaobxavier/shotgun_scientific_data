close all; clear

% workdir='/Users/yanj2/Documents/GitHub/shotgun_scientific_data';
% cd([workdir '/CARDtSNE'])

% wildcard = readtable('index-for-model-sequences_wildcard_3.0.8.txt', 'FileType', 'text');

% read the sample table
tblASVsamples = readtable('../metagenome_data/tblASVsamplesUpdatedWithShotgunWithReadcounts_final.xlsx');

% pick the patient with most shotgun sequenced samples
idxHasShotgun = cellfun(@(x) ~isempty(x), tblASVsamples.AccessionShotgun);
tblShotgun = tblASVsamples(idxHasShotgun,:);

processPatricOutput = 2;
%% read card output table and extract accession number
% groupcounts(card, 'DrugClass')
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
    % fn = 'kma.res';
    % fn2='kma.txt';
    cardTbl = [];
    % f=folders;
    % f=cellfun(@(X) strip(X,'left','.'), f,'UniformOutput', false);
    % f=cellfun(@(X) strsplit(X, '_'), f, 'UniformOutput', false);
    % f=cellfun(@(X) X{1}, f, 'UniformOutput', false);
    for i=1:height(tblShotgun)    
        s = tblShotgun.SampleID{i};
        cardfile =sprintf('../PATRIC_output/CARD/.%s_CARD/kma.txt', s);
        cardfile2 =sprintf('../PATRIC_output/CARD/.%s_CARD/kma.res', s);
        if isfile(cardfile2)
            system(['mv ' cardfile2 ' ' cardfile])
        end


    %     fn  = sprintf('../PATRIC_output/CARD/.%s_CARD/report.txt', s);
    %     srr = tblShotgun.AccessionShotgun{strcmp(tblShotgun.SampleID, s)};

        if isfile(cardfile)
            gz  = sprintf('../PATRIC_output/CARD/.%s_CARD/kma.frag.gz', s);
            if isfile(gz)
                system(['rm ' gz]); % remove the gz file because it takes too much space and we do not need it
            end
    %         % check if matlab table has been precomputed
    %         if isfile(sprintf('../PATRIC_output/CARD/.%s_CARD/tblKraken2.mat', s))
    %             fprintf('Sample %s has been processed by matlab. Loading...\n', s)
    %             load(sprintf('../PATRIC_output/CARD/.%s_CARD/tblKraken2.mat', s));
    %         else
    %             fprintf('Sample %s has not been processed by matlab yet.\n', s)
    %             fprintf('Doing it now...\n')
            srr = tblShotgun.AccessionShotgun{strcmp(tblShotgun.SampleID, s)};
            log = sprintf('../PATRIC_output/CARD/%s_CARD', s);
            [status,result] = system(['cat ' log ' | grep ' srr]);
            if status == 1 || isempty(result)
    %                 error('Sample %s not confirmed by %s', s, log)
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
                    tbl = tbl(:, {'x_Template', 'Template_length', 'Depth'});
                    tbl.sizefactor = repmat(tblShotgun.readcount(strcmp(tblShotgun.SampleID, s)), height(tbl),1);
                    tbl2 = processCARDtable(tbl);
                    tbl2.sample = repmat({s}, height(tbl2), 1);
                    cardTbl = [cardTbl; tbl2];    
                    clear tbl tbl2
                end

            end
    %             % load the table for this file
    %             tblKraken2 = importKraken2Output(fn);
    %             save(sprintf('../PATRIC_output/kraken2/.%s_kraken2/tblKraken2.mat',...
    %                 s), 'tblKraken2');
    %         end
            % compute relative abundances
            % add all the fragments mapped to a genus and show the result
            % ordered by the most adundant in sample

        else
            warning('Sample %s has no CARD output yet...\n', s)
            fprintf('>>>>>> Go PATRIC and request processing %s\n', s);
            fileID = fopen(recordCARDfile,'a');
            fprintf(fileID,'Sample %s has no CARD output yet...>>>>>> Go PATRIC and request processing %s\n', s, s);
            fclose(fileID);
        end
    end

    writetable(cardTbl, 'cardTbl_2021April29.xlsx');
else
    mycardTbl = 'cardTbl_2021April29.xlsx';
    cardTbl = readtable(mycardTbl);
end
%%
idx1 = find(ismember(cardTbl.Properties.VariableNames, 'N/A'));
idx2 = find(ismember(cardTbl.Properties.VariableNames, 'unknown'));
drugs = readtable('DrugClass.xlsx');
allDrugs = cardTbl.Properties.VariableNames(idx1:idx2);
C = cardTbl{:, idx1:idx2};

%% If a gene has more than one drug class
drugs.levels = ones(height(drugs),1);
for i=1:length(drugs.DrugClass)
    idx = find(ismember(cardTbl.Properties.VariableNames, drugs.DrugClass(i)));
    d1 = cardTbl{:, idx};
    if sum(d1)>0
        for j=1:size(C,2)
            drug1 = drugs.DrugClass{i};
            drug2 = allDrugs{j};
            if ~strcmp(drug1, drug2)
                a = d1 .* C(:, j); 
                if sum(a) >0  % if the two drugs overlap in drug class                   
                    sprintf('%s and %s overlaps', drug1, drug2);
                    if  sum(d1 - C(:, j))==0
                        sprintf('%s and %s are the same', drug2, drug1)
                    elseif sum(a-d1)==0 
                        drugs.levels(ismember(drugs.DrugClass, drug2))=drugs.levels(ismember(drugs.DrugClass, drug2))+1;
                        sprintf('%s contains all %s', drug2, drug1);
                    elseif sum(a-C(:, j))==0
                        drugs.levels(ismember(drugs.DrugClass, drug1))=drugs.levels(ismember(drugs.DrugClass, drug1))+1;
                        sprintf('%s contains all %s', drug1, drug2);
                    end
                end
            end
        end
    end
    
    
end



% d2 = cellfun(@(X) strrep(X, ' ' , ''), d2, 'UniformOutput', false);
% setdiff(d2, drugs.DrugClass)
    
    
  %%  
    
% for i=1:3 
%     if isfile([carddir '/' folders{i} '/' fn])
%         movefile([carddir '/' folders{i} '/' fn], [carddir '/' folders{i} '/' fn2]);
%     end
%     s=dir([carddir '/' folders{i} '/' fn2]);
%     if ~isempty(s)
%         if s.bytes<130
%             sprintf('sample %s has empty output', f{i});
%             continue;
%         else
% 
%             tbl = readtable([carddir '/' folders{i} '/' fn2], 'delimiter', '\t');
%             tbl = tbl(:, {'x_Template', 'Template_length', 'Depth'});
%             tbl.sizefactor = repmat(sizefactor(1), height(tbl),1);
%             tbl2 = processCARDtable(tbl);
%             tbl2.sample = repmat(f(i), height(tbl2), 1);
%             cardTbl = [cardTbl; tbl2];    
%             clear tbl tbl2
%         end
%     else
%         warning('%s has no CARD!!!', folders{i})
%     end

% writetable(cardTbl, 'cardCompOfMetagenomeShortgun.xlsx');

%% remove drugs not found in our data set
Z=table2array(cardTbl(:, idx1:idx2));
% Idx = find(sum(Z) == 0);
% cardTbl(:, Idx+idx1-1)=[];
clear Z Idx
%%
% idx2 = find(ismember(cardTbl.Properties.VariableNames, 'unknown'));
drugs = cardTbl.Properties.VariableNames(idx1:idx2);


%%
% CardRAbd = cell2table(f', 'VariableNames', {'sample'});
% CardRAbd = cardTbl(:, {'sample'});
CardRAbd = groupcounts(cardTbl, 'sample');
for i=1:length(drugs)
% for i=1:5
    C=cardTbl{:, {drugs{i}}};
%     Cg = grpstats(C, 'sample', 'sum', 'RelavantPercentage');
    Idx = find(C==1);
    C2=cardTbl(Idx, {'sample', drugs{i}, 'RelavantPercentage'});
    G = grpstats(C, {'sample' drugs{i}}, 'sum', 'DataVars', 'RelavantPercentage');
    G = G(:, {'sample', 'sum_RelavantPercentage'});
    G.Properties.VariableNames = {'sample', drugs{i}};
    CardRAbd = outerjoin(CardRAbd, G, 'Type', 'left','mergeKeys', true);
    clear C G
%     CardRAbd.("sample_G")=[];
% G(G.(drugs{i})==0, :) = [];
end
%%

C=table2array(CardRAbd(:,2:end));
C(isnan(C))=0;
CardRAbd{:,2:end }= C;
writetable(CardRAbd, 'CardDrugclassRelAbundance.xlsx')

%% create a table for antibiotics coloring
% %% 
% T=readtable('cardCompOfMetagenomeShortgun.xlsx');
% d = T.Properties.VariableNames;
% d = d(16:end-1);
% 
% %%
% D = array2table(d', 'VariableNames', {'DrugClass'});
% writetable(D, 'DrugClass.xlsx')


