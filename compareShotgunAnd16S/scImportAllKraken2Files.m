% imports each Kraken2 file for each sample if if doens't exist already
% and saves the table as a .mat file in the same directory

% read the sample table
tblASVsamples = readtable('../tblASVsamplesUpdatedWithShotgun.csv', 'Format', '%s%s%d%s%s%s%d%s');


%% keep only the samples with shotgun sequenced samples
idxHasShotun = cellfun(@(x) ~isempty(x), tblASVsamples.AccessionShotgun);
tblASVsamples(~idxHasShotun, :) = [];

%% make the sample plot to shotgun
addpath('../');
for i = 1:height(tblASVsamples)
    s = tblASVsamples.SampleID{i};
    fn  = sprintf('../PATRIC_output/kraken2/.%s_kraken2/report.txt', s);
    srr = tblASVsamples.AccessionShotgun{i};
    % check if Kraken2 file exist at all
    if isfile(fn)
        % check if matlab table has been precomputed
        if isfile(sprintf('../PATRIC_output/kraken2/.%s_kraken2/tblKraken2.mat', s))
            fprintf('Sample %s has been processed by matlab. skipping...\n', s)
            %load(sprintf('../PATRIC_output/kraken2/.%s_kraken2/tblKraken2.mat', s));
        else
            fprintf('Sample %s has not been processed by matlab yet.\n', s)
            fprintf('Doing it now...\n')
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
