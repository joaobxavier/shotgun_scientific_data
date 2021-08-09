close all; clear

workdir='/Users/yanj2/Documents/GitHub/shotgun_scientific_data_Figure';
cd([workdir '/Figure4'])
% 
CardRabd = readtable('../metagenome_data/CardDrugclassRelAbundance.xlsx');
CardRabd.Properties.VariableNames{1}='SampleID';
samplesWithCard = CardRabd.SampleID;

calculateTSNE = 0;
%% path of tsne
addpath('../utils');

%% path to data
data_path = '../deidentified_data_tables/';

%% load and unstack the counts table
tblcounts = readtable(strcat(data_path, 'counts/tblcounts_asv_melt.csv'));
tblcounts = unstack(tblcounts, 'Count', 'ASV');

%% normalize ASV counts to relative abundance
matcounts = tblcounts{:, 2:end}; % the first column is 'SampleID'
matcounts(isnan(matcounts)) = 0; % replace missing taxa count with 0
matcounts = matcounts ./ sum(matcounts, 2); % convert to relative abundance
tblcounts{:, 2:end} = matcounts;

%% load the taxonomy table
tbltaxonomy = readtable(strcat(data_path, 'taxonomy/tblASVtaxonomy_silva132_v4v5_filter.csv'));

%%
[~,ia,ib]=intersect(CardRabd.SampleID, tblcounts.SampleID);
T = tblcounts;
ln=length(ib);
newT = outerjoin(T, CardRabd, 'Type', 'left','mergeKeys', true);
TM = newT{:, 2:end};
TM(isnan(TM))=0;
newT{:, 2:end}=TM;
clear TM T
idx = max(find(contains(newT.Properties.VariableNames, 'ASV')));
Tasv = newT(:, 1:idx);
Tcard = newT(:, [1, (idx+1):end]);

noshotgunColor=[.3 .6 .6];

Sum = sum(Tcard{:, 2:end}, 2);

%% tSNE using functions downloaded from https://lvdmaaten.github.io/tsne/code/bh_tsne.tar.gz
% IMPORTANT! Run tSNE on the entire dataset takes ~ 30 minutes on my PC (MacBook Pro, 2018, 6 cpu parallelization)
% note: tSNE coefficients vary from run to run
% The output of tSNE is stored so can be called to plot directly

if calculateTSNE==1
    original_dir = pwd;
    software_path = '../utils/';
    fastTSNE_path = strcat(software_path, 'bh_tsne');
    cd(fastTSNE_path);
    system('g++ sptree.cpp tsne.cpp tsne_main.cpp -o bh_tsne -O2');
    fprintf('running tSN;...\n');
    scoreLin = fast_tsne(tblcounts{:, 2:end}, 2, 20, 30, 0.5);
    rng(100)
    scoreLin = fast_tsne(Tasv{:, 2:end}, 2, 20, 30, 0.5)
    cd(original_dir);
    fprintf('tSNE done.\n');
else
    X= load('../savedMat/tSNE.mat');
    scoreLinT=X.scoreLinT;
    clear X
    scoreLin = scoreLinT{:, {'scoreLin1','scoreLin2'}};
end
%% get the dominant taxa of each sample and its rgb color
% each sample in the tSNE plot will be colored by its dominant taxa
[~, dominant_ASV_idx] = max(Tasv{:, 2:end}, [], 2);
ASV_color = hex2rgb(tbltaxonomy(dominant_ASV_idx, :).HexColor);

%% plot tSNE scores in 2 a reduced 2-dimensional space
bgColor = [220,220,220]/255;
figure();
hold on;
cardIdx = find(Sum>0);
scatter(scoreLin(:, 1), scoreLin(:, 2), 100, ASV_color, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
axis square;
box on;
set(gca,'Color',bgColor); % set background color
xlabel('tSNE 1');
ylabel('tSNE 2');

xtick_lb = (floor(min(scoreLin(:, 1))/10)-0.5)*10;
xtick_ub = (ceil(max(scoreLin(:, 1))/10)+0.5)*10;
ytick_lb = (floor(min(scoreLin(:, 2))/10)-0.5)*10;
ytick_ub = (ceil(max(scoreLin(:, 2))/10)+0.5)*10;

axis([xtick_lb,xtick_ub,ytick_lb,ytick_ub]);
set(gca,'XTick',xtick_lb:10:xtick_ub);
set(gca,'YTick',ytick_lb:10:ytick_ub);


%% plot samples with shotgun sequencing
% 
colorPlot = repmat([.9 0 0], height(Tcard), 1);
shotGun = readtable('../metagenome_data/tblASVsamplesUpdatedWithShotgunWithReadcounts_final.xlsx');
shotGunSample = shotGun.SampleID(cellfun(@(X) ~isempty(X), shotGun.AccessionShotgun));
[~, ~, isg]=intersect(shotGunSample, Tcard.SampleID);
G=zeros(height(Tcard),1);
G(isg)=1;
wshotgunColor = [.9 .9 .2];
mycolor = [noshotgunColor; wshotgunColor];
figure()
uG=unique(G, 'stable');
for i=1:length(uG)
scatter(scoreLin(G==uG(i), 1), scoreLin(G==uG(i), 2), 100, mycolor(i,:), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
hold on
% scatter(scoreLin(G==uG(, 1), scoreLin(isg, 2), 100, wshotgunColor, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
end
axis square;
box on;
set(gca,'Color',bgColor); % set background color
xlabel('tSNE 1');
ylabel('tSNE 2');
legend({'no shotgun' 'shotgun'})
xtick_lb = (floor(min(scoreLin(:, 1))/10)-0.5)*10;
xtick_ub = (ceil(max(scoreLin(:, 1))/10)+0.5)*10;
ytick_lb = (floor(min(scoreLin(:, 2))/10)-0.5)*10;
ytick_ub = (ceil(max(scoreLin(:, 2))/10)+0.5)*10;

axis([xtick_lb,xtick_ub,ytick_lb,ytick_ub]);
set(gca,'XTick',xtick_lb:10:xtick_ub);
set(gca,'YTick',ytick_lb:10:ytick_ub);
title('Samples with shotgun sequencing has diverse taxsonomy composition', 'FontSize', 20)
%%
if calculateTSNE==1
    scoreLinT = array2table(scoreLin);
    scoreLinT.SampleID=Tasv.SampleID;
    save('../savedMat/tSNE.mat', 'scoreLinT')
end

%% stool consistency vs log read counts
grouporder = {'formed', 'semi-formed', 'liquid'};
T1 = shotGun(:, {'Consistency' 'readcount'})
stoolC=cellstr(shotGun.Consistency);
figure
vs = violinplot(log10(T1.readcount), stoolC, ...
                'Bandwidth', 0.2, ...
                'ViolinAlpha', 0.35, ...
                'ShowNotches', false, ...
                'ShowMean', true, ...
                'GroupOrder', grouporder);
vs(1).MedianPlot.MarkerFaceColor = [0 0 0];
vs(2).MedianPlot.MarkerFaceColor = [0 0 0];
vs(3).MedianPlot.MarkerFaceColor = [0 0 0];
vs(1).MeanPlot.LineWidth=3;
vs(2).MeanPlot.LineWidth=3;
vs(3).MeanPlot.LineWidth=3;
% vs(1).WhiskerPlot.MarkerFaceColor = [ 0 0 0];
% vs(1).WhiskerPlot.MarkerFaceColor = [ 0 0 0];
vs(1).ViolinColor = [.2 0.8 .3];
vs(2).ViolinColor = [.5 0.5 1];
vs(3).ViolinColor = [.8 0.3 0.2];
vs(1).WhiskerPlot.Color = [.2 .2 .2];
vs(1).ScatterPlot.SizeData=50;
vs(2).ScatterPlot.SizeData=50;
vs(3).ScatterPlot.SizeData=50;
xlabel('Stool Consistancy');
ylabel('log10(counts in shotgun samples)');
set(gca, 'xticklabel', ...
    {sprintf('formed (n=%i)', sum(ismember(T1.Consistency, {'formed'}))), ...
            sprintf('semi-formed (n=%i)', sum(ismember(T1.Consistency, {'semi-formed'}))), ...
            sprintf('liquid (n=%i)', sum(ismember(T1.Consistency, {'liquid'})))}, ...
    'Fontsize',14);
xlim([0.2, 3.8]);
title('shotgun sequencing reads in different stool consistency')