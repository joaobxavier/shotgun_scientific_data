close all; clear


CardRabd = readtable('CardDrugclassRelAbundance.xlsx');
CardRabd.Properties.VariableNames{1}='SampleID';
samplesWithCard = CardRabd.SampleID;
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
% c=1:100;
% T = tblcounts([ib; c'], 1:200);
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


% %% colormap for antibiotics
% 
% % load drug table
% drugClass = readtable('DrugClass.xlsx');
% [maxCardDominance, dominant_abx_idx] = max(Tcard{:, 2:end}, [], 2);
% [reorder, Idx]=sort(maxCardDominance);
% 
% Tcard = Tcard(Idx, :);
% Tasv = Tasv(Idx, :);
% dominant_abx_idx = dominant_abx_idx(Idx);
Sum = sum(Tcard{:, 2:end}, 2);
% %
% % abx_color = hex2rgb(DrugClass(dominant_abx_idx, :).HexColor);
% 
% abx_color = zeros(height(Tcard), 3);
% Abx = Tcard.Properties.VariableNames(2:end);
% drugForLabel = repmat({''}, height(Tcard),1);
% for i=1:length(dominant_abx_idx)
%     if Sum(i)==0
%         abx_color(i,:) = hex2rgb(drugClass(strcmp(drugClass.DrugClass , 'noShortgun'),:).HexColor);
%         drugForLabel{i} = 'noShortgun';
%     else
%         abx_color(i,:) = hex2rgb(drugClass(strcmp(drugClass.DrugClass , Abx{dominant_abx_idx(i)}),:).HexColor);
%         drugForLabel(i)=drugClass(strcmp(drugClass.DrugClass , Abx{dominant_abx_idx(i)}),:).DrugClass;
%     end    
% end
%% tSNE using functions downloaded from https://lvdmaaten.github.io/tsne/code/bh_tsne.tar.gz
% IMPORTANT! Run tSNE on the entire dataset takes ~ 30 minutes on my PC (MacBook Pro, 2018, 6 cpu parallelization)
% note: tSNE coefficients vary from run to run
original_dir = pwd;
software_path = '../utils/';
fastTSNE_path = strcat(software_path, 'bh_tsne');
cd(fastTSNE_path);
system('g++ sptree.cpp tsne.cpp tsne_main.cpp -o bh_tsne -O2');
fprintf('running tSN;...\n');
% scoreLin = fast_tsne(tblcounts{:, 2:end}, 2, 20, 30, 0.5);
rng(100)
scoreLin = fast_tsne(Tasv{:, 2:end}, 2, 20, 30, 0.5)
cd(original_dir);
fprintf('tSNE done.\n');

%% get the dominant taxa of each sample and its rgb color
% each sample in the tSNE plot will be colored by its dominant taxa
[~, dominant_ASV_idx] = max(Tasv{:, 2:end}, [], 2);
ASV_color = hex2rgb(tbltaxonomy(dominant_ASV_idx, :).HexColor);

%% plot tSNE scores in 2 a reduced 2-dimensional space
bgColor = [220,220,225]/255;
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



% %% plot tSNE scores in 2 a reduced 2-dimensional space with dominant antibiotics color
% figure();
% hold on;
% cardIdx = find(Sum>0);
% % scatter(scoreLin(:, 1), scoreLin(:, 2), 100, ASV_color, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
% scatter(scoreLin(:, 1), scoreLin(:, 2), 100, abx_color, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
% scatter(scoreLin(cardIdx, 1), scoreLin(cardIdx, 2), 100, abx_color(cardIdx, :), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
% axis square;
% box on;
% set(gca,'Color',bgColor); % set background color
% xlabel('tSNE 1');
% ylabel('tSNE 2');
% 
% xtick_lb = (floor(min(scoreLin(:, 1))/10)-0.5)*10;
% xtick_ub = (ceil(max(scoreLin(:, 1))/10)+0.5)*10;
% ytick_lb = (floor(min(scoreLin(:, 2))/10)-0.5)*10;
% ytick_ub = (ceil(max(scoreLin(:, 2))/10)+0.5)*10;
% 
% axis([xtick_lb,xtick_ub,ytick_lb,ytick_ub]);
% set(gca,'XTick',xtick_lb:10:xtick_ub);
% set(gca,'YTick',ytick_lb:10:ytick_ub);
% %% make legend
% [udrug, iu]= unique(drugForLabel, 'stable');
% l = length(udrug);
% figure
% for i=1:l
%     scatter(1, i, 400, abx_color(iu(i),:), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
%     hold on
%     text(1.2, i, udrug{i}, 'fontsize', 20)
%     hold on 
% end
% hold off
% set(gca, 'xlim', [0.8 1.8])
% axis off
% %% plot distribution of single antibiotic
% % 
% % colorPlot = repmat([.9 0 0], height(Tcard), 1);
% noshotgunColor=[.3 .6 .6];
% 
% for i=2:length(udrug)
% % for i=2    
%     cardIdx = find(Tcard{:, ismember(Tcard.Properties.VariableNames, udrug{i})} >0);
%     colorPlotX= 1-[0*Tcard{cardIdx, ismember(Tcard.Properties.VariableNames, udrug{i})}  ...
%                     1*Tcard{cardIdx, ismember(Tcard.Properties.VariableNames, udrug{i})} ...
%                     1*Tcard{cardIdx, ismember(Tcard.Properties.VariableNames, udrug{i})}];
%     figure(i+10)
%     scatter(scoreLin(:, 1), scoreLin(:, 2), 100, noshotgunColor, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
%     hold on
%     scatter(scoreLin(cardIdx, 1), scoreLin(cardIdx, 2), 100, colorPlotX, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
%     axis square;
%     box on;
%     set(gca,'Color',bgColor); % set background color
%     xlabel('tSNE 1');
%     ylabel('tSNE 2');
% 
%     xtick_lb = (floor(min(scoreLin(:, 1))/10)-0.5)*10;
%     xtick_ub = (ceil(max(scoreLin(:, 1))/10)+0.5)*10;
%     ytick_lb = (floor(min(scoreLin(:, 2))/10)-0.5)*10;
%     ytick_ub = (ceil(max(scoreLin(:, 2))/10)+0.5)*10;
% 
%     axis([xtick_lb,xtick_ub,ytick_lb,ytick_ub]);
%     set(gca,'XTick',xtick_lb:10:xtick_ub);
%     set(gca,'YTick',ytick_lb:10:ytick_ub);
%     title(udrug{i}, 'FontSize', 20)
% end

%%
%% plot samples with shotgun sequencing
% 
% colorPlot = repmat([.9 0 0], height(Tcard), 1);
shotGun = readtable('../metagenome_data/tblASVsamplesUpdatedWithShotgunWithReadcounts_final.xlsx');
shotGunSample = shotGun.SampleID(cellfun(@(X) ~isempty(X), shotGun.AccessionShotgun));
[~, ~, isg]=intersect(shotGunSample, Tcard.SampleID);
G=zeros(height(Tcard),1);
G(isg)=1;
noshotgunColor=[.4 .5 .6];
wshotgunColor = [.9 .4 .1];
mycolor = [noshotgunColor; wshotgunColor];
figure(99)
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
%% plot samples with vanA PCR
vanA = readtable('../metagenome_data/tblVanA.csv');
vanA_shotgun_pos = intersect(CardRabd.SampleID, vanA.SampleID(vanA.VanA==1));
vanA_shotgun_neg = intersect(CardRabd.SampleID, vanA.SampleID(vanA.VanA==0));
% G2=zeros(height(Tcard),1);
% [~, ichecked]=intersect(Tcard.SampleID, vanA.SampleID);
[~, ~, ivanA_pos]=intersect(vanA_shotgun_pos, Tcard.SampleID);
[~, ~, ivanA_neg]=intersect(vanA_shotgun_neg, Tcard.SampleID);
inotchecked=setdiff(1:length(G2), [ivanA_pos; ivanA_neg]);
% G2(ivanA_neg) = 1;
% G2(ivanA_pos)=2;

%%
notCheckedcol = [.7 .7 .7];
noVanAcolor=[0.05 .7  1];
vanAcolor = [.95 .1 .7];
dotsize = 100;
fg = figure(100)
scatter(scoreLin(inotchecked, 1), scoreLin(inotchecked, 2), dotsize, notCheckedcol, 'filled', 'MarkerEdgeColor', notCheckedcol*0.4, 'LineWidth', 1);
hold on
scatter(scoreLin(ivanA_neg, 1), scoreLin(ivanA_neg, 2), dotsize, noVanAcolor, 'filled', 'MarkerEdgeColor', noVanAcolor*0.4, 'LineWidth', 1);
hold on
scatter(scoreLin(ivanA_pos, 1), scoreLin(ivanA_pos, 2), dotsize, vanAcolor, 'filled', 'MarkerEdgeColor', vanAcolor*0.4, 'LineWidth', 1);
axis square;
box on;
set(gca,'Color',bgColor); % set background color
xlabel('tSNE 1');
ylabel('tSNE 2');
legend({'not checked' 'vanA PCR (-)' 'vanA PCR (+)'})
xtick_lb = (floor(min(scoreLin(:, 1))/10)-0.5)*10;
xtick_ub = (ceil(max(scoreLin(:, 1))/10)+0.5)*10;
ytick_lb = (floor(min(scoreLin(:, 2))/10)-0.5)*10;
ytick_ub = (ceil(max(scoreLin(:, 2))/10)+0.5)*10;

axis([xtick_lb,xtick_ub,ytick_lb,ytick_ub]);
set(gca,'XTick',xtick_lb:10:xtick_ub);
set(gca,'YTick',ytick_lb:10:ytick_ub);
title('Samples with positive vanA PCR', 'FontSize', 20)
fig.Renderer='Painters';
% %% plot vanA gene PCR only in shotgun samples
% G2(G==0) = -1;
% notOverlap = [.8 .8 .8];
% noVanAcolor=[.2 .3 .8];
% vanAcolor = [.8 .5 .3];
% dotsize = 100;
% figure(100)
% scatter(scoreLin(G2==-1, 1), scoreLin(G2==-1, 2), dotsize, notOverlap, 'filled', 'MarkerEdgeColor', notOverlap*0.9, 'LineWidth', 1);
% hold on
% scatter(scoreLin(G2==0, 1), scoreLin(G2==0, 2), dotsize, noVanAcolor, 'filled', 'MarkerEdgeColor', noVanAcolor*0.4, 'LineWidth', 1);
% hold on
% 
% scatter(scoreLin(G2==1, 1), scoreLin(G2==1, 2), dotsize, vanAcolor, 'filled', 'MarkerEdgeColor', vanAcolor*0.4, 'LineWidth', 1);
% axis square;
% box on;
% set(gca,'Color',bgColor); % set background color
% xlabel('tSNE 1');
% ylabel('tSNE 2');
% legend({'no shotgun' 'vanA PCR (-)' 'vanA PCR (+)'})
% xtick_lb = (floor(min(scoreLin(:, 1))/10)-0.5)*10;
% xtick_ub = (ceil(max(scoreLin(:, 1))/10)+0.5)*10;
% ytick_lb = (floor(min(scoreLin(:, 2))/10)-0.5)*10;
% ytick_ub = (ceil(max(scoreLin(:, 2))/10)+0.5)*10;
% 
% axis([xtick_lb,xtick_ub,ytick_lb,ytick_ub]);
% set(gca,'XTick',xtick_lb:10:xtick_ub);
% set(gca,'YTick',ytick_lb:10:ytick_ub);
% title('Samples with vanA PCR and shotgun sequencing', 'FontSize', 20)

%%
save('2021April28', '-v7.3');

