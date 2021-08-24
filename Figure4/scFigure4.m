close all; clear

data_path = '../deidentified_data_tables/';
addpath('../utils')
tblcounts = readtable(strcat(data_path, 'counts/tblcounts_asv_melt.csv'));
tblcounts = unstack(tblcounts, 'Count', 'ASV');

cardT = readtable('../metagenome_data/cardTbl.csv');
cardT = cardT(:, {'resistGene' ,'RelavantPercentInCARD', 'SampleID'});
cardT.Properties.VariableNames = {'resistGene', 'RelavantPercentage', 'SampleID'};
vanApcr = readtable(strcat(data_path,'/meta_data/tblVanA.csv'));
% X = outerjoin(vanA, cardT, 'Type','left','MergeKeys',true);
%% plot vanA PCR in tSNE
%% plot samples with vanA PCR
% shotGun = readtable('../deidentified_data_tables/samples/tblASVsamples.csv');
shotGun = readtable(strcat(data_path, 'samples/tblASVsamples.csv'), ...
                        'Format', '%s%s%d%s%s%s%d%s');
shotGunSample = shotGun.SampleID(cellfun(@(X) ~isempty(X), shotGun.AccessionShotgun));
X= load('../savedMat/tSNE.mat');
scoreLinT=X.scoreLinT;
clear X
scoreLin = scoreLinT{:, {'scoreLin1','scoreLin2'}};

vanA_pcr_pos = intersect(shotGunSample, vanApcr.SampleID(vanApcr.VanA==1));
vanA_pcr_neg = intersect(shotGunSample, vanApcr.SampleID(vanApcr.VanA==0));

[~, ~, ivanA_pos]=intersect(vanA_pcr_pos, scoreLinT.SampleID);
[~, ~, ivanA_neg]=intersect(vanA_pcr_neg, scoreLinT.SampleID);
inotchecked=setdiff(1:height(scoreLinT), [ivanA_pos; ivanA_neg]);

bgColor = [240,240,240]/255;
noVanAcolor=[0.05 .7  1];
vanAcolor = [.95 .1 .7];
notOverlap = [.8 .8 .8];
dotsize = 100;

figure(100)
scatter(scoreLin(inotchecked, 1), scoreLin(inotchecked, 2), dotsize, notOverlap, 'filled', 'MarkerEdgeColor', notOverlap*0.9, 'LineWidth', 1);
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
title('Samples with vanA PCR and shotgun sequencing', 'FontSize', 20)

%% compare vanA and vanB
vanAcomp = cardT(strcmp(cardT.resistGene, 'vanA'), :);
vanBcomp = cardT(strcmp(cardT.resistGene, 'vanB'), :);

vanAB = outerjoin(vanAcomp(:, 2:3), vanBcomp(:, 2:3), 'Keys', 'SampleID');
vanAB.Properties.VariableNames = {'vanA' 'SampleID' 'vanB' 's'};
vanAB.SampleID(cellfun(@(X) isempty(X), vanAB.SampleID)) = vanAB.s(cellfun(@(X) isempty(X), vanAB.SampleID));
vanAB.s =[];
vanAB.vanA(isnan(vanAB.vanA))=0;
vanAB.vanB(isnan(vanAB.vanB))=0;
[R, p] = corrcoef(vanAB.vanA, vanAB.vanB);
corr(vanAB.vanA, vanAB.vanB, 'Type','Spearman')

%%
g = zeros(height(vanAB),1);
g(vanAB.vanA>0 & vanAB.vanB>0)=1;
cmap = [.5 .4 .45;
    .3 .6 .3];
figure
plot(vanAB.vanA(g==0) * 100, vanAB.vanB(g==0) *100, 'o', ...
    'MarkerEdgeColor', cmap(1, :), 'MarkerFaceColor', cmap(1, :),'MarkerSize', 10)
hold on
plot(vanAB.vanA(g==1) * 100, vanAB.vanB(g==1) *100, 'o', ...
    'MarkerEdgeColor', cmap(2, :), 'MarkerFaceColor', cmap(2, :),'MarkerSize', 10)
hold off
xlabel('vanA abundance (%)', 'fontsize', 14)
ylabel('vanB abundance (%)', 'fontsize', 14)
legend({'only vanA or vanB detected' 'both vanA and vanB detected'}, 'fontsize', 12)
%%
vanAcompCARD=grpstats(vanAcomp, 'SampleID', 'sum', 'DataVar', 'RelavantPercentage');
vanAcompCARD.GroupCount=[];

noVanAsample = setdiff(unique(cardT.SampleID), vanAcompCARD.SampleID);
noVanA = array2table(noVanAsample, 'VariableNames', {'SampleID'});
noVanA.sum_RelavantPercentage = zeros( height(noVanA),1);
vanAcompCARD =[vanAcompCARD; noVanA];

clear vanAcomp noVanAsample noVanA
%% find samples with both CARD mapping and PCR
X2=innerjoin(vanAcompCARD, vanApcr, 'Keys', 'SampleID');
X2.sum_RelavantPercentage = X2.sum_RelavantPercentage*100;
%% ranksum test
[p,h,stats] = ranksum(X2.sum_RelavantPercentage(X2.VanA==0), ...
                         X2.sum_RelavantPercentage(X2.VanA==1))


%% beeswarm plot

figure
n1=sum(X2.VanA==1);
n0=sum(X2.VanA==0);
bs=beeswarm( X2.VanA,X2.sum_RelavantPercentage, ...
            'colormap', [noVanAcolor; vanAcolor], ...
         'sort_style','up','dot_size',4,'corral_style', 'gutter', 'overlay_style','ci');
set(gca, 'ylim', [-1 11], 'ytick', 0:2:10, ...
        'xtick', [0 1], ...
        'xticklabel', {'vanA PCR(-)' ,'vanA PCR(+)'})
ylabel('relative abundance of vanA in CARD', 'fontsize', 14) 
text(-0.1, 10.5, sprintf('n=%i', n0),'fontsize', 14)
text(0.9, 10.5, sprintf('n=%i', n1),'fontsize', 14)
