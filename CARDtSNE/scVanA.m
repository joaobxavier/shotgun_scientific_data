close all; clear

% workdir='/Users/yanj2/Documents/GitHub/shotgun_scientific_data';
% cd('/CARDtSNE')
% cardT = readtable('cardCompOfMetagenomeShortgun.xlsx');
cardT = readtable('cardTbl_2021April29.xlsx');
cardT = cardT(:, {'resistGene' ,'RelavantPercentage', 'sample'});
cardT.Properties.VariableNames = {'resistGene', 'RelavantPercentage', 'SampleID'};
data_path = '../deidentified_data_tables/';
vanApcr = readtable('../metagenome_data/tblVanA.csv');

% X = outerjoin(vanA, cardT, 'Type','left','MergeKeys',true);

%%
vanAcomp = cardT(strcmp(cardT.resistGene, 'vanA'), :);
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

%% boxplot
figure
h=boxplot(X2.sum_RelavantPercentage, X2.VanA,'Labels',{'vanA(-)','vanA(+)'},'Whisker',2)
ylabel('vanA in metagenomic mapping', 'fontsize', 12)
set(h, 'linewidth', 2)

%% beeswarm plot
noVanAcolor=[0.05 .7  1];
vanAcolor = [.95 .1 .7];
figure
n1=sum(X2.VanA==1);
n0=sum(X2.VanA==0);
bs=beeswarm( X2.VanA,X2.sum_RelavantPercentage, ...
            'colormap', [noVanAcolor; vanAcolor], ...
         'sort_style','up','dot_size',4,'corral_style', 'gutter', 'overlay_style','ci');
set(gca, 'ylim', [-1 11], 'ytick', 0:2:10, ...
        'xtick', [0 1], ...
        'xticklabel', {'vanA(-)' ,'vanA(+)'})
ylabel('relative abundance of vanA gene in shotgun metagenome', 'fontsize', 14) 
text(-0.1, 10.5, sprintf('n=%i', n0),'fontsize', 14)
text(0.9, 10.5, sprintf('n=%i', n1),'fontsize', 14)
%% plot histogram

figure
number_of_bins = 20;
md ='lowess';

% distr = 'ev';
% h1fit=histfit(X2.sum_RelavantPercentage(X2.VanA==0), number_of_bins, distr);
% hold on
% h1fit(2).Color = [.2 .2 .9];
% h2fit=histfit(X2.sum_RelavantPercentage(X2.VanA==1), number_of_bins, distr);
binw = 0.3;
h1=histogram(X2.sum_RelavantPercentage(X2.VanA==0));
h1.BinWidth = binw;
h1.FaceColor = noVanAcolor;
hold on
h2=histogram(X2.sum_RelavantPercentage(X2.VanA==1));
h2.BinWidth = binw;
h2.FaceColor =vanAcolor;

[hi1, cx1] = hist(X2.sum_RelavantPercentage(X2.VanA==0), number_of_bins);
y1 = smooth(hi1, md);
[hi2, cx2] = hist(X2.sum_RelavantPercentage(X2.VanA==1),number_of_bins);
y2=smooth(hi2, md);
plot(cx1, y1, 'Color', [.1 .3 .95], 'LineWidth', 2)
hold on
plot(cx2, y2, 'Color',[.8 .1 .1], 'LineWidth', 2)
legend({'vanA PCR (-)'  'vanA PCR (+)' })
title('vanA can be detected using both PCR and metagenomics', 'fontsize', 16)
xlabel('relative percentage in CARD', 'fontsize', 14)
ylabel('counts of samples', 'fontsize', 14)
%%
% figure
% [N,edges] = histcounts(X2.sum_RelavantPercentage(X2.VanA==0), 'Normalization','pdf');
% edges = edges(2:end) - (edges(2)-edges(1))/2;
% plot(edges, N);




%%
% X = cardT(strcmp(cardT.sample, '1042X'), {'RelavantPercentage', 'sample'});