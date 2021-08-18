close all; clear;

vre_tr= phytreeread('1044_VRE_tree_1000gene_treeWithGenomeIds.nwk');

%%
plot(vre_tr)
%%%% reroot using the Enterococcus faecalis 62 genome
ind = getbyname(vre_tr,'936153.3'); 
[sel,sel_leaves, sel_branches] = select(vre_tr,'criteria','distance',...
                          'threshold',0.00001,'reference',ind);
                                           
vre_tr2 = reroot(vre_tr, sel);
plot(vre_tr2)

%%
%%%% remove Enterococcus faecalis 62 because it compresses the other
%%%% isolates too much for clear visulization

ind = getbyname(vre_tr2,'936153.3');
vre_tr3 = prune(vre_tr2,ind);
view(vre_tr3)

%% replace the patric ref ID with genome name
genomeTbl = readtable('VRE1044_genomeTbl.xlsx');
genomeTbl = genomeTbl(:, 1:2);
% get the leave names of vre_tr3
leaveNames = get(vre_tr3,'LeafNames');
tree_gNames = repmat({''}, length(leaveNames),1);
for i=1:length(leaveNames)
    ind=find(ismember(genomeTbl.GenomeID, leaveNames{i}));
    tree_gNames{i} = genomeTbl.GenomeName{ind};    
end

%%
tree_gNames=cellfun(@(X) strrep(X, 'Enterococcus faecium ', ''), tree_gNames, 'UniformOutput', false);
ind = find(~contains(tree_gNames, 'VRE'));
% tree_gNames(ind) = strcat('MAG_', tree_gNames(ind));
%%
fill_cmap = [ 255     153    255  
    178    102    255 
    51    153    255 
    0    204    102 
    230    140    0        
    230    90    110
    255 215 0] / 255;
outline_cmap = [.7 .8 .8;
                .1 .1 .1 ];

types={'VRE', 'MAG'};
outline_color = zeros(length(leaveNames), 3);
for i=1:length(types)
    ind = find(contains(tree_gNames, types{i})); 
    outline_color(ind, :) =  repmat(outline_cmap(i, :),length(ind), 1);
end

samples = cellfun(@(X) strrep(X, 'Enterococcus faecium ', ''), genomeTbl.GenomeName, 'UniformOutput', false);
samples = samples(contains(samples, 'VRE'));
samples = cellfun(@(X) X(1:5), samples, 'UniformOutput', false);
samples = unique(samples);
fill_color= zeros(length(leaveNames), 3);
samples_for_each_leaf = repmat({''}, length(tree_gNames),1);
for i=1:length(samples)
    ind = find(contains(tree_gNames, samples{i})); 
    fill_color(ind, :) =  repmat(fill_cmap(i, :),length(ind), 1);
    samples_for_each_leaf(ind)=repmat(samples(i),length(ind), 1);
end
%
[l, ia, ~]=unique(samples_for_each_leaf);
figure
for i=1:length(ia)
    plot(1, 0.6*i, 'o', 'MarkerFaceColor', fill_color(ia(i),:), ...
        'MarkerEdgeColor', fill_color(ia(i),:), 'MarkerSize', 15)
    text(1.1, 0.6*i, samples_for_each_leaf{ia(i)}, 'fontsize', 14)
    hold on
end
plot(1, 0.6*(i+1), 'o', 'MarkerFaceColor', [1 1 1], ...
        'MarkerEdgeColor', outline_cmap(1,:), 'MarkerSize', 15, 'linewidth', 2)
text(1.1, 0.6*(i+1), 'MAG', 'fontsize', 14)
hold on
plot(1, 0.6*(i+2), 'o', 'MarkerFaceColor', [1 1 1], ...
        'MarkerEdgeColor', outline_cmap(2,:), 'MarkerSize', 15, 'linewidth', 2)
text(1.1, 0.6*(i+2), 'isolates', 'fontsize', 14)


hold off
set(gca, 'Ylim', [-1 length(ia)+1])
title('legend', 'fontsize', 12)
box off
%
h = plot(vre_tr3, 'TerminalLabels', false, 'Type','square');
h.LeafDots.Color = [1 1 1];
h.BranchDots.Color = [.2 .2 .2];
h.BranchDots.MarkerSize = 1;
h.LeafDots.MarkerEdgeColor = [1 1 1];
hold on
for i=1:length(tree_gNames)
	plot(h.LeafDots.XData(i),h.LeafDots.YData(i),'o', ...
        'MarkerSize', 9, ...
        'MarkerFaceColor', fill_color(i,:),...
        'MarkerEdgeColor', outline_color(i,:), ...
        'linewidth', 2)
end
hold on
text(h.LeafDots.XData+.00005,h.LeafDots.YData, tree_gNames, 'fontSize', 12, 'Interpreter', 'none')
set(gca, 'xlim', [-0.0003 .0013])

set(gca,'YTick',[]);
set(gcf,'renderer','painters')
