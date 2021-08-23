%% 
close all; clear;
addpath('../utils');
X=load('../savedMat/Tcorr_shannon');
Tcorr_shannon = X.Tcorr_shannon;

Y=load('../savedMat/corrT');
corrT = Y.corrT;
clear X Y

Taxa = {'Phylum' 'Class' 'Order' 'Family' 'Genus'};

%% scatter plot of read counts vs correlation
figure()
cmap = [.9 .6 .2; ...
        .2 .6 .8; ...
        .8, .2 .2;...
        .8 .8 .2];
dotS = 15;
for i=1:length(Taxa)
    subplot(5, 1, i)
    scatter(corrT{corrT.([Taxa{i} '_Pvalue'])<=0.05,[Taxa{i} '_Cor']}, ...
            log10(corrT.readcount(corrT.([Taxa{i} '_Pvalue'])<=0.05)), ...
            dotS, 'MarkerFaceColor', cmap(1, :),'MarkerEdgeColor', cmap(1,:)*0.6)
    n1=length(corrT{corrT.([Taxa{i} '_Pvalue'])<=0.05,[Taxa{i} '_Cor']});
    hold on
    scatter(corrT{corrT.([Taxa{i} '_Pvalue'])>0.05,[Taxa{i} '_Cor']}, ...
            log10(corrT.readcount(corrT.([Taxa{i} '_Pvalue'])>0.05)), ...
            dotS, 'MarkerFaceColor', cmap(2, :), 'MarkerEdgeColor', cmap(2, :)*0.6)
    n2=length(corrT{corrT.([Taxa{i} '_Pvalue'])>0.05,[Taxa{i} '_Cor']});
    legend({sprintf('significant (n=%i)', n1), ...
                sprintf('nonsignificant (n=%i)', n2)}, ...
                'Location','eastOutside')
    title(Taxa{i}, 'FontSize', 14)    
    ylabel('log10(readscount)')
    if i==length(Taxa)
        xlabel('Correlation', 'FontSize', 14)
    end
    set(gca, 'xlim', [-0.2 1])

end


%% violin plot of correlation in each stool consistency
cmap = [.1 .3 .8;
        .6 .7 .1;
        .6 0.1 .4 ];

alphaValue = 0.5;
binwidth = 0.1;
stoolC=cellstr(Tcorr_shannon.Consistency);
stoolG = {'formed' 'semi-formed' 'liquid'};
diverseOrder ={'low diversity', 'middle diversity', 'high diversity'};
dotSz = 50;
fig=figure
for i=1:length(Taxa)

    A = Tcorr_shannon{:, Taxa{i}};
    A_cutoff = prctile(A, [33 67]);
    Diversity =repmat({''}, height(Tcorr_shannon),1);

    Diversity(A <= A_cutoff(1)) = {'low diversity'};
    Diversity(A > A_cutoff(1) & A <= A_cutoff(2)) = {'middle diversity'};
    Diversity(A > A_cutoff(2)) = {'high diversity'};
    
    for j=1:3
        subplot(length(Taxa), 3, (i-1)*3+j)
        idx = ismember(Tcorr_shannon.Consistency, stoolG{j});

        vs = violinplot(Tcorr_shannon{idx , [Taxa{i} '_Cor']}, Diversity(idx), ...
                    'Bandwidth', 0.1, ...
                    'ViolinAlpha', 1, ...
                    'ShowNotches', false, ...
                    'ShowMean', false, ...
                    'GroupOrder', diverseOrder);
        d =Diversity(idx);

        vs(1).MedianPlot.MarkerFaceColor = [0 0 0];
        vs(2).MedianPlot.MarkerFaceColor = [0 0 0];
        vs(3).MedianPlot.MarkerFaceColor = [0 0 0];
        vs(1).MeanPlot.LineWidth=3;
        vs(2).MeanPlot.LineWidth=3;
        vs(3).MeanPlot.LineWidth=3;
        % vs(1).WhiskerPlot.MarkerFaceColor = [ 0 0 0];
        % vs(1).WhiskerPlot.MarkerFaceColor = [ 0 0 0];
        vs(1).ViolinColor = cmap(1,:);
        vs(2).ViolinColor = cmap(2,:);
        vs(3).ViolinColor = cmap(3,:);
        vs(1).WhiskerPlot.Color = [.2 .2 .2];
        vs(1).ScatterPlot.SizeData=dotSz;
        vs(2).ScatterPlot.SizeData=dotSz;
        vs(3).ScatterPlot.SizeData=dotSz;
        set(gca, 'ylim', [-0.25 1.08], 'ytick', 0:0.5:1, ...
                'xticklabel',  {sprintf('n=%i', length(find(ismember(d, 'low diversity')))), ...
                                sprintf('n=%i', length(find(ismember(d, 'middle diversity')))), ...
                                sprintf('n=%i', length(find(ismember(d, 'high diversity'))))}, ...
                 'xlim', [0.5 3.5]) 
        xlabel(stoolG{j}, 'fontsize', 12)
        if j==1
            ylabel('Correlation', 'fontsize', 12)
        else
            set(gca, 'yticklabel', [])
            box off
        end
        if i==1 & j==1
            legend([vs(1).ScatterPlot, vs(2).ScatterPlot, vs(3).ScatterPlot], ...
                    {'low diversity' ,...
                    'middle diversity', ...
                    'high diversity'}, 'fontsize', 10)
        end
                                                                  
    end
end
set(fig,'renderer','painters');
