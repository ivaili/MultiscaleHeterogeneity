function [out, outPC, outNorm] = plotClassifiedEdges2(mat, ids, plotFig, labels)

% Given a classification of nodes into networks (e.g. DMN, FPN, CON),
% count the number of edges within and between each network pair, and
% optionally plot the result as a matrix. Useful for summarising NBS output.
% Counts are on binary topology only.
%
% -------
% INPUTS:
% -------
% mat       - binary N*N adjacency matrix (N = number of nodes). Can be the
%             output of an NBS analysis, e.g. nbs.NBS.con_mat{1}.
%
% ids       - N*1 vector of network ids. Each network is a unique number;
%             each row is a node.
%
% plotFig   - 1 to plot `out`, 2 for `outPC`, 3 for `outNorm`, 0 otherwise.
%             Default is 0.
%
% labels    - M*1 cell of network names used for axis labels when plotting.
%
% -------
% OUTPUTS:
% -------
% out       - M*M matrix of edge counts between each network pair. Diagonal
%             entries are within-network edges.
%             NOTE: sum(sum(triu(out))) = sum(sum(adj))/2.
%
% outPC     - `out` normalised by the total number of unique edges, giving
%             the proportion of edges in each category.
%
% outNorm   - `out` normalised separately for each network pair by the
%             number of possible edges between them. The values represent
%             connection density per subgraph, which corrects for
%             differences in module size that can bias `outPC`.
%
% Alex Fornito, Monash University, Oct 2016
% Updates: Iva Ilioska

%==========================================================================
% Preliminaries
%==========================================================================

if nargin < 3
    plotFig = 0;
    labels = {};
elseif nargin < 4
    labels = {};
end

% Binarise and zero the diagonal
adj = double(logical(mat));
d = logical(eye(size(adj,1)));
adj(d) = 0;

%==========================================================================
% Compute outputs
%==========================================================================

unq = unique(ids);                  % unique network ids
M = length(unq);
out     = zeros(M);
outNorm = zeros(M);

for i = 1:M
    for j = i:M

        % Index nodes by their actual network id, not loop position,
        % so gaps in `ids` (e.g. [1 2 3 5]) are handled correctly.
        inds_i = find(ids == unq(i));
        inds_j = find(ids == unq(j));

        edgeCount = sum(sum(adj(inds_i, inds_j)));

        if i == j
            edgeCount = edgeCount / 2;              % diagonal counted twice
        end

        out(i,j) = edgeCount;
        out(j,i) = edgeCount;

        normFactor = length(inds_i) * length(inds_j);
        if i == j
            normFactor = (length(inds_i)^2 - length(inds_i)) / 2;
        end

        outNorm(i,j) = edgeCount / normFactor;
        outNorm(j,i) = edgeCount / normFactor;
    end
end

% Proportion of all unique edges falling in each block
outPC = out ./ (nnz(adj) / 2);

%==========================================================================
% Plot output
%==========================================================================

if plotFig > 0

    if plotFig == 1
        plotMat = out;
    elseif plotFig == 2
        plotMat = outPC;
    elseif plotFig == 3
        plotMat = outNorm;
    end

    figure
    imagesc(plotMat, [0, max(plotMat(:))]);
    colormap(hot);
    hold on
    set(gca,'XTick',1:1:length(labels));
    set(gca,'XTickLabel',labels,'FontSize',16,'FontWeight','Bold','XTickLabelRotation',45);
    set(gca,'YTick',1:1:length(labels));
    set(gca,'YTickLabel',labels,'FontSize',16,'FontWeight','Bold');
    colormap jet
    colorbar;
    set(gcf,'Color','w');

    % Print values in lower triangle
    for i = 1:size(plotMat,1)
        for j = i:size(plotMat,2)
            text(i, j, num2str(plotMat(i,j), 2), ...
                'HorizontalAlignment','center', ...
                'Color','w','FontSize',16,'FontWeight','Bold');
        end
    end

end
