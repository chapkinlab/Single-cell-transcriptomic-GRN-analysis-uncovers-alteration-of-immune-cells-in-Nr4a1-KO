% =========================================================================
% build_grn_compare.m
% Reproduces the scGEAToolbox "Compare Gene Network" workflow
% (callback_CompareGeneNetwork.m) with the PASTE-genes option and the
% PCR (PC Regression) method, for two batch conditions (WT vs Nr4a1 KO).
% Exports the two adjacency matrices + gene list for plotting in Python.
%
% EXACT ORDER (verified against callback_CompareGeneNetwork.m + i_transformx.m):
%   1. Normalize the FULL sce.X -- all genes, ALL cells (both batches
%      together) -- with transform option (c) = library size + log1p:
%          Xt = sc_norm(sce.X, 'type','libsize');  Xt = log1p(Xt);
%   2. Subset to the pasted panel genes (rows) AND split by batch (cols):
%          x1 = Xt(genePanel, WTcells);  x2 = Xt(genePanel, KOcells);
%   3. PC regression per batch on the subset (no subsetting afterward):
%          A = sc_grn(x, 'pcrnet');  % = net.pcrnet(x,3,false,true,false,false,gpu)
%
% Run from this folder with scGEAToolbox on the path:
%     addpath(genpath('..\scGEAToolbox'));  run('build_grn_compare.m')
% =========================================================================

%% ---- 0. INPUTS -----------------------------------------------------------
% ----- Stem cells (Fig 4B, Canonical Wnt, 25 genes): -----
% matfile  = 'stem_cells_only_exonic.mat';
% genefile = 'genes_stem_wnt_fig4b.txt';
% outtag   = 'stem';
% ----- CD8+ T cells (Fig 6B, 20 genes) -> swap to these three lines: -----
matfile  = 'cd8tcells_cells_only_exonic.mat';
genefile = 'genes_cd8t_fig6b.txt';
outtag   = 'cd8t';

% ---- folder holding the .mat, the gene .txt, and where outputs are written.
% Set to '' to use MATLAB's current folder. Uses full paths so the script
% runs no matter what the current folder is.
datadir = 'C:\Users\selim\Escritorio\claude_coworks\scgeatoolbox_stuff\grns_paper';
if isempty(datadir), datadir = pwd; end

% Fraction of strongest links kept for the layout/plot threshold, per cell
% type: stem = top 25%, cd8t = top 30% (30% keeps Nr4a1's links).
switch outtag
    case 'stem', top_frac = 0.25;
    case 'cd8t', top_frac = 0.30;
    otherwise,   top_frac = 0.25;
end
pctl = 100 * (1 - top_frac);   % percentile cutoff: 75 (stem) / 70 (cd8t)

wtLabel = "WT";
koLabel = "Nr4a1 KO";

S = load(fullfile(datadir, matfile));                     % file stores the object as `sce`
sce = S.sce;

Xall  = sce.X;                         % genes x cells (raw counts, ALL genes/cells)
gall  = string(sce.g(:));
batch = string(sce.c_batch_id(:));     % one label per cell

fprintf('Loaded %s : %d genes x %d cells\n', matfile, numel(gall), numel(batch));
u = unique(batch);
fprintf('Batch ids present: %s\n', strjoin(compose('%s (%d cells)', u, ...
        arrayfun(@(x) sum(batch==x), u)), ', '));
assert(any(batch == wtLabel), "No cells with batch id '%s'.", wtLabel);
assert(any(batch == koLabel), "No cells with batch id '%s'.", koLabel);

%% ---- 1. NORMALIZE FULL MATRIX (all genes, ALL cells together) -----------
% This is i_transformx option (c): library size + log1p, applied to sce.X.
% Library size = per-cell total over ALL genes (not just the panel).
Xt = sc_norm(Xall, 'type', 'libsize');   % (a)
Xt = log1p(Xt);                          % (b)

%% ---- 2. RESOLVE PANEL GENES, THEN SUBSET (rows) + SPLIT BATCH (cols) -----
wanted = strtrim(string(readlines(fullfile(datadir, genefile))));
wanted = wanted(strlength(wanted) > 0);
[tf, loc] = ismember(upper(wanted), upper(gall));   % as the GUI matches (upper)
if any(~tf)
    warning('%d gene(s) not found in %s and skipped: %s', ...
        sum(~tf), matfile, strjoin(wanted(~tf), ', '));
end
loc = loc(tf);
g   = gall(loc);                        % final panel symbols, in list order
fprintf('Using %d / %d panel genes.\n', numel(g), numel(wanted));

x_WT = Xt(loc, batch == wtLabel);       % pasted genes x WT cells
x_KO = Xt(loc, batch == koLabel);       % pasted genes x KO cells

%% ---- 3. PC REGRESSION PER BATCH (on the panel subset) -------------------
% Inlined so it works whether you run the whole file or section-by-section.
% Matches sc_grn(x,'pcrnet') = net.pcrnet(x,3,false,true,false,false,gpu).
% NOTE: the toolbox leaves A directed (sc_grnview2 draws a digraph); your
% plotting uses graph(), which needs symmetry, so we symmetrize A = (A+A')/2.
% Delete the two symmetrize lines to keep the directed networks.
fprintf('Building WT network...\n');
A_WT = net.pcrnet(x_WT, 3, false, true, false, false, pkg.i_usegpu(x_WT));
A_WT = (A_WT + A_WT.') / 2;
fprintf('Building KO network...\n');
A_KO = net.pcrnet(x_KO, 3, false, true, false, false, pkg.i_usegpu(x_KO));
A_KO = (A_KO + A_KO.') / 2;

%% ---- 4. MAP TO A1 / A2 ---------------------------------------------------
% Plot titles and CSV names both use A1 = KO, A2 = WT. To flip, swap RHS.
A1 = A_KO;      % 'Nr4a1 KO' -> A1_KO.csv , plotted as "{\itNr4a1} KO"
A2 = A_WT;      % 'WT'       -> A2_WT.csv , plotted as "WT"
g1 = g;         % shared panel gene list (matches rows/cols of A1, A2)

%% ---- 5. EXPORT (tagged per cell type) -----------------------------------
writematrix(A1, fullfile(datadir, sprintf('A1_KO_%s.csv',  outtag)));
writematrix(A2, fullfile(datadir, sprintf('A2_WT_%s.csv',  outtag)));
writecell(cellstr(g1), fullfile(datadir, sprintf('gene_list_%s.csv', outtag)));
save(fullfile(datadir, sprintf('grn_compare_%s.mat', outtag)), 'A1','A2','g1','A_WT','A_KO','-v7.3');

% Optional shared MATLAB "master" force layout (NetworkX computes its own).
all_weights = [abs(A1(A1 ~= 0)); abs(A2(A2 ~= 0))];
weight_threshold = prctile(all_weights, pctl);   % top_frac per cell type
A1_filtered = A1; A1_filtered(abs(A1_filtered) < weight_threshold) = 0;
A2_filtered = A2; A2_filtered(abs(A2_filtered) < weight_threshold) = 0;
A_master = max(abs(A1_filtered), abs(A2_filtered));
G_master = graph(A_master, string(g1), 'omitselfloops');
figure('Visible','off'); h_master = plot(G_master, 'Layout', 'force');
coords = [h_master.XData', h_master.YData']; close;
writematrix(coords, fullfile(datadir, sprintf('master_coords_%s.csv', outtag)));

fprintf('Exported A1_KO_%s.csv, A2_WT_%s.csv, gene_list_%s.csv, master_coords_%s.csv, grn_compare_%s.mat\n', ...
        outtag, outtag, outtag, outtag, outtag);


%% ---- 6. MATLAB SELF-PLOT (side-by-side check) ---------------------------
% Customization Settings
gene_font_size   = 10;         % font size for gene names
node_marker_size = 6;          % node size
graph_layout     = 'force';  % 'force', 'circle', 'subspace', 'layered'

% 1. Combine non-zero absolute weights & threshold (top links per top_frac)
all_weights = [abs(A1(A1 ~= 0)); abs(A2(A2 ~= 0))];
weight_threshold = prctile(all_weights, pctl);   % top_frac per cell type
A1_filtered = A1; A1_filtered(abs(A1_filtered) < weight_threshold) = 0;
A2_filtered = A2; A2_filtered(abs(A2_filtered) < weight_threshold) = 0;

% 2. Create individual graph objects using full gene lists
gene_labels = string(g1);
G1 = graph(A1_filtered, gene_labels, 'omitselfloops');
G2 = graph(A2_filtered, gene_labels, 'omitselfloops');

% 3. Create Master Graph for synchronized layout
A_master = max(abs(A1_filtered), abs(A2_filtered));
G_master = graph(A_master, gene_labels, 'omitselfloops');
figure('Visible', 'off');
h_master = plot(G_master, 'Layout', graph_layout);
master_X = h_master.XData;
master_Y = h_master.YData;
close;
% Overwrite exported coords with the SAME layout actually plotted below,
% so `python plot_grn_networkx.py --layout matlab` matches this figure.
writematrix([master_X', master_Y'], fullfile(datadir, sprintf('master_coords_%s.csv', outtag)));

% 4. Compute Edge Colors (Blue = Positive, Orange = Negative)
E1_colors = zeros(G1.numedges, 3);
for i = 1:G1.numedges
    if G1.Edges.Weight(i) > 0
        E1_colors(i, :) = [0 0.4470 0.7410];   % Blue for positive
    else
        E1_colors(i, :) = [0.8500 0.3250 0.0980]; % Orange for negative
    end
end
E2_colors = zeros(G2.numedges, 3);
for i = 1:G2.numedges
    if G2.Edges.Weight(i) > 0
        E2_colors(i, :) = [0 0.4470 0.7410];   % Blue for positive
    else
        E2_colors(i, :) = [0.8500 0.3250 0.0980]; % Orange for negative
    end
end

% 5. Render Plot with Outer Box Enclosures
fig = figure('Color', 'w', 'Position', [100 100 1100 500]);
% Left Plot: KO
ax1 = subplot(1,2,1);
h1 = plot(G1, 'XData', master_X, 'YData', master_Y, ...
    'NodeColor', [0 0.4470 0.7410], ...
    'MarkerSize', node_marker_size, ...
    'LineWidth', 1.2, ...
    'NodeLabelColor', [0 0 0], ...
    'NodeFontSize', gene_font_size);
if G1.numedges > 0
    h1.EdgeColor = E1_colors;
end
title('{\itNr4a1} KO', 'FontSize', 14, 'Color', 'k');
set(ax1, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, 'Color', 'w');
box(ax1, 'on');
% Right Plot: WT
ax2 = subplot(1,2,2);
h2 = plot(G2, 'XData', master_X, 'YData', master_Y, ...
    'NodeColor', [0 0.4470 0.7410], ...
    'MarkerSize', node_marker_size, ...
    'LineWidth', 1.2, ...
    'NodeLabelColor', [0 0 0], ...
    'NodeFontSize', gene_font_size);
if G2.numedges > 0
    h2.EdgeColor = E2_colors;
end
title('WT', 'FontSize', 14, 'Color', 'k');
set(ax2, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, 'Color', 'w');
box(ax2, 'on');

% Save the MATLAB figure too
%exportgraphics(fig, sprintf('grn_compare_matlab_%s.png', outtag), 'Resolution', 300);
exportgraphics(fig, fullfile(datadir, sprintf('grn_compare_matlab_%s.svg', outtag)), 'ContentType', 'vector');
