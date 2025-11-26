% canonical Wnt signaling pathway genes
genelist = ["NCOA3" "CSNK1A1" "RYK" "JUP" "WNT5A" "CAV1" "WNT5B" "AXIN2" ...
    "WNT9A" "WNT9B" "WNT16" "SFRP4" "SFRP5" "BCL9" "SFRP2" "CPE" "CD24" ...
    "LEF1" "DDX3X" "GSK3B" "USP34" "LRRK2" "LRP6" "WNT8B" "LRP5" "FZD10" ...
    "WNT8A" "CELSR2" "SFRP1" "WNT6" "FRAT1" "FRAT2" "FRZB" "WNT11" "TMEM67" ...
    "DVL1" "DVL2" "DVL3" "GPC4" "MYOC" "WNT1" "WNT2" "WNT3" "WNT4" "FZD2" ...
    "SIAH1" "FZD1" "FZD4" "TCF7L2" "WNT10B" "TCF7L1" "WNT10A" "FZD6" "WNT3A" ...
    "FZD5" "FZD8" "FZD7" "SIAH2" "WNT7A" "FZD9" "WNT7B" "CSNK1E" "CCDC88C" ...
    "FERMT2" "NR4A2" "TCF7" "PRKD2" "SNAI1" "STRN" "RAB5A" "NXN" "PKD1" ...
    "WNT2B" "ARL6" "RARG" "CTNNB1" "DIXDC1" "KLHL12" "GATA3" "PORCN" ...
    "PLCG2" "VPS35" "BCL9L" "NR4A1"];

% Load data (assuming this block runs successfully)
A = load('stem_cells_only.mat');
sce = copy(A.sce);
clear A

% Normalization should come from the whole matrix, but for now we can
% normalize for the current celltyp/portion of cells
Xsub = pkg.norm_libsize(full(sce.X), 1e4);
% For computing DR 
idg = ismember(upper(sce.g), genelist);
% Find common genes and extract the matrices
genes = sce.g(idg);
Xsub = Xsub(idg, :);

% Define Gray Color
light_gray = [0.9 0.9 0.9];

% --- 1. KO cells Heatmap ---
figure; % Use default figure background
ax_ko = gca;
ax_ko.Color = light_gray; % Set AXES background to gray
idx = contains(sce.c_batch_id, "KO");
Xko = Xsub(:,idx);
Xko_norm = normr(Xko);
Xko2norm = Xko_norm*Xko_norm'; 
imagesc(Xko2norm);
title('KO Co-expression Matrix');
colorbar;

% --- 2. WT cells Heatmap ---
figure; % Use default figure background
ax_wt = gca;
ax_wt.Color = light_gray; % Set AXES background to gray
idx = contains(sce.c_batch_id, "WT");
Xwt = Xsub(:,idx);
Xwt_norm = normr(Xwt);
Xwt2norm = Xwt_norm*Xwt_norm';
imagesc(Xwt2norm);
title('WT Co-expression Matrix');
colorbar;

% --- 3. Differential Co-expression Matrix (Xdiff) Heatmap ---
Xdiff = Xko2norm - Xwt2norm;
figure; % Use default figure background
ax_diff = gca;
ax_diff.Color = light_gray; % Set AXES background to gray
imagesc(Xdiff); % Plot the filtered Xdiff matrix
colorbar; % Add a color bar

% Use a standard divergent colormap and center the color axis
colormap(ax_diff, 'jet'); % Using 'jet' as a standard replacement for 'coolwarm'
max_abs_val = max(abs(Xdiff(:)));
caxis([-max_abs_val, max_abs_val]); % Center colormap around zero

title('Filtered Xdiff Heatmap (KO vs WT Co-expression)');
xlabel('Genes'); 
ylabel('Genes'); 

% Set both X and Y axis tick labels using 'genes'
num_genes = size(Xdiff, 1);
set(gca, 'XTick', 1:num_genes);        
set(gca, 'XTickLabel', genes);            
set(gca, 'XTickLabelRotation', 90);   
set(gca, 'YTick', 1:num_genes);        
set(gca, 'YTickLabel', genes);            
set(gca, 'YTickLabelRotation', 0);    

% Adjust font properties
ax_diff.FontSize = 8;
ax_diff.FontWeight = 'bold';
ax_diff.XLabel.FontSize = 10;
ax_diff.YLabel.FontSize = 10;
ax_diff.Title.FontSize = 12;

% --- 4. Network Plot ---
% Filter for edges (This part remains unchanged and correct)
idx = ~isnan(Xdiff) & Xdiff ~= 0;
non_nan_abs_diffs = abs(Xdiff(idx));

if isempty(non_nan_abs_diffs)
    warning('No significant differences to plot in the network after filtering.');
    disp('Consider checking your Xdiff matrix or adjusting the initial filtering.');
    return;
end

% Use the 95th percentile of absolute differences as the threshold
edge_display_threshold = prctile(non_nan_abs_diffs, 95);
disp(['Edge display threshold (95th percentile): ', num2str(edge_display_threshold)]);

[rows, cols, values] = find(triu(Xdiff, 1));
valid_indices = abs(values) >= edge_display_threshold;
source_nodes = rows(valid_indices);
target_nodes = cols(valid_indices);
edge_weights = values(valid_indices);

if isempty(source_nodes)
    warning('No edges meet the specified display threshold for network plotting.');
    return;
end

% 3. Create a graph object
G = graph(source_nodes, target_nodes, edge_weights, genes);

% 4. Plot the network and capture the handle 'h'
figure; % Create new figure
ax_graph = gca;
ax_graph.Box = 'off'; % Remove the box/border

h = plot(G, ...
         'Layout', 'force', ...          
         'NodeLabel', G.Nodes.Name, ...  
         'NodeLabelColor', [0 0 0], ...  
         'NodeColor', [0.7 0.7 0.7], ... 
         'MarkerSize', 8);      

ax_graph.Color = 'w'; % Set AXES background to WHITE

title('Gene Co-expression Difference Network (Top 5% | KO vs. WT)');

% --- Apply general properties to the VALID object 'h' immediately ---
h.NodeFontSize = 10; 
h.NodeFontWeight = 'bold';
h.EdgeAlpha = 1.0; % Set alpha to 1.0 
% NO general h.EdgeColor setting here, to preserve highlight colors

% --- Customize edge appearance using highlight() ---
pos_edges_idx = find(edge_weights > 0); % KO gain (increased co-expression)
neg_edges_idx = find(edge_weights < 0); % KO loss (decreased co-expression)
fixed_linewidth = 2; 

if ~isempty(pos_edges_idx)
    highlight(h, source_nodes(pos_edges_idx), target_nodes(pos_edges_idx), ...
              'EdgeColor', 'r', ... % Red for increased co-expression in KO
              'LineWidth', fixed_linewidth);
end

if ~isempty(neg_edges_idx)
    highlight(h, source_nodes(neg_edges_idx), target_nodes(neg_edges_idx), ...
              'EdgeColor', 'b', ... % Blue for decreased co-expression in KO
              'LineWidth', fixed_linewidth);
end
