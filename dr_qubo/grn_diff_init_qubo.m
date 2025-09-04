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

%genelist = ["LGR5" "WNT2B"  "CTNNB1"];
% For plotting difference scale matrix and extract genes
X = pkg.norm_libsize(full(sce.X), 1e4);
g = sce.g;
ng = length(g);
K = 200;
K0 = K;
cell_state = sce.list_cell_attributes{10};
cust_idx = ismember(upper(g), genelist);

% Extract sample/condition
idx = sce.c_batch_id == "Nr4a1 KO";
Xko = X(:,idx);
% Row normalize gene vector ensure unitary vector this way
Xko_norm = normr(Xko);
% Co-variance matrix g x g
Xko2norm = Xko_norm*Xko_norm'; % Similarity/adjacency matrix
% Extract desired cell state to include
cs_ko = cell_state(idx);
% Normalize cell state and obtain gene vector
cs_ko = cs_ko / norm(cs_ko);
Xko_cs = Xko_norm*cs_ko; 
% Compute the nearest neighboars of target genes only
Xko_mnn = adjX_mat_construct_sparse_cust_idx(Xko, 'mnn', K0, cust_idx); % large K gets more mnn


% Extract sample/condition
idx = sce.c_batch_id == "WT";
Xwt = X(:,idx);
% Row normalize gene vector ensure unitary vector this way
Xwt_norm = normr(Xwt);
% Co-variance matrix g x g
Xwt2norm = Xwt_norm*Xwt_norm';% Similarity/adjacency matrix
% Extract desired cell state to include
cs_wt = cell_state(idx);
% Normalize cell state and obtain gene vector
cs_wt = cs_wt/norm(cs_wt);
Xwt_cs = Xwt_norm*cs_wt; 
% Compute the nearest neighboars of target genes only
Xwt_mnn = adjX_mat_construct_sparse_cust_idx(Xwt, 'mnn', K0, cust_idx); % large K gets more mnn

% We want negative values being optimized, so Xko as negative
Xdiff = Xwt2norm - Xko2norm; % This gets rid of the 1 diagonal by definition

% Cell state information in genes difference
Vdiff = Xwt_cs - Xko_cs;

% Create interaction network from desired selected genes or pathway genes
idx = ismember(upper(g), genelist);
Xnet_target = zeros(ng, ng);
for ig = 1:ng
    if idx(ig) == true
        Xnet_target(ig, idx) = -1;
        Xnet_target(idx, ig) = -1;
    end
end

Q0 = Xdiff;
% This is first level change
Q = Xnet_target + Xdiff + diag(Vdiff);
% Favor any mutual nearest neigbors
Q = Q - full(Xwt_mnn + Xko_mnn); 
% % Favor any mutual nearest neigbors overlapped
% sym_mnn = full(Xwt_mnn*Xko_mnn);
% sym_mnn = normalize(sym_mnn, 'range', [0 1]);
% sym_mnn = max(sym_mnn, sym_mnn');
% Q = Q - m*sym_mnn;
Q1 = Q;

% Diagonals become zero, that's a good sign
% figure;
% imagesc(Q); % Plot the filtered Xdiff matrix
% colorbar; % Add a color bar
% colormap("parula"); % Consider using a divergent colormap like 'coolwarm' for differences

P = 10*max(abs(Q1(:))); % A sufficiently large penalty. A common heuristic is a multiple of the max
                         % absolute value in your original Q. Adjust as needed.
N = size(Q1, 1); % Number of genes (variables)

% Add penalty for off-diagonal terms
Q_penalty = P * ones(N, N); % Matrix of all P's
Q_penalty = Q_penalty - diag(diag(Q_penalty)); % Set diagonal to zero for off-diagonal penalty
Q = Q1 + Q_penalty;


% Add penalty for diagonal terms
for i = 1:N
    Q(i, i) = Q(i, i) + P * (1 - 2*K);
end

qprob  = qubo(Q);
result = solve(qprob);
idx = result.BestX == 1;
nsols = sum(idx)
sub_g_net = g(idx);
%sub_Q_net = -Q0(idx, idx);
sub_Q_net = -Q1(idx, idx);

figure;
imagesc(sub_Q_net); % Plot the filtered Xdiff matrix
colorbar; % Add a color bar
colormap("parula"); % Consider using a divergent colormap like 'coolwarm' for differences

title('Filtered Xdiff Heatmap (KO vs WT Co-expression)');
xlabel('Genes'); % X-axis represents genes (columns)
ylabel('Genes'); % Y-axis represents genes (rows)

% Set both X and Y axis tick labels using 'g' (your gene names)
num_genes = size(sub_Q_net, 1); % Get the number of genes

set(gca, 'XTick', 1:num_genes);        % Set X-ticks for each column
set(gca, 'XTickLabel', sub_g_net);            % Apply gene names to X-axis
set(gca, 'XTickLabelRotation', 90);   % Rotate X-labels for better readability

set(gca, 'YTick', 1:num_genes);        % Set Y-ticks for each row
set(gca, 'YTickLabel', sub_g_net);            % Apply gene names to Y-axis
set(gca, 'YTickLabelRotation', 0);    % Keep Y-labels horizontal

% --- Set the font weight of the tick labels to 'bold' ---
ax = gca; % Get the current axes handle
ax.FontSize = 8; % Adjust font size as needed (can be larger for bold)
ax.FontWeight = 'bold'; % Make all axis labels (including ticks) bold

ax.XLabel.FontSize = 10;
ax.YLabel.FontSize = 10;
ax.Title.FontSize = 12;


writetable(table(sub_g_net),'qubo_genes_sol_dr.txt');


sub_Q_net = -Q0(idx, idx);
% NEtwork plot
idx = ~isnan(sub_Q_net) & sub_Q_net ~= 0;
non_nan_abs_diffs = abs(sub_Q_net(idx));

if isempty(non_nan_abs_diffs)
    warning('No significant differences to plot in the network after filtering.');
    disp('Consider checking your sub_Q_net matrix or adjusting the initial filtering.');
    return; % Exit the script if no valid edges exist
end

% Example: Use the 80th percentile of absolute differences as the threshold
edge_display_threshold = prctile(non_nan_abs_diffs, 95); % Only show top 20% of abs diffs
disp(edge_display_threshold)


% 2. Extract source, target nodes, and edge weights from Xdiff_filtered
% and self-loops (A-A connections).
[rows, cols, values] = find(triu(sub_Q_net, 1));

% Filter for edges that meet the display threshold
valid_indices = abs(values) >= edge_display_threshold;
source_nodes = rows(valid_indices);
target_nodes = cols(valid_indices);
edge_weights = values(valid_indices);

if isempty(source_nodes)
    warning('No edges meet the specified display threshold for network plotting.');
    return;
end

% 3. Create a graph object
% 'g' (your gene names) will be used as the Node Names.
G = graph(source_nodes, target_nodes, edge_weights, sub_g_net);

% 4. Plot the network
figure;
h = plot(G, ...
         'Layout', 'force', ...          % 'force' layout often provides good spatial separation
         'NodeLabel', G.Nodes.Name, ...  % Use gene names as node labels
         'NodeLabelColor', [0 0 0], ...  % Black labels
         'NodeColor', [0.5 0.5 0.5], ... % Light grey nodes
         'MarkerSize', 5);               % Adjust node marker size

title('Gene Co-expression Difference Network (KO vs. WT)');


% --- Customize edge appearance based on positive/negative differences ---
% Positive edges (increased co-expression in KO vs WT)
pos_edges_idx = find(edge_weights > 0);
if ~isempty(pos_edges_idx)
    % Use a fixed, larger LineWidth for all positive edges
    fixed_positive_linewidth = 2; % Adjust this value as desired
    highlight(h, source_nodes(pos_edges_idx), target_nodes(pos_edges_idx), ...
              'EdgeColor', 'b', ... % Red for positive differences
              'LineWidth', fixed_positive_linewidth);
end

% Negative edges (decreased co-expression in KO vs WT)
neg_edges_idx = find(edge_weights < 0);
if ~isempty(neg_edges_idx)
    % Use a fixed, larger LineWidth for all negative edges
    fixed_negative_linewidth = 2; % Adjust this value as desired
    highlight(h, source_nodes(neg_edges_idx), target_nodes(neg_edges_idx), ...
              'EdgeColor', 'r', ... % Blue for negative differences
              'LineWidth', fixed_negative_linewidth);
end

% Optional: Adjust node label properties for better readability
h.NodeFontSize = 15; % Or any other desired size
h.NodeFontWeight = 'bold'; % If you still want them bold
