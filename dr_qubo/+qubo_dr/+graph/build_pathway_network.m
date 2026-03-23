function Xnet_target = build_pathway_network(g, genelist)
% BUILD_PATHWAY_NETWORK Create pathway network bias matrix for QUBO.
%
%   Xnet_target = BUILD_PATHWAY_NETWORK(g, genelist)
%
% Builds a bias matrix that encourages selection of pathway-focused subnetworks.
% For each gene in the target pathway, sets its row and column to -1,
% biasing the QUBO solver to include genes from the pathway.
%
% INPUT:
%   g           - Gene names of current subset (string array, length ng)
%   genelist    - Target pathway gene names (string array)
%
% OUTPUT:
%   Xnet_target - Pathway bias matrix (ng × ng)
%                 Diagonal and off-diagonals = -1 for pathway genes
%
% AUTHOR: Selim Romero, Texas A&M University

    g = string(g(:));
    genelist = string(genelist(:));

    ng = length(g);
    Xnet_target = zeros(ng, ng);

    % Find indices of genelist genes in g
    [is_in_pathway, idx_in_g] = ismember(genelist, g);
    idx_in_g = idx_in_g(is_in_pathway);

    % For each pathway gene, set its row and column to -1
    % This biases the QUBO solver to favor edges involving pathway genes
    for i = idx_in_g'
        Xnet_target(i, :) = -1;
        Xnet_target(:, i) = -1;
    end

end
