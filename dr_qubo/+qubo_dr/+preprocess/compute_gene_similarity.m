function [Xnorm2, Xnorm] = compute_gene_similarity(X)
% COMPUTE_GENE_SIMILARITY Compute cosine similarity (Gram) matrix for genes.
%
%   [Xnorm2, Xnorm] = COMPUTE_GENE_SIMILARITY(X)
%
% Computes row-normalized expression vectors (L2 norm) and their Gram matrix
% (cosine similarity between genes).
%
% INPUT:
%   X       - Count/expression matrix (genes × cells)
%
% OUTPUT:
%   Xnorm2  - Gram matrix (genes × genes), cosine similarity
%   Xnorm   - Row-normalized matrix (genes × cells), each row is unit vector
%
% AUTHOR: Selim Romero, Texas A&M University

    % Row-normalize using normr (normalize rows to unit L2 norm)
    % normr is from MATLAB's Neural Network Toolbox
    % If not available, use manual normalization
    try
        Xnorm = normr(X);
    catch
        % Manual row normalization if normr unavailable
        row_norms = sqrt(sum(X.^2, 2));
        row_norms(row_norms == 0) = 1;  % Avoid division by zero
        Xnorm = X ./ row_norms;
    end

    % Compute Gram matrix: cosine similarity between genes
    Xnorm2 = Xnorm * Xnorm';

    % Ensure diagonal is zero (no self-similarity)
    Xnorm2(1:size(Xnorm2, 1) + 1:end) = 0;

end
