function  [X, genelist] = matrix_sort(X,genelist)
    %{
    matrix_sort re-sorts the matrix X according to sequencing depth along
    cell direction and gene expression along gene direction.
    INPUT: 
    X --------> Count matrix in gene x cell fashion (g,c)
    genelist -> Gene list belonging yo X count matrix
    OUTPUT:
    X --------> Count matrix (g,c) with most important elements firts
    genelist -> Gene list sorted by gene expression (largerst-smallest)
    %}
    num_cells0 = size(X,2);
    num_genes0 = size(X,1);
    
    % Computing mass and sorting rows
    MM = zeros(num_cells0,1);
    for i = 1:num_cells0
        % Mass must be the sum of all genes in a cell (cell mass)
        MM(i) = sum(X(:,i)); % Alternatively can be 1 each
    end
    [~,jdx] = sortrows(MM,'descend');
    X = X(:,jdx);

    %Sort by expressed genes in cells 
    gexpression = zeros(num_genes0,1);
    for i = 1:num_genes0
        gexpression(i) = sum(X(i,:));
    end
    [~,idx] = sortrows(gexpression,'descend');
    X = X(idx,:);
    genelist =  genelist(idx);

end