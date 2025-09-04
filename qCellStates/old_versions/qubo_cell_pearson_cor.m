function T = qubo_cell_pearson_cor(Xc, genelist, n, k, maxn)
    %{ 
    qubo_cell_cm computes the center of mass cost function with
    a constraint of "radius" of 10 or 10 top genes for each section.
    INPUT:
        Xc ---> scRNAseq count matrix (sce.X) (sparse)
        n ----> number of genes to look per section e.g. 1000 out of 
                total genes
        k ----> number of top k genes to retrieve from cost function
    OUTPUT:
        T-----> Table containing 10 top genes of the cost function 
                per 1000 gene section in the sample. 
    NOTE: Count matrix is sorted by gene expression and cell "mass"
          computed from total number of genes in each cell.
          Xc can contain a subset of cells but must contain all the genes,
          to let decide the algorithm what is important or not.
    %}
    
    X = full(Xc);
    %X = normalize(X,"scale");
    % Initial size
    num_cells0 = size(X,2);
    num_genes0 = size(X,1);
    
    % Sorting matrix
    MM = zeros(num_cells0,1);
    for i = 1:num_cells0
        MM(i) = sum(X(:,i)); 
    end
    [MM,jdx] = sortrows(MM,'descend');
    X = X(:,jdx);

    %Sort by expressed genes in cells 
    gexpression = zeros(num_genes0,1);
    for i = 1:num_genes0
        gexpression(i) = sum(X(i,:));
    end
    [~,idx] = sortrows(gexpression,'descend');
    clear gexpression MM jdx;
    X = X(idx,:);
    genelist =  genelist(idx);


    % Compute gene-gene Pearson correlation
    %X = pearson_residuals(X);

    nslice = ceil(num_genes0/n);
    success = -1*ones(nslice,1);
    relative_err = 10000*ones(nslice,1);
    fval = zeros(nslice,1);
    X_bak = X;
    genes_sol = strings(nslice, k);
    

    kslice = 0;
    % Loop for analyzing all cells and genes provided
    for ii = 1:nslice
        fprintf("Computing %d section out of %d\n", ii,nslice);
        % retrieve the top n and look for top 10 genes with this model
        ibeg = 1 + (ii-1)*n;
        iend = n + (ii-1)*n;
        iend = min( iend, num_genes0 - n*(ii-1)); % For last chunck
        if maxn < iend
            fprintf("Max iteration reached %d in block %d \n",iend,ii);
            break;
        end
        kslice = kslice + 1;
        X = X_bak(ibeg:iend,:);
        % Compute gene-gene Pearson correlation
        X = pearson_residuals(X);

        %fprintf("ibeg %d and iend %d \n",ibeg,iend);

        genelist2 = genelist(ibeg:iend);
        num_cells = size(X,2);
        num_genes = size(X,1);

        % Centering the data X from other data type such as ATAC-seq 
        % may help to provide a better description
        % -----Consider centering and/or normalizing this------
        Q = X*X';% gene interaction 


        % Constraint definition (note: we can have more constraints)
        % M( sum(x) - max)^2 where max is the max simultaneous genes/cells
        % Scaling 1
        avg = mean(Q);
        med = median(avg);
        % Scaling 2
        %scale = mean( min(Q) + max(Q) )/2;
        max_constraint = k;
        n = num_genes;
        % (sum(x) - max)^2 constraint in A c d
        A = ones(n);
        c = -2*max_constraint*ones(n,1);
        d = max_constraint;
        %M = med^2; % This seems to enforce condition
        M = 100;
        % QUBO execution
        qprob1 = qubo(Q + M*A, M*c, M*d);
        sol1 = solve(qprob1);
        % QUBO results
        x = sol1.BestX;
        energy = [x'*Q*x, x'*M*A*x, M*c'*x, x'*M*d*x];
        fval(ii) = sol1.BestFunctionValue;
        kmin = sum(sol1.BestX==1);
        genes_sol(ii,1:kmin) = genelist2(sol1.BestX==1);
        fprintf("Gene #1 in %d: %s %s %s %s %s \n",ii,genes_sol(ii,1:5));

        % Validating constraint
        constrs = qubo(M*A, M*c, M*d); % QUBO problem just for constraints 
        test_constraint = evaluateObjective(constrs,sol1.BestX);
        sol2 = solve(constrs);
        relative_err(ii) = abs( test_constraint - sol2.BestFunctionValue);
        relative_err(ii) = relative_err(ii)/test_constraint;
        if relative_err(ii) > 10^-6
            fprintf("Constrain violated with relative error: %f \n", ...
                    relative_err(ii));
            success(ii) = false;
        else
            fprintf("Condition satisfied!\n")   
            success(ii) = true;
        end
    end 
    fval = fval(1:kslice);
    genes_sol = genes_sol(1:kslice,:);
    relative_err = relative_err(1:kslice);
    success = success(1:kslice);
    T = table( genes_sol, fval, relative_err, success );
end