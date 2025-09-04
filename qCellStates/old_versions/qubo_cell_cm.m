function T = qubo_cell_cm(Xc, genelist, n, k, maxn)
    %{ 
    qubo_cell_cm computes the center of mass cost function with
    a constraint of "radius" of 10 or 10 top genes for each section.
    INPUT:
        Xc ---> scRNAseq count matrix (sce.X) (sparse) g x c
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
    
    % Sort matrix by gene expression and sequencing depth
    [X, genelist] = matrix_sort(X,genelist);

    num_genes0 = size(X,1);
    nslice = ceil(num_genes0/n);

    success = -1*ones(nslice,1);
    relative_err = 10000*ones(nslice,1);
    test_constraint = zeros(nslice,1);
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

        X = X_bak(ibeg:iend,:);
        genelist2 = genelist(ibeg:iend);
        num_genes = size(X,1);

        [CM, mass_zero ] = center_mass(X);

        if mass_zero
            fprintf("Total mass is zero, skipping %d iteration \n",ii);
            continue;
        else 
            kslice = kslice + 1;
        end

        % Centering the data and CM from other data type such as ATAC-seq 
        % may help to provide a better description
        % -----Consider centering and/or normalizing this------
        Q = CM*CM';% gene interaction by center of mass -1*Q will bring HVG

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
        M = med^2; % This seems to enforce condition

        % QUBO execution
        qprob1 = qubo(Q + M*A, M*c, M*d);
        sol1 = solve(qprob1);
        % QUBO results
        %x = sol1.BestX;
        fval(ii) = sol1.BestFunctionValue;
        genes = genelist2(sol1.BestX == 1);
        glen = length(genes);
        %fprintf("Length of jdx %d \n", glen);
        if glen > k
            genes = genes(1:k);
            kg = k;
        elseif glen < k
            kg = glen;
        else
            kg = k;
        end 
        genes_sol(ii,1:kg) = genes;
        fprintf("Gene #1 in %d: %s %s %s %s %s \n",ii,genes_sol(ii,1:5));

        % Validating constraint
        constrs = qubo(M*A, M*c, M*d); % QUBO problem just for constraints 
        test_constraint(ii) = evaluateObjective(constrs,sol1.BestX);
        sol2 = solve(constrs    );
        relative_err(ii) = abs( test_constraint(ii) - sol2.BestFunctionValue);
        relative_err(ii) = relative_err(ii)/test_constraint(ii);
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
    test_constraint = test_constraint(1:kslice);
    success = success(1:kslice);
    T = table( genes_sol, fval, test_constraint, relative_err, success );
end