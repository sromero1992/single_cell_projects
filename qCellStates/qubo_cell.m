function T = qubo_cell(Xc, genelist, n, k, maxn, method, scale)
    %{ 
    qubo_cell_cm computes the center of mass cost function with
    a constraint of "radius" of 10 or 10 top genes for each section.
    INPUT:
        Xc -----> scRNAseq count matrix (sce.X) (sparse) g x c
        n ------> Number of genes to look per section e.g. 1000 out of 
                  total genes
        k ------> Number of top k genes to retrieve from cost function
        method -> Method to retrive matrix:
                    pc = pearson correlation 
                    cm = center of mass
                    jac = jaccard indexing 
        nhvg --> If present work with n Highly variable genes
    OUTPUT:
        T-----> Table containing 10 top genes of the cost function 
                per 1000 gene section in the sample. 
    NOTE: Count matrix is sorted by gene expression and cell "mass"
          computed from total number of genes in each cell.
          Xc can contain a subset of cells but must contain all the genes,
          to let decide the algorithm what is important or not.
    %}
     
    % n arg in check here for nhvg
    fprintf("Qubo cells mode: %s \n",method)
    X = full(Xc);
    
    % Sort matrix by gene expression and sequencing depth
    %[X, genelist] =matrix_sort(X, genelist);

    num_genes0 = size(X,1);
    nslice = ceil(num_genes0/n);

    % Initializations for writting files
    fval = zeros(nslice,1);
    success = -1*ones(nslice,1);
    genes_sol = strings(nslice, k);
    rel_pct_err = 10000*ones(nslice,1);
    test_constraint = zeros(nslice,1);
    constraint_val = zeros(nslice,1);

    X_bak = X;

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
               
        % Stop if matrix is practically zero
        zero_stop = sum(sum(X)) == 0;
        if zero_stop
            fprintf("Total matrix is zero, skipping %d iteration \n",ii);
            continue;
        else 
            kslice = kslice + 1;
        end

        offset = 10;
        % Generate matrix X for QUBO method
        switch method
            case 'cm' 
                [X, zero_stop ] = center_mass(X);
                Q = X*X';
            case 'pc'
                % gene interaction by center of mass -1*Q will bring HVG 
                % depending on the matrix input
                X = pearson_residuals(X);
                % If HVG matrix is input, this gets genes belonging to
                % different cell type
                Q  = X*X';
                % If HVG matrix is input, this gets genes belonging to
                % cell type
                %Q  = -X*X';
            case 'jac'
                X = jaccard_all_pairs(X, false); 
                Q = X;
            case 'pcnet'
                [X] = log( sc_norm(X) + 1); 
                [Q] = sc_pcnet(X); 
                fprintf("PCNET matrix norm: %f \n", norm(Q,2))
                offset = 1;
            otherwise
                X = pearson_residuals(X);
                Q  = X*X';
        end

        if scale
            Q = minmax_scaling(Q);
            Q = Q;
        end
        fprintf("Q matrix norm: %f \n", norm(Q,2))

        % Constraint definition (note: we can have more constraints)
        % M( sum(x) - max)^2 where max is the max simultaneous genes/cells
        % Scaling 1
        avg = mean(Q);
        med = abs(median(avg));
        % Scaling 2
        %scale = mean( min(Q) + max(Q) )/2;
        max_constraint = k;
        n = num_genes;
        % (sum(x) - max)^2 constraint in A c d
        A = ones(n);
        c = -2*max_constraint*ones(n,1);
        d = max_constraint;
        if ~scale 
            M = ( med + offset)^2; % This enforce condition
        else
            M = med + offset;
        end
        fprintf("Annealing constant M = %f \n",M);

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
        % Evaluate sol1 X values
        test_constraint(ii) = evaluateObjective(constrs,sol1.BestX);
        sol2 = solve(constrs);
        % Get constraint solution
        constraint_val(ii) = sol2.BestFunctionValue;
        rel_pct_err(ii) = test_constraint(ii) - constraint_val(ii);
        rel_pct_err(ii) = abs(rel_pct_err(ii)/constraint_val(ii))*100;
        if rel_pct_err(ii) > 5.0
            fprintf("Constrain violated with relative percentage error: %f \n", ...
                    rel_pct_err(ii));
            success(ii) = false;
        else
            fprintf("Condition satisfied!\n")   
            success(ii) = true;
        end
    end 
    fval = fval(1:kslice);
    genes_sol = genes_sol(1:kslice,:);
    rel_pct_err = rel_pct_err(1:kslice);
    test_constraint = test_constraint(1:kslice);
    constraint_val = constraint_val(1:kslice);
    success = success(1:kslice);
    T = table( genes_sol, fval, test_constraint, constraint_val, rel_pct_err, success );
end