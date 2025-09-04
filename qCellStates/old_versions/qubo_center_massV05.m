Xc = sce.X; % cell description
genelist = sce.g;
 maxn = 5000;
n = 1000;
k = 100;

    X = full(Xc);
    %X = normalize(X,"scale");
    % Initial size
    num_cells0 = size(X,2);
    num_genes0 = size(X,1);
    
    % Computing mass and sorting rows
    MM = zeros(num_cells0,1);
    for i = 1:num_cells0
        % Mass must be the sum of all genes in a cell (cell mass)
        MM(i) = sum(X(:,i)); % Alternatively can be 1 each
    end
    [MM,jdx] = sortrows(MM,'descend');
    X = X(:,jdx);

    %Sort by expressed genes in cells 
    gexpression = zeros(num_genes0,1);
    for i = 1:num_genes0
        gexpression(i) = sum(X(i,:));
    end
    [~,idx] = sortrows(gexpression,'descend');
    clear gexpression;
    X = X(idx,:);
    genelist =  genelist(idx);

    nslice = ceil(num_genes0/n);

    
    success = -1*ones(nslice,1);
    relative_err = 10000*ones(nslice,1);
    test_constraint = zeros(nslice,1);
    fval = zeros(nslice,1);
    %x    = zeros(nslice,k);
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
        num_cells = size(X,2);
        num_genes = size(X,1);

        MM = zeros(num_cells,1);
        for i = 1:num_cells
            % Mass must be the sum of all genes in each cell (cell mass)
            MM(i) = sum(X(:,i)); % Alternatively can be 1 each
        end
        % Getting center of mass for a subsample set of cells? (here are all cells)
        RI = zeros(num_genes,1);
        CM = zeros(num_genes,num_cells); 
        M = sum(MM); %  mass for gene's center 
        if M == 0
            fprintf("Total mass is zero, skipping %d iteration \n",ii);
            continue;
        end

        kslice = kslice + 1;
        for i = 1:num_genes
            % Center of mass (vector) RI = sum_j^ncells mj * x_ij) for i-th gene 
            RI(i) = dot( MM(:), X(i,:) )/M; 
            % Center of mass equation must be sum_j( RI - mj * x_ij)
            % Center of mass equation matrix form is RI_i -mj*x_ij
            % where x_ij is the i-th gene for the j-th cell
            CM(i,:) = MM(:).*( RI(i) - X(i,:)'); 
        end
        % Centering the data and CM from other data type such as ATAC-seq 
        % may help to provide a better description
        % -----Consider centering and/or normalizing this------
        Q = -1*CM*CM';% gene interaction by center of mass
        
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
        genes_sol(ii,1:k) = genelist2(sol1.BestX==1);
        fprintf("Gene #1 in %d: %s %s %s %s %s \n",ii,genes_sol(ii,1:5));

        % Validating constraint
        constrs = qubo(M*A, M*c, M*d); % QUBO problem just for constraints 
        test_constraint(ii) = evaluateObjective(constrs,sol1.BestX);
        sol2 = solve(constrs );
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