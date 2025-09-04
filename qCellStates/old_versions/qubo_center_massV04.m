n = 1000;
Xc = sce.X(1:n,:); % cell description


X = Xc;
X = full(X);
%X = normalize(X,"scale");
num_cells = size(X,2);
num_genes = size(X,1);

%Sort by expressed genes in cells 
gnorm = zeros(num_genes,1);
for i = 1:num_genes
    gnorm(i) = norm(X(i,:)); 
end
[gnorm,idx] = sortrows(gnorm,'descend');
X = X(idx,:);
genelist =  sce.g(idx);

% Computing mass
MM = zeros(num_cells,1);
for i = 1:num_cells
    % Mass must be the sum of all genes in a cell (cell mass)
    %MM(i) = sum(X(:,i)); 
    MM(i) = 1;
end
[MM,jdx] = sortrows(MM,'descend');
X = X(:,jdx);
cells = sce.c_cell_cycle_tx(jdx);

%X(1:5,1:5)

% Getting center of mass for a subsample set of cells? (here are all cells)
RI = zeros(num_genes,1);
CM = zeros(num_genes,num_cells); 
cmtmp = zeros(num_genes,num_cells); 
M = sum(MM.^2); %  mass for gene's center 
for i = 1:num_genes
    % Center of mass or centroid centered by summing cells
    RI(i) = dot( MM(:),X(i,:) )/M; 
    cmtmp(i,:) = - MM(:).*X(i,:)';
    CM(i,:) = RI(i) - MM(:).*X(i,:)'; % fist matrix of center of mass
end
% Centering data and ATAC-seq info may help to provide a better description
% Consider centering this
Q = CM*CM';% gene interaction by center of mass

% Create and solve a new problem with a linear term, a constant term, 
% and the constraint multiplier M set to 1.
% M( sum(x) - max)^2 where max is the max simultaneous genes/cells
avg = mean(Q);
med = median(avg);
max_constrain = 100;
% (sum(x) - max)^2 constraint in A c d
A = ones(n);
c = -2*max_constrain*ones(n,1);
d = max_constrain;
M = med; % Very sensitive to this...
qprob1 = qubo(Q + M*A, M*c, M*d);
sol1 = solve(qprob1);
x = sol1.BestX
fval = sol1.BestFunctionValue
sum(sol1.BestX)
genes_sol = genelist(sol1.BestX==1);

%Test1
%constrs = qubo(M*A, M*(c-c2), M*d); % QUBO problem just for constraints 
constrs = qubo(M*A, M*c, M*d); % QUBO problem just for constraints 
test_constraint = evaluateObjective(constrs,sol1.BestX);
if abs(test_constraint) > 10^-4
    fprintf("Constrain violated with : %f \n", test_constraint);
end
sol_best = sol1.AlgorithmResult.AllX;
sol_best_vals = sol1.AlgorithmResult.AllFunctionValues;
