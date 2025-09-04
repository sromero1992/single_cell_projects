n = 1000;
Xc = sce.X(1:n,:); % cell description

X = Xc;
X = full(X);
%X = normalize(X,"scale");
num_cells = size(X,2);
num_genes = size(X,1);

%Sort by more genes in cells 
gnorm = zeros(num_genes,1);
for i = 1:num_genes
    gnorm(i) = norm(X(i,:)); 
end
%[gnorm,idx] = sortrows(gnorm,'descend');
%X = X(idx,:);

% Getting center of mass for a subsample set of cells? (here are all cells)
mi = 1;
RI = zeros(num_genes,1);
M = 0;
for i = 1:num_genes
    RI(i) = sum( mi*X(i,1:num_cells) ); % gene center
end
M = mi*num_cells; %  mass for gene's center
% Center of mass or centroid
RI = RI/M;

% Centering data and ATAC-seq info may help to provide a better description
% Consider centering this
Q = X*X';% gene interaction

% Create and solve a new problem with a linear term, a constant term, 
% and the constraint multiplier M set to 1.
% M( sum(x) - max)^2 where max is the max simultaneous genes/cells
avg = mean(Q);
med = median(avg);
%max_genes = med;
%max_constrain = max_genes*med;
max_constrain = 10;
A = ones(n);
c = -2*max_constrain*ones(n,1);
d = max_constrain^2;
%M = 1*max_constrain^2;
M = med^1;
% this works
qprob1 = qubo(Q + M*A, M*c, M*d);
sol1 = solve(qprob1);
sol1.BestX

genes_sol = sce.g(sol1.BestX==1);

constrs = qubo(M*A, M*c, M*d); % QUBO problem just for constraints 
test_constrain = evaluateObjective(constrs,sol1.BestX);
if abs(test_constrain) > 10^-4
    fprintf("Constrain violated with : %f", test_constrain);
end
sol_best = sol1.AlgorithmResult.AllX;
sol_best_vals = sol1.AlgorithmResult.AllFunctionValues;

qprob2 = qubo(Q + M*A, M*(c-RI), M*d);
sol2 = solve(qprob2);
sol2.BestX
constrs = qubo(M*A, M*(c-RI), M*d); % QUBO problem just for constraints 
evaluateObjective(constrs,sol2.BestX)
genes_sol2 = sce.g(sol2.BestX==1);

