% Qubo center of mass
% Fist sample a set of cells  or all cells that contain 
% some genes, e.g. 10 to 100

n = 100;
Xc = sce.X(1:n,:); % cell description
%Xg = sce.X(:,1:n); % gene description

X = Xc;
X  = full(X);

% Getting center of mass for a subsample set of cells? (here are all cells)
num_cells = size(X,2);
num_genes = size(X,1);
mi = 1;
RI = zeros(n,1);
M = 0;
for i = 1:n
    %RI(i) = sum( mi*X(:,i) ); % Cell center
    RI(i) = sum( mi*X(i,:) ); % gene center
end
% M = mi*num_cells; %  mass for gene's center
M = mi*num_genes; % mass for cell's center
RI = RI/M;

% Consider centering this
Qc = X'*X;% Cell interaction
qprob = qubo(Qc);
resultc = solve(qprob)
resultc.BestX

% Non-sense
Qg = X*X';% Gene interaction
qprob = qubo(Qg);
resultg = solve(qprob)
resultg.BestX

% Create and solve a new problem with a linear term, a constant term, 
% and the constraint multiplier M set to 1.
% M( sum(x) - max_cells)^2 

max_cells = 10;
A = ones(n);
c = -2*max_cells*ones(n,1);
d = max_cells^2;
M = 1*max_cells^3;
M2 = M*max_cells;
qprob2 = qubo(Q + M*A, M*c -M2*RI, M*d);
sol2 = solve(qprob2)
sol2.BestX

evaluateObjective(constrs,sol2.BestX)

% Check whether the solution satisfies the constraints.
constrs = qubo(M*A, M*c, M*d); % QUBO problem just for constraints 
evaluateObjective(constrs,sol2.BestX)

% is not a solution because constraints do not evaluate to zero
M = 10;
qprob2 = qubo(Q + M*A, M*c, M*d);
sol3 = solve(qprob2)

% Check the result for feasibility.
constrs = qubo(M*A, M*c, M*d); % QUBO problem just for constraints 
evaluateObjective(constrs,sol3.BestX)