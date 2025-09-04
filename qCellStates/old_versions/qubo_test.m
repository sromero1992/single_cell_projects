
% Website tutorial
Q = [0 -1 2;...
    -1 0 4;...
    2 4 0];
c = [-5 6 -4];
d = 12;
qprob = qubo(Q,c,d);
result = solve(qprob);

sol=result.BestX;

f1 = @(x1,x2,x3) -2*x1*x2 + 4*x1*x3 + 8*x2*x3 -5*x1 +6*x2 -4*x3 +12;

x1g = 1:10;
x2g = 1:10;
x3g = 1:10;
[X1,X2,X3] = meshgrid(x1g,x2g,x3g);
F = -2*X1.*X2 + 4*X1.*X3 + 8*X2.*X3 -5*X1 +6*X2 -4*X3 +12;
gridsize = size(F);

[XX1,XX2,XX3] = ngrid(x1g,x2g,x3g);

scatter3(X1,X2,X3)
plot(F)
% solving minimum of 
% f(x) = x'Q x + c x + d  
% Q is square matrix
% c is a vector with "weights" for binary vector x
% d constant
% x is binary vector that solves the problem
% Unconstrained qubo

n = 2;
Q = -2.0*eye(n);
cv = ones(n,1);
%cv = [-2 ; -1];
d = 1;
qprob = qubo(Q,cv,d); 
result = solve(qprob);


% Constrained example
% To ensure that constraints are satisfied at a solution for a QUBO problem, 
% add a positive multiplier M times a quadratic function 
% (Const = 0)^2, e.g.    M (sum(x)-2)^2, with penalty M
% quandratic constraing comes as 
% x'Qx + c x + M( sum(x)-2)^2
% Pentalty term is defined as (sum(x))^2 - 4sum(x) +4
% A = ones(n) for x'Ax = sum(x)^2

% Original Qubo problem
Q = [0 -5 2 -6
    -5 0 -1 3
    2 -1 0 -4
    -6 3 -4 0];

qprob = qubo(Q);
result = solve(qprob)
result.BestX

% Create and solve a new problem with a linear term, a constant term, 
% and the constraint multiplier M set to 1.
A = ones(4);
c = -4*ones(4,1);
d = 4;
M = 1;
qprob2 = qubo(Q + M*A, M*c, M*d);
sol2 = solve(qprob2)
sol2.BestX

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


% Examine other returned solutions for feasibility by looking at the returned 
% quboResult.AlgorithmResult object. The AllX property has potential solutions 
% that you can examine for feasibility.
% sol2.AlgorithmResult.AllX

