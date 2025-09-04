fname = 'pbmc3k_Bcells';
load(fname)
my_sce = sce;

%ns = size(my_sce.X);
%rate = 0.1; % can be 0.1
%n = ceil(rate*ns(1));
%maxn = n;
%rate2 = 0.1;
%k = ceil(rate2*n);

[~, X, g] = sc_hvg(my_sce.X, my_sce.g);

% Max number of genes to look at
n = 5000;
% maxn for each slice, if requiered
maxn = 5000;
% k features to extract
k = 100;

X = X(1:n,:);
g = g(1:n);

scale = true;
T = qubo_cell(X, g, n, k, maxn,'pcnet', scale);

ftable_name = strcat(fname,".csv");
writetable(T,ftable_name);
ftext_name = strcat(fname,'Transpose.txt');
writematrix(T.genes_sol',ftext_name);

