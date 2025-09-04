fname = 'pbmc3k_Bcells';
my_sce = sce;
ns = size(my_sce.X);

rate = 0.1; % can be 0.1
n = ceil(rate*ns(1));
maxn = n;
rate2 = 0.1;
k = ceil(rate2*n);
Tcm = qubo_cell_cm(my_sce.X, my_sce.g, n, k, maxn);
Tpc = qubo_cell_pearson_cor(my_sce.X, my_sce.g, n, k, maxn);

ftable_name = strcat(fname,".csv");
writetable(Tcm,ftable_name)
ftext_name = strcat(fname,'cmTranspose.txt');
writematrix(Tcm.genes_sol',ftext_name)

ftable_name = strcat(fname,".csv");
writetable(Tpc,ftable_name);
ftext_name = strcat(fname,'pcTranspose.txt');
writematrix(Tpc.genes_sol',ftext_name);
