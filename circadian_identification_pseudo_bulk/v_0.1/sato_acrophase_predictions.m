
path="macrophages_only_batcha.mat";
path = "fibroblasts_only_batcha_05.13.2024.mat";
data  = load(path);
sce = data.sce;
clear data;

cell_type_list = unique(sce.c_cell_type_tx);
cell_type = cell_type_list(3);
idx = find(sce.c_cell_type_tx == cell_type);

X = sce.X(:,idx);
% Consider pearson transformation
X = full(X);
g = sce.g;
cell_batch = sce.c_batch_id(idx);

X = sc_norm(X);
% This is not best for Relative error since gets 0 fit values
%X = sc_transform(X,"type","PearsonResiduals");

sce_sub = SingleCellExperiment(X,g);
sce_sub.c_batch_id = cell_batch;
sce_sub.c_cell_type_tx = sce.c_cell_type_tx(idx);
sce_sub.qcfilter;
%clear sce;

%gene_list = ["Clock","Arntl","Per1","Per2","Per3","Cry1","Cry2","Nr1d1","Nr1d2","Rora",...
%            "Rorc","Sirt1","Bhlhe40","Bhlhe41","Timeless", "Xbp1", "Atf4", "Atf6", "Hif1a"];
gene_list = sce_sub.g;
gene_list = "Clock";

m = length(gene_list);

% Batch name for each time point 
time_batch = unique(sce_sub.c_batch_id);
% Number of time points
nzts = length(time_batch);

% Analyze pseudo-bulk wise
ngenes = length(gene_list);
ncell_small = 1e50;
ncell_large = 0;
ncells_by_time_batch = size(nzts,ngenes);
for ign = 1:ngenes
    ig = find(sce_sub.g==gene_list(ign));
    for it = 1:nzts
        ic = find(sce_sub.c_batch_id == time_batch(it));
        Xg_tmp = sce_sub.X;
        Xg_tmp = full(Xg_tmp(ig,ic));
        ncells = size(Xg_tmp,2);
        ncells_by_time_batch(it,ign) = ncells;
        fprintf("Cells in %s for gene %s is %d \n",time_batch(it), gene_list(ign), ncells );
        if ncells < ncell_small
            ncell_small = ncells;
        elseif ncells < 50
            fprintf("Warning: cell partition for pseudobulk is smaller than 50 cells")
            ncell_small = ncells;
        end
        if  ncells > ncell_large
            ncell_large = ncells;
        end
    end
end

% Pseudo bulks
cell_chunk_base = ceil(ncell_small/4);
cell_chunks = ceil(ncells_by_time_batch/cell_chunk_base);

%rand_cols = randperm(size(X,2));

genename = strings(m,1);
acrov = zeros(m,1);
ci_all = zeros(m,1);
amp_all = zeros(m,1);
T_all = zeros(m,1);
mesor_all = zeros(m,1);
ZTs_all = zeros(m,nzts);
mae = zeros(m,1);
rmse = zeros(m,1);
mae_rel = zeros(m,1);
rmse_rel = zeros(m,1);

k = 0; 
for j = 1:1
    for i = 1:m
        my_gene = gene_list(i);
        %fprintf("Processing batch %s with gene %s \n", batch_custm, my_gene);
        time_cycle = 21; % Based on Sato data
        time_step = 3; % Based on Sato data
        [acro, amp, T, mesor, ci, gof, R, failed]= estimate_phase_simpler(sce_sub, my_gene, ...
                                                time_cycle, time_step, false);
        if failed
            fprintf("Gene %s failed %d \n", my_gene, failed);
            continue;
        end
        k = k + 1;

        genename(k) = my_gene;

        % Fitting function values
        acrov(k) = acro;
        amp_all(k) = amp;
        T_all(k) = T; % 24
        mesor_all(k) = mesor;

        % Evaluate cos function
        t = 0:3:21;
        fval = amp * cos( 2*pi*( t - acro)/T ) + mesor; 

        % Measure error
        lsize = length(t);
        abs_err = abs(R - fval);
        mae(k) = sum(abs_err)/lsize;
        rmse(k) = sqrt( sum(abs_err.^2)/lsize);

        % Measure relative error
        abs_err_rel = abs((R - fval)./fval);
        mae_rel(k) = sum(abs_err_rel)/lsize;
        rmse_rel(k) = sqrt( sum(abs_err_rel.^2)/lsize );

        % Bulk-like cell percentage or mean expression
        ZTs_all(k,1:nzts) = R(1:nzts);
        ci_all(k) = ci;

        % if isstruct(gof)
        %     errs(k,:) = [gof.sse gof.rmse];
        % end
    end
end

gene = genename(1:k);
% Tempo output
prior_acrophase_loc = acrov(1:k);
prior_acrophase_95_interval = ci_all(1:k);
% Circadian fit values for cosine function
acrov = acrov(1:k);
amp_all = amp_all(1:k);
T_all = T_all(1:k);
mesor_all = mesor_all(1:k);
% Predicted mean expression/cell percentage per time point
ZTs_all = ZTs_all(1:k,:);
% Errors
mae = mae(1:k);
rmse = rmse(1:k);
mae_rel = mae_rel(1:k);
rmse_rel = rmse_rel(1:k);

T3 = table(gene, amp_all, mesor_all,  acrov, T_all, mae, ...
           rmse, mae_rel, rmse_rel  );
T3.Properties.VariableNames = ["Genes","Amp","Mesor","Acrophase","Period", ...
                               "MAE","RMSE",...
                               "MAE_rel","RMSE_rel"];

%T3 = sortrows(T3,[9,7],{'ascend','ascend'});
%T3 = sortrows(T3,[9,8,7],{'ascend','ascend','ascend'});
T3 = sortrows(T3,[9,7],{'ascend','ascend'});

T = table(gene, prior_acrophase_loc, prior_acrophase_95_interval);
% fname = strcat(cell_type,"_all"); % Run this block if you want a new file
% ftable_name = strcat(fname,"_format.csv");
% writetable(T,ftable_name);


T2 = table(gene, ZTs_all(:,1), ZTs_all(:,2), ZTs_all(:,3), ZTs_all(:,4), ...
           ZTs_all(:,5), ZTs_all(:,6), ZTs_all(:,7), ZTs_all(:,8));
T2.Properties.VariableNames = ["Genes","ZT0","ZT3","ZT6","ZT9","ZT12","ZT15","ZT18","ZT21"];
% fname = strcat(cell_type,"_all"); % Run this block if you want a new file
% ftable_name = strcat(fname,"_macro_circadian.csv");
% writetable(T2,ftable_name);




