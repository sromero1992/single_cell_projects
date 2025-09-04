path = "T-cells_Treg_k_means_120_6.20.24_ASR_ZT.mat";
data  = load(path);
sce = data.sce;
clear data;

% All cell types available
cell_type_list = unique(sce.c_cell_type_tx);
ncell_types = length(cell_type_list);

% Merge a and b replicates for this analysis
%sce.c_batch_id = strrep(sce.c_batch_id,"b","a");


% Rename batches 
batches = unique(sce.c_batch_id);
for ib = 1:length(batches)
    str_test = batches(ib);
    if strlength(str_test) < 5
        % Get batches to re-name
        idx = find(sce.c_batch_id == batches(ib));
        % Get last two characters
        str_test = erase(str_test,"ZT");
        % Renaming to same lenght, so it's naturally ordered
        batches(ib) = strcat( "ZT0", str_test);
        sce.c_batch_id(idx) = batches(ib); 
    end
end

% Reduce 16 batches to 8 batches
batches = unique(sce.c_batch_id);
for ib = 1:length(batches)
    str_test = batches(ib);
    if contains(str_test, 'b') 
        % Get batches to re-name b batch to a
        idx = find(sce.c_batch_id == batches(ib));
        % Renaming to same lenght, so it's naturally ordered
        batches(ib) = strrep(str_test,'b','a');
        sce.c_batch_id(idx) = batches(ib); 
    end
end
batches = unique(sce.c_batch_id);

% Compute R with cell percentage (true)
% Compute R with gene expression (false)
cell_pct = false;
% Predict period
predictT = false;
% N pseudo-bulk chunks
n_bulk = 500;

for icell_type = 1:ncell_types
    % Count matrix only for a cell type
    cell_type = cell_type_list(icell_type);
    idx = find(sce.c_cell_type_tx == cell_type);
    X = sce.X(:,idx);
    X = full(X);
    fprintf("Processing cell type %s \n", cell_type);
    % Normalizing count data for that cell type / (pearson residuals?)
    X = sc_norm(X);
    %X = sc_transform(X,"type","PearsonResiduals"); % This is not best for Relative error since gets 0 fit values
    
    % Gene genes and cells information
    g = sce.g;
    cell_batch = sce.c_batch_id(idx);
    
    sce_sub = SingleCellExperiment(X,g);
    sce_sub.c_batch_id = cell_batch;
    sce_sub.c_cell_type_tx = sce.c_cell_type_tx(idx);
    %sce_sub.qcfilter;
    
    % gene_list = ["Clock","Arntl","Per1","Per2","Per3","Cry1","Cry2","Nr1d1","Nr1d2","Rora",...
    %              "Rorc","Sirt1","Bhlhe40","Bhlhe41","Timeless", "Xbp1", "Atf4", "Atf6", "Hif1a"];
    %gene_list = "Sirt1";
    gene_list = sce_sub.g;            
    ngene = length(gene_list);

    % Gene names stored 
    gene_names = strings(ngene,1);
    % Batch name for each time point 
    batch_time = unique(sce_sub.c_batch_id);
    % Number of time points
    nzts = length(batch_time);

    R0 = zeros(ngene,1);
    R = zeros(n_bulk, nzts);
    acro = zeros(ngene,1);
    amp =zeros(ngene,1);
    T = zeros(ngene,1);
    mesor = zeros(ngene,1);
    failed = zeros(ngene,1);
    
    mae = zeros(ngene,1);
    rmse = zeros(ngene,1);
    mae_rel = zeros(ngene,1);
    rmse_rel = zeros(ngene,1);
    
    %textprogressbar('Processing genes: ');
    for igene = 1:ngene
        prog = floor(igene/ngene*100);
        textprogressbar(prog);
        % Gene index to work on from list
        gene_names(igene) = gene_list(igene);
        ig = find(sce_sub.g== gene_names(igene));
        %fprintf("Processing gene: %s \n", sce_sub.g(ig))
        % Set the same gene name for pseudo-bulks
        range = 1:n_bulk;
        range = range';
        % Pre-screening and computing the Gene exp for all time points
        % with all the cells
        for it = 1:nzts
            ic = find(sce_sub.c_batch_id == batch_time(it));
            % Gene expression for a time point
            Xg_tp = full(sce_sub.X(ig,ic));
            R0(igene, it) = compute_pseudoB_R( Xg_tp, 1, size(Xg_tp,2), cell_pct);
        end
        %if ~all(R0(igene,:)) % skip execution when there are at least one or more zero elements in the row.
        if ~any(R0(igene,:) > 0.1) % skip execution when all elements are below a certain threshold.
            failed(igene) = true;
            continue;
        end
        % Compute Pseudo-bulk gene expression by time point in batches of
        % 100 if possible with Bootstrap method
        for it = 1:nzts
            ic = find(sce_sub.c_batch_id == batch_time(it));
            % Gene expression for a time point
            Xg_tp = full(sce_sub.X(ig,ic));
            % Bootstrap subsampling in batches of 20 cells
            ncellb = min(20, ceil(size(Xg_tp,2)/2));
            %fprintf()
            R(range,it) = compute_pseudoB_R( Xg_tp, n_bulk, ncellb, cell_pct);
        end
        % if failed skip this gene
        time_cycle = 21; % Based on Sato data
        time_step = 3; % Based on Sato data
    
        % Sine function allocations R is nbulks * nzts
        [acro(igene), amp(igene), T(igene), mesor(igene)] = ...
                estimate_phaseR( R, time_cycle, time_step, predictT);
    
        % Allocation for errors
        mae(igene) = 10e10;
        rmse(igene) = 10e10;
        mae_rel(igene) = 10e10;
        rmse_rel(igene) = 10e10;
    
        t = 0:3:21;
        % Evaluate cos function for every time point cell percentage
        fval = zeros(1, nzts);
        fval(:) = amp(igene).*cos( 2*pi*( t - acro(igene))./T(igene) ) + mesor(igene);
    
        % Measure error
        abs_err = abs(R - fval);
        abs_err = sum(abs_err,1)/n_bulk;
        mae(igene) = sum(abs_err,2)/nzts;
        rmse(igene) = sqrt( sum(abs_err.^2, 2)/nzts);
    
        % Measure relative error
        abs_err_rel = abs_err./abs(fval);
        abs_err_rel = sum(abs_err_rel,1)/n_bulk;
        mae_rel(igene) = sum(abs_err_rel, 2)/nzts;
        rmse_rel(igene) = sqrt( sum(abs_err_rel.^2, 2)/nzts );

        % Compare the derivative direction of fval and R0 if no match, true
        diff_bool = diff_test(fval,R0(igene, :));
        if diff_bool
            failed(igene) = diff_bool;
            continue;
        end
    
    end
    %textprogressbar('done!');
    
    T1 = table(gene_names, amp, abs(amp), mesor, acro, T, mae, ...
               rmse, mae_rel, rmse_rel, failed );
    
    clear amp mesor acro T mae rmse mae_rel rmse_rel failed;
    
    T1.Properties.VariableNames = ["Genes","Amp", "Abs_Amp", "Mesor","Acrophase","Period", ...
                                   "MAE","RMSE",...
                                   "MAE_rel","RMSE_rel", "Failed"];
    
    T2 = table(gene_names, R0(:,1), R0(:,2), R0(:,3), R0(:,4), ...
               R0(:,5), R0(:,6), R0(:,7), R0(:,8));
    T2.Properties.VariableNames = ["Genes","ZT0","ZT3","ZT6","ZT9","ZT12","ZT15","ZT18","ZT21"];
    
    % Remove failed
    idx = T1.Failed == 0;
    T1 = T1(idx,:);
    T2 = T2(idx,:);
    
    % % Expression cutoff lower bound 
    idx = abs(T1.Mesor) - abs(T1.Amp) > 0;
    T1 = T1(idx,:);
    T2 = T2(idx,:);

    % % Remove Mesor + amplitude larger than 1 (100% cells) (upper bound)
    % idx = abs(T1.Mesor) + abs(T1.Amp) < 1;
    % T1 = T1(idx,:);
    % T2 = T2(idx,:);
    % 
    % % Remove Mesor + amplitude larger than 1 (100% cells) (lower bound)
    % idx = abs(T1.Mesor) - abs(T1.Amp) > 0;
    % T1 = T1(idx,:);
    % T2 = T2(idx,:);
    
    % pd_mae = fitdist(T1.MAE,'Normal'); % Lognormal looks better but don't understand it
    % pd_rmse = fitdist(T1.RMSE,'Normal');
    % pd_mae_rel = fitdist(T1.MAE_rel,'Normal');
    % pd_rmse_rel = fitdist(T1.RMSE_rel,'Normal');
    % 
    % ci95_mae = paramci(pd_mae);
    % ci95_rmse = paramci(pd_rmse);
    % ci95_mae_rel = paramci(pd_mae_rel);
    % ci95_rmse_rel = paramci(pd_rmse_rel);
    % 
    % % hist(T1.RMSE_rel,200)
    % % y = pdf(pd4,sort(T1.MAE_rel));
    % % plot(sort(T1.MAE_rel),y)
    % 
    % % Upper 95% confidence mean of pdf + 1 standard deviation 
    % mae_cutoff = ci95_mae(2,1) + ci95_mae(2,2);
    % rmse_cutoff = ci95_rmse(2,1) + ci95_rmse(2,2);
    % mae_rel_cutoff = ci95_mae_rel(2,1) + ci95_mae_rel(2,2);
    % rmse_rel_cutoff = ci95_rmse_rel(2,1) + ci95_rmse_rel(2,2);
    % 
    % idx = or( T1.MAE_rel <= mae_rel_cutoff, T1.MAE <= mae_cutoff);
    % T1 = T1(idx,:);
    % T2 = T2(idx,:);
    % 
    % idx = or( T1.RMSE_rel <= rmse_rel_cutoff, T1.RMSE <= rmse_cutoff);
    % T1 = T1(idx,:);
    % T2 = T2(idx,:);

    % % Max-min data compared to fit function filter
    % nrows = size(T1,1);
    % failed = zeros(nrows,1);
    % t0 = 0:0.1:21;
    % t = 0:3:21;
    % for i = 1:nrows
    %     fval = T1.Amp(i).*cos( 2*pi*( t0 - T1.Acrophase(i))./T1.Period(i) ) + T1.Mesor(i);
    %     Rzts = table2array( T2(i,2:end) );
    %     max_diff = abs( max(fval) - max(Rzts) );
    %     min_diff = abs( min(fval) - min(Rzts) );
    %     diff1 = abs( (T1.Abs_Amp(i) - max_diff)/T1.Abs_Amp(i) );
    %     diff2 = abs( (T1.Abs_Amp(i) - min_diff)/T1.Abs_Amp(i) );
    %     if diff1 > 1 || diff2 > 1 
    %         failed(i) = true;
    %     end
    % end
    % idx = ~failed;
    % T1 = T1(idx,:);
    % T2 = T2(idx,:);

    % Sort by acrophase and amplitude to classify better
    [T1, idx] = sortrows(T1,[5,3],{'ascend','descend'});
    T2 = T2(idx,:);
    
    fname = strcat(cell_type,"_all"); % Run this block if you want a new file
    ftable_name = strcat(fname,"_macro_circadian_analysis.csv");
    writetable(T1,ftable_name);
    
    ftable_name = strcat(fname,"_macro_circadian_ZTs.csv");
    writetable(T2,ftable_name);
end

% t0 = 0:0.1:21;
% t = 0:3:21;
% igene = find(T1.Genes=="Sirt1");
% for i = igene
%     fval = T1.Amp(i).*cos( 2*pi*( t0 - T1.Acrophase(i))./T1.Period(i) ) + T1.Mesor(i);
%     Rzts = table2array( T2(i,2:end) );
%     data_range = abs( max(Rzts) - min(Rzts) );
%     plot(t0, fval);
%     hold on;
%     plot(t,Rzts);
%     hold on;
%     % scatter(t,R(1:100,:));
%     % hold on;
% end

% inv_genes = readtable( "invasion_genes.txt", 'readvariablenames', false);
% inv_genes = string(table2array(inv_genes));


% for i = 1:length(inv_genes)
%     targets = find(T1.Genes==inv_genes(i));
%     targets
% end

