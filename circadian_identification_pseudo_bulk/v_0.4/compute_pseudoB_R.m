function  R = compute_pseudoB_R( Xg , n_pseudo_bulk, cell_max, cell_pct)

    % Expression for n-pseudo bulk and nzts time points

    % Shuffle cells for this time point 
    ncells = size(Xg,2);
    R = zeros(n_pseudo_bulk, 1);
    if cell_pct
        for kb = 1:n_pseudo_bulk
            ic_new =  randperm(ncells);
            Xg0 = Xg(ic_new(1:cell_max));
            R(kb) = sum( Xg0 > 0)./length(Xg0);
        end
    else
        for kb = 1:n_pseudo_bulk
            ic_new =  randperm(ncells);
            Xg0 = Xg(ic_new(1:cell_max));
            R(kb) = mean(Xg0);
        end
    end
end