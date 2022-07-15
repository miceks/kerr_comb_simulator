function batch = remove_system(batch, idx_range)
% Remove systems in idx_range

    if all(idx_range == 0)
        idx_range = length(batch.systems);
    end

    if all(idx_range < 0)
        idx_range = (length(batch.systems) + idx_range + 1):length(batch.systems);
    end

    for id = sort(idx_range, 'descend')
        batch.systems(id) = [];
    end

end

    
