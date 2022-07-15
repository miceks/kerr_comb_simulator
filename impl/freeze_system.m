function batch = freeze_system(batch, idx_range)
% Freeze systems in idx_range

    if all(idx_range == 0)
        idx_range = length(batch.systems);
    end

    for id = idx_range
        batch.systems(id).frozen = true;
        batch.systems(id).freeze_in = 0;
    end

end