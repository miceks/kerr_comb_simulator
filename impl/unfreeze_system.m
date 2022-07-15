function batch = unfreeze_system(batch, idx_range)
% Unfreeze systems in idx_range

    if all(idx_range == 0)
        idx_range = length(batch.systems);
    end

    for id = idx_range
        batch.systems(id).frozen = false;
        batch.systems(id).freeze_in = batch.systems(id).freeze_lim;
    end

end