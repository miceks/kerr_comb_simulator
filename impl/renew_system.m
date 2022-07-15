function batch = renew_system(batch, idx_range, budget)
% Increment budget of systems in idx_range

    if ~isnumeric(idx_range)
        idx_range = find(contains({batch.systems.state}, idx_range));
    end

    if all(idx_range == 0)
        idx_range = length(batch.systems);
    end

    if all(idx_range < 0)
        variances = arrayfun(@(sys) sys.variance(end)*~sys.frozen, batch.systems);
        [~, idx] = sort(variances, 'descend', 'MissingPlacement', 'last');
        idx_range = idx(1:min(abs(idx_range), length(batch.systems)));
    end

    for id = idx_range
        batch.systems(id).budget = batch.systems(id).budget + budget;
    end

end