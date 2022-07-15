function batch = fan_system(batch, idx_range, R, N, T, nT, detun_delta, power_delta, ...
                            sweep_series, budget, buffer_size, noise, freeze_lim, shrink)
% Autofork systems which lie on border of systems in idx_range in direction
% detun delta or detun power, where the border tightness is controlled by
% shrink (0 - 1).

    if all(idx_range == 0)
        idx_range = 1:length(batch.systems);
    end
    
    assert(xor(detun_delta, power_delta), "Cannot fan in mixed or no direction")
    
    
    points = [[batch.systems.detun_end]; ...
              [batch.systems.power_end]];
    
    
    bound = boundary(points', shrink);
    if isempty(bound)
        bound = idx_range;
    end
    
    fork_candidates = intersect(bound, idx_range);
    
    for idx = sort(fork_candidates', 'descend')
    
        [in, on] = inpolygon(points(1, idx) + detun_delta, ...
                             points(2, idx) + power_delta, ...
                             points(1, bound), ...
                             points(2, bound));
    
        if ~in && ~on
            batch = fork_system(batch, idx, R, N, T, nT, ...
                                detun_delta, power_delta, ...
                                sweep_series, budget, buffer_size, noise, freeze_lim);
        end
    
    end

end