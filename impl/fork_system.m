function batch = fork_system(batch, idx_range, R, N, T, nT, detun_delta, power_delta, ...
                             sweep_series, budget, buffer_size, noise, freeze_lim)
% Create new systems using fields of systems in idx_range with new final sweep
% values determined by detun_delta and power_delta.

    if all(idx_range < 0)
        idx_range = (length(batch.systems) + idx_range + 1):length(batch.systems);
    end

    states = arrayfun(@(sys) sys.state(end), batch.systems(idx_range));
    assert(~any(states == 'sweep'), 'Cannot fork system in sweep state');

    for idx = sort(idx_range, 'descend')

        sys = batch.systems(idx);
        detun_range = [0 detun_delta] + sys.detun_end;
        power_range = [0 power_delta] + sys.power_end;
        parent_field = sys.field_buffer(mod(sys.series - 1, ...
                                            sys.buffer_size) + 1, :);

        % Up- or downsample
        parent_field = interp1(linspace(-sys.T/2, sys.T/2, sys.nT), parent_field, linspace(-T/2, T/2, nT), 'linear', 0);

        
        batch = add_system(batch, parent_field, R, N, T, nT, ...
                           detun_range, power_range, ...
                           sweep_series, budget, buffer_size, noise, freeze_lim);
    end

end