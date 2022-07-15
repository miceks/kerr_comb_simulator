function batch = add_system(batch, field, R, N, T, nT, ...
                            detun_range, power_range, ...
                            sweep_series, budget, buffer_size, noise, freeze_lim)
% Create and add a new system to batch. 

    sys.id = batch.id_cnt;
    batch.id_cnt = batch.id_cnt + 1;
    
    sys.series = 1;
    sys.budget = budget;
    sys.frozen = false;
    sys.alpha  = batch.alpha;
    sys.s      = batch.s;

    sys.R  = R;     % Roundtrips per series
    sys.N  = N;     % Steps per roundtrip
    sys.T  = T;     % Time window size
    sys.nT = nT;    % Time window points

    
    % Initialize sweep parameters
    sys.sweep_series = sweep_series;

    sys.detun_start = detun_range(1);
    sys.detun_end = detun_range(end);
    sys.detun = detun_range(1);

    sys.power_start = power_range(1);
    sys.power_end = power_range(end);
    sys.power = power_range(1);


    % Initialize circular field buffer
    sys.buffer_size = buffer_size;
    sys.field_buffer = NaN(buffer_size, nT);
    sys.field_buffer(1, :) = field;


    % Generate CW reference for sweep values
    cw_upper = NaN(1, sweep_series);
    cw_lower = NaN(1, sweep_series);

    for n = 1:sweep_series

        detun = interp1([1 sys.sweep_series], [detun_range(1) detun_range(end)], n);
        power = interp1([1 sys.sweep_series], [power_range(1) power_range(end)], n);

        switch batch.type
            case 'Ikeda'
                [upper, ~, lower] = cw_Ikeda(batch.alpha, power, ...
                                             batch.cw_newton_tol, ...
                                             batch.cw_curve_tol);
            case 'LLE'
                [upper, ~, lower] = cw_LLE(batch.alpha, power, ...
                                           batch.cw_power_levels);
        end
              
        cw_upper(n) = interp1(upper(1, :), upper(2, :), max(detun, -pi));
        cw_lower(n) = interp1(lower(1, :), lower(2, :), min(detun, pi));
    end

    sys.cw_upper_power = @(N) interp1(1:sweep_series, cw_upper, ...
                                      N, 'linear', cw_upper(end));
    sys.cw_lower_power = @(N) interp1(1:sweep_series, cw_lower, ...-
                                      N, 'linear', cw_lower(end));


    %Initialize statistics
    sys.freeze_lim = freeze_lim;
    sys.freeze_in  = freeze_lim;

    sys.variance      = NaN;
    sys.cw_lower_diff = NaN;
    sys.cw_upper_diff = NaN;

    states = {'sweep', 'fluctuating',  'cw', 'stationary'};
    sys.state = categorical({'sweep'}, states, 'protected', true);


    % Initialize noise generator
    gen = RandStream('simdTwister', 'Seed', sum(abs(field)));
    sys.noise = noise;
    sys.noise_gen = @(noise) noise*2*((rand(gen, 1, nT) - 0.5) ...
                                 + 1i*(rand(gen, 1, nT) - 0.5));


    % Add system to batch
    batch.systems(end + 1) = sys;

end