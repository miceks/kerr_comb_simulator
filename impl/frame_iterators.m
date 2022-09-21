function [time_it, freq_it, res_it] = frame_iterators(batch, idx, time_ax, freq_ax, res_ax, conv_ax)
% Produce iterators used to plot the time (power), freqeuncy (spectrum) and
% CW resonance reference at each value in the field buffer of the system at
% idx. The iterators accept the number of frames to advance, and advancing
% one advances all by the same amount.

    sys = batch.systems(idx);

    t = (sys.T/sys.nT)*(-sys.nT/2:sys.nT/2-1);              % Time grid
    omega = fftshift((2*pi/sys.T)*(-sys.nT/2:sys.nT/2-1));  % Frequency grid
    rho = sqrt(1 - sys.alpha)*exp(-sys.alpha/2);
 
    switch batch.type
        case 'Ikeda'
            [cw_upper, cw_middle, cw_lower] = cw_Ikeda(sys.alpha, sys.power_end, ...
                                                       batch.cw_newton_tol, ...
                                                       batch.cw_curve_tol);
        case 'LLE'
            [cw_upper, cw_middle, cw_lower] = cw_LLE(sys.alpha, sys.power_end, ...
                                                     batch.cw_power_levels);
    end
    
    % Determine plot limits
    cw_power = [cw_upper(2, :) cw_middle(2, :) cw_lower(2, :)];
    cw_detun = [cw_upper(1, :) cw_middle(1, :) cw_lower(1, :)];

    field_peak = max(abs(sys.field_buffer).^2, [], 'all');
    cw_peak    = sys.alpha/(1 - rho)^2*sys.power_end;   % Detuning and power
                                                        % peak of CW reference.
    power_peak = max(field_peak, cw_peak);

    detun_lim(1) = min([cw_detun(find(cw_power > 0.01*cw_peak, 1, 'first')), ...
                        sys.detun_start, sys.detun_end]);
    detun_lim(2) = max([cw_detun(find(cw_power > 0.01*cw_peak, 1, 'last')), ...
                        cw_peak, sys.detun_start, sys.detun_end]);
    power_lim(1) = 0;
    power_lim(2) = 1.2*power_peak;


    % Initialize series varibales and map states to colors
    series_start = sys.series - sum(~isnan(sys.field_buffer(:, 1))) + 1;
    series_end   = sys.series;
    it_series    = series_start;

    states = {'sweep',   'fluctuating', 'cw',      'stationary'};
    colors = {'#4DBEEE', '#0070ff',     '#7E2F8E', '#77AC30'   };

    field_col = categorical(sys.state(series_start:series_end), states, colors);


    % Plot variance in segements with color of each corresponding to state
    hold(conv_ax, 'on')

    state_switch_idx = find([1 diff(double(sys.state)) 1]);
    for idx = 1:(length(state_switch_idx) - 1)
        first = state_switch_idx(idx);
        last = state_switch_idx(idx + 1) - 1;
        plot(conv_ax, first:last, sys.variance(first:last), ...
             'Color', colors{double(sys.state(first))}, 'LineWidth', 1);
    end

    text(conv_ax, 0.03, 0.95, ...
         sprintf("ID: %d", sys.id), ...
         'Units', 'normalized', ...
         'FontSize', 10, ...
         'FontWeight', 'bold')

    set(conv_ax, 'YScale', 'log')


    % Iterator producing frames for the time view of the field power
    function done = time_view(frames)

        done = (it_series > series_end);
        if done
            return
        end

        series_mod = mod(it_series - 1, sys.buffer_size) + 1;
        field = sys.field_buffer(series_mod, :);

        plot(time_ax, t, abs(field).^2, ...
            'Color', string(field_col(series_mod)), ...
            'LineWidth', 1)
        axis(time_ax, [t(1) t(sys.nT) power_lim(1) power_lim(2)]);

        text(time_ax, 0.03, 0.95, ...
             sprintf("ID: %d   Series: %d/%d", sys.id, it_series, series_end), ...
             'Units', 'normalized', ...
             'FontSize', 10, ...
             'FontWeight', 'bold')

        it_series = it_series + frames;
    end


    % Iterator producing frames for the field spectrum
    function done = freq_view(frames)
        
        done = (it_series > series_end);
        if done
            return
        end

        series_mod = mod(it_series - 1, sys.buffer_size) + 1;
        field = sys.field_buffer(series_mod, :);
        field_hat = fft(field);

        stem(freq_ax, omega, 10*log10(abs(field_hat).^2), ...
             'Marker', 'none', ...
             'Color', string(field_col(series_mod)), ...
             'BaseValue', -100);
        axis(freq_ax, [min(omega)/2 max(omega)/2 -100 10]);

        text(freq_ax, 0.03, 0.95, ...
             sprintf("ID: %d   Series: %d/%d", sys.id, it_series, series_end), ...
             'Units', 'normalized', ...
             'FontSize', 10, ...
             'FontWeight', 'bold')

        it_series = it_series + frames;
    end


    % Iterator producing frames for CW resonance curve
    function done = res_view(frames)
    
        done = (it_series > series_end);
        if done
            return
        end

        % Recalculate resonance curve if needed
        if sys.state(it_series) == "sweep"

            power = interp1([1 sys.sweep_series], [sys.power_start sys.power_end], ...
                            it_series, 'linear', sys.power_end);

            switch batch.type
                case 'Ikeda'
                    [branch{1}, branch{2}, branch{3}] = ...
                        cw_Ikeda(sys.alpha, power, batch.cw_newton_tol, batch.cw_curve_tol);
                case 'LLE'
                    [branch{1}, branch{2}, branch{3}] = ...
                        cw_LLE(sys.alpha, power, batch.cw_power_levels);
            end

        else
            branch{1} = cw_upper;
            branch{2} = cw_middle;
            branch{3} = cw_lower;
        end


        cla(res_ax)
        hold(res_ax, 'on')

        mi_upper = branch{1}(2, :) > batch.alpha & branch{1}(1, :) <= 2*branch{1}(2, :) & batch.s < 0;
        mi_lower = branch{3}(2, :) > batch.alpha & branch{3}(1, :) >= 2*branch{3}(2, :) & batch.s > 0;

        upper_A = find(mi_upper, 1, 'first');
        upper_B = find(mi_upper, 1, 'last');

        if isempty(upper_A) && isempty(upper_B)
            upper_A = length(mi_upper);
            upper_B = length(mi_upper);
        end

        % Stable upper A
        plot(res_ax, branch{1}(1, 1:upper_A), branch{1}(2, 1:upper_A), ...
             'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560])

        % Unstable upper
        plot(res_ax, branch{1}(1, upper_A:upper_B), branch{1}(2, upper_A:upper_B), ... 
             'LineWidth', 1.5, 'Color', [0 0.4470 0.7410])

        % Stable Upper B
        plot(res_ax, branch{1}(1, upper_B:end), branch{1}(2, upper_B:end), ... 
             'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560])


        if ~isequal(branch{1}, branch{3})
            % Middle
            plot(res_ax, branch{2}(1, :), branch{2}(2, :), ...      
                 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.5)

            % Unstable lower
            plot(res_ax, branch{3}(1,  mi_lower), branch{3}(2,  mi_lower), ...  
                 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410])

            % Stable lower
            plot(res_ax, branch{3}(1, ~mi_lower), branch{3}(2, ~mi_lower), ...  
                 'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560]) 
        end

        hold(res_ax, 'off')

        
        detun = interp1([1 sys.sweep_series], [sys.detun_start sys.detun_end], ...
                        it_series, 'linear', sys.detun_end);
        xline(res_ax, detun);
        yline(res_ax, max(0, sys.cw_upper_power(it_series)));
        yline(res_ax, max(0, sys.cw_lower_power(it_series)));

        axis(res_ax, [detun_lim power_lim]);

        
        text(res_ax, 0.03, 0.95, ...
             sprintf("ID: %d   Series: %d/%d", sys.id, it_series, series_end), ...
             'Units', 'normalized', ...
             'FontSize', 10, ...
             'FontWeight', 'bold')

        it_series = it_series + frames;
    end


    % Return iterators
    time_it = @time_view;
    freq_it = @freq_view;
    res_it  = @res_view;

    % Plot first frame
    time_it(0);
    freq_it(0);
    res_it(0);

end