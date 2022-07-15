function sys = update_state(sys, series, cw_tol, var_tol)
% Update state of system at given series with CW tolerance and variance
% tolerance

    if series < sys.sweep_series
        sys.state(series) = categorical({'sweep'});

    elseif sys.cw_upper_diff(series) < cw_tol || sys.cw_lower_diff(series) < cw_tol
        sys.state(series) = categorical({'cw'});

    elseif sys.variance(series) < var_tol
        sys.state(series) = categorical({'stationary'});

    else
        sys.state(series) = categorical({'fluctuating'});
        
    end

end