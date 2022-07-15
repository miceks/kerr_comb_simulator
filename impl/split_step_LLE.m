function sys = split_step_LLE(sys)
% Execute one series of LLE roundtrips using the split step scheme.

    dR = 1/sys.N;           % Step size
    Q = round(sys.R/dR);    % Number of steps in series
    omega = fftshift((2*pi/sys.T)*(-sys.nT/2:sys.nT/2-1)); % Frequency grid

    % Retrieve current field
    psi = sys.field_buffer(mod(sys.series - 1, sys.buffer_size) + 1, :);

    % Update parameter sweeps
    sys.detun = interp1([1 sys.sweep_series], ...
                        [sys.detun_start sys.detun_end], ...
                        sys.series, 'linear', sys.detun_end);
    sys.power = interp1([1 sys.sweep_series], ...
                        [sys.power_start sys.power_end], ...
                        sys.series, 'linear', sys.power_end);

    P_hat = fft(sqrt(sys.power)*ones(1, sys.nT)); % Frequency domain power
    D1 = -sys.alpha - 1i*sys.detun + 1i*sys.s*omega.^2;
    D2 = exp(D1*dR);
    
    % Add noise
    psi = psi + sys.noise_gen(sys.noise(end));

    % Initial linear half-step
    psi = ifft(fft(psi).*sqrt(D2) + ...
               (sqrt(sys.alpha)*P_hat./D1).*(sqrt(D2) - 1));

    for q = 1:Q
        % Nonlinear step
        psi = psi.*exp(1i*dR*abs(psi).^2);

        % Linear step
        psi = ifft(fft(psi).*D2 + (sqrt(sys.alpha)*P_hat./D1).*(D2 - 1));
    end

    % Inverse linear half-step
    psi = ifft(fft(psi).*1./sqrt(D2) + ...
               (sqrt(sys.alpha)*P_hat./D1).*(1./sqrt(D2) - 1));
    
    % Store field
    sys.series = sys.series + 1;
    sys.field_buffer(mod(sys.series - 1, sys.buffer_size) + 1, :) = psi;

    % Update statistics
    P2 = abs(psi).^2;
    P1 = abs(sys.field_buffer(mod(sys.series - 2, sys.buffer_size) + 1, :)).^2;
    power_norm = sum(P2);

    sys.variance(sys.series)      = sum(abs(P2 - P1))/power_norm;
    sys.noise(sys.series)         = sys.noise(end);
    sys.cw_upper_diff(sys.series) = sum(abs(P2 - sys.cw_upper_power(sys.series))/power_norm);
    sys.cw_lower_diff(sys.series) = sum(abs(P2 - sys.cw_lower_power(sys.series))/power_norm);
    
end