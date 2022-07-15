function sys = split_step_Ikeda(sys)
% Execute one series of Ikeda map roundtrips using the split step scheme.

    dz = 1/sys.N;  % Step size
    omega = fftshift((2*pi/sys.T)*(-sys.nT/2:sys.nT/2-1));  % Frequency grid

    % Retrieve current field
    psi = sys.field_buffer(mod(sys.series - 1, sys.buffer_size) + 1, :);

    % Update parameter sweeps
    sys.detun = interp1([1 sys.sweep_series], ...
                        [sys.detun_start sys.detun_end], ...
                        sys.series, 'linear', sys.detun_end);
    sys.power = interp1([1 sys.sweep_series], ...
                        [sys.power_start sys.power_end], ...
                        sys.series, 'linear', sys.power_end);

    D = exp((-sys.alpha/2 + 1i*sys.s*omega.^2)*dz);

    % Add noise
    psi = psi + sys.noise_gen(sys.noise(end));
    
    % Linear half-step
    psi = ifft(fft(psi).*sqrt(D));

    for r = 1:sys.R
        for n = 1:sys.N
            % Nonlinear step
            psi = psi.*exp(1i*abs(psi).^2*dz);

            % Linear step
            psi = ifft(fft(psi).*D);
        end

        % Boundary condition
        psi = sqrt(1 - sys.alpha)*exp(-1i*sys.detun)*psi + ...
              sqrt(sys.alpha*sys.power);
    end

    % Inverse linear half-step
    psi = ifft(fft(psi).*1./sqrt(D));

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