function [cw_upper, cw_middle, cw_lower] = cw_LLE(alpha, P_in, N)
% Calculate points on CW resonance curve for LLE with cavity loss alpha 
% and pump power P_in. The points are partitioned into an upper, middle and 
% lower branch (equal if no bistability occurs) with detuning values in the top 
% row and intracavity power in the second row. The number of intracavity
% power slices N used to find corresponding detuning values governs the
% accuracy.

    detun_bound = max(pi, P_in/alpha);
    
    r = roots([1, 2*detun_bound, (alpha^2 + detun_bound^2), -alpha*P_in]);
    power_min = min(abs(real(r)));
    power_max = P_in/alpha;

    power_inc = linspace(power_min, power_max, N);
    power_dec = linspace(power_max, power_min, N);
    
    % Solve for corresponding detunings at each power level
    detun_inc = power_inc - sqrt(abs(P_in./power_inc - alpha)*alpha);
    detun_dec = power_dec + sqrt(abs(P_in./power_dec - alpha)*alpha);

    points = [detun_inc detun_dec;
              power_inc power_dec];
    points(:, N) = [];
    
    % Partition curve according to possible branch switches
    idx_switch_1 = find(islocalmax(points(1, :)), 1);
    idx_switch_2 = find(islocalmin(points(1, :)), 1);
    
    if ~isempty([idx_switch_1 idx_switch_2])
        cw_upper  = points(:, 1:idx_switch_1);
        cw_middle = points(:, idx_switch_1:idx_switch_2);
        cw_lower  = points(:, idx_switch_2:end);

    else
        cw_upper  = points;
        cw_middle = points;
        cw_lower  = points;

    end

end