function [cw_upper, cw_middle, cw_lower] = cw_Ikeda(alpha, P_in, ntol, ctol)
% Calculate points on CW resonance curve for Ikeda map with cavity loss alpha 
% and pump power P_in. The points are partitioned into an upper, middle and 
% lower branch (equal if no bistability occurs) with detuning values in the top 
% row and intracavity power in the second row. The newton iteration uses ntol 
% as tolerance for each point, and the recursive addition of more points is 
% determined by the curve tolerance ctol.

    rho = sqrt(1 - alpha)*exp(-alpha/2);
    phi = @(delta, P) -delta + P*(1 - exp(-alpha))/alpha;

    % cw solutions satisfy F = 0
    F  = @(delta, P) ((1 - rho)^2 + 4*rho*sin(phi(delta, P)/2).^2).*P - alpha*P_in;
    dF = @(delta, P) [-2*rho*sin(phi(delta, P))*P;
                      (1 - rho)^2 + 4*rho*sin(phi(delta, P)/2)^2 + 2*rho*sin(phi(delta, P))*(1 - exp(-alpha))/alpha*P];
     
    detun_peak = (1 - exp(-alpha))/(1 - rho)^2*P_in;
    power_peak = alpha/(1 - rho)^2*P_in;

    % Use resonance peak and points at +-pi as solution anchor points
    p_peak  = [detun_peak; power_peak];
    p_min = newton([-pi; 0], F, dF, [0; 1], ntol);
    p_max = newton([ pi; 0], F, dF, [0; 1], ntol);

    % Produce interior points
    ctol = ctol*sqrt(P_in); % Scale with pump power to normalize number of points.
    points = [p_min, ...
              interior(p_min, p_peak, F, dF, ntol, ctol), ...
              p_peak, ...
              interior(p_peak, p_max, F, dF, ntol, ctol), ...
              p_max];

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



function pts = interior(pA, pB, F, dF, n_tol, c_tol)
% Recursively return points between pA and pB based on local curvature.

    % Produce new middle point by searching in direction normal to the vector
    % connecting pA and pB
    pC = (pA + pB)/2;
    v = pB - pA;       
    v = [v(2); -v(1)];   
    v = v/norm(v);  
    [pC, h] = newton(pC, F, dF, v, n_tol);
    
    % Terminate recursion when step size in normal direction becomes small,
    % otherwise produce more interior points on either side
    if abs(h) < c_tol
        pts = pC;
        
    else
        pts = [interior(pA, pC, F, dF, n_tol, c_tol), ...
               pC, ...
               interior(pC, pB, F, dF, n_tol, c_tol)];
    end

end


function [p, h] = newton(p0, F, dF, dir, n_tol)
% Newton-Raphson using directional derivative along dir

    h = 0;
    p = p0;
    step = -F(p0(1), p0(2))/(dF(p0(1), p0(2))'*dir);

    while abs(step) > n_tol
        h = h + step;
        p = p0 + h*dir;
        step = -F(p(1), p(2))/(dF(p(1), p(2))'*dir);
    end

end