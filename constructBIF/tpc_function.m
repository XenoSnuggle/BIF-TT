function result = tpc_function(ds, param)

    wc = param.wc;        % Cutoff frequency
    L = param.L;          % Number of modes (length of w and c arrays)
    % xi = param.xi;        % Parameter xi
    beta = param.beta;    % Inverse temperature (1/kT)


    wmax = 4 * wc;
    

    w = zeros(L, 1);      % Array for frequency modes
    c = zeros(L, 1);      % Array for coupling coefficients
    

    for i = 1:L
        w(i) = -wc * log(1 - (1 - exp(-wmax / wc)) * i / L);
        % c(i) = w(i) * sqrt(xi * wc / L * (1 - exp(-wmax / wc)));
        c(i) = w(i) * sqrt(wc / L * (1 - exp(-wmax / wc)));
    end
    

    result = 0 + 1i * 0;
    

    for i = 1:L
        real_part = cos(w(i) * ds) / tanh(beta * w(i) / 2);
        imag_part = -1i * sin(w(i) * ds);
        result = result + (c(i) * c(i) / w(i)) * (real_part + imag_part);
    end
    
    result = 0.5 * result;
end