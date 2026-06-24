function [stmu1, stv1, stmu2, stv2, stmu3, stv3, stPI, lik, it, resp] = mmfit3(x, maxiters, tol, MM)
    % MM=2: Gaussian + 2 Gammas
    % MM=3: Gaussian + 2 Inverse-Gammas
    
    % Initialize parameters
    init_mu1 = 0;
    init_v1 = 1;
    init_mu2 = 3;
    init_v2 = 1;
    init_mu3 = -3;
    init_v3 = 1;
    
    init_PI = [1/3; 1/3; 1/3];
    
    % Initialize Gamma parameters based on MM
    if MM == 2
        init_a1 = alphaGm(init_mu2, init_v2);
        init_b1 = betaGm(init_mu2, init_v2);
        init_a2 = alphaGm(-init_mu3, init_v3);
        init_b2 = betaGm(-init_mu3, init_v3);
    elseif MM == 3
        init_a1 = alphaIG(init_mu2, init_v2);
        init_b1 = betaIG(init_mu2, init_v2);
        init_a2 = alphaIG(-init_mu3, init_v3);
        init_b2 = betaIG(-init_mu3, init_v3);
    end
    
    % Copy to tmp variables
    tmp_mu1 = init_mu1;
    tmp_v1 = init_v1;
    tmp_a1 = init_a1;
    tmp_b1 = init_b1;
    tmp_a2 = init_a2;
    tmp_b2 = init_b2;
    tmp_PI = init_PI;
    
    % Identify positive and negative indices
    xneg = find(x < 1e-14);
    xpos = find(x > -1e-14);
    eps = 1e-14;
    
    % Initialize likelihood
    real_lik = zeros(maxiters + 2, 1);
    it = 0;
    
    % First E-step
    pGa = normpdf(x, tmp_mu1, sqrt(tmp_v1));
    pGa(pGa == 0) = eps;
    
    if MM == 2
        dum2 = gampdf(x, tmp_a1, tmp_b1);
        dum3 = gampdf(-x, tmp_a2, tmp_b2);
    elseif MM == 3
        dum2 = invgampdf(x, tmp_a1, tmp_b1);
        dum3 = invgampdf(-x, tmp_a2, tmp_b2);
    end
    
    % Mask tails
    dum2(xneg) = 0;
    dum3(xpos) = 0;
    
    % Calculate responsibilities
    D1 = tmp_PI(1) * pGa;
    D1(D1 < eps) = eps;
    D2 = tmp_PI(2) * dum2;
    D2(D2 < eps) = eps;
    D3 = tmp_PI(3) * dum3;
    D3(D3 < eps) = eps;
    
    D = D1 + D2 + D3;
    R1 = D1 ./ D;
    R2 = D2 ./ D;
    R3 = D3 ./ D;
    
    % M-step
    N = [sum(R1); sum(R2); sum(R3)];
    tmp_PI = N / sum(N);
    
    % Calculate likelihood
    real_lik(it + 1) = sum(log(tmp_PI(1)*pGa + tmp_PI(2)*dum2 + tmp_PI(3)*dum3));
    
    % EM iterations
    flag = 0;
    while flag == 0
        it = it + 1;
        
        % Update Gaussian parameters
        tmp_mu1 = sum(R1 .* x) / N(1);
        tmp_v1 = sum(R1 .* (x - tmp_mu1).^2) / N(1);
        
        % Update component 2 parameters
        tmp_mu2 = sum(R2 .* x) / N(2);
        tmp_v2 = sum(R2 .* (x - tmp_mu2).^2) / N(2);
        
        % Update component 3 parameters
        tmp_mu3 = sum(R3 .* x) / N(3);
        tmp_v3 = sum(R3 .* (x - tmp_mu3).^2) / N(3);
        
        % Update Gamma parameters
        if MM == 2
            tmp_a1 = alphaGm(tmp_mu2, tmp_v2);
            tmp_b1 = betaGm(tmp_mu2, tmp_v2);
            tmp_a2 = alphaGm(-tmp_mu3, tmp_v3);
            tmp_b2 = betaGm(-tmp_mu3, tmp_v3);
        elseif MM == 3
            tmp_a1 = alphaIG(tmp_mu2, tmp_v2);
            tmp_b1 = betaIG(tmp_mu2, tmp_v2);
            tmp_a2 = alphaIG(-tmp_mu3, tmp_v3);
            tmp_b2 = betaIG(-tmp_mu3, tmp_v3);
        end
        
        % E-step
        pGa = normpdf(x, tmp_mu1, sqrt(tmp_v1));
        pGa(pGa == 0) = eps;
        
        if MM == 2
            dum2 = gampdf(x, tmp_a1, tmp_b1);
            dum3 = gampdf(-x, tmp_a2, tmp_b2);
        elseif MM == 3
            dum2 = invgampdf(x, tmp_a1, tmp_b1);
            dum3 = invgampdf(-x, tmp_a2, tmp_b2);
        end
        
        dum2(xneg) = 0;
        dum3(xpos) = 0;
        dum2(isnan(dum2) | isinf(dum2)) = eps;
        dum3(isnan(dum3) | isinf(dum3)) = eps;
        
        D1 = tmp_PI(1) * pGa;
        D1(D1 < eps) = eps;
        D2 = tmp_PI(2) * dum2;
        D2(D2 < eps) = eps;
        D3 = tmp_PI(3) * dum3;
        D3(D3 < eps) = eps;
        
        D = D1 + D2 + D3;
        R1 = D1 ./ D;
        R2 = D2 ./ D;
        R3 = D3 ./ D;
        
        % M-step
        N = [sum(R1); sum(R2); sum(R3)];
        tmp_PI = N / sum(N);
        tmp_PI(tmp_PI < eps) = eps;
        
        % Calculate likelihood
        real_lik(it + 1) = sum(log(tmp_PI(1)*pGa + tmp_PI(2)*dum2 + tmp_PI(3)*dum3));
        
        % Check convergence
        if abs((real_lik(it + 1) - real_lik(it)) / real_lik(it)) < tol || it > maxiters
            flag = 1;
        end
    end
    
    % Return final parameters
    stmu1 = tmp_mu1;
    stv1 = tmp_v1;
    stmu2 = tmp_mu2;
    stv2 = tmp_v2;
    stmu3 = tmp_mu3;
    stv3 = tmp_v3;
    stPI = tmp_PI;
    lik = real_lik(1:it+1);
    resp = [R1, R2, R3];
end

% Helper functions
function a = alphaGm(mu, var)
    a = mu^2 / var;
end

function b = betaGm(mu, var)
    b = var / mu;
end

function a = alphaIG(mu, var)
    a = (mu^2 / var) + 2;
end

function b = betaIG(mu, var)
    b = mu * ((mu^2 / var) + 1);
end

function p = invgampdf(x, a, b)
    p = (b^a / gamma(a)) .* (1./x).^(a+1) .* exp(-b./x);
    p(x <= 0) = 0;
end
