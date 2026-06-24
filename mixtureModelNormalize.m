function [CorrM_norm] = mixtureModelNormalize(CorrM)
% mixtureModelNormalize  Gaussian-Gamma mixture model normalization of a
%                        connectivity matrix.
%
%   [CorrM_norm, params] = mixtureModelNormalize(CorrM)
%
%   Input:
%       CorrM - square connectivity matrix (e.g. Fisher z-transformed
%               correlation matrix). Assumed symmetric.
%
%   Output:
%       CorrM_norm - matrix normalized by the fitted Gaussian (null)
%                    component: edges expressed as how far they sit from
%                    the noise distribution.
%       params     - struct with the fitted Gaussian parameters and
%                    z-scoring constants, for inspection/QC.

    N = size(CorrM, 1);
    if size(CorrM, 2) ~= N
        error('mixtureModelNormalize:notSquare', ...
              'Input must be a square matrix.');
    end

    % Flatten to the edge vector the MM operates on
    CorrMf = flattenMat(CorrM);

    % Fitting parameters
    maxiters = 100;
    tol      = 1e-6;
    MM       = 2;   % 1 Gaussian + 2 Gammas

    % Step 1: z-score the edge values
    prov2_mean = mean(CorrMf);
    prov2_std  = std(CorrMf);
    prov2_norm = (CorrMf - prov2_mean) / prov2_std;

    % Step 2: fit Gaussian-Gamma mixture
    [stmu1, stv1, stmu2, stv2, stmu3, stv3, stPI, lik, it, resp] = ...
        mmfit3(prov2_norm, maxiters, tol, MM);

    % Step 3: Gaussian component params (normalized space)
    m_gaus_norm   = stmu1;
    var_gaus_norm = stv1;

    % Step 4: denormalize back to original scale
    m_gaus   = prov2_std * m_gaus_norm + prov2_mean;
    var_gaus = (prov2_std^2) * var_gaus_norm;

    % Step 5: normalize edges by the Gaussian component only
    CorrMf_normalized = (CorrMf - m_gaus) / sqrt(var_gaus);

    % Rebuild the square matrix at the input's own size
    CorrM_norm = recon3dMat(CorrMf_normalized, N);

    % Optional diagnostics
    params = struct( ...
        'm_gaus', m_gaus, 'var_gaus', var_gaus, ...
        'm_gaus_norm', m_gaus_norm, 'var_gaus_norm', var_gaus_norm, ...
        'zscore_mean', prov2_mean, 'zscore_std', prov2_std, ...
        'PI', stPI, 'loglik', lik, 'iters', it, ...
        'gamma1', [stmu2, stv2], 'gamma2', [stmu3, stv3]);
end
