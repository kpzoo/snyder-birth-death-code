% Function to test if the empirical distributions from MCMC (samples) and
% the Snyder solution are not from the same distribution (h = 1)
function [h, p, sampSny] = compareSnyMCMC3(samples, qmarg, xset, numRV, sigLevel, interpMeth, treatDiscrete)

% Assumptions and modifications
% - allows sampling from a discrete distribution or a continuous equivalent
% - works with trapz based qmarg from marginalise2
% - allows for setting interp method
% - test for each parameter marginals
% - can choose significance level of test with sigLevel
% - h = 0 indicates the null hypothesis cannot be rejected at sigLevel
% - p is the probability of observing a more extreme test statistic value
% - the test statistic is calculated inside kstest2


% Check formatting consistent
if size(samples, 2) ~= numRV && length(qmarg) ~= numRV
    error('Marginals and samples do not match numRV');
else
    % No. of samples per parameter and initialise Snyder samples
    sampSny = zeros(size(samples));
end

if ~treatDiscrete
    % Interpolate the Snyder marginals to more points on parameter space
    nPts = 100;
    xm = cell(1, numRV);
    cdfm = cell(1, numRV);
    for i = 1:numRV
        xm{i} = linspace(min(xset{i}), max(xset{i}), nPts);
        cdfm{i} = interp1(xset{i}, cumtrapz(xset{i}, qmarg{i}), xm{i}, interpMeth);
    end
    
    % Generate nSamp samples from the Snyder marginal cdfs using uniform rvs
    for i = 1:numRV
        % This takes first unique value of any set
        [cdfm{i}, iduniq] = unique(cdfm{i});
        xm{i} = xm{i}(iduniq);
        if length(iduniq) <= 2
            error(['CDF has too few unique values for variable i = ' num2str(i)]);
        end
    end
    
    u = rand(size(samples));
    for i = 1:numRV
        % Interpolate from cdf the xset value corresponding to u values
        sampSny(:, i) = interp1(cdfm{i}, xm{i}, u(:, i), 'linear', 'extrap');
    end
    % Check samples sensible
    if sum(any(isnan(sampSny))) ~= 0
        error('Interpolation of cdf gave NaN values');
    end
else
    % In the discrete case the probabilities can take frequentist approach
    for i = 1:numRV
        nSamps = size(samples, 1);
        % Frequentist approximation to probabilities
        sampSny(:, i) = datasample(xset{i}, nSamps, 'Weights', qmarg{i});
    end
    
end

% Two sample Kolmogorov-Smirnoff test, null hypothesis that both
% distributions are the same
h = -1*ones(1, numRV);
p = -1*ones(1, numRV);
for i = 1:numRV
    % Run K-S test for each parameter (marginal)
    [h(i), p(i)] = kstest2(samples(:, i), sampSny(:, i), 'Alpha', sigLevel);
end