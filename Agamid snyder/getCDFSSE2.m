% Function to calculate theoretical CDFs and compare to empirical
function [cdf, ss_err] = getCDFSSE2(numFiles, t, nInt, T, xhatSet, pinterp, data)

% Assumptions and modifications
% - xhat now used as a cell input

% Store theoretical cdf and square errors
cdf = zeros(numFiles, length(t));
ss_err = zeros(1, numFiles);

% Calculate theoretical CDFs from each model according to Paradis
for j = 1:numFiles
    % Use estimate to get best cdf
    pi = zeros(size(t));
    for i = 1:length(t)
        % Define the decreasing time span over which integrals are performed
        tstart = t(i);
        tspan = tstart + (0:nInt-1).*(T - tstart)/nInt;
        
        % Rate functions r(t, s) = e^(int mut-lamt) and integrand of PtT
        mut = data{j}.mut;
        lamt = data{j}.lamt;
        rhotT = data{j}.rhotT;
        rst = @(x, t1, t2) exp(rhotT(x, t1, t2));
        inttT = @(x, t1, t2) mut(x, t2).*rst(x, t1, t2);
        
        % Snyder estimates from each dataset
        xhat = xhatSet{j};
        
        % Function values with time for integration
        intSet = inttT(xhat, tstart, tspan);
        rho0t = rhotT(xhat, 0, tstart);
        
        % Inline version of trapz for increased speed, assumes row vectors
        trapzIn = diff(tspan)*(intSet(1:nInt-1) + intSet(2:nInt))'/2;
        PtT2 = (1 + trapzIn)^(-2);
        % CDF point value pi(t(i))
        pi(i) = PtT2*exp(-rho0t)*lamt(xhat, tstart);
    end
    % Get complete CDF for this parameter set
    cdf_p = cumtrapz(t, pi, 2)/trapz(t, pi, 2);
    cdf(j, :) = cdf_p;
    
    % Get sum of square errors
    ss_err(j) = sum((cdf_p - pinterp).^2);
end