% Function to calculate time-varying Snyder BDP log likelihood
function logf = likSnyLogTimeVary(ts, xin, nLin, lamt, inttT)

% Assumptions and modifications
% - calculates for a given parameter set xin
% - uses log vs log10
% - nLin must match ts speciation times 
% - starting from 2 lineages so tspec(1) cannot be 0


% Length of speciation times sequence
lenS = length(ts);
if lenS ~= length(nLin)
    error('Times do not match lineages in length');
end

% Max speciation time - first time with n lineages
T = max(ts);

% Likelihood variable, as it is between events it is 1 less in length
H = zeros(1, lenS-1);

% Number of integration points
nInt = 500;

% Minimise inputs to function
PtT = @(t1, t2) calcPtT(t1, t2, inttT, xin);

% Loop across events and calculate log likelihood
for i = 1:lenS-1
    % Inputs to likelihood
    nlin = nLin(i);
    t1 = ts(i);
    t2 = ts(i+1);
    % Time span for integration
    dt = t1 + (0:nInt-1).*(t2 - t1)/nInt;
    
    % Pre-input start time and then evaluate at all end times
    PtT1 = @(tst) PtT(tst, T);
    
    % Rate function across time steps without lineage dependence
    rate = lamt(xin, dt).*arrayfun(PtT1, dt')';
    
    % Continous likelihood component
    Hc = nlin*trapz(dt, rate);
    
    % Update likelihood component
    Hu = log(nlin*lamt(xin, t2)*PtT(t2, T));
    % Overall likelihood exponent
    H(i) = Hu - Hc;
end

% Get desired log likelihood
logf = sum(H); 
