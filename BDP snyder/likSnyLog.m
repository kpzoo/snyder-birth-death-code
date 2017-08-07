% Function to calculate Snyder constant BDP log likelihood
function logf = likSnyLog(ts, sig, rho, nLin)

% Assumptions and modifications
% - nLin must match ts speciation times e.g. if starting from 2 lineages
% then tspec(1) cannot be 0

% Length of speciation times sequence
lenS = length(ts);
if lenS ~= length(nLin)
    error('Times do not match lineages in length');
end

% Max speciation time - first time with n lineages
T = max(ts);

% Likelihood variable, as it is between events it is 1 less in length
H = zeros(1, lenS-1);

% Loop across events and calculate log likelihood
for i = 1:lenS-1
    % Inputs to likelihood
    nlin = nLin(i);
    t1 = ts(i);
    t2 = ts(i+1);
    
    % Continous likelihood component
    Hc = nlin*sig*(t2 - t1) + nlin*log((1 - rho*exp(sig*(t1 - T)))./(1- rho*exp(sig*(t2 - T))));
    % Update likelihood component
    Hu = log(nlin*sig./(1 - rho*exp(sig*(t2 - T))));
    % Overall likelihood exponent
    H(i) = Hu - Hc;
end

% Get desired log10 likelihood
L = exp(sum(H)); % raw likelihood
logf = log10(L);
