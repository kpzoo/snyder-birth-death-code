% Function to calculate time-varying Hohna BDP log likelihood
function logf = likHohnaLogTimeVary(ts, xin, nLin, lamt, inttT, rst)

% Assumptions and modifications
% - uses likelihood fro Hohna TESS tutorial
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

% Minimise inputs to function
PtT = @(t1, t2) calcPtT(t1, t2, inttT, xin);

% Prob of 1 descendant at T
p1 = @(t1, t2) (PtT(t1, t2).^2).*rst(xin, t1, t2);

% Hohna pre-factor (not in sum) with lineages
f0 = log((p1(0, T).^2)./(PtT(0, T).^2)) +log(prod(nLin(1:end-1)));

% Hohna lineage based factors
fprod = zeros(1, lenS-1);
for i = 1:lenS-1
    t1 = ts(i);
    t2 = ts(i+1);
    % Non-lineage dependent part of product in raw likelihood
    fprod(i) = log(lamt(xin, t2)) + log(p1(t2, T));
end

% Get desired log likelihood
logf = f0 + sum(fprod); 
