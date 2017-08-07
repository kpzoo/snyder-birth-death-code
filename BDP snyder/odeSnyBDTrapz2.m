% Snyer ODE set for time varying birth-death process

% Assumptions
% - uses inline version of trapz in getNeeRatesTrapz2
% - uses trapz integral so faster
% - uses function as arguments
% - assumes time varying birth-death from Nee 1994
% - uses equivalent birth rate for rate matrix from Nee 1994
% - y is a column vector, ts a vector just included for ode113

function dy = odeSnyBDTrapz2(ts, y, xsetMx, numRV, Tsp, nLinCurr, rhotT, lamt, mut)

% Ensure input dimension correct
sz = size(y);
if sz(2) ~= 1
    error('y is not a column vector');
end

% Diagonal of time dependent rate matrix for reconstructed birth (Nee 1994)
rateDiag = getNeeRatesTrapz2(ts, xsetMx, numRV, Tsp, nLinCurr, rhotT, lamt, mut);

% Solve non-linear differential equation set - RV filtering
nonLinDiag = rateDiag*y;
dy = y'.*(-rateDiag + nonLinDiag);

% Ensure output is column vector like the input
dy = dy';