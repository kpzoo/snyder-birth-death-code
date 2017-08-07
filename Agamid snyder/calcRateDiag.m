% Function allows cellfun calculation of rateDiag in getNeeRatesTrapz2
function rdiag = calcRateDiag(xk, inttT, tstart, tspan, lam1, nInt, ncurr)

% Get integrand vector for integration
intSet = inttT(xk, tstart, tspan);

% Trapz version of integral for PtT and rate diagonal
trapzIn = diff(tspan)*(intSet(1:nInt-1) + intSet(2:nInt))'/2;
PtT = (1 + trapzIn)^(-1);
rdiag = ncurr*lam1(xk)*PtT;