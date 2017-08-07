% Calculate the prob a lineage at t1 survives to t2 (>= 1 descendant)
function PtT = calcPtT(t1, t2, inttT, xin)

% Assumptions and modifications
% - expects function inttT - the integrand, contains BDP rate functions

% Define time increments for integral
nInt = 500;
tspan = t1 + (0:nInt-1).*(t2 - t1)/nInt;

% Get integrand vector for integration
intSet = inttT(xin, t1, tspan);

% Trapz version of integral for PtT and rate diagonal
trapzIn = diff(tspan)*(intSet(1:nInt-1) + intSet(2:nInt))'/2;
PtT = (1 + trapzIn)^(-1);


