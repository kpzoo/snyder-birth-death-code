% Function to sample from a constant birth-death process according to
% Hartmann et al 2008
function [T, tbranch] = genBirDeaConst(n, sig, rho)

% Assumptions and modifications
% - inputs - no. extant lineages, n
% - outputs - speciation times
% - only for constant birth and death rates

% Sample n uniform variables in [0, 1] for (n-1) times and one origin
r = rand(1, n-1);
r0 = rand;

% Draw a tree age
T = (1/sig)*log((1 - rho*(r0^(1/n)))/(1 - r0^(1/n)));

% Calculate branching times
num = 1 - rho*exp(-sig*T) - rho*(1 - exp(-sig*T))*r;
den = 1 - rho*exp(-sig*T) - (1 - exp(-sig*T))*r;
tbranch = (1/sig)*log(num./den);
