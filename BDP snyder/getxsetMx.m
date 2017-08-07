% Function to calculate xsetMx, a Kronecker space across xset values that
% combinatorically explores it
function [xset, m, xsetMx, IDMx] = getxsetMx(numRV, xmin, xmax, mi, logSP)

% Assumptions and modifications
% - allowed for log10 spacing between xmin and xmax with vector of bools

% Calculate parameter space
m = prod(mi);
xset = cell(1, 1);
for i = 1:numRV
    if ~logSP(i)
        % Linear spacing 
        xset{i} = linspace(xmin(i), xmax(i), mi(i));
    else
        % Logarithmic spacing
        xset{i} = logspace(log10(xmin(i)), log10(xmax(i)), mi(i));
    end
end

% Create a matrix of identifiers to tell which xset{i} values are used for
% each entry lam(t) calculations
IDMx = zeros(numRV, m);
% Initialise with first variable which has no element repetitions
idxset = 1:mi(1);
IDMx(1, :) = repmat(idxset, 1, m/mi(1));
for i = 2:numRV
    % For further variables numReps gives the number of set repetitions
    % while kronVec gives the number of element repetitions
    idxset = 1:mi(i);
    numReps = m/prod(mi(1:i));
    kronVec = ones(1, prod(mi(1:i-1)));
    IDMx(i, :) = repmat(kron(idxset, kronVec), 1, numReps);
end

% Get the values corresponding to the matrix
xsetMx = zeros(numRV, m);
for i = 1:numRV
    xsetMx(i, :) = xset{i}(IDMx(i, :));
end