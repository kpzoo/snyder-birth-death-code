% Function to marginalise a multivariate joint distribution
function qmarg = marginalise2(numRV, IDMx, q, mi, xset)

% Assumptions and modifications
% - removed probSums as want trapz version if distr truly continuous
% - the IDMx is properly formatted for id sums

% Marginal probability cell
qmarg = cell(1, 1);

% Loop across RVs and get ids to sum for each
for i = 1:numRV
    qcomps = zeros(1, mi(i));
    for j = 1:mi(i)
        % All the entries of the ith variable of id value j
        qcomps(j) = sum(q(IDMx(i, :) == j));
    end
    % Normalise each marginal based on xset values
    Z = trapz(xset{i}, qcomps);
    qcomps = qcomps/Z;
    % Assign marginal
    qmarg{i} = qcomps;
end

