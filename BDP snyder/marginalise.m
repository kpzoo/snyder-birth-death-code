% Function to marginalise a multivariate joint distribution
function [qmarg, probSums] = marginalise(numRV, IDMx, q, mi)

% Assumptions and modifications
% - the IDMx is properly formatted for id sums

% Marginal probability cell
qmarg = cell(1, 1);
probSums = zeros(1, numRV);

% Loop across RVs and get ids to sum for each
for i = 1:numRV
    qcomps = zeros(1, mi(i));
    for j = 1:mi(i)
        % All the entries of the ith variable of id value j
        qcomps(j) = sum(q(IDMx(i, :) == j));
    end
    qmarg{i} = qcomps;
    % Get sums of probability vectors as test feature
    probSums(i) = sum(qmarg{i});
end