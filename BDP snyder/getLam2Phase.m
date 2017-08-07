% Simple function for integral of mut-lamt of 2 phase case
function lam = getLam2Phase(x, tx)

% Check length, even though x(4) never used as death parameter
if length(x) ~= 4
    error('Incorrect parameter input size');
end

% Assume row vector input but will make output consistent in dimension
tsz = size(tx);
if tsz(1) > 1 && tsz(2) == 1
    % Make column vector row vector
    wasCol = 1;
    tx = tx';
elseif tsz(1) > 1 && tsz(2) > 1
    error('Inputs must be vectors');
else
    wasCol = 0;
end

% Use vector of times
tph1 = tx < x(3);
tph2 = tx >= x(3);

% Check timing errors
if all(~tph1) && all(~tph2)
    error('Incorrect time partition');
end

% Calculate birth rate
if all(~tph1)
    % Phase 2 only
    lam = x(2)*ones(size(tx));
elseif all(~tph2)
    % Phase 1 only
    lam = x(1)*ones(size(tx));
else
    lam = x(1)*tph1 + x(2)*tph2;
end

if wasCol
    lam = lam';
end
