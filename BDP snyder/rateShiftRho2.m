% Simple function to calculate rhotT for the rateshift BDP
function rhotT = rateShiftRho2(x, t1, t2, changeTime)

% Assumptions and modifications
% - removed global variable
% - the conditional must be specified here and in tessSnyder <---------

% Calculation based on shift at tcond
if t2 < changeTime
    % No shift
    rhotT = -x(1)*(t2 - t1);
elseif t1 >= changeTime
    % Only in shifted state
    rhotT = -x(3)*(t2 - t1);
else
    % Part of integral before and part after shift
    rhotT = -x(1)*(changeTime - t1) -x(2)*(t2 - changeTime);
end