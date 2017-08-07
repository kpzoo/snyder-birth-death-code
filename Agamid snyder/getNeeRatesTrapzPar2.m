% Function to calculate time dependent rate matrix diagonals for the
% time varying birth-death process according to Nee 1994
function rateDiag = getNeeRatesTrapzPar2(tstart, xsetMx, Tsp, ncurr, rhotT, lamt, mut)

% Assumptions and modifications
% - uses cellfun to speed up assignment of xsetCol
% - uses mat2cell and a pre-assign to lamt as lam1 
% - uses parfor to speed up assignment
% - added inline version of trapz and linspace
% - version of getNeeRates which uses trapz
% - uses function handle input
% - simplified to use integral function
% - Tsp is max speciation time
% - time input must be a scalar

% Attempt to speed up intT by defining here
inttT = @(x, t1, t2) mut(x, t2).*exp(rhotT(x, t1, t2));

% Decompose xsetMx into column and row cell arrays
gridSz = size(xsetMx);
xsetCol  = mat2cell(xsetMx, gridSz(1), ones(1, gridSz(2)));

% Calculate P(t, T) from Nee 1994 using exp(r(t, T)) = rst across parameter 
% space for one time set from t to Tsp
rateDiag = zeros(1, gridSz(2));

% Inline alternative to linspace
nInt = 200;
tspan = tstart + (0:nInt-1).*(Tsp - tstart)/nInt;

% Integrand to be used with trapz to get PtT
%intSet = zeros(size(tspan));

% Can pre-assign 1 input to lamt 
lam1 = @(x)lamt(x, tstart);
% Pre-assign times in inttT so only function of params to estimate
intX = @(x) inttT(x, tstart, tspan);

% Replace loop with cellfun
calcX = @(x) calcRateDiag(x, inttT, tstart, tspan, lam1, nInt, ncurr);
rateDiag = cellfun(calcX, xsetCol);


% % Loop across parameter space and calculate rate diagonal
% %delete(gcp('nocreate')); % remove active session
% %parpool(3); % use 3 cores
% parfor k = 1:gridSz(2)
%     % Parameter value from space
%     %xk = xsetMx(1:end, k);
%     xk = xsetCol{k};
%     
%     % Get integrand vector for integration
%     %intSet = inttT(xk, tstart, tspan);
%     intSet = intX(xk);
%     
%     % Trapz version of integral for PtT and rate diagonal
%     %PtT = (1 + trapz(tspan, intSet, 2))^(-1);
%     trapzIn = diff(tspan)*(intSet(1:nInt-1) + intSet(2:nInt))'/2;
%     PtT = (1 + trapzIn)^(-1);
%     rateDiag(k) = ncurr*lam1(xk)*PtT;
% end

% Remove any negative values of rateDiag (accounts for entries where mut >
% lamt in its parameter space)
rateDiag(rateDiag < 0) = 0;