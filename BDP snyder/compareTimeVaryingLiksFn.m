% Calculate deterministically time-varying BDP likelihoods and compare
function [Lsny, Lhoh, LmargS, LmargH, outPm] = compareTimeVaryingLiksFn(rateID, mstar)

% Assumptions and modifications
% - compares Snyder and Hohna based likelihoods
% - relevant tree data available from TESS

%% Initialisation of parameters and rate functions

% Choose function want to evaluate

% Birth and death rate types also match subfolders, rateName used for saves
rateSet = {'hohna', 'rateshift', 'logistic', 'constant'};
rateName = rateSet{rateID}; % don't use rateSet(rateID) else rateName is a cell
datapath = strcat(rateName);
disp('--------------------------------------------------------------------');
disp(['Calculating likelihood for: ' rateName]);
disp('--------------------------------------------------------------------');

% Get path of code
thisDir = cd;

% Read the true value and range limits from TESS csv
cd(datapath);
xdata = dlmread('init.csv');
cd(thisDir);
% The initialisation vector xdata must have at least 3 rows
x = xdata(1, :); xmin = xdata(2, :); xmax = xdata(3, :);

% Check for boolean about MRCA located in 5th row
if size(xdata, 1) >= 5
    useMRCA = xdata(5, 1);
    if ~useMRCA
        error('Tree did not start from MRCA');
    end
end

% Set parameters based on function choices
switch(rateID)
    case 1
        % Exponentially decreasing speciation rate from Hohna 2014
        
        % Set parameter names from Hohna 2014
        xdefs = {'\delta', '\lambda', '\alpha'};
        % Points for discretised inference space and min/max values
        numRV = length(xdefs);
        mi = mstar*ones(size(xdefs));
        
        % Define rate functions and integrals, lamt and mut return vectors
        lamt = @(x, tx) x(1) + x(2)*exp(-x(3)*tx);
        mut = @(x, tx) x(1)*ones(size(tx));
        
        % Integral of netRate is rhotT
        rhotT = @(x, t1, t2) (x(2)/x(3))*(exp(-x(3)*t2) - exp(-x(3)*t1));
        
    case 2
        % Simultaneous rateshift on both birth and death rates but
        % parametrised in terms of sig = lam - mu and mu, tshift known
        
        % Set parameter names from Hohna tutorial
        xdefs = {'\sigma_1', '\mu_1', '\sigma_2', '\mu_2'};
        % Points for discretised inference space and min/max values
        numRV = length(xdefs);
        mi = mstar*ones(size(xdefs));
        
        % Tree initialisation data and change times
        changeTime = xdata(4, 1);
        
        % Define rate functions and integrals, lamt and mut return vectors
        lamt1 = @(x, tx, tch) (x(1) + x(2)).*(tx < tch) + (x(3)+ x(4)).*(tx >= tch);
        mut1 = @(x, tx, tch) x(2).*(tx < tch) + x(4).*(tx >= tch);
        % Add the changeTime variable and create funcs in x and tx only
        lamt = @(x, tx)lamt1(x, tx, changeTime);
        mut = @(x, tx)mut1(x, tx, changeTime);
        
        % Integral of netRate is rhotT
        rhotT = @(x, t1, t2)rateShiftRho2(x, t1, t2, changeTime);
        
    case 3
        % Logistic birth and death rates from Paradis 2010
        
        % Set parameter names
        xdefs = {'\beta_\lambda', '\alpha_\lambda', '\beta_\mu', '\alpha_\mu'};
        % Points for discretised inference space and min/max values
        numRV = length(xdefs);
        mi = mstar*ones(size(xdefs));
        
        % Define rate functions and integrals, lamt and mut return vectors
        lamt = @(x, tx) 1./(1 + exp(-x(1)*tx + x(2)));
        mut = @(x, tx) 1./(1 + exp(-x(3)*tx + x(4)));
        
        % Integral of netRate is rhotT
        rhotT = @(x, t1, t2) (1/x(3))*log((1 + exp(x(3)*t2 - x(4)))./(1 + exp(x(3)*t1 - x(4)))) -...
            (1/x(1))*log((1 + exp(x(1)*t2 - x(2)))./(1 + exp(x(1)*t1 - x(2))));
    otherwise
        disp('The specified rate id is not supported');
end


%% Read tree data and branching times

% Branching times, tspec for single tree with n tips and and max T
cd(datapath);
tspec = dlmread('branch.csv');
cd(thisDir)
% Ensure times are sorted so each row has increasing times
if tspec ~= sort(tspec, 2);
    disp('Imported speciation times were not sorted');
    tspec = sort(tspec, 2);
end
% Get all the maximum times
T = max(tspec);
% Data from TESS usually has time referenced to tips vs root
tspec = T - tspec;
tspec = sort(tspec);
if tspec(1) ~= 0
    % Add necessary 0 for 2 lineages (see matlab tree)
    tspec = [0, tspec];
end
% Get number of lineages with account for crown
n = length(tspec);
n = n + 1; % because the length doesn't count starting at 2

% Lineages input to Snyder filter
nLin = 2:n;


%% Calculate Stadler and Snyder likelihoods

% Get space for evaluating likelihood function
[xset, m, xsetMx, IDMx] = getxsetMx(numRV, xmin, xmax, mi, zeros(1, numRV));
% Decompose xsetMx into column and row cell arrays
gridSz = size(xsetMx);
xsetCol  = mat2cell(xsetMx, gridSz(1), ones(1, gridSz(2)));

% Rate function r(t, s) = e^(int mut-lamt) and integrand of PtT
netRate = @(x, tx) mut(x, tx) - lamt(x, tx);
rst = @(x, t1, t2) exp(rhotT(x, t1, t2));
inttT = @(x, t1, t2) mut(x, t2).*rst(x, t1, t2);

% Loop storage of log likelihoods
Lsny = zeros(1, m);
Lhoh = zeros(1, m);

% Snyder correction
correctSny = mi(1)*log(factorial(n-1));

% Loop across parameter space and evaluate likelihood functions
parfor i = 1:m
    % Snyder log likelihood
    Lsny(i) = likSnyLogTimeVary(tspec, xsetCol{i}, nLin, lamt, inttT);
    % Hohna log likelihood
    Lhoh(i) = likHohnaLogTimeVary(tspec, xsetCol{i}, nLin, lamt, inttT, rst);
    % Progress
    disp(['Completed: ' num2str(i) ' of ' num2str(m)]);
end

%% Plotting and post analysis

% Marginalise the log likelihoods
[LmargS, ~] = marginalise(numRV, IDMx, Lsny, mi);
[LmargH, ~] = marginalise(numRV, IDMx, Lhoh, mi);

% Ouput structure of parameters
outPm.x = x;
outPm.xset = xset;
outPm.nLin = nLin;
outPm.tspec = tspec;
outPm.correctSny = correctSny;

save([rateName '.mat']);

% % Legend names
%leg = {'Snyder', 'Hohna'};

% % Plot the likelihoods over the parameter set (unordered)
% figure;
% plot(1:m, Lsny, 1:m, Lhoh);
% xlabel('points');
% ylabel('log likelihood');
% legend(leg, 'location', 'best');
% title(['Likelihoods for ' rateName ' model']);
% 
% % Plot the marginals
% figure;
% for j = 1:numRV
%     subplot(ceil(numRV/2), 2, j);
%     hold on
%     plot(xset{j}, LmargS{j}, xset{j}, LmargH{j});
%     h = gca;
%     plot([x(j) x(j)], h.YLim, 'k:'); 
%     hold off
%     xlabel(['x_' num2str(j)]);
%     legend(leg, 'true', 'location', 'best');
%     ylabel(['log L(data | x_' num2str(j) ')']);
% end
