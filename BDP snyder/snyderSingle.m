% Script to Snyder filter single functions (no batch) to compare to TESS

% Assumptions and modifications
% - uses parfor to improve performance
% - can get any prior by assigning uniform probs to different space
% increments - like inverse probability transform
% - modified to account for continuous distributions more sensibly
% - crown conditioning matches Nee likelihood and Tanja's likelihood (5)
% - allows process to start at crown age so N(tcrown) = 2
% - naming convention on files: samples, tesstree, init and 
% - specify folder name to determine where to import and export data
% - functions taken from those set in TESS
% - allow branching times to be imported as well as priors and x
% - fixes T to max(tspec) which makes more sense esp for Paradis
% - uses algorithms from Hohna 2013
% - complete isochronous sampling assumed

clearvars
clc
close all

% Start timing
tic;

% Booleans of importance delimiting the functions we want
rateID = 2;
reverseData = 1; % convert branch times to T - times (always on)
wantTree = 0;
crownStart = 1; % condition on starting with 2 lineages and survival
empiricalData = 0; % if data real or simulated
interpMeth = 'previous'; % interpolation method for sampling Snyder
nosave = 0; % option to not save - for testing
treatDiscrete = 1; % decide if to sample from discrete probs or cont pdf

%% Initial parameters and definitions

% Birth and death rate types also match subfolders, rateName used for saves
rateSet = {'hohna', 'logistic'};
rateName = rateSet{rateID}; % don't use rateSet(rateID) else rateName is a cell
datapath = strcat(rateName);
disp('--------------------------------------------------------------------');
disp(['Simulation of function: ' rateName]);
disp(['[Logspace grid, starting at crown] = ' num2str(0) ', ' num2str(crownStart)]);
disp('--------------------------------------------------------------------');

% Get path of code
thisDir = cd;

% Read the true value and range limits from TESS csv
cd(datapath);
xdata = dlmread('init.csv');
cd(thisDir);
% The initialisation vector xdata must have at least 3 rows
if ~empiricalData
    x = xdata(1, :);
else
    % If data empirical, no true values so put dummy
    x = -ones(1, numRV);
    disp('Empirical dataset being processed');
end
xmin = xdata(2, :);
xmax = xdata(3, :);

% Check for boolean about MRCA and nIters from MCMC, located in 5th row,
% 4th row specific to rate functions
if size(xdata, 1) >= 5
    useMRCA = xdata(5, 1);
    nIters = xdata(5, 2);
    disp(['MCMC simulation done for useMRCA = ' num2str(useMRCA)]);
end

% Set parameters based on function choices
switch(rateID)
        
    case 1
        % Exponentially decreasing speciation rate from Hohna 2014
        
        % Set parameter names and estimates from Hohna 2014
        xdefs = {'\delta', '\lambda', '\alpha'};
        % Points for discretised inference space and min/max values
        numRV = length(xdefs);
        mi = 20*ones(size(xdefs));
        
        % Extra initial conditions to read if available
        if size(xdata, 1) > 3
            TsetR = xdata(4, 1);
            nBatchR = xdata(4, 2);
            nTax = xdata(4, 3);
        end
        
        % Define rate functions and integrals, lamt and mut return vectors
        lamt = @(x, tx) x(1) + x(2)*exp(-x(3)*tx);
        mut = @(x, tx) x(1)*ones(size(tx));
        
        % Integral of netRate is rhotT
        rhotT = @(x, t1, t2) (x(2)/x(3))*(exp(-x(3)*t2) - exp(-x(3)*t1));
        
    case 2
        % Logistic birth and death rates from Paradis 2010
        
        % Set parameter names and estimates from Paradis 2010
        xdefs = {'\beta_\lambda', '\alpha_\lambda', '\beta_\mu', '\alpha_\mu'};
        % Points for discretised inference space and min/max values
        numRV = length(xdefs);
        mi = 15*ones(size(xdefs));
        
        % Define rate functions and integrals, lamt and mut return vectors
        lamt = @(x, tx) 1./(1 + exp(-x(1)*tx + x(2)));
        mut = @(x, tx) 1./(1 + exp(-x(3)*tx + x(4)));
        %lamt = @(x, tx) 1./(1 + x(2)*(x(1).^tx));
        %mut = @(x, tx) 1./(1 + x(4)*(x(3).^tx));
        
        % Integral of netRate is rhotT
        rhotT = @(x, t1, t2) (1/x(3))*log((1 + exp(x(3)*t2 - x(4)))./(1 + exp(x(3)*t1 - x(4)))) -...
            (1/x(1))*log((1 + exp(x(1)*t2 - x(2)))./(1 + exp(x(1)*t1 - x(2))));
        %rhotT = @(x, t1, t2) (-1/log(x(3)))*log((1 + exp(-log(x(3))*t2 - log(x(4))))./(1 + exp(-log(x(3))*t1 - log(x(4))))) -...
        %    (-1/log(x(1)))*log((1 + exp(-log(x(1))*t2 - log(x(2))))./(1 + exp(-log(x(1))*t1 - log(x(2)))));
        
    otherwise
        disp('The specified rate id is not supported');
end

% Rate function r(t, s) = e^(int mut-lamt) and integrand of PtT
netRate = @(x, tx) mut(x, tx) - lamt(x, tx);
rst = @(x, t1, t2) exp(rhotT(x, t1, t2));
inttT = @(x, t1, t2) mut(x, t2).*rst(x, t1, t2);


%% Read tree data and branching times

% If want to investigate the input tree
if wantTree
    % Read newick string of tree
    cd(datapath);
    fid = fopen('tesstree.txt');
    datastr = textscan(fid, '%c');
    cd(thisDir);
    % This will be a cell so get in char row vector format to use regexp
    datastr = datastr{1}';
    % Starting index of newick part
    st = regexp(datastr, '(');
    treeStr = datastr(st(1):length(datastr));

    % Get tree object evaluate and plot
    tree = phytreeread(treeStr);
    % Get the branch to leaf distances and the no. lineages
    [M, ID, D] = getmatrix(tree);
    n0 = get(tree, 'NumLeaves');
    % Plot the tree
    h = plot(tree, 'orient', 'top');
    ylabel('Time')
    set(h.terminalNodeLabels,'Rotation',65)
end

% Branching times, tspec for single tree with n tips and and max T
cd(datapath);
tspec = dlmread('branch.csv');
cd(thisDir);
% Ensure tspec is a row vector
if ~isrow(tspec)
    tspec = tspec';
end

% Data from TESS usually has time referenced to tips vs root
if reverseData
    tspec = max(tspec) - tspec;
    tspec = sort(tspec);
end

if ~crownStart
    % Assume starting with 1 lineage at time 0
    [nBatch, n] = size(tspec);
    if ~any(tspec == 0)
        error('Speciation times do no include 0');
    end
else
    % Assume start with 2 lineages at time tcrown
    tspec = sort(tspec);
    if tspec(1) ~= 0
        % Add necessary 0 for 2 lineages (see matlab tree)
        tspec = [0, tspec];
    end
    % Get number of lineages with account for crown
    [nBatch, n] = size(tspec);
    n = n + 1; % because the length doesn't count starting at 2
end

% No. batches, TMRCA (T) and no. lineages and events (nData)
if nBatch ~= 1
    error('Script only works for single tree analysis');
end
T = max(tspec, [], 2)';
% Lineages input to Snyder filter
if ~crownStart
    nLin = 1:n;
else
    % See tspec creation for why this works
    nLin = 2:n;
end
% No. of iterations, equivalent to counting events
nData = length(nLin)-1;

%% Snyder filtering of the reconstructed process

% Get space for rate function that combines all the parameter spaces
[xset, m, xsetMx, IDMx] = getxsetMx(numRV, xmin, xmax, mi, zeros(1, numRV));
% Write xsetMx to evaluate likelihood surface in TESS
csvwrite('xset.csv', xsetMx');

% Variables for results
xhat = zeros(nBatch, numRV);
perc = zeros(nBatch, numRV);

% Set uniform joint prior - assumes uniformly spaced xset values
q0 = ones(1, m)/m;
% Marginal priors that keep integral to 1
q0marg = marginalise2(numRV, IDMx, q0, mi, xset);

% Calculate initial estimate of mean from priors
xest.priors = zeros(1, numRV);
for i = 1:numRV
    % Integral definition of mean
    xest.priors(i) = trapz(xset{i}, q0marg{i}.*xset{i});
end

% Main Snyder filtering code
[~, tn, qev] = snyderFilterPar(m, tspec, nData, xsetMx, rhotT, lamt, mut, nLin, q0);
%[~, tn, qev] = snyderFilterCrown(m, tspec, nData, xsetMx, numRV, rhotT, lamt, mut, nLin);

% Last posterior and marginals from event times
qnlast = qev(end, :); % joint sums to 1 
qmarg = marginalise2(numRV, IDMx, qnlast, mi, xset); % integrates to 1
qmargD = marginalise(numRV, IDMx, qnlast, mi); % marginals sum to 1

% Conditional estimates at last event time
for i = 1:numRV
    % Integral definition of mean
    xhat(i) = trapz(xset{i}, qmarg{i}.*xset{i});
end
perc(1:numRV) = 100*((1 - xhat(1:numRV)./x).^2);

% Simulation time and clear iterative posteriors as too large
tsimSny = toc/60;
disp(['Simulation time = ' num2str(tsimSny) ' mins']);
% Largest sized variables, remove to save space
clear qn 


%% Read in MCMC TESS results and compare to Snyder

% Read csv with MCMC samples
cd(datapath);
samples = csvread('samples.csv');
cd(thisDir);

% Data structure comparing estimates
xest.name = xdefs;
xest.true = x;
xest.sny = xhat;
xest.mcmc = mean(samples);
disp(xest);

% Combined norm on parameters (attempt to avoid Stein's paradox)
eSny = x - xhat;
ePrior = x - xest.priors;
eMCMC = x - xest.mcmc;
normSet.name = {'priors', 'mcmc', 'snyder'};
normSet.val = [norm(ePrior), norm(eMCMC), norm(eSny)];
disp(normSet);

% Test the possiblility that the marginals are statistically the same
try
    % Get discrete marginals
    if treatDiscrete
        [h, p, sampSny] = compareSnyMCMC3(samples, qmargD, xset, numRV, 0.01, interpMeth, treatDiscrete);
    else
        [h, p, sampSny] = compareSnyMCMC3(samples, qmarg, xset, numRV, 0.01, interpMeth, treatDiscrete);
    end
    sampleFail = 0;
    % Write Snyder samples for use in R's ggplot
    if ~nosave
        cd(datapath);
        % The 0, 0 arguments mean write first row and column, without it
        % sampSny will miss its first row of data
        csvwrite('sampSny.csv', sampSny, 0, 0);
        % Check size of written file
        tempSamp = csvread('sampSny.csv');
        if ~all(size(tempSamp) == size(samples))
            error('The written Snyder samples do not have the same dimensions as the MCMC');
        end
        cd(thisDir);
    end
catch failSamp
    % Get error (the exception thrown)
    disp('Sampling from Snyder failed');
    disp(failSamp);
    h = 0;
    p = 0;
    sampSny = -1;
    sampleFail = 1;
end

% Get true birth and death rates if available and estimated rates always
tset = linspace(tn(1), tn(end), 1000);
if ~exist('changeTime', 'var')
    changeTime = [];
end
[lamhat, muhat, lamn, mun] = getLamMuhatTESS(xsetMx, tset, qnlast, rateID,...
    mut, lamt, empiricalData, x, changeTime);

% Do plots of posteriors, autocorr and other diagnostics
plotTESS2(numRV, samples, qmarg, xset, x, datapath, nosave, sampSny, treatDiscrete,...
    xdefs, tset, lamn, mun, lamhat, muhat, sampleFail, empiricalData, thisDir, qev, IDMx, mi);

% Just for check of convegence print discrete based estimate across events
xhatEv = qev*xsetMx';

%% Storage of data and saving of dependencies

% Get all m files called in creating this simulated data
currName = dbstack; % gives a struct with current m file name
currName = currName.file;
[fList,pList] = matlab.codetools.requiredFilesAndProducts(currName);

% Store all data
if ~nosave
    clear samples qev
    cd(datapath);
    save([rateName '_' num2str(rateID) '_' num2str(nBatch)]);
    cd(thisDir);
end