% Script to Snyder filter time-varying birth-death models and compare the
% results to TESS R package MCMC samples

% Assumptions and modifications
% - uses norms and separates bias and var of estimator, plotBatch2
% - assumes fixed BDP with known parameters that are also fixed
% - naming convention on files: samples, tesstree, init and branch
% - specify folder name to determine where to import and export data
% - functions taken from those set in TESS package by Hohna 2015
% - allow branching times to be imported as well as priors and x
% - uses algorithms from Hohna 2013
% - assumes complete and isochronous sampling

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
datapath = strcat(rateName, '/batch');
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
    xest.true = x;
else
    % If data empirical, no true values so put dummy
    x = -ones(1, numRV);
    disp('Empirical dataset being processed');
end
xmin = xdata(2, :);
xmax = xdata(3, :);
% Tree initialisation data
nBatch = xdata(4, 2);
nTax = xdata(4, 3);
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
        % Speciation-decay model from Hohna 2014
        
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
        % Logistic model from Paradis 2010
        
        % Set parameter names and estimates from Paradis 2010
        xdefs = {'\beta_\lambda', '\alpha_\lambda', '\beta_\mu', '\alpha_\mu'};
        % Points for discretised inference space and min/max values
        numRV = length(xdefs);
        mi = 15*ones(size(xdefs));
        
        % Define rate functions and integrals, lamt and mut return vectors
        lamt = @(x, tx) 1./(1 + exp(-x(1)*tx + x(2)));
        mut = @(x, tx) 1./(1 + exp(-x(3)*tx + x(4)));
        
        % Integral of netRate is rhotT
        rhotT = @(x, t1, t2) (1/x(3))*log((1 + exp(x(3)*t2 - x(4)))./(1 + exp(x(3)*t1 - x(4)))) -...
            (1/x(1))*log((1 + exp(x(1)*t2 - x(2)))./(1 + exp(x(1)*t1 - x(2))));
        
    otherwise
        disp('The specified rate id is not supported');
end

% Rate function r(t, s) = e^(int mut-lamt) and integrand of PtT
netRate = @(x, tx) mut(x, tx) - lamt(x, tx);
rst = @(x, t1, t2) exp(rhotT(x, t1, t2));
inttT = @(x, t1, t2) mut(x, t2).*rst(x, t1, t2);


%% Read tree data and branching times

% Branching times, tspec for single tree with n tips and and max T
cd(datapath);
tspec = dlmread('branchBatch.csv');
cd(thisDir)
% Ensure times are sorted so each row has increasing times
if tspec ~= sort(tspec, 2);
    disp('Imported speciation times were not sorted');
    tspec = sort(tspec, 2);
end
% Get all the maximum times
T = max(tspec, [], 2)';
% Data from TESS usually has time referenced to tips vs root
if reverseData
    for i = 1:nBatch
        tspec(i, :) = T(i) - tspec(i, :);
        tspec(i, :) = sort(tspec(i, :));
    end
end

if ~crownStart
    % Assume starting with 1 lineage at time 0
    [nBatch, n] = size(tspec);
    if ~all(tspec(:, 1) == 0)
        % Should start with a 0 at which there is 1 lineage
        error('No root time of 0 in tspec');
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

% Read ess values to find any runs we should discard
cd(datapath);
ess = csvread('ess.csv');
cd(thisDir);

% Check for NaN values in case and replace with 0
if any(any(isnan(ess)))
    ess(isnan(ess)) = 0;
end

% Any row possessing at least 1 ess < 500 mark as unsuitable
[rowFail, ~] = find(ess < 500);
rowFail = unique(rowFail);
nFail = length(rowFail);
disp(['No. non-converged MCMC runs = ' num2str(nFail) ' out of ' num2str(nBatch)]);
% Particularly check for MCMC that did not even really start
efail = ess(rowFail, :);
[row0, ~] = find(efail == 0);
disp(['No. failed MCMC runs = ' num2str(length(unique(row0))) ' out of ' num2str(nBatch)]);

% Only Snyder filter the supposedly converged cases
nGood = nBatch - nFail;
iGood = setdiff(1:nBatch, rowFail);
if nGood < 1
    error('All MCMC simulations failed');
end

% Get space for rate function that combines all the parameter spaces
[xset, m, xsetMx, IDMx] = getxsetMx(numRV, xmin, xmax, mi, zeros(1, numRV));

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

% Variables for results
xhat = zeros(nBatch, numRV);
perc = zeros(nBatch, numRV);
qnlast = cell(1, nBatch);
qmarg = cell(1, nBatch);
qmargD = cell(1, nBatch);
xhatEv = cell(1, nBatch);
tn = cell(1, nBatch);

% Main Snyder filtering code, use qev for final posterior
for i = iGood
    [~, tn{i}, qev] = snyderFilterPar(m, tspec(i, :), nData, xsetMx, rhotT, lamt, mut, nLin, q0);
    % Last posterior and marginals and times
    qnlast{i} = qev(end, :);
    [qmargD{i}, ~] = marginalise(numRV, IDMx, qnlast{i}, mi);
    qmarg{i} = marginalise2(numRV, IDMx, qnlast{i}, mi, xset);
    % Conditional estimates, xhatEv gives estimates at event times
    xhatEv{i} = qev*xsetMx';
    % Conditional estimates at last event time
    for j = 1:numRV
        % Integral definition of mean
        xhat(i, j) = trapz(xset{j}, qmarg{i}{j}.*xset{j});
    end
    perc(i, 1:numRV) = 100*((1 - xhat(i, 1:numRV)./x).^2);
    disp('**********************************************************');
    disp(['Finished Snyder batch: ' num2str(i) ' of ' num2str(nBatch)]);
    disp(['Square percent errors: ' num2str(perc(i, 1:numRV))]);
    disp('**********************************************************');
    clear qn qev
end

% Save Snyder results in case of later errors
cd(datapath);
save([rateName '_' num2str(rateID) '_' num2str(nBatch) '_sny']);
cd(thisDir);

%% Read in MCMC TESS results and compare to Snyder

% Error arrays across runs
ePrior = x - xest.priors;
eSny = zeros(nBatch, numRV);
eMCMC = zeros(nBatch, numRV);
xMCMC = zeros(nBatch, numRV);

% Variables for sampling for Snyder marginals
samples = cell(1, nBatch);
sampSny = cell(1, nBatch);
sampleFail = zeros(1, nBatch);

% Variables for rate function estimates
lamhat = cell(1, nBatch);
muhat = cell(1, nBatch);

% Fixed time frame for comparing rate estimates based on mean T
Tmean = mean(T);
tset = linspace(0, Tmean, 1000);

% True rates across tset
lamn = lamt(x, tset);
mun = mut(x, tset);

% Loop through Snyder estimates and relevant MCMC samples
%j = 0; % counter in case xhat only on iGood
for i = iGood
    % Read relevant MCMC samples csv file
    cd(datapath);
    samples{i} = csvread(strcat(['samples_' num2str(i) '.csv']));
    cd(thisDir);
    
    % Get MCMC estimates
    xMCMC(i, 1:numRV) = mean(samples{i});
    
    % Get errors and norms
    eSny(i, 1:numRV) = x - xhat(i, :);
    eMCMC(i, 1:numRV) = x - xMCMC(i, :);

    % Sample from each Snyder marginal set
    try
        % Sample from staircase CDF then backsample with first unique id
        if treatDiscrete
            [~, ~, sampSny{i}] = compareSnyMCMC3(samples{i}, qmargD{i}, xset, numRV, 0.01, interpMeth, treatDiscrete);
        else
            [~, ~, sampSny{i}] = compareSnyMCMC3(samples{i}, qmarg{i}, xset, numRV, 0.01, interpMeth, treatDiscrete);
        end
        % Write Snyder samples to csv
        cd(datapath);
        % The 0, 0 arguments mean write first row and column, without it
        % sampSny will miss its first row of data
        csvwrite(['sampSny' num2str(i) '.csv'], sampSny{i}, 0, 0);
        cd(thisDir);
    catch failSamp
        % Get error (the exception thrown)
        disp(['Sampling from Snyder failed at i = ' num2str(i)]);
        disp(failSamp);
        sampSny{i} = -1;
        sampleFail(i) = 1;
    end
    
    % Get true birth and death rates if available and estimated rates always
    [lamhat{i}, muhat{i}, ~, ~] = getLamMuhatTESS(xsetMx, tset, qnlast{i}, rateID,...
        mut, lamt, empiricalData, x, changeTime);
end

% Calculate sampled Snyder stats
xhatSamp = zeros(nBatch, numRV);
eSnySamp = zeros(nBatch, numRV);
if ~any(sampleFail(iGood))
    sampleTrue = 1;
    for i = iGood
        % Sampled estimate and error
        xhatSamp(i, 1:numRV) = mean(sampSny{i});
        eSnySamp(i, 1:numRV) = x - xhatSamp(i, 1:numRV);
    end
    % Remove all 0s for failed MCMCs
    xhatSamp = xhatSamp(iGood, :);
    eSnySamp = eSnySamp(iGood, :);
else
    sampleTrue = 0;
end

% Check the sample means
figure;
plot(xhat, 'b');
hold on
plot(xhatSamp, 'rs');
hold off
xlabel('runs');
ylabel('mean estimates');
legend('posterior based', 'sample based', 'location', 'best');


% Remove all the 0s from the failed MCMCs
eSny = eSny(iGood, :);
eMCMC = eMCMC(iGood, :);
xhat = xhat(iGood, :);
xMCMC = xMCMC(iGood, :);

% Assign mean across runs to overall estimates
xest.sny = mean(xhat);
if sampleTrue
    xest.snySamp = mean(xhatSamp);
end
xest.mcmc = mean(xMCMC);
xest.name = xdefs;

% Bias and variance of estimators
biasEst.mcmc = mean(-eMCMC);
biasEst.sny = mean(-eSny);
varEst.mcmc = var(-eMCMC);
varEst.sny = var(-eSny);
disp(biasEst);
disp(varEst);
% The MSE is what the filter must minimise
mseEst.mcmc = mean(eMCMC.^2);
mseEst.sny = mean(eSny.^2);
disp(mseEst);

% Variables to store metrics per run
normSny = -ones(nGood, 1);
normSnySamp = -ones(nGood, 1);
normMCMC = -ones(nGood, 1);
ssSny = zeros(nGood, 1);
ssSnySamp = zeros(nGood, 1);
ssMCMC = zeros(nGood, 1);

% Get statistics and metrics across good runs
for i = 1:nGood
    % 2-norm values
    normSny(i) = norm(eSny(i, :));
    normMCMC(i) = norm(eMCMC(i, :));
    % Relative percent SSE
    ssSny(i) = sum((100*eSny(i, :)./x).^2);
    ssMCMC(i) = sum((100*eMCMC(i, :)./x).^2);
    % Same stats but for sampled Snyder
    if sampleTrue
        normSnySamp(i) = norm(eSnySamp(i, :));
        ssSnySamp(i) = sum((100*eSnySamp(i, :)./x).^2);
    end
end

% Calculate metrics and error indices
normSet.name = {'priors', 'mcmc', 'snyder', 'snyder sampled'};
normSet.val = [norm(ePrior), mean(normMCMC), mean(normSny), mean(normSnySamp)];
disp(normSet);
ssSet.name = {'priors', 'mcmc', 'snyder', 'snyder sampled'};
ssSet.val = [sum((100*ePrior./x).^2), mean(ssMCMC), mean(ssSny), mean(ssSnySamp)];
disp(ssSet);

% Main plotting function
plotTESSPub(numRV, samples, xhat, xMCMC, x, datapath, iGood, nGood, sampSny, ssSny, ssMCMC, treatDiscrete,...
    normSny, normMCMC, tset, lamhat, muhat, sampleTrue, empiricalData, T, nTax, thisDir, xdefs, lamn, mun);


%% Final storage of data and saving of dependencies

% Get all m files called in creating this simulated data
currName = dbstack; % gives a struct with current m file name
currName = currName.file;
[fList,pList] = matlab.codetools.requiredFilesAndProducts(currName);

% Store all data
clear samples sampSny
cd(datapath);
save([rateName '_' num2str(rateID) '_' num2str(nBatch)]);
cd .. 
cd ..