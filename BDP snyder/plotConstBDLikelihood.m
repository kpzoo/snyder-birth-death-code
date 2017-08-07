% Script to plot the likelihood functions of Stadler and Snyder
clearvars
clc
close all

% Assumptions and modifications
% - does 2 Snyder likelihoods idx = 2 and 5
% - does all 7 Stadler likelihoods in a plot
% - looking at constant rate BDPs only
% - generate a tree from Stadler based on conditionings
% - plot likelihood across xsetMx grid

%% Simulate the constant rate BDP

% Initialise BDP tree variables, n tips
n = 100;
mi = 150*ones(1, 2);
numRV = length(mi);
nData = n-1;
nLin = 1:n;

% Parameterisation of BDP and space limits
xmin = [0.01 0.01];
xmax = [0.99 100];
xdefs = {'\rho', '\sigma'};
% Get parameter space joint grid
[xset, m, ~, ~] = getxsetMx(numRV, xmin, xmax, mi, zeros(1, numRV));

% Set true values of process as mean or other function
x = cellfun(@mean, xset);
% Set parameters for true constant birth-death
rho0 = x(1);
sig0 = x(2);
lam0 = sig0/(1 - rho0);
mu0 = lam0 - sig0;

% Simulate the constant BD process according to Hartmann 2008
[T, tbranch] = genBirDeaConst(n, sig0, rho0);
% Convert branch times into coalescents (0 is present)
if max(tbranch) >= T
    error('The branching times are not proper');
end
tcoal = sort([0 tbranch]);

% Convert times forward from T
tspec = sort(T - tcoal);
tspec = tspec - tspec(1); % set zero as first speciation
t0 = max(tspec);

%% Every Stadler likelihood

% Shortened version of function as tree fixed - Stadler likelihoods
likpar = @(idx, lam, mu, sig)likStadlerLog(idx, tcoal(1:end-1), t0, lam, mu, sig, n);
% Span all the idx values
idx = 1:7;
numStad = length(idx);

% Param space corresponding to all values in L and all other L's
xvals = zeros(mi);

% Loop across parameter space
L = zeros([length(idx), mi]);
for k = 1:numStad
    for i = 1:mi(1)
        rho = xset{1}(i);
        for j = 1:mi(2)
            % Re-parametrise for likelihood inputs
            sig = xset{2}(j);
            lam = sig/(1 - rho);
            mu = lam - sig;
            % Log likelihood at desired parameter set
            L(k, i, j) = likpar(idx(k), lam, mu, sig);
        end
    end
end

%% Snyder likelihood in terms of speciation times also from 2 lineages

% Calculate across parameter space
Lsny = zeros([2, mi]);
idxSny = [2 5]; % the 2 Stadler likelihoods that should match Sny

for k = 1:2
    % Adjust inputs based on if starting on 2 lineages
    idSny = idxSny(k);
    switch(idSny)
        case 2
            % Starting from 1 lineage at 0
            ts = tspec(1:end);
            nLin2 = nLin(1:end);
        case 5
            % Starting from 2 lineages so tspec(2) is first element
            ts = tspec(2:end);
            nLin2 = nLin(2:end);
    end
    % Get log likelihoods
    for i = 1:mi(1)
        rho = xset{1}(i);
        for j = 1:mi(2)
            % Re-parametrise for likelihood inputs
            sig = xset{2}(j);
            % Log likelihood at desired parameter set
            Lsny(k, i, j) = likSnyLog(ts, sig, rho, nLin2);
        end
    end
end


%% Marginalise and process the likelihoods

% Get all likelihoods in total (using 3d arrays) and combine
numLik = size(Lsny, 1) + size(L, 1);
Ltot = [L; Lnee; Lsny];
snyID = numLik - size(Lsny, 1) + 1; % starting index for Snyder lik
neeID = numStad+1; % index of Nee lik
nameLik = cell(1, numLik);
for i = 1:numLik
    if i <= numStad
        nameLik{i} = ['Stadler_' num2str(i)];
    end
    if i == numStad+1
        nameLik{i} = 'Nee';
    end
    if i > numStad+1
        nameLik{i} = ['Snyder_' num2str(i - numStad - 1)];
    end
end

% Check for inf and NaN
if any(any(any(isnan(Ltot)))) || any(any(any(isinf(Ltot))))
    error('The likelihood values are inf or NaN');
end

% Marginalise  all the likelihoods for both parameters
Lmarg = cell(numLik, 2);
for i = 1:2
    for k = 1:numLik
        Ltemp = Ltot(k, :, :);
        % Transpose Ltemp so indices of summation match i
        Ltemp = squeeze(Ltemp)'; % squeeze forces 2d array
        Lmarg{k, i} = sum(Ltemp, i);
        if i == 2
            % Transpose to keep as row vectors
            Lmarg{k, i} = Lmarg{k, i}';
        end
    end
end

% Get Snyder corrections to Stadler ones at idx 2 and 5
correctSny = mi(1)*log10(factorial(n-1));
kid = 0;
dsnystad = cell(2, 2); % hold difference in likelihoods
for i = idxSny
    % Starting index of relevant Snyder liks
    startSny = snyID + kid;
    for j = 1:2
        dsnystad{kid+1, j} = (Lmarg{startSny, j} - Lmarg{i, j})/correctSny;
    end
    % Count to move along Snyder indices
    kid = kid + 1;
end

% Check corrections now between Stadler and Snyder
dsize = size(dsnystad);
countCorr = 0;
for i = 1:dsize(1)
    for j = 1:dsize(2)
        if max(abs(dsnystad{i, j} - 1)) < 10^-10
            % Update counter for which correction works
            countCorr = countCorr + 1;
        end
    end
end
if countCorr == sum(dsize)
    disp(['Snyder correction to log likelihood value works: ' num2str(countCorr)]);
else
    disp('Correction did not work');
end

%% Get ML values from indices in grid

% Variable for ML parameter values
x_ml = zeros(numLik, 2);
for k = 1:numLik
    % Convert to 2D array
    Ltemp = squeeze(Ltot(k, :, :));
    % Get max value
    mL = max(max(Ltemp));
    % Get argument corresponding to this
    [arg1, arg2] = find(Ltemp == mL);
    % Get parameter grid values corresponding to args
    x_ml(k, 1) = xset{1}(arg1);
    x_ml(k, 2) = xset{2}(arg2);
end

%% Plot the log likelihood across parameter space

% Plot all the Stadler likelihoods
figure;
idnames = cell(1, length(idx));
for i = 1:numRV
    subplot(ceil(numRV/2), 2, i);
    hold on
    for j = idx
        plot(xset{i}, Lmarg{j, i}, 'linewidth', 2);
        idnames{j} = num2str(idx(j));
    end 
    hold off
    xlabel(xdefs{i});
    ylabel('log_{10} L');
    title('Stadler marginal likelihood');
    grid;
    legend(idnames, 'location', 'best');
end

% Compare Snyder likelihoods to Stadler ones
figure;
idnames = {'Stadler root', 'corrected Snyder root', 'Stadler mrca', 'corrected Snyder mrca', 'true'};
for i = 1:numRV
    subplot(ceil(numRV/2), 2, i);
    hold on
    for j = [idxSny(1), snyID, idxSny(2), snyID+1]
        % Correct the Snyder marginals only
        if ismember(j, idxSny)
            plot(xset{i}, Lmarg{j, i}, ':', 'linewidth', 2);
        else
            plot(xset{i}, Lmarg{j, i} - correctSny, 's');
        end
    end 
    h = gca;
    plot([x(i) x(i)], h.YLim, 'k', 'linewidth', 2);
    hold off
    xlabel(xdefs{i});
    ylabel('log_{10} L');
    title('Marginal likelihoods');
    grid;
    legend(idnames, 'location', 'best');
end

% Plot approximate ML values for each method
figure;
plot([x(1) x(1)], [xmin(2) xmax(2)], 'k');
hold on
plot([xmin(1) xmax(1)], [x(2) x(2)], 'k');
for i = 1:numLik
    plot(x_ml(i, 1), x_ml(i, 2), 'o', 'linewidth', 2);
end
hold off
xlabel(xdefs{1});
ylabel(xdefs{2});
grid;
legend(['true \rho', 'true \sigma', nameLik], 'location', 'best');