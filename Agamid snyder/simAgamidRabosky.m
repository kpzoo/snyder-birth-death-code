% Script to simulate and estimate time varying birth-death models from the
% Agamid dataset

% Assumptions and modifications
% - properly handles branching times from getBtimes in R
% - looks at all Rabosky functions for Agamid data
% - allowed for duplicates
% - data obtained by extracting from the R package laser
% - discretised time and used T, mu, a, b (= k in paper) from Nee 1994

clearvars
clc
close all

tic;

% Booleans of interest
crownStart = 1; % start at crown (MRCA)

% Define rate function type
rateID = 1;

% Displaty function to be evaluated and settings
rateSet = {'const', 'spvar', 'exvar', 'bothvar'};
rateName = rateSet{rateID}; % don't use rateSet(rateID) else rateName is a cell
datapath = strcat(rateName);
disp('--------------------------------------------------------------------');
disp(['Simulation of function: ' rateName]);
disp(['Starting at crown = ' num2str(crownStart)]);
disp('--------------------------------------------------------------------');


%% Process branch times and other data from R laser package

% Read CSV branching times file
agfile = 'branch.csv';
tspec = dlmread(agfile);
% Ensure tspec is a row vector
if ~isrow(tspec)
    tspec = tspec';
end

% Data from R is referenced to tips vs root
tspec = max(tspec) - tspec;
tspec = sort(tspec);

if ~crownStart
    % Assume starting with 1 lineage at time 0
    n = length(tspec);
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
    n = length(tspec);
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

% Get Rabosky simulation values
rabos{1} = dlmread('const.csv');
rabos{2} = dlmread('spvar.csv');
rabos{3} = dlmread('exvar.csv');
rabos{4} = dlmread('bothvar.csv');

% Take mean as estimate
xRab = median(rabos{rateID});

%% Initial parameters and definitions

% Birth and death rate types
xEst.type = rateName;

% Set parameters based on function choices, now manually set xmin and xmax
% and include the estimates from Paradis and Rabosky
switch(rateID)
    case 1
        % Constant rate birth-death (null model Rabosky 2008)
        
        % Set parameter names and estimates
        xdefs = {'\lambda', '\mu'};
        xEst.x = [6.80747 0];
        
        % Points for discretised inference space and min/max values
        numRV = length(xdefs);
        mi = 30*ones(size(xdefs));
        xmin = [0 0];
        xmax = [10 1];
        
        % Define rate functions and integrals
        lamt = @(x, tx) x(1)*ones(size(tx));
        mut = @(x, tx) x(2)*ones(size(tx));
        
        % Integral of netRate is rhotT
        netRate = @(x, tx) mut(x, tx) - lamt(x, tx);
        rhotT = @(x, t1, t2) -(x(1) - x(2))*(t2 - t1);
        
    case 2
        % Best fitting function from Rabosky 2008 - SPVAR: exponentially
        % decreasing birth rate and constant death rate (3 parameters)
     
        % Set parameter names and estimates from Rabosky 2008
        xdefs = {'\lambda_0', 'k', '\mu_0'};
        
        % Points for discretised inference space and minimum search values
        numRV = length(xdefs);
        mi = 30*ones(size(xdefs));
        xmin = [1 1 0.0001];
        xmax = [100 25 0.01];
        
        % Define rate functions and integrals, lamt and mut return vectors
        lamt = @(x, tx) x(1)*exp(-x(2)*tx);
        mut = @(x, tx) x(3)*ones(size(tx));
        
        % Integral of netRate is rhotT
        netRate = @(x, tx) mut(x, tx) - lamt(x, tx);
        rhotT = @(x, t1, t2) x(3)*(t2 - t1) + (x(1)/x(2))*(exp(-x(2)*t2) - exp(-x(2)*t1));
        
    case 3
        % EXVAR is constant birth rate and exponentially increasing death
        % (3 parameters)

        % Set parameter names and estimates from Rabosky 2008
        xdefs = {'\lambda_0', 'z', '\mu_0'};
        
        % Points for discretised inference space and min/max values
        numRV = length(xdefs);
        mi = 30*ones(size(xdefs));
        xmin = [0.01 0.01 0.0001];
        xmax = [10 10 0.1];
        
        % Define rate functions and integrals, lamt and mut return vectors
        lamt = @(x, tx) x(1)*ones(size(tx));
        mut = @(x, tx) x(3)*(1 - exp(-x(2)*tx));
        
        % Integral of netRate is rhotT
        netRate = @(x, tx) mut(x, tx) - lamt(x, tx);
        rhotT = @(x, t1, t2) (x(3)- x(1))*(t2 - t1) + (x(3)/x(2))*(exp(-x(2)*t2) - exp(-x(2)*t1));
        
    case 4
        % BOTHVAR is exponentially decreasing birth rate and exponentially 
        % increasing death rate (4 parameters)

        % Set parameter names and estimates from Rabosky 2008
        xdefs = {'\lambda_0', 'k', 'z', '\mu_0'};
        
        % Points for discretised inference space and min/max values
        numRV = length(xdefs);
        mi = 20*ones(size(xdefs));
        xmin = [1 1 0.01 0.0001];
        xmax = [100 25 1 0.01];
        
        % Define rate functions and integrals, lamt and mut return vectors
        lamt = @(x, tx) x(1)*exp(-x(2)*tx);
        mut = @(x, tx) x(4)*(1 - exp(-x(3)*tx));
        
        % Integral of netRate is rhotT
        netRate = @(x, tx) mut(x, tx) - lamt(x, tx);
        rhotT = @(x, t1, t2) x(4)*(t2 - t1) + (x(4)/x(3))*(exp(-x(3)*t2) - exp(-x(3)*t1))...
            + (x(1)/x(2))*(exp(-x(2)*t2) - exp(-x(2)*t1));
        
        
    otherwise
        disp('The specified rate id is not supported');
end

% Rate function r(t, s) = e^(int mut-lamt) and integrand of PtT
rst = @(x, t1, t2) exp(rhotT(x, t1, t2));
inttT = @(x, t1, t2) mut(x, t2).*rst(x, t1, t2);


%% Snyder filtering of the reconstructed process

% Get space for rate function that combines all the parameter spaces
[xset, m, xsetMx, IDMx] = getxsetMx(numRV, xmin, xmax, mi, zeros(1, numRV));

% Set uniform joint prior - assumes uniformly spaced xset values
q0 = ones(1, m)/m;

% Main Snyder filtering code
[~, tn, qev] = snyderFilterPar2(m, tspec, nData, xsetMx, rhotT, lamt, mut, nLin, q0);

% Last posterior and marginals from event times
qnlast = qev(end, :); % joint sums to 1 
qmarg = marginalise2(numRV, IDMx, qnlast, mi, xset); % integrates to 1

% Conditional estimates at last event time
xhat = zeros(1, numRV);
for i = 1:numRV
    % Integral definition of mean
    xhat(i) = trapz(xset{i}, qmarg{i}.*xset{i});
end

% Get estimated rates based on function type
tset = linspace(tn(1), tn(end), 1000);
[lamhat, muhat] = plotLamMuhatRabosky(xsetMx, xhat, tset, qnlast, rateID, mut, lamt);

% Simulation time and data store
tsimSny = toc/60;
disp(['Simulation time = ' num2str(tsimSny) ' mins']);
clear qev
save(['ag_' num2str(rateID) '_' num2str(mi)]);

%% Analyse and plot results

% Plot birth and death rates
figure;
[hAx, hLine1, hLine2] = plotyy(tset, lamhat, tset, muhat);
xlabel('time');
hAx(1).XLim = [tn(1) tn(end)];
hAx(2).XLim = [tn(1) tn(end)];
ylabel(hAx(1),'birth rate') % left y-axis
ylabel(hAx(2),'death rate') % right y-axis
% Birth rates
hLine1(1).LineStyle = '-';
hLine1(1).Color = 'b';
% Death rates
hLine2(1).LineStyle = '-';
hLine2(1).Color = 'r';
legend('\lambda est', '\mu est', 'location', 'best');
title('Estimated rates');
grid;

% Net diversification rate
figure;
plot(tset, lamhat - muhat);
xlabel('time');
ylabel('net rate');


% Plot the marginal posteriors with inclusion of comparative values
figure;
for i = 1:numRV
    % Subplot each parameter with its posterior mean
    subplot(ceil((numRV)/2), 2, i);
    plot(xset{i}, qmarg{i}, 'b', [xhat(end, i) xhat(end, i)], [0 max(qmarg{i})], 'r');
    hold on
    % Add values from Rabosky
    plot([xRab(i) xRab(i)], [0 max(qmarg{i})], 'ks-')
    hold off
    xlabel(['x_' num2str(i)]);
    ylabel(['P(x_' num2str(i) ')']);
    if exist('xEst', 'var')
        legend('marginal', 'cond mean', xEst.type, 'location', 'best');
    else
        legend('marginal', 'cond mean', 'location', 'best');
    end
end
saveas(gcf, ['pdfag_' rateName]);