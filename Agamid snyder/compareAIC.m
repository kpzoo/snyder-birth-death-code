% Evaluate Rabosky and Snyder AICs and other model selection statistics, as
% well as plot comparisons of estimates
clearvars
clc
close all

% Assumptions and modifications
% - aic and lik data files from R in current directory
% - files from Matlab runs of each model also present

%% Data extraction

% Estimated ML Rabosky parameters from laser package
rabos{1} = dlmread('const.csv');
rabos{2} = dlmread('spvar.csv');
rabos{3} = dlmread('exvar.csv');
rabos{4} = dlmread('bothvar.csv');

% AIC and max log likelihood from laser package
aicRab = dlmread('aic.csv');
likRab = dlmread('lik.csv');

% Get all the Snyder run files, assumes ag_ identifier
files = dir(['ag_' '*.mat']);
nF = length(files);
if nF ~= 4
    error('Expect 4 files, 1 for each model');
end
% Variables to store posteriors and estimates
Xhat = cell(1, nF);
Qmarg = cell(1, nF);
Qjnt = cell(1, nF);
Xset = cell(1, nF);
XsetMx = cell(1, nF);
rateNames = cell(1, nF);
Tspec = cell(1, nF);
nParam = zeros(1, nF);
% Extract data from each file
for i = 1:nF
    % Load relevant data
    temp = load(files(i).name, 'xhat', 'qmarg', 'qnlast', 'xset', 'xsetMx', 'rateName', 'tspec');
    % Assign into cells
    Xhat{i} = temp.xhat;
    Qmarg{i} = temp.qmarg;
    Qjnt{i} = temp.qnlast;
    Xset{i} = temp.xset;
    XsetMx{i} = temp.xsetMx;
    rateNames{i} = temp.rateName;
    nParam(i) = length(Xset{i});
    Tspec{i} = temp.tspec;
end

%% Akaike based metrics

% Take median AIC, lik and estimates for Rabosky as differences small and
% due to different initial conditions
aic1 = median(aicRab);
lik1 = median(likRab);
est1 = cell(1, nF);
for i = 1:nF
    est1{i} = median(rabos{i});
end

% Calculate an AIC (sort of) for the Snyder results
aic2temp = zeros(1, nF);
lik2temp = zeros(1, nF);
for i = 1:nF
    qtemp = Qjnt{i};
    [maxq, idq] = max(qtemp);
    % MAP value from grid (proportional to likelihood)
    lik2temp(i) = maxq;
    % AIC with this MAP value
    aic2temp(i) = 2*(nParam(i) - log(maxq));
end

%% Paradis based metrics

% Compare with Paradis metric, first ensure tspec same for each file
tspec = Tspec{1};
for i = 2:nF
    if all(tspec ~= Tspec{i})
        error('Branching data not consistent');
    end
end

% Time grid over which integrals are performed 
nInt = 2000;
nt = 5000;
T = max(tspec);
t = (0:nt-1).*T/nt;

% Get branching data into form of an empirical cdf (matches cdfplot matlab)
uspec = unique(tspec);
cdfu = (1:length(uspec))/length(uspec);
% Interpolate points from unique speciation times with extrap and zoh
pinterp = interp1(uspec, cdfu, t, 'previous', 'extrap');

% Construct dataset for Paradis evaluation of Snyder and Rabosky
data = cell(1, nF);
for i = 1:nF
   data{i} = load(files(i).name, 'lamt', 'mut', 'rhotT');
end

% Series of SSE from various methods
[cdfsny, ss_sny] = getCDFSSE2(nF, t, nInt, T, Xhat, pinterp, data);
[cdfrab, ss_rab] = getCDFSSE2(nF, t, nInt, T, est1, pinterp, data);
r_sny = ss_sny/max(ss_sny);
r_rab = ss_rab/max(ss_rab);

%% Plotting and visualisation

% Plot marginals for Snyder against Rabosky estimates
for j = 1:nF
    % Specific data
    xhat = Xhat{j};
    xset = Xset{j};
    qmarg = Qmarg{j};
    xRab = est1{j};
    % All marginals for a given model
    figure;
    for i = 1:nParam(j)
        % Subplot each parameter with its posterior mean
        subplot(ceil(nParam(j)/2), 2, i);
        plot(xset{i}, qmarg{i}, 'b', [xhat(end, i) xhat(end, i)], [0 max(qmarg{i})], 'r');
        hold on
        % Add values from Rabosky
        plot([xRab(i) xRab(i)], [0 max(qmarg{i})], 'ks-')
        hold off
        xlabel(['x_' num2str(i)]);
        ylabel(['P(x_' num2str(i) ')']);
    end
end

% Plot the SSE for each method normalised to maximum
figure;
plot(r_sny, 'bo-.', 'MarkerSize', 10)
hold on
plot(r_rab, 'kx-.', 'MarkerSize', 10)
hold off
ylabel('sum of square error ratio');
legend('Snyder', 'Rabosky', 'location', 'best');
% Label with models
h = gca;
set(h, 'XTick', 1:5)
set(h, 'XTickLabel', rateNames);
title('Comparison of model fits across estimators');

% Plot theoretical cdfs with square errors against data
figure;
stairs(uspec, cdfu, 'k');
hold all
for i = 1:nF
    plot(t, cdfsny(i, :));
    legName{i} = [rateNames{i} ', ss = ' num2str(ss_sny(i))];
end
hold off
grid;
xlabel('time');
ylabel('cdfs');
title('Theretical CDF from Snyder estimates across models');
legend(legName, 'location', 'best');
xlim([0 T]);
%axes('Position', [.2 .5 .2 .2]) % for inset of other plot
