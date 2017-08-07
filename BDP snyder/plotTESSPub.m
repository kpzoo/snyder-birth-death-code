% Function to do main plotting of single run tessSnyderBatch (publication)
function plotTESSPub(numRV, samples, xhat, xMCMC, x, datapath, iGood, nGood, sampSny, ssSny, ssMCMC, treatDiscrete,...
    normSny, normMCMC, tset, lamhat, muhat, sampleTrue, empiricalData, T, nTax, thisDir, xdefs, lamn, mun)

% Assumptions and modifications
% - replaces raw trajectories with quantiles and mean
% - matched bandwidths in combined posteriors from all samples
% - iGood is index to specify converged MCMC runs with high enough ess
% - plotBatchTESS had boxplots and assumed iGood was needed here
% - assuming samples and sampSny etc are cells to be unpacked

%% Simple smoothed distribution of maximum speciation times

% The smoother is optimised for normal densities, bw is bandwith of moving
% average window, fT is probability and Ti the value
[fT, Ti, bw] = ksdensity(T, 'support', 'positive');
figure;
plot(Ti, fT);
xlabel('max speciation time');
ylabel('smoothed pdf');
grid;
title(['Density of T for [nTax bandwidth] = [' num2str(nTax) ' ' num2str(bw) ']']);

%% Density of estimates across runs from MCMC and Snyder

% Smoothed density on the parameter estimates
figure;
for i = 1:numRV
    subplot(ceil(numRV/2), 2, i);
    % Main density of each parameter for both methods
    [fmcmc, xmcmc] = ksdensity(xMCMC(:, i));
    [fsny, xsny] = ksdensity(xhat(:, i));
    plot(xmcmc, fmcmc, '-', xsny, fsny, ':', 'linewidth', 2);
    hold on
    % Get current axes and plot true value if not empirical data
    h = gca;
    if ~empiricalData
        plot(x(i)*ones(1, 2), h.YLim, 'k', 'linewidth', 2);
    end
    grid;
    xlabel('parameter space');
    ylabel('pdf of estimates');
    hold off
    title(['Estimates across runs of ' xdefs{i} ' = ' num2str(x(i))]);
    legend('mcmc', 'snyder', 'true', 'location', 'best');
end
cd(datapath);
saveas(gcf, 'paramEst');
cd(thisDir);

% Smoothed density of 2 norm across parameters
figure;
[fmcmc, xmcmc] = ksdensity(normMCMC);
[fsny, xsny] = ksdensity(normSny);
plot(xmcmc, fmcmc, xsny, fsny, 'linewidth', 2);
grid;
xlabel('norm across parameters');
ylabel('pdf of 2 norm');
hold off
title(['Square metrics across ' num2str(nGood) ' runs']);
legend('mcmc', 'snyder', 'location', 'best');
cd(datapath);
saveas(gcf, 'normBatch');
cd(thisDir);


% Smoothed density on the relative SSE across parameters
figure;
% Smoothed density of %SSE across parameters
[fmcmc, xmcmc] = ksdensity(ssMCMC);
[fsny, xsny] = ksdensity(ssSny);
plot(xmcmc, fmcmc, xsny, fsny, 'linewidth', 2);
grid;
xlabel('%SSE across parameters');
ylabel('pdf of relative %SSE');
hold off
title(['Relative SSE across ' num2str(nGood) ' runs']);
legend('mcmc', 'snyder', 'location', 'best');


%% Birth and death rate comparisons

% Get rates in concatenated form for easy statistics
lamh = cell2mat(lamhat);
muh = cell2mat(muhat);
% Assume each cell has the same number of elements
lamh = reshape(lamh, [length(lamhat{2}) nGood]);
muh = reshape(muh, [length(muhat{2}) nGood]);
% Get statistics
lammean = mean(lamh, 2);
lamlb = quantile(lamh', 0.025);
lamub = quantile(lamh', 1-0.025);
mumean = mean(muh, 2);
mulb = quantile(muh', 0.025);
muub = quantile(muh', 1-0.025);

% Plot trajectory statistics of lam(t), mu(t) from batch runs
figure;
subplot(2, 1, 1);
if ~empiricalData
    plot(tset, lamn, 'b', 'linewidth', 2);
end
xlabel('time');
xlim([tset(1) tset(end)]);
h = subplot(2, 1, 1);
ylim([0 1.1*h.YLim(2)]);
ylabel('birth rate');
hold on
% Mean and 95% trajectory bounds
plot(tset, lammean, 'r--', 'linewidth', 2);
plot(tset, lamlb,'g:', tset, lamub, 'g:', 'linewidth', 2);
grid;
if ~empiricalData
    legend('true', 'mean', '2.5%', '97.5%', 'location', 'best');
else
    legend('mean', '2.5%', '97.5%', 'location', 'best');
end
hold off
subplot(2, 1, 2);
if ~empiricalData
    plot(tset, mun, 'b', 'linewidth', 2);
end
xlabel('time');
xlim([tset(1) tset(end)]);
h = subplot(2, 1, 2);
ylim([0 1.1*h.YLim(2)]);
ylabel('death rate');
hold on
% Mean and 95% trajectory bounds
plot(tset, mumean, 'r--', 'linewidth', 2);
plot(tset, mulb,'g:', tset, muub, 'g:', 'linewidth', 2);
grid;
if ~empiricalData
    legend('true', 'mean', '2.5%', '97.5%', 'location', 'best');
else
    legend('mean', '2.5%', '97.5%', 'location', 'best');
end
hold off
cd(datapath);
saveas(gcf, 'rateBatch');
cd(thisDir);


%% Combined smoothed posteriors on same parameter space

if sampleTrue
    % Autocorrelation functions of both sample sets
    % Look at autocorrelation function of samples of parameters
    figure;
    subplot(2, 1, 1);
    hold all
    for i = 1:numRV
        for j = iGood
            [acf,lags,~] = autocorr(samples{j}(:, i));
            plot(lags, acf, 'o:');
        end
    end
    hold off
    grid;
    xlabel('lag');
    ylabel('autocorrelation');
    legend(xdefs, 'location', 'best');
    title('Autcorrelation of MCMC samples');
    subplot(2, 1, 2);
    hold all
    for i = 1:numRV
        for j = iGood
            [acf,lags,~] = autocorr(sampSny{j}(:, i));
            plot(lags, acf, 'o:');
        end
    end
    hold off
    grid;
    xlabel('lag');
    ylabel('autocorrelation');
    legend(xdefs, 'location', 'best');
    title('Autcorrelation of Snyder samples');
    
    % Averaged density across all samples from all files at once
    combMCMC = samples{iGood(1)};
    combSny = sampSny{iGood(1)};
    % Stack all the samples
    for i = iGood(2:end)
        combMCMC = [combMCMC; samples{i}];
        combSny = [combSny; sampSny{i}];
    end
    % Obtain full length of sample set
    lenSny = size(combSny, 1);
    lenMCMC = size(combMCMC, 1);
    % Check on consistency
    if lenMCMC ~= lenSny
        error('Combined sample lengths not equal');
    end
    
    % Plot overall smoothed density
     figure;
    for i = 1:numRV
        subplot(ceil(numRV/2), 2, i);
        hold on
        % Kernel smoothing, 100 pts with normal density
        [pSny, xSny, bwSny] = ksdensity(combSny(:, i), 'kernel', 'normal');
        if treatDiscrete
            % Increase the bandwidth since discrete marginals
            newBw = 3*bwSny;
            [pSny, xSny] = ksdensity(combSny(:, i), 'width', newBw);
        else
            % Keep bdwidth same across sample plots even if ~treatDiscrete
            newBw = bwSny;
        end
        [pMCMC, xMCMC] = ksdensity(combMCMC(:, i), 'kernel', 'normal', 'width', newBw);
        plot(xMCMC, pMCMC, xSny, pSny, 'linewidth', 2);
        % True values
        h1 = gca;
        yl = h1.YLim;
        plot([x(i) x(i)], yl, 'k', 'linewidth', 2);
        hold off
        box on
        xlabel(xdefs{i});
        ylabel(['P(' xdefs{i} '| data)']);
        legend('Hohna', 'Snyder', 'true', 'location', 'best');
        title('Marginal smoothed posteriors');
        grid;
    end
    cd(datapath);
    saveas(gcf, 'combpdf');
    cd(thisDir);
    
else
    % MCMC autocorrelation
    % Look at autocorrelation function of samples of parameters
    figure;
    hold all
    for j = iGood
        [acf,lags,~] = autocorr(samples{j}(:, i));
        plot(lags, acf, 'o:');
    end
    hold off
    grid;
    xlabel('lag');
    ylabel('autocorrelation');
    legend(xdefs, 'location', 'best');
    title('Autcorrelation of MCMC samples');
    cd(datapath);
    saveas(gcf, 'autocorr');
    cd ..
    cd ..
end


