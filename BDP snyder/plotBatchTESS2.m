% Function to do main plotting of single run tessSnyderBatch
function plotBatchTESS2(numRV, samples, xhat, xMCMC, x, datapath, iGood, nGood, sampSny, ssSny, ssMCMC, treatDiscrete,...
    normSny, normMCMC, tset, lamhat, muhat, sampleTrue, empiricalData, T, nTax, thisDir, xdefs, lamn, mun)

% Assumptions and modifications
% - matched bandwidths in combined posteriors from all samples
% - iGood still relevant for lamhat and other cells
% - now iGood already dealt with at input stage and do densities
% - plotBatchTESS had boxplots and assumed iGood was needed here
% - assumes removed the failed MCMC trees from the metrics and estimates
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

%% Boxplot of estimates across runs from MCMC and Snyder

% Smoothed density on the parameter estimates
figure;
for i = 1:numRV
    subplot(ceil(numRV/2), 2, i);
    % Main density of each parameter for both methods
    [fmcmc, xmcmc] = ksdensity(xMCMC(:, i));
    [fsny, xsny] = ksdensity(xhat(:, i));
    plot(xmcmc, fmcmc, xsny, fsny, 'linewidth', 2);
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


% Smoothed density on the norm across parameters
figure;
% Smoothed density of 2 norm across parameters
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

% Plot multiple reconstructions of lam(t), mu(t) from batch runs
figure;
subplot(2, 1, 1);
if ~empiricalData
    plot(tset, lamn, 'r', 'linewidth', 2);
end
xlabel('time');
xlim([tset(1) tset(end)]);
ylabel('birth rate');
hold on
for j = iGood
    % Estimates for this run
    ha = plot(tset, lamhat{j}, 'g:');
    set(ha, 'Color', [0.5 0.5 0.5]);
end
grid;
if ~empiricalData
    legend('true trajectory', 'batch estimates', 'location', 'best');
else
    legend('batch estimates across runs', 'location', 'best');
end
hold off
subplot(2, 1, 2);
if ~empiricalData
    plot(tset, mun, 'r', 'linewidth', 2);
end
xlabel('time');
xlim([tset(1) tset(end)]);
ylabel('death rate');
hold on
for j = iGood
    % Estimates for this run
    hb = plot(tset, muhat{j}, 'g:');
    set(hb, 'Color', [0.5 0.5 0.5]);
end
grid;
if ~empiricalData
    legend('true trajectory', 'batch estimates', 'location', 'best');
else
    legend('batch estimates across runs', 'location', 'best');
end
hold off
cd(datapath);
saveas(gcf, 'rateBatch');
cd(thisDir);


%% Combined smoothed posteriors on same parameter space

if sampleTrue
    % Plot comparative smoothed densities
    figure;
    for i = 1:numRV
        subplot(ceil(numRV/2), 2, i);
        hold on
        for j = iGood
            % Kernel smoothing, 100 pts with normal density
            [pSny, xSny] = ksdensity(sampSny{j}(:, i), 'kernel', 'normal');
            [pMCMC, xMCMC] = ksdensity(samples{j}(:, i), 'kernel', 'normal');
            % ISSUE: 'support', [xmin(i) xmax(i)], get limits <----------
            plot(xMCMC, pMCMC, 'b-');
            plot(xSny, pSny, 'r-');
        end
        % True values
        h1 = gca;
        yl = h1.YLim;
        plot([x(i) x(i)], yl, 'k', 'linewidth', 2);
        hold off
        box on
        xlabel(xdefs{i});
        ylabel(['P(' xdefs{i} '| data)']);
        legend('Hohna', 'Snyder', 'location', 'best');
        title('Marginal posteriors around true (black) value');
    end
    
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


