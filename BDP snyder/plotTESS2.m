% Function to do main plotting of single run tessSnyder3
function plotTESS2(numRV, samples, qmarg, xset, x, datapath, nosave, sampSny, treatDiscrete,...
    xdefs, tn, lamn, mun, lamhat, muhat, sampleFail, empiricalData, thisDir, qev, IDMx, mi)

% Assumptions and modifications
% - added plot of net diversification and turnover
% - removed ratename, added nosave
% - added smoothed density plot
% - removed histogram comparison
% - stopped saving raw data to rateName


% Direct marginal cdf comparison
figure;
for i = 1:numRV
    subplot(ceil(numRV/2), 2, i);
    h1 = histogram(samples(:, i), 'Normalization', 'probability');
    % Values of histogram at midpoints of bins
    edg = h1.BinEdges;
    eav = 0.5*(edg(2:end) + edg(1:end-1));
    pav = h1.Values;
    plot(eav, cumsum(pav), 'r');
    hold on
    % Plot all posteriors Snyder
    qval = qmarg{i};
    if sum(qval) ~= 1
        % Marginals that integrate to 1
        plot(xset{i}, cumtrapz(xset{i}, qval), 'b-.');
    else
        % Marginals that sum to 1
        plot(xset{i}, qval, 'b-.');
    end
    % True values
    plot([x(i) x(i)], [0 1], 'r', 'linewidth', 2);
    hold off
    xlabel(['x_' num2str(i)]);
    ylabel(['P(x_' num2str(i) ')']);
    legend('Hohna', 'Snyder', 'true', 'location', 'best');
    grid;
    title('Marginal cdfs');
end

if ~sampleFail
%     % Plot comparative histograms
%     figure;
%     for i = 1:numRV
%         subplot(ceil(numRV/2), 2, i);
%         hbin = histogram(samples(:, i), 'Normalization', 'probability');
%         binWdth = hbin.BinWidth;
%         hold on
%         % Plot all posterior histograms from Snyder
%         hbin = histogram(sampSny(:, i), 'Normalization', 'probability');
%         hbin.BinWidth = binWdth; % to keep bin widths same for comparison
%         % True values
%         h1 = gca;
%         yl = h1.YLim;
%         plot([x(i) x(i)], yl, 'k', 'linewidth', 2);
%         hold off
%         xlabel(['x_' num2str(i)]);
%         ylabel(['P(x_' num2str(i) ')']);
%         legend('Hohna', 'Snyder', 'true', 'location', 'best');
%         grid;
%         title('Marginal posterior histograms');
%     end
%     cd(datapath);
%     saveas(gcf, 'histoData');
%     cd ..
%     cd ..
    
    % Plot comparative smoothed densities
    figure;
    for i = 1:numRV
        subplot(ceil(numRV/2), 2, i);
        % Kernel smoothing, 100 pts with normal density
        [pSny, xSny, bwSny] = ksdensity(sampSny(:, i));
        if treatDiscrete
            % Increase the bandwidth since discrete marginals
            [pSny, xSny] = ksdensity(sampSny(:, i), 'width', 5*bwSny);
        end
        [pMCMC, xMCMC] = ksdensity(samples(:, i));
        plot(xMCMC, pMCMC, xSny, pSny, 'linewidth', 2);
        hold on
        % True values
        h1 = gca;
        yl = h1.YLim;
        plot([x(i) x(i)], yl, 'k', 'linewidth', 2);
        hold off
        xlabel(['x_' num2str(i)]);
        ylabel(['P(x_' num2str(i) ' | data)']);
        legend('Hohna', 'Snyder', 'true', 'location', 'best');
        grid;
        title('Smoothed marginal posteriors');
    end
    if ~nosave
        cd(datapath);
        saveas(gcf, 'pdf');
        cd(thisDir);
    end
    
    % Autocorrelation functions of both sample sets
    % Look at autocorrelation function of samples of parameters
    figure;
    subplot(2, 1, 1);
    hold all
    for i = 1:numRV
        [acf,lags,~] = autocorr(samples(:, i));
        plot(lags, acf, 'o:');
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
        [acf,lags,~] = autocorr(sampSny(:, i));
        plot(lags, acf, 'o:');
    end
    hold off
    grid;
    xlabel('lag');
    ylabel('autocorrelation');
    legend(xdefs, 'location', 'best');
    title('Autcorrelation of Snyder samples');
    if ~nosave
        cd(datapath);
        saveas(gcf, 'autocorr');
        cd(thisDir);
    end
else
    % MCMC autocorrelation
    % Look at autocorrelation function of samples of parameters
    figure;
    hold all
    for i = 1:numRV
        [acf,lags,~] = autocorr(samples(:, i));
        plot(lags, acf, 'o:');
    end
    hold off
    grid;
    xlabel('lag');
    ylabel('autocorrelation');
    legend(xdefs, 'location', 'best');
    title('Autcorrelation of MCMC samples');
    if ~nosave
        cd(datapath);
        saveas(gcf, 'autocorr');
        cd(thisDir);
    end
end


% Plot birth and death rate estimates for single run with true values if
% known on subplots else on 1 plot
if ~empiricalData
    figure;
    subplot(2, 1, 1)
    plot(tn, lamn, tn, lamhat, 'linewidth', 2);
    xlabel('time');
    ylabel('birth rates');
    legend('true', 'estimated');
    grid;
    subplot(2, 1, 2)
    plot(tn, mun, tn, muhat, 'linewidth', 2);
    xlabel('time');
    ylabel('death rates');
    legend('true', 'estimated');
    grid;
    if ~nosave
        cd(datapath);
        saveas(gcf, 'rates');
        cd(thisDir);
    end
else
    % No true rates known so plot estimates
    figure;
    [hAx, hLine1, hLine2] = plotyy(tn, lamhat, tn, muhat);
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
end

% Plot net diversif and turnover estimates if true values
if ~empiricalData
    figure;
    subplot(2, 1, 1)
    plot(tn, lamn-mun, tn, lamhat-muhat, 'linewidth', 2);
    xlabel('time');
    ylabel('net diversification');
    legend('true', 'estimated');
    grid;
    subplot(2, 1, 2)
    plot(tn, mun./lamn, tn, muhat./lamhat, 'linewidth', 2);
    xlabel('time');
    ylabel('turnover ratio');
    legend('true', 'estimated');
    grid;
    if ~nosave
        cd(datapath);
        saveas(gcf, 'ratesTrans');
        cd(thisDir);
    end
end

% Plot the event changes of the marginals (that sum to 1) with time
[nIt, ~] = size(qev);
figure;
for i = 1:nIt
    % Get event marginals which sum to 1
    qmargEv = marginalise(numRV, IDMx, qev(i, :), mi);
    for j = 1:numRV
        subplot(ceil(numRV/2), 2, j);
        hold on
        plot(xset{j}, qmargEv{j});
        hold off
        xlabel(['x_' num2str(j)]);
        ylabel(['P(x_' num2str(j) ' | data)']);
    end
end
        
        

