% Compares the likelihood functions from Snyder and Hohna (TESS)
% simulations by numerically calculating across the parameter space
clearvars
clc
close all

% Simulation parameters
%mstar = [20 15];
%rateID = [1 3];
mstar = 15;
rateID = 2;
nSim = length(rateID);

% Set of possible functions
rateSet = {'hohna', 'rateshift', 'logistic', 'constant'};

% Variables to store output
Lsny = cell(1, nSim);
Lhoh = cell(1, nSim);
LmargS = cell(1, nSim);
LmargH = cell(1, nSim);
outPm = cell(1, nSim);

% Loop across rate functions and calculate log likelihoods
for i = 1:nSim
    % Main code
    [Lsny{i}, Lhoh{i}, LmargS{i}, LmargH{i}, outPm{i}] = compareTimeVaryingLiksFn(rateID(i), mstar(i));
    % Save data sequentially and update progress
    save(['lik_' num2str(i) '.mat']);
    disp(['Finished: ' num2str(i) 'of ' num2str(nSim)]);
end

% Legend names
leg = {'Snyder', 'Hohna'};

% Plot the likelihoods over the parameter set (unordered)
for i = 1:nSim
    % Extract some parameter data
    numRV = length(outPm{i}.x);
    m = mstar(i)^numRV;
    xset = outPm{i}.xset;
    rateName = rateSet{rateID(i)};
    
    figure;
    plot(1:m, Lsny{i}, 1:m, Lhoh{i});
    xlabel('points');
    ylabel('log likelihood');
    legend(leg, 'location', 'best');
    title(['Likelihoods for ' rateName ' model']);
    
    % Plot the marginals
    figure;
    for j = 1:numRV
        subplot(ceil(numRV/2), 2, j);
        hold on
        plot(xset{j}, LmargS{i}{j}, xset{j}, LmargH{i}{j});
        hold off
        xlabel(['x_' num2str(j)]);
        legend(leg, 'location', 'best');
        ylabel(['log L(data | x_' num2str(j) ')']);
    end
end
