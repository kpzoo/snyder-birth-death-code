% Function to perform time varying BDP Snyder filter
function [qn, tn, qev] = snyderFilterPar(m, tspec, nData, xsetMx, rhotT, lamt, mut, nLin, q0)

% Assumptions and modifications
% - accounts for identical speciation times
% - uses parfor in getNeeRatesTrapzPar or cellfun for getNeeRatesTrapz2
% - took prior out of this function
% - does standard Snyder analysis outputs time and posterior
% - added output of event time posteriors for visualisation

% Set maximum speciation time
Tsp = max(tspec);

% Posterior vectors on events
qev = zeros(nData+1, m);
qev(1, :) = q0;
tev = zeros(nData, 1);

% Cell to save output of ODE solver and set options
qset = cell(1, 1);
tset = cell(1, 1);
options = odeset('NonNegative', 1:m);
elemLen = zeros(1, nData);

% Loop across coalescent events and perform filtering
for i = 1:nData
    % Current no. lineages
    nLinCurr = nLin(i);
    
    if tspec(i) ~= tspec(i+1)
        % Solve linear ODEs continuously with setting of options, no Q
        [tsol, qsol] = ode113(@(ts, y) odeSnyBDTrapzPar(ts, y, xsetMx, Tsp, nLinCurr, rhotT, lamt, mut),...
            [tspec(i) tspec(i+1)], qev(i, :)', options);
    else
        % Identical coalescent times - just assign last value
        tsol = tspec(i+1);
        qsol = qev(i, :);
        disp(['Duplicate event time at ' num2str(i)]);
    end
    
    % Assign the output values of time and posteriors
    qset{i} = qsol;
    tset{i} = tsol;
    elemLen(i) = length(tsol);
    
    % Perturb the q posterior for the new event
    pert = getNeeRatesTrapzPar2(tsol(end), xsetMx, Tsp, nLinCurr, rhotT, lamt, mut);
    qev(i+1, :) = qsol(end, :);
    tev(i+1) = tsol(end);
    qev(i+1, :) = qev(i+1, :).*pert./(qev(i+1, :)*pert');
    
    disp(['Finished: ' num2str(i) ' of ' num2str(nData)]);
    
end

% Get full length of ODE solution data and assign appending vectors
lenFull = sum(elemLen);
stop = 0;
qn = -ones(lenFull, m);
tn = -ones(lenFull, 1);

% Append the cell based posterior and time data into a single structure
for i = 1:nData
    % Loop calculates the start and end points along the array successively
    % and then assigns the appropriate cell element
    start = stop + 1;
    stop = elemLen(i) + start - 1;
    tn(start:stop) = tset{i};
    qn(start:stop, :) = qset{i};
end
