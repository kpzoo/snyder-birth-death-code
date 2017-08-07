% Function calculates Stadler likelihoods from Stadler 2013
function logf = likStadlerLog(idx, ts, t0, lam, mu, sig, n)

% Assumptions and modifications
% - everything done in log form for numerical stability
% - ts: vector of times in the coalescent direction, does not include t0
% - sampling prob of 1
% - entire function refers to reversed times, hence starting at t0 really
% means the process ends at t0


%% Likelihood calculations

% Check size of ts - cannot include t0
if length(ts) ~= n-1
    error('The coalescent times should not have t0');
end

% Shortened version of pn function, with lam, mu, sig already included
p = @(nid, t)pn(nid, t, lam, mu, sig);
% Further brief function of only t, 1 lineage only
p1 = @(t)pn(1, t, lam, mu, sig);

% Depending on if conditioning on crown
if(idx > 3 && idx < 7)
    % Separate the crown age in reverse time, last values as ts sorted
    t1 = ts(end);
    % Branching times now n-2 in length
    ts = ts(1:end-1); 
end

% All likelihoods require this probability calculation on branching times
pset = zeros(size(ts));
for i = 1:length(ts)
    pset(i) = log10(p1(ts(i)));
end
prodSet = sum(pset); % it is a sum of log probs

% Get likelihood functions based on different conditionings
switch(idx)
    case 1
        % Condition on process starting at t0
        %f = p1(t0)*(lam^(n-1))*prodSet;
        logf = log10(p1(t0)) + (n-1)*log10(lam) + prodSet;
    case 2
        % Condition on process between t0 and 0 surviving
        %f = p1(t0)*(lam^(n-1))*prodSet/(1 - p(0, t0));
        logf = log10(p1(t0)/(1 - p(0, t0))) + (n-1)*log10(lam) + prodSet;
    case 3
        % Condition on observing n species
        q0 = (lam/mu)*p(0, t0);
        %f = ((lam/q0)^(n-1))*prodSet;
        logf = (n-1)*log10(lam/q0) + prodSet;
        
    case 4
        % Condition on 2 BDPs starting at t1
        %f = (p1(t1)^2)*(lam^(n-2))*prodSet;
        logf = 2*log10(p1(t1)) + (n-2)*log10(lam) + prodSet;
        
    case 5
        % Same as 4 but also conditioned on both BDPs surviving
        %f = (p1(t1)^2)*(lam^(n-2))*prodSet/((1 - p(0, t1))^2);
        logf = 2*log10(p1(t1)/(1 - p(0, t1))) + (n-2)*log10(lam) + prodSet;
    case 6
        % Condition on n species and crown t1
        q0 = (lam/mu)*p(0, t0);
        %f = ((lam/q0)^(n-2))*prodSet/(n-1);
        logf = (n-2)*log10(lam/q0) - log10(n-1) + prodSet;
        
    case 7
        % Condition on n integrated on all possible t0 (uniform, improper)
        t1 = ts(end);
        %f = n*p1(t1)*(lam^(n-1))*prodSet/(1 - p(0, t1));
        logf = log10(n*p1(t1)/(1 - p(0, t1))) + (n-1)*log10(lam) + prodSet;
end

%% Sub-functions to get probabilities of interest

% Prob of observing n lineages after sampling in a birth-death tree of age t
function prob = pn(nid, t, lam, mu, sig)

switch(nid)
    case 0
        % Zero lineages sampled
        prob = 1 - sig/(lam - mu*exp(-sig*t));
    case 1
        % Single lineage
        prob = ((sig^2)*exp(-sig*t))/((lam - mu*exp(-sig*t))^2);
    otherwise
        % More than 1 lineage
        p0_l = 1 - sig/(lam - mu*exp(-sig*t));
        p1_l = ((sig^2)*exp(-sig*t))/((lam - mu*exp(-sig*t))^2);
        q = (lam/mu)*p0_l;
        prob = p1_l*q^(nid-1);
end