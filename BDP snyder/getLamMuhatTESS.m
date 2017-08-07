% Function for plotting estimated rates for use with tessSnyder
function [lamhat, muhat, lamn, mun] = getLamMuhatTESS(xsetMx, tn, qnlast, rateID, mut, lamt, empiricalData, x, changeTime)

% Obtain true values if available
if ~empiricalData
    lamn = lamt(x, tn);
    mun = mut(x, tn);
end

% Obtain estimates of parametric birth/death rates from final posterior
switch(rateID) 
    case 1
        % TESS function from Hohna 2014
        muhat = zeros(size(tn));
        lamhat = zeros(size(tn));
        for i = 1:length(tn)
            lamhat(i) = qnlast*(xsetMx(1, :) + xsetMx(2, :).*exp(-xsetMx(3, :)*tn(i)))';
            muhat(i) = qnlast*(xsetMx(1, :))';
        end
        
    case 2
        % Rateshift function from Hohna TESS tutorial
        muhat = zeros(size(tn));
        lamhat = zeros(size(tn));
        % Possible values to switch between pre-computed outside loop
        lamval1 = qnlast*(xsetMx(1, :) + xsetMx(2, :))';
        lamval2 = qnlast*(xsetMx(3, :) + xsetMx(4, :))';
        muval1 = qnlast*xsetMx(2, :)';
        muval2 = qnlast*xsetMx(4, :)'; 
        % Choose rates based on change time, which is same for lam and mu
        for i = 1:length(tn)
            if tn(i) < changeTime
                lamhat(i) = lamval1;
                muhat(i) = muval1;
            else
                lamhat(i) = lamval2;
                muhat(i) = muval2;
            end
        end  
        
    case 3
        % Logistic birth and death functions from Paradis 2010
        muhat = zeros(size(tn));
        lamhat = zeros(size(tn));
        for i = 1:length(tn)
            lamhat(i) = qnlast*(1./(1 + exp(-xsetMx(1, :)*tn(i) + xsetMx(2, :))))';
            muhat(i) = qnlast*(1./(1 + exp(-xsetMx(3, :)*tn(i) + xsetMx(4, :))))';
        end
        
        
    otherwise
        error('Rate id input not supported');
end