% Function for plotting estimated rates in Agamid data
function [lamhat, muhat] = plotLamMuhatRabosky(xsetMx, xhat, tn, qnlast, rateID, mut, lamt)

% Obtain estimates of parametric birth/death reates from posterior
switch(rateID)
    case 1
        % Constant rate birth-death (can directly subst cond means)
        muhat = mut(xhat, tn);
        lamhat = lamt(xhat, tn);
        
    case 2
        % Rabosky 2008 SPVAR model
        muhat = mut(xhat, tn);
        lamhat = zeros(size(tn));
        for i = 1:length(tn)
            lamhat(i) =  qnlast*(xsetMx(1, :).*exp(-xsetMx(2, :)*tn(i)))';
        end
        
    case 3
        % Rabosky 2008 EXVAR model
        lamhat = lamt(xhat, tn);
        muhat = zeros(size(tn));
        for i = 1:length(tn)
            muhat(i) =  qnlast*(xsetMx(3, :).*(1 - exp(-xsetMx(2, :)*tn(i))))';
        end
        
    case 4
        % Rabosky 2008 BOTHVAR model
        muhat = zeros(size(tn));
        lamhat = zeros(size(tn));
        for i = 1:length(tn)
            lamhat(i) =  qnlast*(xsetMx(1, :).*exp(-xsetMx(2, :)*tn(i)))';
            muhat(i) =  qnlast*(xsetMx(4, :).*(1 - exp(-xsetMx(2, :)*tn(i))))';
        end
        
    case 5
        % Modified SPVAR model with no deaths
        muhat = zeros(size(tn));
        lamhat = zeros(size(tn));
        for i = 1:length(tn)
            lamhat(i) =  qnlast*(xsetMx(1, :).*exp(-xsetMx(2, :)*tn(i)))';
        end
end