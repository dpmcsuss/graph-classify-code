function [Lhat Lsem] = discriminant_classifiers(preds,ys,params,discrim)
% performs various machine learning type classifiers on feature spaces
% INPUT
%   preds:  a matrix (k x s) of k feature vectors for s samples
%   ys:  a vector (k x 1) of the binary class identity for each of s samples
%   params: parameters for classifiers
%   discrim:specifies which discriminant classifiers to use
% OUTPUT:
%   Lhat:   a structure containing Lhat for each of the classifiers used


siz = size(preds);
S   = siz(end);                         % # of samples (robust to preds being an arbitrary array)


if isfield(discrim,'LDA'), LDA=1; else LDA=0; end
if isfield(discrim,'QDA'), QDA=1; else QDA=0; end
if isfield(discrim,'dLDA'), dLDA=1; else dLDA=0; end
if isfield(discrim,'dQDA'), dQDA=1; else dQDA=0; end

for i=1:S
    if LDA
        e0=(preds(:,i)-params.mu0);
        e1=(preds(:,i)-params.mu1);
        yhat.LDA(i)=e0'*(params.Sig\e0) - params.lnprior0 > e1'*(params.Sig\e1) - params.lnprior1;
    end
    
    if dLDA
        e0=(preds(:,i)-params.mu0);
        e1=(preds(:,i)-params.mu1);
        yhat.dLDA(i)=e0'*params.InvdSig*e0 - params.lnprior0 > e1'*params.InvdSig*e1 - params.lnprior1;
    end
    
    if QDA
        e0=(preds(:,i)-params.mu0);
        e1=(preds(:,i)-params.mu1);
        yhat.QDA(i)=e0'*(params.Sig0\e0) - params.lnprior0 > e1'*(params.Sig1\e1) - params.lnprior1;
    end
    
    if dQDA
        e0=(preds(:,i)-params.mu0);
        e1=(preds(:,i)-params.mu1);
        yhat.dQDA(i)=e0'*params.InvdSig0*e0 - params.lnprior0 > e1'*params.InvdSig1*e1 - params.lnprior1;
    end
    
end


fn = fieldnames(yhat);                  % names of subspaces
n_subspace = length(fn);                % # of different subspace

for j=1:n_subspace                      % for each, compute stats
    correct_vect = abs(yhat.(fn{j})-ys);
    Lhat.(fn{j}) = mean(correct_vect);  % percent correct
    Lsem.(fn{j}) = std(correct_vect)/sqrt(S);   % variance of correct
end

end