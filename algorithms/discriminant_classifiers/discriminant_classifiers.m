function Lhat = discriminant_classifiers(preds,targs,params,alg)
% performs various machine learning type classifiers on feature spaces
% INPUT
%   preds:  a matrix (k x s) of k feature vectors for s samples
%   targs:  a vector (k x 1) of the binary class identity for each of s samples
%   params: parameters for classifiers
%   alg:    specifies which algorithms to use
% OUTPUT:
%   Lhat:   a structure containing Lhat for each of the classifiers used


siz = size(preds);
S   = siz(end);                         % # of samples (robust to preds being an arbitrary array)


if isfield(alg,'LDA'), LDA=1; else LDA=0; end
if isfield(alg,'QDA'), QDA=1; else QDA=0; end

for i=1:S
    if LDA
        e0=(preds(:,i)-params.mu0);
        e1=(preds(:,i)-params.mu1);       
        yhat.LDA(i)=e0'*params.InvSig*e0 - params.lnprior0 > e1'*params.InvSig*e1 - params.lnprior1;
    end
    
    if QDA
        e0=(preds(:,i)-params.mu0);
        e1=(preds(:,i)-params.mu1);       
        yhat.QDA(i)=e0'*params.InvSig0*e0 - params.lnprior0 > e1'*params.InvSig1*e1 - params.lnprior1;
    end        
end


fn = fieldnames(yhat);                  % names of subspaces
n_subspace = length(fn);                % # of different subspace

for j=1:n_subspace                      % for each, compute stats
    correct_vect = abs(yhat.(fn{j})-targs);
    Lhat.(fn{j}) = mean(correct_vect);  % percent correct
    Lvar.(fn{j}) = var(correct_vect);   % variance of correct
end

end