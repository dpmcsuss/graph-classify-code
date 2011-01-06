function [Lstats yhat] = plugin_classify(As,classes,params,subspace)
% this script implements the bayes plugin classifier, using a (sub-)space
% specified by user
%
% INPUT:
%   As:         an array (currently typically n x n x s) of predictor semiables
%   classes:    a list of classes (typically 1 x s)
%   params:     a structure containing parameters for the algorithm to use
%   subspace:   structure containing name, indices, and size of each subspace
%
% OUTPUT:
%   Lstats:     structure containing mean, s.e.m for each subspace
%   yhat:       list of estimated class identity for each graph

siz     = size(As);
s       = siz(end);                         % # of samples (robust to As being an arbitrary array)
poiss   = any(As(:))>1;          % whether data is Poisson or Bernoulli
n_subs  = length(subspace);
yhat    = nan(s,n_subs);

for i=1:s
    datum = As(:,:,i);               % this line only makes sense when data are graphs    
    for j=1:n_subs
        ind = subspace(j).indices;         % for code legibility
        
        data_tmp=datum(ind);
        
        if poiss                        % if poisson
            post0=sum(sum(data_tmp.*params.lnE0(ind) - params.E0(ind)))+params.lnprior0;
            post1=sum(sum(data_tmp.*params.lnE1(ind) - params.E1(ind)))+params.lnprior1;
        else                            % if bernoulli
            post0=sum(data_tmp.*params.lnE0(ind)+(1-data_tmp).*params.ln1E0(ind))+params.lnprior0;
            post1=sum(data_tmp.*params.lnE1(ind)+(1-data_tmp).*params.ln1E1(ind))+params.lnprior1;
        end
        
        [~, bar] = sort([post0, post1]); % find the bigger one
        yhat(i,j)=bar(2)-1;
    end
end

siz=size(classes);
if siz(1)<siz(2), classes=classes'; end

s_sqrt = sqrt(s);
for j=1:n_subs                                  % for each, compute stats
    incorrect_vect = yhat(:,j)~=classes;
    Lstats.avg(j) = mean(incorrect_vect);         % percent correct
    Lstats.sem(j) = std(incorrect_vect)/s_sqrt;   % s.e.m of correct
end