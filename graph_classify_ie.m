function [Lhat Lvar ind P yhat] = graph_classify_ie(Atrn,Gtrn,alg,Atst,Gtst)
% this script classifies using a number of different approaches
% INPUT:
% Atrn:     n x n x s array, where |V|=n, and s is the number of samples
% Gtrn:     a structure containing constants for training data (see get_constants.m for details)
% alg:      a structure containing fields specifying algorithm parameters.  list of fields:
%   signal_subgraph_ind:    if exists, do naive bayes classification using these indices only
%   nb_ind:                 if exists, do naive bayes classification using all indices
%   num_inc_edges:          if exists, do incoherent signal classification, using this number of edges
%   num_coh_edges:          if exists, do incoherent signal classification, using this number of edges
% Atst:     same as Atrn but for test data
% Gtst:     same as Gtrn but for test data
% 
% OUTPUT:
%   Lhat:   misclassification rate for each algorithm implemented
%   Lvar:   misclassification variance for each algorithm implemented
%   ind:    index of edges used for each algorithm implemented
%   P:      parameter estimates
%   yhat:   list of estimated class identity for each graph

if nargin==3                        % in-sample classifier for debugging purposes
    P = get_params(Atrn,Gtrn); 
    Atst = Atrn;
    Gtst = Gtrn;
else
    P = get_params(Atrn(:,:,1:Gtrn.s),Gtrn); % update parameters using only other samples    
end 

% initialize stuff for the different classifiers

if isfield(alg,'signal_subgraph_ind'),  tru=true;                                       else tru= false;  end
if isfield(alg,'nb_ind'),               nb=true;                                        else nb = false;  end
if isfield(alg,'num_inc_edges'),        inc=true;  ind(Gtst.s).inc = zeros(1,Gtrn.n);   else inc = false; end
if isfield(alg,'num_coh_edges'),        coh=true;  ind(Gtst.s).coh = zeros(1,Gtrn.n);   else coh = false; end

for i=1:Gtst.s

    if tru                          % classify using only true signal edges
        yhat.tru(i)  = ie_classify(Atst(:,:,i),P,alg.signal_subgraph_ind);
    end

    if nb                           % naive bayes classifier
        yhat.nb(i)  = ie_classify(Atst(:,:,i),P,alg.nb_ind); % estimated class identities
    end

    if inc                          % incoherent edge classifier
        ind(i).inc = get_inc_edges(P.d_pos,alg.num_inc_edges);
        yhat.inc(i)= ie_classify(Atst(:,:,i),P,ind(i).inc);
    end

    if coh                          % coherent edge classifier
        ind(i).coh  = get_max_edges(P.d_pos);
        yhat.coh(i) = ie_classify(Atst(:,:,i),P,ind(i).coh);
    end

end

fn=fieldnames(yhat);                % names of classifiers
for i=1:length(fn)                  % for each, compute Lhat
    correct_vect = abs(yhat.(fn{i})-Gtst.ys);
    Lhat.(fn{i}) = mean(correct_vect);  % percent correct
    Lvar.(fn{i}) = var(correct_vect);   % variance of correct
end

end

function y = ie_classify(datum,P,ind)
% this function classifies using independent edge assumption.  each edge
% could be distributed according to a poisson distribution, or a
% bernoulli. the class conditional posterior is computed as appropriate.
% 
% INPUT
% datum:    the graph to be classified
% P:        structure of parameter estimates to use for classification
% ind:      indices to use in classifier
% 
% OUTPUT:
% y:        estimated class

if any(datum(ind))>1            % if poisson
    post0=sum(sum(datum(ind).*P.lnE0(ind) - P.E0(ind)))+P.lnprior0;
    post1=sum(sum(datum(ind).*P.lnE1(ind) - P.E1(ind)))+P.lnprior0;
else                            % if bernoulli
    post0=sum(datum(ind).*P.lnE0(ind)+(1-datum(ind)).*P.ln1E0(ind))+P.lnprior0;
    post1=sum(datum(ind).*P.lnE1(ind)+(1-datum(ind)).*P.ln1E1(ind))+P.lnprior0;
end

[foo bar] = sort([post0, post1]);
y=bar(2)-1;

end



