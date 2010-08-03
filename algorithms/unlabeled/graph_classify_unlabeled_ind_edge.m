function [Lhat ind Lvar P yhat] = graph_classify_unlabeled_ind_edge(Atrn,Gtrn,alg,Atst0,Atst1,Gtst,P)
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
%   ind:    index of edges used for each algorithm implemented
%   Lvar:   misclassification variance for each algorithm implemented
%   P:      parameter estimates
%   yhat:   list of estimated class identity for each graph

if nargin==6                            % in-sample classifier for debugging purposes
    P = get_params_mle(Atrn,Gtrn); 
%     Atst = Atrn;
%     Gtst = Gtrn;
else                                    % update parameters using only training data   
%     P = get_params_mle(Atrn(:,:,1:Gtrn.s),Gtrn);  
end 

% initialize stuff for the different classifiers
if isfield(alg,'signal_subgraph_ind'),  tru=true;                                       
else tru= false;  end

if isfield(alg,'nb_ind'),               nb=true;                                        
else nb = false;  end

if isfield(alg,'num_inc_edges'),        inc=true;  ind(1).inc = zeros(1,alg.num_inc_edges);    
else  inc = false; end

if isfield(alg,'num_coh_vertices'),     coh=true;  ind(1).coh = zeros(1,alg.num_coh_vertices); 
else coh = false; end

if isfield(alg,'num_inc_edges') || isfield(alg,'num_coh_vertices')
    ind(2:Gtst.s)=ind(1);               % pre-allocate memory to speed stuff up
else
    ind=[];
end

for i=1:Gtst.s

    if tru                          % classify using only true signal edges
        yhat.tru(i)  = ie_classify(Atst0(:,:,i),Atst1(:,:,i),P,alg.signal_subgraph_ind);
    end

    if nb                           % naive bayes classifier
        yhat.nb(i)  = ie_classify(Atst0(:,:,i),Atst1(:,:,i),P,alg.nb_ind); % estimated class identities
    end

    if inc                          % incoherent edge classifier
        ind(i).inc = get_inc_edges(P.d_opt,alg.num_inc_edges);
        yhat.inc(i)= ie_classify(Atst0(:,:,i),Atst1(:,:,i),P,ind(i).inc);
    end

    if coh                          % coherent edge classifier
        ind(i).coh  = get_max_edges(P.d_opt,alg.num_coh_vertices);
        yhat.coh(i) = ie_classify(Atst0(:,:,i),Atst1(:,:,i),P,ind(i).coh);
    end

end

fn=fieldnames(yhat);                % names of classifiers
for i=1:length(fn)                  % for each, compute Lhat
    correct_vect = abs(yhat.(fn{i})-Gtst.ys);
    Lhat.(fn{i}) = mean(correct_vect);  % percent correct
    Lvar.(fn{i}) = var(correct_vect);   % variance of correct
end

end

function y = ie_classify(datum0,datum1,P,ind)
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

if any(datum0(ind))>1 || any(datum1(ind))>1            % if poisson
    post0=sum(sum(datum(ind).*P.lnE0(ind) - P.E0(ind)))+P.lnprior0;
    post1=sum(sum(datum(ind).*P.lnE1(ind) - P.E1(ind)))+P.lnprior1;
else                            % if bernoulli
    post0=sum(datum0(ind).*P.lnE0(ind)+(1-datum0(ind)).*P.ln1E0(ind))+P.lnprior0;
    post1=sum(datum1(ind).*P.lnE1(ind)+(1-datum1(ind)).*P.ln1E1(ind))+P.lnprior1;
end

[foo bar] = sort([post0, post1]);
y=bar(2)-1;

end



