function [Lhat ind Lvar P yhat] = graph_classify_ER(Atrn,Gtrn,alg,Atst,Gtst)

if nargin==3                            % in-sample classifier for debugging purposes
    Atst = Atrn;
    Gtst = Gtrn;
end 

P.eta0 = mean(mean(mean(Atrn(:,:,Gtrn.y0))));
P.eta1 = mean(mean(mean(Atrn(:,:,Gtrn.y1))));

for i=1:Gtst.s

        yhat.tru(i)  = ie_classify(Atst(:,:,i),P,alg.signal_subgraph_ind);

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
    post1=sum(sum(datum(ind).*P.lnE1(ind) - P.E1(ind)))+P.lnprior1;
else                            % if bernoulli
    post0=sum(datum(ind).*P.lnE0(ind)+(1-datum(ind)).*P.ln1E0(ind))+P.lnprior0;
    post1=sum(datum(ind).*P.lnE1(ind)+(1-datum(ind)).*P.ln1E1(ind))+P.lnprior1;
end

[foo bar] = sort([post0, post1]);
y=bar(2)-1;

end



