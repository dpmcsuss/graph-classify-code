function [Lhat Lvar P yhat] = graph_classify_ER(Atrn,Gtrn,Atst,Gtst)

if nargin==3                            % in-sample classifier for debugging purposes
    Atst = Atrn;
    Gtst = Gtrn;
end 

P.eta0 = mean(mean(mean(Atrn(:,:,Gtrn.y0))));
P.eta1 = mean(mean(mean(Atrn(:,:,Gtrn.y1))));

P.pi0=Gtrn.s0/Gtrn.s;
P.pi1=Gtrn.s1/Gtrn.s;

for i=1:Gtst.s
        yhat.er(i)  = er_classify(Atst(:,:,i),P);
end

fn=fieldnames(yhat);                % names of classifiers
for i=1:length(fn)                  % for each, compute Lhat
    correct_vect = abs(yhat.(fn{i})-Gtst.ys);
    Lhat.(fn{i}) = mean(correct_vect);  % percent correct
    Lvar.(fn{i}) = var(correct_vect);   % variance of correct
end

end

function y = er_classify(datum,P)

sum0=sum(datum(:));
sum1=sum(1-datum(:));

y= (log(P.eta1)*sum0+log(1-P.eta1)*sum1) + log(P.pi1) ...
    > log(P.eta0)*sum0+log(1-P.eta0)*sum1 + log(P.pi0);

end



