function [Lhat P yhat] = graph_classify_bayes_plugin(Atrn,Gtrn,alg,Atst,Gtst)

P.rho1_mle = zeros(Gtrn.p,1);
P.rho0_mle = P.rho1_mle;
graphs=dec2bin(0:Gtrn.p-1);

if nargin==3                            % in-sample classifier for debugging purposes
    Atst = Atrn;
    Gtst = Gtrn;
end


% get mle params
for l=Gtrn.y0
    A=Atrn(:,:,l);
    A=num2str(A(Gtrn.inds))';
    for m=1:Gtrn.p
        if strcmp(A,graphs(m,:)),
            P.rho0_mle(m)=P.rho0_mle(m)+1;
        end
    end
end

for l=Gtrn.y1
    A=Atrn(:,:,l);
    A=num2str(A(Gtrn.inds))';
    for m=1:Gtrn.p
        if strcmp(A,graphs(m,:)),
            P.rho1_mle(m)=P.rho1_mle(m)+1;
        end
    end
end
P.rho0_mle=P.rho0_mle+1/(Gtrn.s*2);
P.rho0_mle=P.rho0_mle/sum(P.rho0_mle);

P.rho1_mle=P.rho1_mle+1/(Gtrn.s*2);
P.rho1_mle=P.rho1_mle/sum(P.rho1_mle);

if isfield(alg,'tru'), tru=true; else tru = false;  end
if isfield(alg,'mle'), mle=true; else mle = false;  end

for i=1:Gtst.s

    if tru                          % classify using only true signal edges
        yhat.tru(i)  = bayes_plugin(Atst(:,:,i),alg.rho_tru);
    end

    if mle                           % naive bayes classifier
        yhat.mle(i)  = bayes_plugin(Atst(:,:,i),P.rho1_mle, P.rho0_mle); % estimated class identities
    end

end

fn=fieldnames(yhat);                % names of classifiers
for i=1:length(fn)                  % for each, compute Lhat
    correct_vect = abs(yhat.(fn{i})-Gtst.ys);
    Lhat.(fn{i}) = mean(correct_vect);  % percent correct
    Lvar.(fn{i}) = var(correct_vect);   % variance of correct
end


    function y = bayes_plugin(datum,rho1,rho0)
        A=num2str(datum(Gtrn.inds))';
        for m=1:Gtrn.p
            if strcmp(A,graphs(m,:)),
                break
            end
        end
        y=rho1(m)>rho0(m);
    end

end


