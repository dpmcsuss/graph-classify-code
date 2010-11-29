function Lhat = graph_classify_spectral(Atrn,Gtrn,alg,Atst,Gtst)

if nargin==3                            % in-sample classifier for debugging purposes
    P = get_params_mle(Atrn,Gtrn); 
    Atst = Atrn;
    Gtst = Gtrn;
else                                    % update parameters using only training data   
    P = get_params_mle(Atrn(:,:,1:Gtrn.s),Gtrn);  
end

lap0 = get_laplacian(P.E0);
lap1 = get_laplacian(P.E1);

keyboard

Lhat=0.5;

end
