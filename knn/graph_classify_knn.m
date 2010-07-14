function Lhat = graph_classify_knn(Atrn,Gtrn,alg,Atst,Gtst)

d=Gtrn.n^2;

xTr = double(reshape(Atrn,d,Gtrn.s0+Gtrn.s1));
yTr = zeros(1,Gtrn.s0+Gtrn.s1);
yTr(Gtrn.y0) = 0;
yTr(Gtrn.y1) = 1;

xTe = double(reshape(Atst,d,Gtst.s0+Gtst.s1));
yTe = zeros(1,Gtst.s0+Gtst.s1);
yTe(Gtst.y0) = 0;
yTe(Gtst.y1) = 1;

k=round(sqrt(min(Gtrn.s0,Gtrn.s1)));

if ~isfield(alg,'knn_vanilla'),  alg.knn_vanilla= true; end
if ~isfield(alg,'knn_lmnn'),     alg.knn_lmnn   = false; end

if alg.knn_vanilla
    knnerrI     = knnclassify(eye(d),xTr,yTr,xTe,yTe,k);  
    Lhat.knn    = knnerrI(2); 
end

if alg.knn_lmnn
    L           = lmnn(xTr,yTr,'quiet',1);
    knnerrL     = knnclassify(L,xTr,yTr,xTe,yTe,k);       
    Lhat.lmnn   = knnerrL(2); 
    % enerr=energyclassify(L,xTr,yTr,xTe,yTe,k);
    % Lhats.energy = enerr(2);
end

