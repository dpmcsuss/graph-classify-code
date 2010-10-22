function [LAP QAP] = classify_unlabeled_graphs(adjacency_matrices,class_labels,alg,params)
%%
constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code

P       = params;
P.lnE0  = log(P.E0);
P.ln1E0 = log(1-P.E0);
P.lnE1  = log(P.E1);
P.ln1E1 = log(1-P.E1);

trn0= randperm(constants.s0);
trn1= randperm(constants.s1)+constants.s0;
tst = randperm(constants.s);
I   = eye(constants.n);

if strcmp(alg.names(1),'LAP'),
    LAP.do      = true;
    LAP.yhat    = zeros(constants.s/2,1);
    LAP.correct = zeros(constants.s/2,1);
    LAP.work0   = zeros(constants.s/2,1);
    LAP.work1   = zeros(constants.s/2,1);
else LAP.do = false;
end

if strcmp(alg.names(2),'QAP'),
    QAP.do      = true;
    QAP.yhat    = NaN(constants.s/2,alg.QAP_max_iters);
    QAP.correct = NaN(constants.s/2,alg.QAP_max_iters);
    QAP.work0   = NaN(constants.s/2,alg.QAP_max_iters);
    QAP.work1   = NaN(constants.s/2,alg.QAP_max_iters);
    QAP.obj0    = NaN(constants.s/2,alg.QAP_max_iters);
    QAP.obj1    = NaN(constants.s/2,alg.QAP_max_iters);
else QAP.do = false;
end

for j=1:constants.s/2
    
    Atrn(:,:,1) = adjacency_matrices(:,:,trn0(j));
    Atrn(:,:,2) = adjacency_matrices(:,:,trn1(j));
    Atst        = adjacency_matrices(:,:,tst(j));
    
    % do LAP
    if LAP.do
        LAP_ind0 = munkres(-Atrn(:,:,1)*Atst');     %Note that sum(sum(-A(LAP,:).*Atst))==f
        LAP.work0(j)=norm(Atrn(:,:,1)-Atst(LAP_ind0,:)) < norm(Atrn(:,:,1)-Atst); % check that LAP work
        
        LAP_ind1 = munkres(-Atrn(:,:,2)*Atst');     %Note that sum(sum(-A(LAP,:).*Atst))==f
        LAP.work1(j)=norm(Atrn(:,:,2)-Atst(LAP_ind1,:)) < norm(Atrn(:,:,2)-Atst); % check that LAP work
        
        LAP_lik0=sum(sum(Atst(LAP_ind0,:).*P.lnE0+(1-Atst(LAP_ind0,:)).*P.ln1E0));
        LAP_lik1=sum(sum(Atst(LAP_ind1,:).*P.lnE1+(1-Atst(LAP_ind1,:)).*P.ln1E1));
        
        [~, bar] = sort([LAP_lik0, LAP_lik1]);
        LAP.yhat(j)=bar(2)-1;
        LAP.correct(j)=(LAP.yhat(j)==class_labels(tst(j)));
    end
    
    % do QAP
    if QAP.do
        
        [~,~,~,iter0,~,QAP_inds0]=sfw(Atrn(:,:,1),-Atst,alg.QAP_max_iters,I);
        for ii=1:iter0
            QAP.obj0(j,ii) = norm(Atrn(:,:,1) - Atst(QAP_inds0{ii},QAP_inds0{ii}));
            QAP.work0(j,ii)= QAP.obj0(j,ii) < norm(Atrn(:,:,1)-Atst); % check that QAP works
        end
        
        [~,~,~,iter1,~,QAP_inds1]=sfw(Atrn(:,:,2),-Atst,alg.QAP_max_iters,I);
        for ii=1:iter1
            QAP.obj1(j,ii) = norm(Atrn(:,:,2) - Atst(QAP_inds1{ii},QAP_inds1{ii}));
            QAP.work1(j,ii)= QAP.obj1(j,ii) < norm(Atrn(:,:,1)-Atst); % check that QAP works
        end
        
        for ii=1:min(iter0,iter1)
            QAP_lik0=sum(sum(Atst(QAP_inds0{ii},QAP_inds0{ii}).*P.lnE0+(1-Atst(QAP_inds0{ii},QAP_inds0{ii})).*P.ln1E0));
            QAP_lik1=sum(sum(Atst(QAP_inds1{ii},QAP_inds1{ii}).*P.lnE1+(1-Atst(QAP_inds1{ii},QAP_inds1{ii})).*P.ln1E1));
            
            [~, bar] = sort([QAP_lik0, QAP_lik1]);
            QAP.yhat(j,ii)=bar(2)-1;
            QAP.correct(j,ii)=(QAP.yhat(j)==class_labels(tst(j)));
        end
        
    end
end

if LAP.do
    LAP.Lhat = 1-mean(LAP.correct);
    LAP.Lvar = var(LAP.correct);
end

if QAP.do
    for ii=1:alg.QAP_max_iters
        corrects = QAP.correct(:,ii);
        keeper   = ~isnan(corrects);
        corrects = corrects(keeper);
        
        QAP.Lhat(ii)= 1-mean(corrects);
        QAP.Lvar(ii)= var(corrects);
        QAP.Lstd(ii)= std(corrects);
        QAP.num(ii) = length(corrects);
        QAP.obj0_avg(ii) = mean(QAP.obj0(keeper,ii));
        QAP.obj1_avg(ii) = mean(QAP.obj1(keeper,ii));        
        QAP.obj0_var(ii) = var(QAP.obj0(keeper,ii));
        QAP.obj1_var(ii) = var(QAP.obj1(keeper,ii));        
    end
end
