function [Lhats inds num_iters] = classify_unlabeled_graphs(adjacency_matrices,class_labels,alg,params)

constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code

ytrn = [0 1];

P       = params;
P.lnE0  = log(P.E0);
P.ln1E0 = log(1-P.E0);
P.lnE1  = log(P.E1);
P.ln1E1 = log(1-P.E1);

yhat    = zeros(constants.s/2,1);
correct = zeros(constants.s/2,1);

%% do LAP

if strcmp(alg.names(1),'LAP')
    
    trn0=randperm(constants.s0);
    trn1=randperm(constants.s1)+constants.s0;
    tst =randperm(constants.s);
    
    for j=1:constants.s/2
        
        Atrn(:,:,1)=adjacency_matrices(:,:,trn0(j));
        Atrn(:,:,2)=adjacency_matrices(:,:,trn1(j));
        Atst=adjacency_matrices(:,:,tst(j));
        ytst=class_labels(tst(j));
        
        B=Atst;
        [myp0,f]=munkres(-Atrn(:,:,1)*B');     %Note that sum(sum(-A(myp,:).*B))==f
        %         Atst1=B(myp1,:);
        
        worked0(j)=norm(Atrn(:,:,1)-B(myp0,:)) < norm(Atrn(:,:,1)-B);
        
        B=Atst;
        [myp1,f]=munkres(-Atrn(:,:,2)*B');     %Note that sum(sum(-A(myp,:).*B))==f
        %         Atst2=B(myp2,:);

        worked1(1)=norm(Atrn(:,:,2)-B(myp1,:)) < norm(Atrn(:,:,2)-B);

        lik0=sum(sum(Atst(myp0,:).*P.lnE0+(1-Atst(myp0,:)).*P.ln1E0));
        lik1=sum(sum(Atst(myp1,:).*P.lnE1+(1-Atst(myp1,:)).*P.ln1E1));
        
        [~, bar] = sort([lik0, lik1]);
        yhat(j)=bar(2)-1;
        correct(j)=(yhat(j)==class_labels(tst(j)));
                
    end
   
    Lhats.LAP = 1-mean(correct);
 
end

%% do QAP

Lhats.QAP    = zeros(alg.QAP_max_iters,1);

Atst=zeros(constants.n,constants.n,constants.s-2);
num_iters=zeros(1,alg.QAP_max_iters);

inds{alg.QAP_max_iters} = [];


k=0;
for j=tst_ind
    k=k+1;
    if j < constants.s0,
        B=Atrn(:,:,1);
    else
        B=Atrn(:,:,2);
    end
    A=adjacency_matrices(:,:,j);
    %FW approximate solution of QAP
    [f,myp,x,iter,fs,myps]=sfw(A,-B,alg.QAP_max_iters,A)
    %         [f,myp,x,iter]=sfw(A,-B,ii); % Note that sum(sum(-A.*B(myp,myp)))==f
    Atst(myp,myp,k)=A;
    num_iters(ii,j)=iter;
    % myp=1:constants.n;
    
    Gtrn = get_constants(Atrn,ytrn);
    Gtst = get_constants(Atst,ytst);
    
    LQAP = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst,params);
    
    Lhats.QAP(ii)   = LQAP.nb;
    
    
end


if alg.save, save([alg.datadir alg.fname '_results'],'Lhats','inds'); end
