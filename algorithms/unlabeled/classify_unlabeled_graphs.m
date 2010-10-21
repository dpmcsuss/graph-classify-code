function [Lhats inds num_iters] = classify_unlabeled_graphs(adjacency_matrices,class_labels,alg,params)

constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code

% get indices for test data
tst_ind = [2:constants.s0, constants.s0+2:constants.s];

ytrn = [0 1];
ytst = class_labels(tst_ind);

Atrn(:,:,1) = adjacency_matrices(:,:,1);
Atrn(:,:,2) = adjacency_matrices(:,:,constants.s0+1);
Atst        = zeros(constants.n,constants.n,constants.s-2); % pre-allocate memory for Atst


%% do LAP

if strcmp(alg.names(1),'LAP')
    k=0;
    for j=tst_ind
        k=k+1;
        if j < constants.s0,
            B=Atrn(:,:,1);
        else
            B=Atrn(:,:,2);
        end
        A=adjacency_matrices(:,:,j);
        [myp,f]=munkres(-B*A');     %Note that sum(sum(-A(myp,:).*B))==f
        Atst(:,:,k)=A(myp,:);
    end
    
    Gtrn = get_constants(Atrn,ytrn);
    Gtst = get_constants(Atst,ytst);
    
    LLAP = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst,params);
    
    Lhats.LAP = LLAP.nb;
end

%% do QAP

Lhats.QAP    = zeros(alg.QAP_max_iters,1);

Atst=zeros(constants.n,constants.n,constants.s-2);
num_iters=zeros(1,alg.QAP_max_iters);

inds{alg.QAP_max_iters} = [];

for ii=1:alg.QAP_max_iters
    
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
        [f,myp,x,iter]=sfw(A,-B,ii); % Note that sum(sum(-A.*B(myp,myp)))==f
        Atst(myp,myp,k)=A;
        num_iters(ii,j)=iter;
        % myp=1:constants.n;
    end    
    
    Gtrn = get_constants(Atrn,ytrn);
    Gtst = get_constants(Atst,ytst);
    
    LQAP = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst,params);
    
    Lhats.QAP(ii)   = LQAP.nb;
end

if alg.save, save([alg.datadir alg.fname '_results'],'Lhats','inds'); end
