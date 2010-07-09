function [Lhats inds num_iters] = wrapper_hold_out_unbalanced_unlabeled_training_data(adjacency_matrices,class_labels,alg)

constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code
inds{alg.max_fw_iters} = [];

if isfield(alg,'signal_subgraph_ind'), Lhats.tru = zeros(alg.max_fw_iters,1); end
Lhats.nb    = zeros(alg.max_fw_iters,1);
Lhats.inc   = zeros(alg.max_fw_iters,1);
Lhats.coh   = zeros(alg.max_fw_iters,1);

Atrn(:,:,1) = adjacency_matrices(:,:,1);
Atrn(:,:,2) = adjacency_matrices(:,:,constants.s0+1);
Atst=zeros(constants.n,constants.n,constants.s-2);
num_iters=zeros(1,alg.max_fw_iters);

for ii=1:alg.max_fw_iters

    tst_ind = 1:constants.s; tst_ind(constants.s0)=[]; tst_ind(1)=[];
    k=0;
    for j=tst_ind
        k=k+1;
        if j < constants.s0, B=Atrn(:,:,1); else B=Atrn(:,:,2); end
        A=adjacency_matrices(:,:,j);
        [f,myp,x,iter]=sfw(A,-B,ii);
        num_iters(ii)=iter;
        Atst(:,:,k)=A(myp,myp);
    end

    ytrn = [0 1];
    ytst = class_labels(tst_ind);

    Gtrn = get_constants(Atrn,ytrn);
    Gtst = get_constants(Atst,ytst);

    [Lhat Lvar ind] = graph_classify_ie(Atrn,Gtrn,alg,Atst,Gtst);

    if isfield(alg,'signal_subgraph_ind'),
        Lhats.tru(ii)  = Lhat.tru;
        if j==1, inds{i} = ind; end % just store indices from first run, ignoring nb (as that is all), needed for plotting
    end
    Lhats.nb(ii)   = Lhat.nb;
    Lhats.inc(ii)  = Lhat.inc;
    Lhats.coh(ii)  = Lhat.coh;
end

if alg.save, save([alg.datadir alg.fname '_results'],'Lhats','inds'); end