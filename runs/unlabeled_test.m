% this code generalizes unlabeled_test2, by making class conditional
% differences a block model.

clear; clc

n = 10;     % # of vertices
s = 102;    % # of samples

which_sim = 'block';

switch which_sim
    case 'point_mass'
        E0=zeros(n);
        E0([1 12 23])=1;

        E1=zeros(n);
        E1(1:3,1:3)=1;
    case 'diag_block'
        E0=0.5*ones(n);
        E0([1 12 23])=0.25;

        E1=0.5*ones(n);
        E1(1:3,1:3)=0.25;
    case 'block'
        m=3;
        E0=0.5*ones(n);
        E0(1:m,1:m)=0.25;

        E1=0.5*ones(n);
        E1(1:m,1:m)=0.75;
    case 'dense'
        E0=rand(n);
        E1=rand(n);
end


for ii=1:100

    A0 = repmat(E0,[1 1 s/2]) > rand(n,n,s/2);         % class 0 training samples
    A1 = repmat(E1,[1 1 s/2]) > rand(n,n,s/2);         % class 1

    siz0=size(A0);                                          % size of array storing class 0 adjacency matrices
    siz1=size(A1);                                          % size of array storing class 1 adjacency matrices
    adjacency_matrices(:,:,1:siz0(3))=A0;                   % create adjacency_matrices to store all sampled graphs
    adjacency_matrices(:,:,siz0(3)+1:siz0(3)+siz1(3))=A1;   % add class 1 samples
    class_labels=[zeros(1,siz0(3)) ones(1,siz1(3))];        % vector of class labels
    constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code

    %% save stuff
    params.n = n;
    params.s = s;

    params.E0=E0;
    params.E1=E1;

    params.lnE0  = log(params.E0);
    params.ln1E0 = log(1-params.E0);
    params.lnE1  = log(params.E1);
    params.ln1E1 = log(1-params.E1);

    % log-priors
    params.lnprior0 = log(constants.s0/constants.s);
    params.lnprior1 = log(constants.s1/constants.s);

    % various measures to compute difference matrix
    params.d_pos = abs(params.E0-params.E1);           % position difference
    params.d_opt = abs(params.E0./sqrt(params.E0.*(1-params.E0)) - params.E1./sqrt(params.E1.*(1-params.E1))); % optimal difference


    alg.datadir = '~/Research/data/sims/unlabeled/';
    alg.figdir  = '~/Research/figs/sims/unlabeled/';
    alg.fname   = which_sim;
    alg.save = 1;

    save([alg.datadir alg.fname],'adjacency_matrices','class_labels','params','alg','constants')

    %% setup algorithmic parameters

    alg.nb_ind              = 1:n^2;            % use independent edge classifier for directed loopy graphs
    % alg.num_inc_edges       = 3^2;              % use incoherent classifier with num_signal_vertices^2 edges
    % alg.num_coh_vertices    = 3;                % use coherent classifier with num_signal_vertices^2 edges
    % alg.signal_subgraph_ind = find(E0~=E1);

    %% test using in-sample training data, using vertex labels
    [Lhatin ind Pin yhatin] = graph_classify_ind_edge(adjacency_matrices,constants,alg); % compute in-sample classification accuracy
    disp(Lhatin)

    %% test using out-of-sample data
    % but use vertex labels

    tst_ind=1:constants.s; tst_ind(constants.s0+1)=[]; tst_ind(1)=[];
    ytrn=[0 1];
    ytst=constants.ys(tst_ind);

    Atrn(:,:,1)=adjacency_matrices(:,:,1);
    Atrn(:,:,2)=adjacency_matrices(:,:,constants.s0+1);

    Atst=adjacency_matrices(:,:,tst_ind);

    Gtrn=get_constants(Atrn,ytrn);
    Gtst=get_constants(Atst,ytst);

    Lhat_labeled = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst)

    %% permute testing data

    As=adjacency_matrices;
    for i=1:constants.s
        q=randperm(constants.n);
        A=As(:,:,i);
        As(:,:,i)=A(q,q);
    end
    Atst=As(:,:,tst_ind);

    % test classification performance when vertex labels are permuted
    % (should be just less than 1/2)
    % (in other words, we didn't try to solve the isomorphism problem first)
    Lhat_permuted = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst)

    %% approximately solve isomorphism problem

    k=0;
    alg.fw_max_iter=30;
    B=Atrn(:,:,1);
    for j=tst_ind
        k=k+1;
        A=As(:,:,j);
        [f,myp,x,iter,fs,myps{k}]=sfw(B,-A,alg.fw_max_iter);
    end

    B=Atrn(:,:,2);
    for j=tst_ind
        k=k+1;
        A=As(:,:,j);
        [f,myp,x,iter,fs,myps{k}]=sfw(B,-A,alg.fw_max_iter);
    end

    %% compute

    Atst0=zeros(n,n,length(tst_ind));
    Atst1=zeros(n,n,length(tst_ind));

    for j=1:alg.fw_max_iter
        k=0;
        for l=tst_ind
            k=k+1;
            A=As(:,:,l);

            len0=length(myps{k});
            if j>len0, j0=len0; else j0=j; end
            Atst0(:,:,k)=A(myps{k}{j0},myps{k}{j0});

            len1=length(myps{k+Gtst.s});
            if j>len1, j1=len1; else j1=j; end
            Atst1(:,:,k)=A(myps{k+Gtst.s}{j1},myps{k+Gtst.s}{j1});

        end
        Lhats{j} = graph_classify_unlabeled_ind_edge(Atrn,Gtrn,alg,Atst0,Atst1,Gtst,params);
    end


    rates.nb(ii,1)=Lhat_permuted.nb;
    for j=1:alg.fw_max_iter
        rates.nb(ii,j+1)=Lhats{j}.nb;
    end
    rates.nb(ii,j+2)=Lhat_labeled.nb;

end



%% make plots

est_params  = get_params(adjacency_matrices,constants);         % estimate parameters from data

plot_params(est_params,alg,params)                              % plot params and estimated params
% plot_recovered_subspaces(constants,est_params,alg)              % plot recovered subspaces
xxx = plot_unlabeled_rates(Lhat_permuted,Lhats,Lhat_labeled,alg);