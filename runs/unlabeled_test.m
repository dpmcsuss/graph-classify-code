% this run generates F_{GY}, where F_0 and F_1 are both point masses.
% this code is for debugging purposes.

clear; clc

n = 10;     % # of vertices
s = 102;    % # of samples

E0=zeros(n);
E0([1 12 23])=1;

E1=zeros(n);
E1(1:3,1:3)=1;

A0 = repmat(E0,[1 1 s/2]) > rand(n,n,s/2);         % class 0 training samples
A1 = repmat(E1,[1 1 s/2]) > rand(n,n,s/2);         % class 1

siz0=size(A0);                                          % size of array storing class 0 adjacency matrices
siz1=size(A1);                                          % size of array storing class 1 adjacency matrices
adjacency_matrices(:,:,1:siz0(3))=A0;                   % create adjacency_matrices to store all sampled graphs
adjacency_matrices(:,:,siz0(3)+1:siz0(3)+siz1(3))=A1;   % add class 1 samples
class_labels=[zeros(1,siz0(3)) ones(1,siz1(3))];        % vector of class labels

%% save stuff
params.n = n;
params.s = s;

params.E0=E0;
params.E1=E1;

alg.datadir = '~/Research/data/sims/';
alg.figdir  = '~/Research/figs/sims/';
alg.fname   = 'unlabeled_test';
alg.save = 1;

save([alg.datadir alg.fname],'adjacency_matrices','class_labels','params','alg')

%% setup algorithmic parameters

alg.nb_ind              = 1:n^2;            % use independent edge classifier for directed loopy graphs
alg.num_inc_edges       = 3^2;              % use incoherent classifier with num_signal_vertices^2 edges
alg.num_coh_vertices    = 3;                % use coherent classifier with num_signal_vertices^2 edges
alg.signal_subgraph_ind = find(E0~=E1);

%% test using in-sample training data, using vertex labels
constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code
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

Lhatin2 = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst)

%% permute testing data

As=adjacency_matrices;
for i=1:constants.s
    q=randperm(constants.n);
    A=As(:,:,i);
    As(:,:,i)=A(q,q);
end
Atst=As(:,:,tst_ind);

% test classification performance when data is permuted (should be just
% less than 1/2)
Lhatin3 = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst)

%% approximately solve isomorphism problem

k=0;
alg.fw_max_iter=30;
for j=tst_ind
    k=k+1;
    if j <= constants.s0, B=Atrn(:,:,1); else B=Atrn(:,:,2); end
    A=As(:,:,j);
    [f,myp,x,iter,fs,myps{k}]=sfw(B,-A,alg.fw_max_iter);
end


%%

Atst=zeros(n,n,length(tst_ind));
for j=1:alg.fw_max_iter
    k=0;
    for l=tst_ind
        k=k+1;
        len=length(myps{k});
        if j>len, jj=len; else jj=j; end
        A=As(:,:,l);
        Atst(:,:,k)=A(myps{k}{jj},myps{k}{jj});
    end
    Lhats{j} = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst);
end


%% make plots

est_params  = get_params(adjacency_matrices,constants);         % estimate parameters from data

plot_params(est_params,alg,params)                              % plot params and estimated params
plot_recovered_subspaces(constants,est_params,alg)              % plot recovered subspaces

%%

figure(3), clf, hold all
for j=1:alg.fw_max_iter
    xxx.nb(j)=Lhats{j}.nb;
    xxx.tru(j)=Lhats{j}.tru;
    xxx.inc(j)=Lhats{j}.inc;
    xxx.coh(j)=Lhats{j}.coh;
end


plot(xxx.tru)
plot(xxx.nb)
plot(xxx.inc)
plot(xxx.coh)

axis([0 alg.fw_max_iter 0 .6])
legend('tru','nb','inc','coh','location','best')
