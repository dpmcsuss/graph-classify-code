%% simulate kidney-egg problem
clear; clc

alg.datadir = '~/Research/data/graph_sims/unlabeled/';
alg.figdir  = '~/Research/figs/graph_sims/unlabeled/';
alg.fname   = 'hetero';
alg.save    = 1;

n   = 10;   % # of vertices
S   = 102;  % # of samples

E0=rand(n);   % params in class 0
E1=rand(n);   % params in class 1

A0 = repmat(E0,[1 1 S/2]) > rand(n,n,S/2); % class 0 samples
A1 = repmat(E1,[1 1 S/2]) > rand(n,n,S/2); % class 1 samples

siz0=size(A0);                                          % size of array storing class 0 adjacency matrices
siz1=size(A1);                                          % size of array storing class 1 adjacency matrices
adjacency_matrices(:,:,1:siz0(3))=A0;                   % create adjacency_matrices to store all sampled graphs
adjacency_matrices(:,:,siz0(3)+1:siz0(3)+siz1(3))=A1;   % add class 1 samples
class_labels=[zeros(1,siz0(3)) ones(1,siz1(3))];        % vector of class labels

params.n=n; params.S=S; params.E0=E0; params.E1=E1;

save([alg.datadir alg.fname],'adjacency_matrices','class_labels','params','alg')

%% knobs to tweak in classifiers

alg.names = [{'LAP'}; {'QAP'}];
% alg.names = [{'LAP'}];

alg.QAP_max_iters = 10;

% provide indices to classifiers
alg.nb_ind = 1:n^2;                 % use all edges
% signal_subgraph=zeros(n);
% signal_subgraph(egg,egg)=1;
% alg.signal_subgraph_ind = find(signal_subgraph); % find indices of signal subgraph

%% performance when not trying to solve GIP


constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code
tst_ind = [2:constants.s0, constants.s0+2:constants.s];

adj_mat=0*adjacency_matrices;
% make data unlabeled
for i=constants.s
    q=randperm(constants.n);
    A=adjacency_matrices(:,:,i);
    adj_mat(:,:,i)=A(q,q);
end

Atrn=adj_mat(:,:,[1 constants.s0+1]);
Atst=adj_mat(:,:,tst_ind);

ytrn=[0 1];
ytst=constants.ys(tst_ind);

Gtrn=get_constants(Atrn,ytrn);
Gtst=get_constants(Atst,ytst);

Lhat = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst,params);
Lhats.rand = Lhat.nb;


%% performance using true parameters and labels

constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code
Atrn=adjacency_matrices(:,:,[1 constants.s0+1]);
tst_ind = [2:constants.s0, constants.s0+2:constants.s];
Atst=adjacency_matrices(:,:,tst_ind);
ytrn=[0 1];
ytst=constants.ys(tst_ind);

Gtrn=get_constants(Atrn,ytrn);
Gtst=get_constants(Atst,ytst);

Lhat = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst,params);
Lhats.star = Lhat.nb;

%% test using hold-out training data
[Lhat inds num_iters] = classify_unlabeled_graphs(adjacency_matrices,class_labels,alg,params);
Lhats.LAP = Lhat.LAP;
Lhats.QAP = Lhat.QAP;

%% make plots


plot_model(params,alg)                              % plot params and estimated params

figure(3), clf, hold all
plot(1:alg.QAP_max_iters,Lhats.rand*ones(alg.QAP_max_iters,1),'g','linewidth',2)
plot(1:alg.QAP_max_iters,Lhats.LAP*ones(alg.QAP_max_iters,1),'k','linewidth',2)
plot(max(num_iters'),Lhats.QAP,'b','linewidth',2)
plot(1:alg.QAP_max_iters,Lhats.star*ones(alg.QAP_max_iters,1),'r','linewidth',2)
legend('rand','LAP','QAP','L^*')
axis([0 alg.QAP_max_iters 0 .5])
% legend('tru','nb','inc','coh','location','best')
