%% simulate kidney-egg problem
clear; clc

n = 10;     % # of vertices
p = 0.5;    % prob of connection for kidney
q0 = 0.25;   % prob of connection for egg
q1 = 0.75;   % prob of connection for egg
s = 102;    % # of samples

num_signal_vertices  = 3;
E0=p*ones(n);
E0(1:num_signal_vertices,1:num_signal_vertices)=q0;

E1=p*ones(n);
E1(1:num_signal_vertices,1:num_signal_vertices)=q1;

A0      = repmat(E0,[1 1 s/2]) > rand(n,n,s/2);         % class 0 training samples
A1      = repmat(E1,[1 1 s/2]) > rand(n,n,s/2);         % class 1

siz0=size(A0);                                          % size of array storing class 0 adjacency matrices
siz1=size(A1);                                          % size of array storing class 1 adjacency matrices
adjacency_matrices(:,:,1:siz0(3))=A0;                   % create adjacency_matrices to store all sampled graphs
adjacency_matrices(:,:,siz0(3)+1:siz0(3)+siz1(3))=A1;   % add class 1 samples
class_labels=[zeros(1,siz0(3)) ones(1,siz1(3))];        % vector of class labels

%% save stuff
params.n = n;
params.p = p;
params.q0 = q0; 
params.q1 = q1;
params.s = s;

params.num_signal_edges = num_signal_vertices^2;                   % # of edges in signal subgraph
params.num_signal_vertices = sqrt(num_signal_edges);                                                 % # of vertices containing signal
params.signal_subgraph=abs(E0-E1)>1e-10;                    % signal subgraph
params.signal_subgraph_ind = find(params.signal_subgraph);                      % find indices of signal subgraph

params.E0=E0;                    
params.E1=E1;                    

alg.datadir = '~/Research/necog/data/sims/';
alg.figdir  = '~/Research/necog/figs/sims/';
alg.fname   = 'cep_unlabeled';
alg.save = 1;

save([alg.datadir alg.fname],'adjacency_matrices','class_labels','params','alg')

%% setup algorithmic parameters

alg.signal_subgraph_ind = params.signal_subgraph_ind;               % use true signal subgraph
alg.nb_ind              = 1:params.n^2;                             % use naive bayes classifier
alg.num_inc_edges       = params.num_signal_edges;                  % use incoherent classifier with num_signal_vertices^2 edges
alg.num_coh_vertices    = params.num_signal_vertices;               % use coherent classifier with num_signal_vertices^2 edges
alg.num_signal_edges    = params.num_signal_edges;                  % # of signal edges
alg.max_fw_iters        = 30;

%% test using in-sample training data
constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code
[Lhatin Lvarin ind Pin yhatin] = graph_classify_ind_edge(adjacency_matrices,constants,alg); % compute in-sample classification accuracy
disp(Lhatin)

%% the results from this should be the performance given the same training data as the unlabeled classifier gets

Atrn=adjacency_matrices(:,:,[1 constants.s0+1]);
tst_ind=1:constants.s; tst_ind(constants.s0)=[]; tst_ind(1)=[];
Atst=adjacency_matrices(:,:,tst_ind);
ytrn=[0 1];
ytst=constants.ys(tst_ind);

Gtrn=get_constants(Atrn,ytrn);
Gtst=get_constants(Atst,ytst);

Lhatin2 = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst)

%% same as above, but permute data and don't tell anybody (shhhhh)
% make data unlabeled
As=adjacency_matrices;
for i=1:constants.s
    q=randperm(constants.n);
    A=As(:,:,i);
    As(:,:,i)=A(q,q);
end


Atrn=As(:,:,[1 constants.s0+1]);
tst_ind=1:constants.s; tst_ind(constants.s0)=[]; tst_ind(1)=[];
Atst=As(:,:,tst_ind);
ytrn=[0 1];
ytst=constants.ys(tst_ind);

Gtrn=get_constants(Atrn,ytrn);
Gtst=get_constants(Atst,ytst);

Lhatin3 = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst)

%% test using hold-out training data
[Lhats inds num_iters] = get_Lhat_unlabeled_graphs(adjacency_matrices,class_labels,alg);

%% make plots

est_params  = get_params(adjacency_matrices,constants);         % estimate parameters from data

plot_params(est_params,alg,params)                              % plot params and estimated params
plot_recovered_subspaces(constants,est_params,alg)              % plot recovered subspaces 

%%

figure(3), clf, hold all
plot(max(num_iters'),Lhats.tru)
plot(max(num_iters'),Lhats.nb)
plot(max(num_iters'),Lhats.inc)
plot(max(num_iters'),Lhats.coh)
axis([0 alg.max_fw_iters 0 .6])
legend('tru','nb','inc','coh','location','best')
