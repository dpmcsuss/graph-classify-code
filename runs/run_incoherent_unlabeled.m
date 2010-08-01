%% simulate kidney-egg problem
clear; clc

n = 10;     % # of vertices
p = 0.5;    % prob of connection for kidney
q0 = 0.05;   % prob of connection for egg
q1 = 0.95;   % prob of connection for egg
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
params.num_signal_vertices = num_signal_vertices;                                                 % # of vertices containing signal
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

%% test using in-sample training data, using vertex labels
constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code
[Lhatin Lvarin ind Pin yhatin] = graph_classify_ind_edge(adjacency_matrices,constants,alg); % compute in-sample classification accuracy
disp(Lhatin)

%% the results from this should be the performance given the same training data as the unlabeled classifier gets

%% make training data 

ATRN=zeros(n);
ind=[1:2:10];
for j=1:n-1
    ind=[ind (j*10+rem(j,2)+1):2:(j+1)*n];
end
% ind=[1:2:10 12:2:20 21:2:30 32:2:40];
ATRN(ind)=1;
ATRN0=ATRN;
ATRN0(1:num_signal_vertices,1:num_signal_vertices)=0;

ATRN1=ATRN;
ATRN1(1:num_signal_vertices,1:num_signal_vertices)=1;
Atrn(:,:,1)=ATRN0;
Atrn(:,:,2)=ATRN1;

% Atrn=adjacency_matrices(:,:,[1 constants.s0+1]);
tst_ind=1:constants.s; tst_ind(constants.s0)=[]; tst_ind(1)=[];
Atst=adjacency_matrices(:,:,tst_ind);
ytrn=[0 1];
ytst=constants.ys(tst_ind);

Gtrn=get_constants(Atrn,ytrn);
Gtst=get_constants(Atst,ytst);


Lhatin2 = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst)

%% get classification performance when *not* using vertex labels, and using
%% only a single training datum (should be about 0.5)

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



%% approximately solve isomorphism problem

k=0;
max_iter=30;
for j=tst_ind
    k=k+1;
    if j < constants.s0, B=Atrn(:,:,1); else B=Atrn(:,:,2); end
    A=adjacency_matrices(:,:,j);
    [f,myp,x,iter,fs,myps{k}]=sfw(A,-B,max_iter);
end


%%

Atst=zeros(n,n,length(tst_ind));
for j=1:max_iter
    k=0;
    for l=tst_ind
        k=k+1;
        len=length(myps{k});
        if j>len, jj=len; else jj=j; end
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
for j=1:max_iter
    xxx.nb(j)=Lhats{j}.nb;
    xxx.tru(j)=Lhats{j}.tru;
    xxx.inc(j)=Lhats{j}.inc;
    xxx.coh(j)=Lhats{j}.coh;
end


plot(xxx.tru)
plot(xxx.nb)
plot(xxx.inc)
plot(xxx.coh)

axis([0 alg.max_fw_iters 0 .6])
legend('tru','nb','inc','coh','location','best')
