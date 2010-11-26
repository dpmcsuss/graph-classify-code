clear; clc

alg.datadir = '~/Research/data/sims/labeled/';
alg.figdir  = '~/Research/figs/sims/labeled/';
alg.fname   = 'xor';
alg.save    = 0;

n = 10;
d = choose(n,2);
p = 0.5;
E = p*ones(n);
s = 2000;
alg.loopy = 0;
alg.directed = 0;
adjacency_matrices = repmat(E,[1 1 s]) > rand(n,n,s);         % class 0 training samples

class_labels=zeros(1,s);
num_edges=zeros(s,1);
for l=1:s
    adjacency_matrices(:,:,l)=triu(adjacency_matrices(:,:,l),+1);
    num_edges(l)=sum(sum(adjacency_matrices(:,:,l)));
    if num_edges(l)>d/2, class_labels(l)=1; 
    elseif  num_edges(l)<d/2, class_labels(l)=0; 
    else class_labels(l)=double(rand>0.5);
    end
end
sum(class_labels)
constants = get_constants(adjacency_matrices,class_labels,alg.directed,alg.loopy);     % get constants to ease classification code
params.n=n; params.d=d; params.p=p; params.s=s;

save([alg.datadir alg.fname],'adjacency_matrices','class_labels','params','constants')
%% load([alg.datadir alg.fname])                                      % loads adjacency_matrices, class_labels, algs, params, and misc.


alg.ind_edge           = true;
% alg.signal_subgraph_ind= params.signal_subgraph_ind;               % use true signal subgraph
alg.nb_ind             = constants.inds;                             % use naive bayes classifier
alg.num_inc_edges      = n;                  % use incoherent classifier with num_signal_vertices^2 edges
alg.num_coh_vertices   = round(n/10);               % use coherent classifier with num_signal_vertices^2 edges
% alg.num_signal_edges   = params.num_signal_edges;                  % # of signal edges

alg.er              = true;
% alg.bayes_plugin    = true;

alg.knn             = true;
alg.knn_vanilla     = true;
alg.knn_lmnn        = false;
alg.knn_mmlmnn      = false;

alg.num_splits      = 3;
alg.num_repeats     = 2;
alg.min_samples     = 2;
alg.max_samples     = min(constants.s0,constants.s1)-4;

%% test using in-sample training data
[Lhatin Lvarin Pin yhatin] = graph_classify_ER(adjacency_matrices,constants,alg); % compute in-sample classification accuracy
disp(Lhatin)

%% test using hold-out training data
[Lhats alg inds] = get_Lhat_hold_out(adjacency_matrices,class_labels,alg);

%% make plots

est_params  = get_params(adjacency_matrices,constants);         % estimate parameters from data
% plot_params(est_params,alg,params)                              % plot params and estimated params
plot_recovered_subspaces(constants,est_params,alg)              % plot recovered subspaces 
% plot_edge_identification_rates(inds,constants,alg)              % plot misclassification rates and edge detection rates
[mean_hat std_hat] = plot_Lhats(Lhats,alg)                                           % plot misclassification rates