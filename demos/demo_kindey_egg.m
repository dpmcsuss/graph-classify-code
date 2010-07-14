clear; clc

alg.datadir = '~/Research/necog/data/sims/';
alg.figdir  = '~/Research/necog/figs/sims/';
alg.fname   = 'kidney_egg';
alg.save = 1;

load([alg.datadir alg.fname])                                      % loads adjacency_matrices, class_labels, algs, params, and misc.

alg.ind_edge           = true;
alg.signal_subgraph_ind= params.signal_subgraph_ind;               % use true signal subgraph
alg.nb_ind             = 1:params.n^2;                             % use naive bayes classifier
alg.num_inc_edges      = params.num_signal_edges;                  % use incoherent classifier with num_signal_vertices^2 edges
alg.num_coh_vertices   = params.num_signal_vertices;               % use coherent classifier with num_signal_vertices^2 edges
alg.num_signal_edges   = params.num_signal_edges;                  % # of signal edges

alg.knn             = true;
alg.knn_vanilla     = true;
alg.knn_lmnn        = false;
alg.knn_mmlmnn      = false;

alg.num_splits      = 2;
alg.num_repeats     = 3;

%% test using in-sample training data
constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code
[Lhatin Lvarin ind Pin yhatin] = graph_classify_ind_edge(adjacency_matrices,constants,alg); % compute in-sample classification accuracy
disp(Lhatin)

%% test using hold-out training data
[Lhats alg inds] = get_Lhat_hold_out(adjacency_matrices,class_labels,alg);

%% make plots

est_params  = get_params(adjacency_matrices,constants);         % estimate parameters from data
plot_params(est_params,alg,params)                              % plot params and estimated params
plot_misclassification(Lhats,inds,constants,alg)                % plot misclassification rates and edge detection rates
plot_recovered_subspaces(constants,est_params,alg)              % plot recovered subspaces 
plot_Lhats(Lhats,alg)                                           % plot misclassification rates