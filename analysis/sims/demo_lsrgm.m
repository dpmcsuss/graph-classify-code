clear; clc

alg.datadir = '~/Research/necog/data/sims/';
alg.figdir  = '~/Research/necog/figs/sims/';
alg.fname   = 'lsrgm';

load([alg.datadir alg.fname])                                      % loads adjacency_matrices, class_labels, algs, params, and misc.

alg             = catstruct(alg,params);
alg.save        = 0;

alg.ind_edge    = true;
alg.knn         = false;
alg.knn_vanilla = true;
alg.knn_lmnn    = true;
% alg.knn_mmlmnn      = false;

alg.num_splits      = 10;
alg.num_repeats     = 10;

%% test using in-sample training data
constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code
[Lhatin Lvarin ind Pin yhatin] = graph_classify_ind_edge(adjacency_matrices,constants,alg); % compute in-sample classification accuracy
disp(Lhatin)

%% test using hold-out training data
[Lhats alg inds] = get_Lhat_hold_out(adjacency_matrices,class_labels,alg);

%% make plots

est_params  = get_params(adjacency_matrices,constants);         % estimate parameters from data
plot_params(est_params,alg,params)                              % plot params and estimated params
plot_recovered_subspaces(constants,est_params,alg)              % plot recovered subspaces 
plot_edge_identification_rates(inds,constants,alg)              % plot misclassification rates and edge detection rates
plot_Lhats(Lhats,alg)                                           % plot misclassification rates