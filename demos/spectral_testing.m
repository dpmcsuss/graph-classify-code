clear; clc

alg.datadir = '~/Research/necog/data/sims/';
alg.figdir  = '~/Research/necog/figs/sims/';
alg.fname   = 'kidney_egg';
alg.save    = 0;

load([alg.datadir alg.fname])                                      % loads adjacency_matrices, class_labels, algs, params, and misc.

alg.ind_edge        = false;
alg.knn             = false;
alg.spectral        = true;


%% test using in-sample training data
constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code

Lhat = graph_classify_spectral(adjacency_matrices,constants,alg);
disp(Lhat)

%% test using hold-out training data

alg.num_splits      = 3;
alg.num_repeats     = 2;

[Lhats alg inds] = get_Lhat_hold_out(adjacency_matrices,class_labels,alg);

%% make plots

est_params  = get_params(adjacency_matrices,constants);         % estimate parameters from data
plot_params(est_params,alg,params)                              % plot params and estimated params
plot_recovered_subspaces(constants,est_params,alg)              % plot recovered subspaces 
plot_edge_identification_rates(inds,constants,alg)              % plot misclassification rates and edge detection rates
plot_Lhats(Lhats,alg)                                           % plot misclassification rates