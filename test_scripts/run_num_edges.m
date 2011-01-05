clear; clc

alg.datadir = '~/Research/data/sims/labeled/';
alg.figdir  = '~/Research/figs/sims/labeled/';
alg.fname   = 'twenty_edges';
alg.save    = 0;

n = 20;
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


%% test using in-sample training data

% alg.bayes_plugin=true;
alg.ind_edge=true;
alg.nb_ind=constants.inds;
% alg.mle=true;
alg.knn=true;
alg.knn_vanilla=true;
% alg.knn_lmnn=true;

% [Lhat_bayes Pin] = graph_classify_bayes_plugin(adjacency_matrices,constants,alg);
% disp(Lhat_bayes)

%% test using hold-out training data

alg.num_splits      = 4;
alg.num_repeats     = 2;

[Lhats alg inds] = get_Lhat_hold_out(adjacency_matrices,class_labels,alg);

%% make plots

% est_params  = get_params(adjacency_matrices,constants);         % estimate parameters from data
% plot_params(est_params,alg,params)                              % plot params and estimated params
% plot_recovered_subspaces(constants,est_params,alg)              % plot recovered subspaces 
% plot_edge_identification_rates(inds,constants,alg)              % plot misclassification rates and edge detection rates
plot_Lhats(Lhats,alg)                                           % plot misclassification rates