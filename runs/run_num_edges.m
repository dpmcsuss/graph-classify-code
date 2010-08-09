clear; clc

alg.datadir = '~/Research/data/sims/labeled/';
alg.figdir  = '~/Research/figs/sims/labeled/';
alg.fname   = 'num_edges';
alg.save    = 0;

n = 5;
d = choose(n,2);
p = 0.5;
E = p*ones(n);
s = 1000;
loopy = 0;
directed = 0;
adjacency_matrices = repmat(E,[1 1 s]) > rand(n,n,s);         % class 0 training samples


class_labels=zeros(s,1);
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
constants = get_constants(adjacency_matrices,class_labels,directed,loopy);     % get constants to ease classification code


%% test using in-sample training data

Lhat_bayes = graph_classify_bayes(adjacency_matrices,constants,alg);

[Lhatin Lvarin ind Pin yhatin] = graph_classify_ind_edge(adjacency_matrices,constants,alg); % compute in-sample classification accuracy
disp(Lhatin)

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