%% simulate kidney-egg problem
clear; clc

params.n = 50;  n = params.n; % # of vertices
params.p = 0.4; p = params.p; % prob of connection for kidney
params.q = 0.1; q = params.q; % prob of connection for egg
params.s = 220; s = params.s; % # of samples

params.num_signal_vertices = 5;                                                 % # of vertices containing signal
params.signal_vertices = 1:params.num_signal_vertices;                          % index of vertices containing signal
params.signal_subgraph = zeros(n);                                              % pre-allocate memory for signal subgraph
params.signal_subgraph(params.signal_vertices,params.signal_vertices) = 1;      % create signal subgraph
params.signal_subgraph_ind = find(params.signal_subgraph);                      % find indices of signal subgraph
params.num_signal_edges = length(params.signal_subgraph_ind);                   % # of edges in signal subgraph

params.E0 = p*ones(n); E0=params.E0;                                            % true class 0 probabilities
params.E1 = p*ones(n);                                                          % true class 1 probabilities for kidney
params.E1(params.signal_subgraph_ind) = q; E1=params.E1;                        % true class 1 probabilities for egg

A0      = repmat(E0,[1 1 s/2]) > rand(n,n,s/2);         % class 0 training samples
A1      = repmat(E1,[1 1 s/2]) > rand(n,n,s/2);         % class 1

siz0=size(A0);                                          % size of array storing class 0 adjacency matrices
siz1=size(A1);                                          % size of array storing class 1 adjacency matrices
adjacency_matrices(:,:,1:siz0(3))=A0;                   % create adjacency_matrices to store all sampled graphs
adjacency_matrices(:,:,siz0(3)+1:siz0(3)+siz1(3))=A1;   % add class 1 samples
class_labels=[zeros(1,siz0(3)) ones(1,siz1(3))];        % vector of class labels

% % make vertices unlabeled
% for i=1:params.s
%     q=randperm(n);
%     A=adjacency_matrices(:,:,i);
%     adjacency_matrices(:,:,i)=A(q,q);
% end

datadir = '~/Research/necog/data/sims/';                % name of directory in which to store data upon saving stuff
figdir  = '~/Research/necog/figs/sims/';                % name of directory in which to store figures upon saving stuff
fname   = 'kidney_egg';                                 % name of files to save (to be appended with more details)

save([datadir fname],'adjacency_matrices','class_labels','params','datadir','figdir','fname')