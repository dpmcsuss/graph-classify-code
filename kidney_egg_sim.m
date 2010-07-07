%% simulate kidney-egg problem
clear; clc

params.n = 50;  n = params.n;
params.p = 0.4; p = params.p;
params.q = 0.1; q = params.q;
params.s = 220; s = params.s;

params.num_signal_vertices = 5;
params.signal_vertices = 1:params.num_signal_vertices;
params.signal_subgraph = zeros(n);
params.signal_subgraph(params.signal_vertices,params.signal_vertices) = 1;
params.signal_subgraph_ind = find(params.signal_subgraph);
params.num_signal_edges = length(params.signal_subgraph_ind);

params.E0 = p*ones(n); E0=params.E0;
params.E1 = p*ones(n); 
params.E1(params.signal_subgraph_ind) = q; E1=params.E1;

A0      = repmat(E0,[1 1 s/2]) > rand(n,n,s/2);     % class 0 training samples
A1      = repmat(E1,[1 1 s/2]) > rand(n,n,s/2);     % class 1

siz0=size(A0);
siz1=size(A1);
adjacency_matrices(:,:,1:siz0(3))=A0;
adjacency_matrices(:,:,siz0(3)+1:siz0(3)+siz1(3))=A1;
class_labels=[zeros(1,siz0(3)) ones(1,siz1(3))];

datadir = '~/Research/necog/data/sims/';
figdir  = '~/Research/necog/figs/sims/';
fname   = 'kidney_egg';

save([datadir fname],'adjacency_matrices','class_labels','params','datadir','figdir','fname')