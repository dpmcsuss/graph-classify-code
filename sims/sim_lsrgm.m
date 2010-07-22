%% simulate kidney-egg problem
clear; clc

n = 10;
m = 3;
s = 250;
p = 0.2;
b = log(p/(1-p));
a = -b;

x0 = zeros(n,1);
x1 = zeros(n,1);
x1(1:m) = 1;

d0 = x0*x0';
d1 = x1*x1';

% this function makes p=.5 in egg, and p=.1 everywhere else
eta0 = 1./(1+exp(-(a*d0+b)));
eta1 = 1./(1+exp(-(a*d1+b)));

A0      = repmat(eta0,[1 1 s/2]) > rand(n,n,s/2);         % class 0 training samples
A1      = repmat(eta1,[1 1 s/2]) > rand(n,n,s/2);         % class 1

siz0=size(A0);                                          % size of array storing class 0 adjacency matrices
siz1=size(A1);                                          % size of array storing class 1 adjacency matrices
adjacency_matrices(:,:,1:siz0(3))=A0;                   % create adjacency_matrices to store all sampled graphs
adjacency_matrices(:,:,siz0(3)+1:siz0(3)+siz1(3))=A1;   % add class 1 samples
class_labels=[zeros(1,siz0(3)) ones(1,siz1(3))];        % vector of class labels

alg.datadir = '~/Research/necog/data/sims/';                % name of directory in which to store data upon saving stuff
alg.figdir  = '~/Research/necog/figs/sims/';                % name of directory in which to store figures upon saving stuff
alg.fname   = 'lsrgm';                                 % name of files to save (to be appended with more details)

params.n = n;
params.m = m;
params.p = p;

params.E0 = eta0;
params.E1 = eta1;

params.signal_subgraph_ind=find(eta0-eta1);
params.nb_ind               = 1:n^2;
params.num_inc_edges        = m^2;
params.num_coh_vertices     = m;
params.num_signal_edges     = m^2;


save([alg.datadir alg.fname],'adjacency_matrices','class_labels','params','alg')