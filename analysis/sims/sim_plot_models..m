clear; clc

alg.datadir = '~/Research/data/graph_sims/unlabeled/';
alg.figdir  = '~/Research/figs/graph_sims/unlabeled/';
alg.fname   = 'models';
alg.save    = 1;

n   = 10;   % # of vertices
p   = 0.5;  % prob of connection for kidney
q0  = 0.25; % prob of connection for egg
q1  = 0.75; % prob of connection for egg
egg = 1:3;  % vertices in egg    

E0=p*ones(n);   % params in class 0
E0(egg,egg)=q0; % egg params in class 0

E1=p*ones(n);   % params in class 1
E1(egg,egg)=q1; % egg params in class 1

params1.E0=E0; params1.E1=E1;


n   = 10;   % # of vertices

E0=rand(n);   % params in class 0
E1=rand(n);   % params in class 1


params2.E0=E0; params2.E1=E1;


plot_models(params1,params2,alg)