% generate simulations for testing purposes
% for each simulation, we have the following structures:
% 
% name:     name of simulation
% model:    model (n,s, and parameters)
% As:       adjacency matrices
% targs:    target variables
% params:   parameters derived from model
% Lstar:    bayes optimal classification performance


clear, clc

%% kidney egg1
i=1;

M.n = 10;
M.m = 3;
M.s = 1000;

p = 0.5;
q = 0.25;
q0 = p-q;
q1 = p+q;

M.E0 = p.*ones(M.n); M.E0(1:M.m,1:M.m) = q0;
M.E1 = p.*ones(M.n); M.E1(1:M.m,1:M.m) = q1;

sim{i} = gen_ind_edge_data(M);
sim{i}.model = M;
sim{i}.name='kidney_egg1';
sim{i}.Lstar = 1-binocdf(floor(M.m^2/2),M.m^2,p-q);


%% kidney egg2
i=2; clear M

M.n = 10;
M.m = 3;
M.s = 1000;

p = 0.5;
q = 0.05;
q0 = p-q;
q1 = p+q;

M.E0 = p.*ones(M.n); M.E0(1:M.m,1:M.m) = q0;
M.E1 = p.*ones(M.n); M.E1(1:M.m,1:M.m) = q1;

sim{i} = gen_ind_edge_data(M);
sim{i}.model = M;
sim{i}.name='kidney_egg2';
sim{i}.Lstar = 1-binocdf(floor(M.m^2/2),M.m^2,p-q);


%% kidney egg3
i=3; clear M

M.n = 100;
M.m = 5;
M.s = 1000;

p = 0.5;
q = 0.1;
q0 = p-q;
q1 = p+q;

M.E0 = p.*ones(M.n); M.E0(1:M.m,1:M.m) = q0;
M.E1 = p.*ones(M.n); M.E1(1:M.m,1:M.m) = q1;

sim{i} = gen_ind_edge_data(M);
sim{i}.model = M;
sim{i}.name='kidney_egg3';
sim{i}.Lstar = 1-binocdf(floor(M.m^2/2),M.m^2,p-q);


%% hetero kidney egg
i=4; clear M

M.n = 10;
M.m = 3;
M.s = 1000;

M.E0 = rand(M.n).*ones(M.n);
M.E1 = M.E0; 
M.E1(1:M.m,1:M.m) = rand(M.m);

sim{i} = gen_ind_edge_data(M);
sim{i}.model = M;
sim{i}.name='hetero_kidney_egg';

%% kidney egg2
i=5; clear M

M.n = 10;
M.m = 3;
M.s = 100;

p = 0.5;
q = 0.05;
q0 = p-q;
q1 = p+q;

M.E0 = p.*ones(M.n); M.E0(1:M.m,1:M.m) = q0;
M.E1 = p.*ones(M.n); M.E1(1:M.m,1:M.m) = q1;

sim{i} = gen_ind_edge_data(M);
sim{i}.model = M;
sim{i}.name='kidney_egg4';
sim{i}.Lstar = 1-binocdf(floor(M.m^2/2),M.m^2,p-q);


%% save
save('~/Research/data/sims/test_data','sim');
