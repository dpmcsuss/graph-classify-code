% this script loads a number of simulated datasets and then tests various
% classification algorithms to ensure that any modifications to the code
% have not broken anything, and that everything works as it should
clc, clear,
% load('~/Research/data/sims/test_data') %generate_sims
generate_sims


%% test Ltrue < Lhat_in < Lhat_holdout

Lhat_true_in_holdout

%% test incoherent model selection

signal_subgraph_model_selection

%% test algs on impossible data

clear
M.n = 10;
M.m = M.n^2;
M.s = 100000;

M.E0 = rand(M.n);
M.E1 = M.E0;

i=6;
sim{i} = gen_ind_edge_data(M);
sim{i}.model = M;
sim{i}.name='impossible_data';
sim{i}.Lstar = 0.5;

constants = get_constants(sim{i}.As,sim{i}.targs);

alg.naive_bayes=1;
subspace.naive_bayes=1:M.n^2;
subspace.fisher_incoherent=10;
subspace.incoherent=10;
subspace.cor_incoherent=10;
subspace.coherent=3;

alg.num_repeats=1;

[Lhat Lsem alg] = classify_holdout(sim{i}.As,constants,alg,subspace);

fprintf('\n\n**********************************')
fprintf(['\nsome impossible data, Lhat should be about 0.5\n'])
fprintf('**********************************\n\n')


fn = fieldnames(Lhat);
for i=1:length(fn)
    fprintf('Lhat (sem)\n')
    fprintf('%1.3g (%1.3g) \n\n', Lhat.(fn{i}),Lsem.(fn{i}))

    if Lhat.(fn{i})-2*Lsem.(fn{i})<=0.5 && Lhat.(fn{i})+2*Lsem.(fn{i})>=0.5
        disp([(fn{i}) ' gets around 50% when it is impossible'])
    else
        error([(fn{i}) ' performance is unacceptable good!'])
    end

end

%% test "manual subspace"

clear,
M.n = 10;
M.m = M.n^2;
M.s = 100000;

M.E0 = rand(M.n);
M.E1 = M.E0;

i=6;
sim{i} = gen_ind_edge_data(M);
sim{i}.model = M;
sim{i}.name='impossible_data';
sim{i}.Lstar = 0.5;

constants = get_constants(sim{i}.As,sim{i}.targs);

alg.naive_bayes=1;
signal_edges=find(sim{i}.model.E0-sim{i}.model.E1);
manual=1:sim{i}.model.n^2;
manual(signal_edges)=[];
subspace.manual=manual;
alg.num_repeats=1;

[Lhat Lsem alg] = classify_holdout(sim{i}.As,constants,alg,subspace);

fprintf('\n\n****************************************************')
fprintf('\nchecking manual subspace, discarding all signal edges,\n')
fprintf('then classifying, Lhat should be about 0.5\n')
fprintf('********************************************************\n\n')


fn = fieldnames(Lhat);
for i=1:length(fn)
    fprintf('Lhat (sem)\n')
    fprintf('%1.3g (%1.3g) \n\n', Lhat.(fn{i}),Lsem.(fn{i}))

    if Lhat.(fn{i})-2*Lsem.(fn{i})<=0.5 && Lhat.(fn{i})+2*Lsem.(fn{i})>=0.5
        disp([(fn{i}) ' gets around 50% when it is impossible'])
    else
        error([(fn{i}) ' performance is unacceptable good!'])
    end

end

%% test invars doesn't work when it shouldn't

clear;

n=10;
s=10000;
p=0.5;

As=rand(n,n,s)<p;
targs=[zeros(1,s/2) ones(1,s/2)];

x = get_graph_invariants(As);

constants = get_constants(As,targs);

[Atrn Gtrn Atst Gtst inds] = crossvalprep(As,constants,s/4,s/4,s/4,s/4);

clear discrim
discrim.dLDA=1;
% discrim.QDA=1;

params = get_discriminant_params(x,inds,discrim);

[Lhat Lsem] = discriminant_classifiers(x(:,inds.ytst),targs(inds.ytst),params,discrim)


fprintf('\n\n****************************************************')
fprintf('\nchecking discrim analysis should get Lhat=0.5 here\n')
fprintf('********************************************************\n\n')


fn = fieldnames(Lhat);
for i=1:length(fn)
    fprintf('Lhat (sem)\n')
    fprintf('%1.3g (%1.3g) \n\n', Lhat.(fn{i}),Lsem.(fn{i}))

    if Lhat.(fn{i})-2*Lsem.(fn{i})<=0.5 && Lhat.(fn{i})+2*Lsem.(fn{i})>=0.5
        disp([(fn{i}) ' gets around 50% when it is impossible'])
    else
        error([(fn{i}) ' performance is unacceptable good!'])
    end

end


%% test invars does work when choosing class based on invars

clear
n=10;
s=1000;
p=0.5;

As=rand(n,n,s)<p;
targs=nan(1,s);

for i=1:s
    if sum(sum(As(:,:,i)))>n^2/2
        targs(i)=1;
    else
        targs(i)=0;
    end
end

x = get_graph_invariants(As);
constants = get_constants(As,targs);
[Atrn Gtrn Atst Gtst inds] = crossvalprep(As,constants,s/4,s/4,s/4,s/4);

clear discrim
discrim.dLDA=1;

params = get_discriminant_params(x,inds,discrim);

[Lhat Lsem] = discriminant_classifiers(x(:,inds.ytst),targs(inds.ytst),params,discrim)



















