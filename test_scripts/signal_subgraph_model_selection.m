

%% for cross-validation
clear 
i=5;

M.n = 100;
M.m = 3;
M.s = 500;

p = 0.5;
q = 0.25;
q0 = p-q;
q1 = p+q;

M.E0 = p.*ones(M.n); M.E0(1:M.m,1:M.m) = q0;
M.E1 = p.*ones(M.n); M.E1(1:M.m,1:M.m) = q1;

sim{i} = gen_ind_edge_data(M);
sim{i}.model = M;
sim{i}.name='kidney_egg4';
sim{i}.Lstar = 1-binocdf(floor(M.m^2/2),M.m^2,p-q);


alg.num_repeats = 1;        % # of times to repeat using the same # of hold-outs
alg.num_splits  = 1;        % # of different hold-outs
alg.naive_bayes = true;
alg.s0_trn      = sim{i}.model.s/4;     % list of #'s of hold-out samples for class 0
alg.s1_trn      = sim{i}.model.s/4;     % list of #'s of hold-out samples for class 1
alg.trn_s       = alg.s0_trn+alg.s1_trn;    % total # of training samples

constants = get_constants(sim{i}.As,sim{i}.targs);     % get constants to ease classification code

% subspace.naive_bayes = find(sum(sim{i}.As,3)>0);

nsteps=20;
inc=round(logspace(0,log10(sim{i}.model.n^2),nsteps));
coh=round(logspace(log10(1),log10(sim{i}.model.n),nsteps));

for j=1:nsteps
    subspace.incoherent = inc(j);
    subspace.coherent   = coh(j);
    [Lhat_tmp Lsem_tmp] = classify_holdout(sim{i}.As,constants,alg,subspace);
    Lhat_holdout.incoherent(j) = Lhat_tmp.incoherent;
    Lsem_holdout.incoherent(j) = Lsem_tmp.incoherent;

    Lhat_holdout.coherent(j) = Lhat_tmp.coherent;
    Lsem_holdout.coherent(j) = Lsem_tmp.coherent;
end

%%
clf
subplot(121)
errorbar([0 inc],[0.5 Lhat_holdout.incoherent],[0 Lsem_holdout.incoherent],...
    'color','k','linestyle','-','linewidth',2);
hold all
plot([sim{i}.model.m^2 sim{i}.model.m^2],[0 0.5])
axis('tight')
set(gca,'XScale','log','XTick',inc)
xlabel('assumed signal subgraph size')
ylabel('Lhat')
title('incoherent signal subgraph estimation')

subplot(122)
errorbar([0 coh],[0.5 Lhat_holdout.coherent],[0 Lsem_holdout.coherent],...
    'color','k','linestyle','-','linewidth',2);
hold all
plot([sim{i}.model.m sim{i}.model.m],[0 0.5])
axis('tight')
set(gca,'XScale','log','XTick',unique(coh))
xlabel('assumed # of nodes in egg')
ylabel('Lhat')
title('coherent signal subgraph estimation')


