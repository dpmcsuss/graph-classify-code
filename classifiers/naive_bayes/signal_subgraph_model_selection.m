

%% for cross-validation
clear
i=5;

M.n = 100;
M.m = 3;
M.s = 200;

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

inc=round(logspace(0,log10(sim{i}.model.n^2),20));

alg.num_repeats = 1;        % # of times to repeat using the same # of hold-outs
alg.num_splits  = 1;        % # of different hold-outs
alg.naive_bayes = true;
alg.s0_trn      = sim{i}.model.s/4;     % list of #'s of hold-out samples for class 0
alg.s1_trn      = sim{i}.model.s/4;     % list of #'s of hold-out samples for class 1
alg.trn_s       = alg.s0_trn+alg.s1_trn;    % total # of training samples

constants = get_constants(sim{i}.As,sim{i}.targs);     % get constants to ease classification code

% subspace.naive_bayes = find(sum(sim{i}.As,3)>0);

for j=1:length(inc)
%     clear subspace Lhat_tmp Lsem_tmp Lhat_holdout Lsem_holdout
    subspace.incoherent  = inc(j);
    [Lhat_tmp Lsem_tmp] = classify_holdout(sim{i}.As,constants,alg,subspace);
    Lhat_holdout.incoherent(j) = Lhat_tmp.incoherent;
    Lsem_holdout.incoherent(j) = Lsem_tmp.incoherent;
end


clf
errorbar([0 inc],[0.5 Lhat_holdout.incoherent],[0 Lsem_holdout.incoherent],...
    'color','k','linestyle','-','linewidth',2);
hold all
plot([sim{i}.model.m^2 sim{i}.model.m^2],[0 0.5])
axis('tight')
set(gca,'XScale','log')
set(gca,'XTick',inc)
xlabel('assumed signal subgraph size')
ylabel('Lhat (sem)')

print('-dpdf','/Users/jovo/Research/figs/sims/trunk_like')
