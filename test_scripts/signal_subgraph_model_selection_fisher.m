

%% for cross-validation
clear 
i=5;

M.n = 100;
M.m = 5;
M.s = 1000;

p = rand(M.n);
q = p+randn(M.n)*.01;
q(1:M.m,1:M.m)=rand(M.m);
q(q<0)=eps;
q(q>1)=1-eps;

M.E0 = p.*ones(M.n); 
M.E1 = q.*ones(M.n); 

sim{i} = gen_ind_edge_data(M);
sim{i}.model = M;
sim{i}.name='kidney_egg4';
sim{i}.Lstar = 1-binocdf(floor(M.m^2/2),M.m^2,p-q);


alg.num_repeats = 10;        % # of times to repeat using the same # of hold-outs
alg.num_splits  = 1;        % # of different hold-outs
alg.naive_bayes = true;
alg.s0_trn      = 20;     % list of #'s of hold-out samples for class 0
alg.s1_trn      = 20;     % list of #'s of hold-out samples for class 1
alg.trn_s       = alg.s0_trn+alg.s1_trn;    % total # of training samples

constants = get_constants(sim{i}.As,sim{i}.targs);     % get constants to ease classification code

subspace.naive_bayes = find(sum(sim{i}.As,3)>0);
subspace.fisher_incoherent=M.m^2;
subspace.incoherent=M.m^2;


%%

params = get_params_mle(sim{i}.As,constants);

if isfield(subspace,'fisher_incoherent')
    A0 = sim{i}.As(:,:,constants.y0);
    A1 = sim{i}.As(:,:,constants.y1);
   params.d_pval = get_Fisher_pvals(A0,A1,constants.d,constants.s0,constants.s1); 
end


% subspace = get_subspace(params,subspace);
% plot(subspace.incoherent, subspace.fisher_incoherent)

[pos indpos] = sort(params.d_pos(:),'descend');
[pval indpval] = sort(params.d_pval(:),'ascend');
[cor indcor] = sort(params.d_opt(:),'ascend');

plot(indpos,indpval)



%%
[Lhat_tmp Lsem_tmp] = classify_holdout(sim{i}.As,constants,alg,subspace);

%%
figure(1), clf
errorbar([Lhat_tmp.naive_bayes Lhat_tmp.incoherent Lhat_tmp.fisher_incoherent],...
    [Lsem_tmp.naive_bayes Lsem_tmp.incoherent Lsem_tmp.fisher_incoherent]) %,...
%     'color','k','linestyle','-','linewidth',2);
axis('tight')
legend('nb','inc','fisher')

figure(2), 
boxplot([Lhat_tmp.naive_bayes Lhat_tmp.incoherent Lhat_tmp.fisher_incoherent])
set(gca,'XTick',1:3,'XTickLabel',[{'nb'}; {'inc'}; {'fisher'}])
