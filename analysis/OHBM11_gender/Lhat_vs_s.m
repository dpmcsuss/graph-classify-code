clear, clc, clf
alg.fname='BLSA50';
datadir = '~/Research/data/MRI/BLSA/';
alg.datadir = datadir;
alg.postdir = [alg.datadir 'results/'];
alg.figdir = '~/Research/figs/MRI/BLSA/';

load([alg.postdir alg.fname '_data'])
s=length(uniq_keeps);
As=nan(70,70,s);
targs=nan(70);
for i=1:s
    AA = uniq_keeps(i).Adj;
    AA(isnan(AA))=0;
    AA(AA<0.2)=0;
    AA(AA>0.2)=1;
    As(:,:,i)=AA;
    targs(i) = uniq_keeps(i).sex;
end

constants = get_constants(As,targs);     % get constants to ease classification code

%% 

alg.s0_trn = 15; %[2:2:18 18 18 18 18];      % list of #'s of hold-out samples for class 0
alg.s1_trn = 15; %2:2:26;                    % list of #'s of hold-out samples for class 1
alg.num_repeats = 1000;                 % # of times to repeat using the same # of hold-outs

% get fisher significant
d_pval = run_get_fisher_pvals(As,constants);
[~, fis_ind]=sort(d_pval(:),'ascend');

deg=sum(d_pval,1)+sum(d_pval,2)';
[foo deg_ind] = sort(deg,'ascend');

% get pos and opt rankings
params = get_params_mle(As,constants);
[~, inc_ind]=sort(params.d_pos(:),'descend');
[~, cor_ind]=sort(params.d_opt(:),'descend');

% choose subspaces
subspace.inc = inc_ind(1:40);
subspace.fis = fis_ind(1:40);
subspace.cor = cor_ind(1:40);
subspace.nb  = find(params.d_pos>0);

ind_max_mat=zeros(70);
ind_max_mat(deg_ind(1:16),deg_ind(1:16))=1;
ind_max_mat=triu(ind_max_mat,+1);
subspace.deg=find(ind_max_mat);

% perform classification
[Lhat_tmp Lsem_tmp] = classify_holdout(As,constants,alg,subspace);



%% plot performance

% clear, load('/Users/jovo/Research/data/MRI/BLSA/results/BLSA50_Lhat_vs_s.mat')

prep_plots;

alg.trn_s=alg.s0_trn+alg.s1_trn;


for i=1:n_subspaces
        errorbar(alg.trn_s,mean(Lhats.(fn{i})),std(Lhats.(fn{i}))./sqrt(alg.trn_s),...
            'color',colors{i},'linestyle',ls{i},'linewidth',2);
end

axis('tight')
legend('norm. inc.','fish. inc.','inc.','n.b.','coh_f_i_s')
xlabel('# of training samples','fontsize',fs)
ylabel('misclassification rate','fontsize',fs)


alg.save=0;
alg.ext='_Lhat_vs_s';
save_plots;