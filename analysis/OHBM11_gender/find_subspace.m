clear, clc, clf
alg.fname='BLSA50';
datadir = '~/Research/data/MRI/BLSA/';
alg.datadir = datadir;
alg.postdir = [alg.datadir 'results/'];
alg.figdir = '~/Research/figs/MRI/BLSA/';

load([alg.postdir alg.fname '_data'])
s=length(uniq_keeps);
As=zeros(70,70,s);
for i=1:s
    AA = uniq_keeps(i).Adj;
    AA(isnan(AA))=0;
    AA(AA<0.2)=0;
    AA(AA>0.2)=1;
    As(:,:,i)=AA;
    targs(i) = uniq_keeps(i).sex;
end

constants = get_constants(As,targs);     % get constants to ease classification code


alg.s0_trn = 15; %constants.s0-5;     % list of #'s of hold-out samples for class 0
alg.s1_trn = 15; %constants.s1-5;     % list of #'s of hold-out samples for class 1

alg.num_repeats = 1000;        % # of times to repeat using the same # of hold-outs
alg.num_splits  = 1;        % # of different hold-outs
alg.naive_bayes = true;
alg.trn_s       = alg.s0_trn+alg.s1_trn;    % total # of training samples


subspace.naive_bayes = find(sum(As,3)>0);
subspace.incoherent = [];
% subspace.fisher_incoherent = [];
% subspace.cor_incoherent = [];

subspace_names  = fieldnames(subspace);             % pre-allocate space for Lhats
n_subspaces=length(subspace_names);



%%
num_coh_vertices=1:70;
num_inc_edges=50; %5:20:2000;

for j=1:length(num_inc_edges)
    disp(num_inc_edges(j))
    subspace.incoherent = num_inc_edges(j);
%     subspace.coherent = num_coh_vertices(j);
    [Lhat_tmp Lsem_tmp] = classify_holdout(As,constants,alg,subspace);
    
    for k=1:n_subspaces
        Lhats.(subspace_names{k})(j)=mean(Lhat_tmp.(subspace_names{k}));
        Lsems.(subspace_names{k})(j)=std(Lsem_tmp.(subspace_names{k}))/sqrt(alg.num_repeats);
    end
    
end

%%
fn=subspace_names;
fs=12;

colors{1}='r';
colors{2}='g';
colors{3}='b';
colors{4}='k';

ls{1}='--';
ls{2}='-';
ls{3}=':';
ls{4}='-.';

for i=1:70, blowup(i)=choose(i,2); end

figure(1), clf, hold on
for i=1:n_subspaces
    if ~strcmp(fn{i},'coherent')
        plot(num_inc_edges,Lhats.(fn{i}),... Lsems.(fn{i}),...
            'color',colors{i},'linestyle',ls{mod(i,3)+1},'linewidth',2);
    else
        plot(blowup(7:end),Lhats.(fn{i})(7:end),... Lsems.(fn{i}),...
            'color',colors{i},'linestyle',ls{mod(i,3)+1},'linewidth',2);
    end
end

axis('tight')
legend('inc','cor','fisher','coh')
set(gca,'XScale','log')
xlabel('size of signal subgraph','fontsize',12)
ylabel('misclassification rate','fontsize',12)


wh=[3 2]*1.5;   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname=[alg.figdir alg.fname '_Lhats_subspace'];
print('-dpdf',figname)
print('-deps',figname)
saveas(gcf,figname)


