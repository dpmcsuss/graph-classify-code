%% load data and prepare stuff for saving

clear, clc, clf

etc.dataset     = 'BLSA50';
etc.dir_data    = '~/Research/data/MRI/BLSA/';
etc.dir_results = [etc.dir_data 'results/'];
etc.dir_figs    = '~/Research/figs/MRI/BLSA/';
etc.fname       = '_Lhats_loo';

load([etc.dir_results etc.dataset '_data'])
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

%% set algorithm stuff

xval.s0_trn = constants.s0-1;   % list of #'s of hold-out samples for class 0
xval.s1_trn = constants.s1-1;   % list of #'s of hold-out samples for class 1

subspace(1).name='naive_bayes';

subspace(2).name='incoherent';
subspace(2).size=40;

subspace(3).name='incoherent_z';
subspace(3).size=40;

subspace(4).name='incoherent_f';
subspace(4).size=40;

%% perform classification
[Lhats xval] = classify_par(As,constants,xval,subspace);
mean(Lhats)


%% plot performance

figure(1), clf, hold on

errorbar(mean(Lhats),std(Lhats)/sqrt(xval.num_iters),'.')

XTickLabel=[];
for i=1:length(subspace)
    XTickLabel=[XTickLabel {subspace(i).name}];
end
set(gca,'XTick',1:length(subspace),'XTickLabel',XTickLabel,'fontsize',12)
ylabel('misclassification rate')


etc.save=1;
if etc.save==1
    figname=[etc.dir_figs etc.dataset etc.fname];
    wh=[6 3];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    print('-dpdf',figname)
    print('-dpng',figname)
    saveas(gcf,figname)
    save([etc.dir_results etc.dataset etc.fname])
end