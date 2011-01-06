%% load data and prepare stuff for saving

clear, clc, clf

etc.dataset     = 'BLSA50';
etc.dir_data    = '~/Research/data/MRI/BLSA/';
etc.dir_results = [etc.dir_data 'results/'];
etc.dir_figs    = '~/Research/figs/MRI/BLSA/';
etc.fname       = '_Lhats_loo_vs_s_out';

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

s0_trn = [1:constants.s0-1, (constants.s0-1)*ones(1,8)];   % list of #'s of hold-out samples for class 0
s1_trn = 1:constants.s1-1;   % list of #'s of hold-out samples for class 1

subspace(1).name='naive_bayes';
subspace(2).name='incoherent';
subspace(3).name='incoherent_z';
subspace(4).name='incoherent_f';

%% perform classification

for i=1:length(s0_trn)
    disp(i)
    
    for j=1:length(subspace)
        subspace(j).indices = [];
        subspace(j).size    = 40; 
    end

    xval.s0_trn = s0_trn(i);
    xval.s1_trn = s1_trn(i);
    
    max_iters=choose(constants.s0,xval.s0_trn)*choose(constants.s1,xval.s1_trn);
    xval.num_iters=min(1000,max_iters);
    
    Lhats{i} = classify_par(As,constants,xval,subspace);
    mean(Lhats{i})
end


%% plot performance

figure(1), clf, hold on

num_iters=21*29;

for i=1:length(s0_trn)
    Lavg(i,:)=mean(Lhats{i});
    Lsem(i,:)=std(Lhats{i})/sqrt(num_iters);
end

h = errorbar(repmat((s1_trn+s0_trn)',1,length(subspace)),Lavg,Lsem);
% h = plot(repmat(s0_trn',1,length(subspace)),Lavg,'linewidth',2);

set(h(1),'Color',0.5*[1 1 1],'linestyle','-');
set(h(2),'Color',0.5*[1 1 1],'linestyle','--');
set(h(3),'Color',0.0*[1 1 1],'linestyle','-');
set(h(4),'Color',0.0*[1 1 1],'linestyle','--');


XTickLabel=[];
for i=1:length(subspace)
    XTickLabel=[XTickLabel {subspace(i).name}];
end
ylabel('misclassification rate')
xlabel('number of training samples')
title('Leave-two-out Cross-Validation')
legend(XTickLabel,'Location','Best')
% set(gca,'Xscale','log')
% axis([1 length(s0_trn) 0 0.5])

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
