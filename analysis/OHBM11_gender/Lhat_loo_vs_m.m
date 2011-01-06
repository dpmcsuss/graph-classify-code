%% load data and prepare stuff for saving

clear, clc, clf

etc.dataset     = 'BLSA50';
etc.dir_data    = '~/Research/data/MRI/BLSA/';
etc.dir_results = [etc.dir_data 'results/'];
etc.dir_figs    = '~/Research/figs/MRI/BLSA/';
etc.fname       = '_Lhats_loo_vs_m';

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

subspace(1).name='incoherent';
subspace(2).name='incoherent_z';
subspace(3).name='incoherent_f';

%% perform classification

params = get_ind_edge_params(As,constants);
max_edges = length(find(abs(params.E0-params.E1)>0));
num_edges=1:max_edges;

for i=1:length(num_edges)
    disp(i)
    
    for j=1:length(subspace)
        subspace(j).indices=[];
        subspace(j).size=num_edges(i); 
    end
    [subspace params] = get_subspace_indices(params,subspace,As,constants);
    
    Lhats{i} = classify_par(As,constants,xval,subspace);
    mean(Lhats{i})
end


%% plot performance

figure(1), clf, hold on

num_iters=21*29;

for i=1:length(num_edges)
    Lavg(i,:)=mean(Lhats{i});
    Lsem(i,:)=std(Lhats{i})/sqrt(num_iters);
end

% errorbar(repmat(num_edges',1,length(subspace)),Lavg,Lsem)
h = plot(repmat(num_edges',1,length(subspace)),Lavg,'linewidth',2);

set(h(1),'Color',0.67*[1 1 1],'linestyle','--');
set(h(2),'Color',0.33*[1 1 1],'linestyle','--');
set(h(3),'Color',0.00*[1 1 1],'linestyle','-');


XTickLabel=[];
for i=1:length(subspace)
    XTickLabel=[XTickLabel {subspace(i).name}];
end
ylabel('misclassification rate')
xlabel('size of signal subgraph')
title('Leave-one-out Incoherent Signal Subgraph Performance estimates')
legend(XTickLabel,'Location','Best')
set(gca,'Xscale','log')
axis([1 length(num_edges) 0 0.5])

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
