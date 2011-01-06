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
    As(:,:,i)=AA;
    targs(i) = uniq_keeps(i).sex;
end

constants = get_constants(As,targs);     % get constants to ease classification code

A0 = As(:,:,constants.y0);
A1 = As(:,:,constants.y1);

IX=ones(70);
IX=find(triu(IX,+1));

vecA0=nan(length(IX),constants.s0);
for i=1:constants.s0
    tmp=A0(:,:,i);    
    vecA0(:,i) = tmp(IX);
end
binA0 = vecA0;
binA0(binA0<=0.2)=0;
binA0(binA0>0.2)=1;


vecA1=nan(length(IX),constants.s1);
for i=1:constants.s1
    tmp=A1(:,:,i);    
    vecA1(:,i) = tmp(IX);
end
binA1 = vecA1;
binA1(binA1<=0.2)=0;
binA1(binA1>0.2)=1;


[h,p,ci,stats] = ttest2(vecA0',vecA1');
[p_ttest ind_ttest] = sort(p);
plot(p_ttest)

pval = get_Fisher_pvals(binA0,binA1,2415,constants.s0,constants.s1);
[p_fisher ind_fisher] = sort(pval);

%%
for i=2415:-1:100
    
    figure(1), clf
    subplot(121)
    hist(vecA0(ind_fisher(i),:))
    set(gca,'XLim',[0 1])
    
    subplot(122)
    hist(vecA1(ind_fisher(i),:))
    set(gca,'XLim',[0 1])

    pause
end

