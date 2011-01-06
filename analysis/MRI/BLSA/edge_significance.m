clear, clc
load('/Users/jovo/Research/data/MRI/BLSA/results/BLSA50_data.mat')

%% get data

n=70;
s=length(uniq_keeps);
Ab=zeros(n,n,s);
Aw=zeros(n,n,s);
t=0.5;
for i=1:s
    A=uniq_keeps(i).Adj;
    A(isnan(A))=0;
    Aw(:,:,i)=A;

%     t=quantile(A(:),0.2);
    A(A<t)=0;
    A(A>t)=1;
    Ab(:,:,i)=A;
    ys(i)=uniq_keeps(i).sex;
end
constants = get_constants(Ab,ys);

%% binary significance

pvals=run_get_fisher_pvals(Ab,constants);
inds=find(pvals<1);

figure(1), clf
subplot(121)
plot(sort(pvals(inds)))
axis('tight')


%% get weighted significance

A0 = Aw(:,:,constants.y0);
A1 = Aw(:,:,constants.y1);

vecA0 = reshape(A0,[constants.d constants.s0]);
vecA1 = reshape(A1,[constants.d constants.s1]);

[h,p,ci,stats] = ttest2(vecA0',vecA1');

inds=find(p<1);

subplot(122)
plot(sort(p(inds)))
axis('tight')