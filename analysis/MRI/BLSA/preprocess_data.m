
% get BLSA79 data

clear all, clc
alg.datadir = '~/Research/data/MRI/BLSA/';
alg.figdir  = '~/Research/figs/MRI/BLSA/';
alg.fname   = 'BLSA79';
alg.save    = 1;


files   = dir([alg.datadir '2010_11_07/quotient_graphs/']); files(1:2)=[];
load([alg.datadir '2010_11_07/IACL-blsa-subjectdata-1110.mat']);
S       = length(files);
n       = 70;
As    = zeros(n,n,S);
targs = nan(S,1);
ids   = char(zeros(S,4));

for s=1:S
    A= importdata([alg.datadir '2010_11_07/quotient_graphs/' files(s).name]);
    A(isnan(A))=0;
    %     A(A>=0.2)=1;
    %     A(A<0.2)=0;
    As(:,:,s)=A(2:n+1,2:n+1);
    targs(s)=blsaSubjects.male(s);
    ids(s,:)=blsaSubjects.id(s,:);
end

save([alg.datadir '/results/BLSA79_data'],'As','targs','alg')

As_raw=As;


%%

for i=1:S
    A=As(:,:,i);
    A(A>=0.2)=1;
    A(A<0.2)=0;
    As(:,:,i)=A;
end

As_bin=As;

%%
    
As=As_bin;

subspace.naive_bayes=find(sum(As,3));

constants = get_constants(As,targs);     % get constants to ease classification code
[Lhat params subspace] = run_naive_bayes(As,constants,As,constants,subspace);
disp(Lhat)

%%
i0=1; i1=1; clear A0 A1
subspace.naive_bayes=find(sum(As,3));
A0=zeros(length(subspace.naive_bayes),constants.s-sum(targs));
A1=zeros(length(subspace.naive_bayes),sum(targs));
for i=1:constants.s
    A=As(:,:,i);
    if targs(i)==0
        A0(:,i0)=A(subspace.naive_bayes);
        i0=i0+1;
    else
        A1(:,i1)=A(subspace.naive_bayes);
        i1=i1+1;
    end
end

%%

s0=constants.s0;
s1=constants.s1;
for i=1:length(subspace.naive_bayes)
    ss0=sum(A0(i,:));
    x(1,1)=ss0;
    x(1,2)=s0-ss0;
    
    ss1=sum(A1(i,:));
    x(2,1)=ss1;
    x(2,2)=s1-ss1;

    pvalue=FisherExactTest22(x);
    pval(i)=pvalue(3);
    
end

% subspace.incoherent=find(pval<0.04);


%%

% As=As_bin;
% constants = get_constants(As,targs);     % get constants to ease classification code

[spvals ixpvals] = sort(pval);
Lhat_inc=nan(1,length(subspace.naive_bayes));

for i=1:length(subspace.naive_bayes)
    subspace.incoherent=subspace.naive_bayes(ixpvals(1:i));
    [Lhat params subspace] = run_naive_bayes(As,constants,As,constants,subspace);
    Lhat_inc(i)=Lhat.incoherent;
end

plot(Lhat_inc)
    
    
    
%% 

As=As_raw;



