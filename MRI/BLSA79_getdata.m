clear, clc
alg.datadir = '~/Research/data/MRI/BLSA/BLSA79/';
alg.figdir  = '~/Research/figs/MRI/BLSA/';
alg.fname   = 'BLSA79';
alg.save    = 1;


files   = dir([alg.datadir 'quotient_graphs/']); files(1:2)=[];
load([alg.datadir 'IACL-blsa-subjectdata-1110.mat']);
S       = length(files);
n       = 70;
As      = zeros(n,n,S);
targs   = nan(1,S);

for s=1:S
    A= importdata([alg.datadir 'quotient_graphs/' files(s).name]);
    A(isnan(A))=0;
    A(A>=0.2)=1;
    A(A<0.2)=0;
    As(:,:,s)=A(2:n+1,2:n+1);
    targs(s)=blsaSubjects.male(s);
end

save([alg.datadir 'BLSA79_data'],'As','targs','alg')

