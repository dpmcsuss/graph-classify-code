% clear, clc
% alg.datadir = '~/Research/data/MRI/BLSA/2010_11_07/';
% alg.figdir  = '~/Research/figs/MRI/BLSA/';
% alg.fname   = 'BLSA79';
% alg.save    = 1;
% 
% 
% files   = dir([alg.datadir 'quotient_graphs/']); files(1:2)=[];
% load([alg.datadir 'IACL-blsa-subjectdata-1110.mat']);
% S       = length(files);
% n       = 70;
% As      = zeros(n,n,S);
% targs   = nan(1,S);
% 
% for s=1:S
%     A= importdata([alg.datadir 'quotient_graphs/' files(s).name]);
%     A(isnan(A))=0;
%     
%     A(A>=0.2)=1;
%     A(A<0.2)=0;
%     As(:,:,s)=A(2:n+1,2:n+1);
%     targs(s)=blsaSubjects.male(s);
% end
% 
% save([alg.datadir 'BLSA79_data'],'As','targs','alg')

%%

clear, clc
alg.datadir = '~/Research/data/MRI/BLSA/';
alg.figdir  = '~/Research/figs/MRI/BLSA/';
subdir      = '2010_11_07/';
alg.fname   = 'BLSA79';
alg.save    = 1;


files   = dir([alg.datadir subdir 'IACL-blsadata-79subjects-1110/']); files(1:2)=[];
csvfile = importdata([alg.datadir subdir 'IACL-blsa-subjectdata-1110.csv']);
S       = length(files);
n       = 70;
As      = nan(n,n,S);
ys   = nan(1,S);

A=ones(71);
A=triu(A,+1);
imagesc(A)
nb=find(A);
j=1;

for s=1:S
    A= importdata([alg.datadir subdir 'IACL-blsadata-79subjects-1110/' files(s).name]);

    num_zeros=length(find(A(nb)<0.2))
    if num_zeros>0
%         imagesc(A);
%         keyboard
    else
        j=j+1;
    end
    A(isnan(A))=0;
    As(:,:,j)=A(2:n+1,2:n+1);
    for k=1:length(csvfile.textdata)
        if strcmpi(files(s).name(24:27),csvfile.textdata{k,1})
            if strcmpi(csvfile.textdata{k,3},'male')
                ys(j)=1;
            else
                ys(j)=0;
            end
        end
    end

end

As=As(:,:,1:j);
ys=ys(1:j);

save([alg.datadir 'BLSA' num2str(j) '_data'],'As','ys','alg')
