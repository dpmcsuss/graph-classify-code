clear, clc,
load('/Users/jovo/Research/data/MRI/BLSA/results/edge_significance.mat')

[foo IX] = sort(d_pval(:));
[I J] = ind2sub([70 70],IX);
CorticalLabels = importdata('/Users/jovo/Research/data/MRI/BLSA/misc/IACL-corticallabeldescriptions-1110.txt',' ');

clear orderedregions
m=find(foo<0.999999999999, 1, 'last' );
maxlen=19;
for i=1:m
    
    if I(i)<36, ii=I(i)+2; else ii=I(i)+3; end
    if J(i)<36, jj=J(i)+2; else jj=J(i)+3; end
    
    if ii<=11, ii_start=7; elseif ii<39, ii_start=8; else ii_start=9; end
    if jj<=11, jj_start=7; elseif jj<39, jj_start=8; else jj_start=9; end
    
    label1=CorticalLabels{ii}(ii_start:end-1);
    if length(label1)>maxlen, label1=label1(1:maxlen); end
    label2=CorticalLabels{jj}(jj_start:end-1);
    if length(label2)>maxlen, label2=label2(1:maxlen); end
    
    orderedregions{i,1}=foo(i);
    orderedregions{i,2}=label1;
    orderedregions{i,3}=label2;
    
end


%%
clf, clc

n=70;
AA=linspace(0,1,n^2);
AA=reshape(AA',[n n]);
AA=(AA+AA')/2;
AA(2:13,38:39)=1;
AA=triu(AA,+1);
AAT=nan(n);

for i=1:n
    for j=i:n
        
        clear I J
        if i<=n/2 && j>n/2
            J=i+n/2; I=j-n/2;
        elseif  i>n/2 && j<=n/2
            I=i-n/2; J=j+n/2;
        elseif i<=n/2 && j<=n/2
            I=i+n/2; J=j+n/2;
        elseif i>n/2 && j>n/2
            I=i-n/2; J=j-n/2;
        end
        
        AAT(I,J)=AA(i,j);
    end
end

figure(1), clf
subplot(221), imagesc(AA); colorbar
subplot(222), imagesc(AAT); colorbar



%%

d_pvalT=ones(n);
for i=1:n
    for j=i:n
        
        clear I J
        if i<=n/2 && j>n/2
            J=i+n/2; I=j-n/2;
        elseif  i>n/2 && j<=n/2
            I=i-n/2; J=j+n/2;
        elseif i<=n/2 && j<=n/2
            I=i+n/2; J=j+n/2;
        elseif i>n/2 && j>n/2
            I=i-n/2; J=j-n/2;
        end
        
        d_pvalT(I,J)=d_pval(i,j);
    end
end

subplot(223), imagesc(d_pval); colorbar
subplot(224), imagesc(d_pvalT); colorbar


%%
ass=ones(n);
inds=find(triu(ass,+1));

figure(2), clf
scatter(d_pval(inds),d_pvalT(inds))


%%
% params = get_ind_edge_params(As,constants);
% d_inc=abs(params.E0-params.E1);
% [foo IX] = sort(d_inc(:),'descend');
% [I J] = ind2sub([70 70],IX);
% 
% m=find(foo>0, 1, 'last' );
% for i=1:m
%     
%     incs(i)=d_inc(I(i),J(i));
%     
%     if I(i)<=35; II=I(i)+35; else II=I(i)-35; end
%     if J(i)<=35; JJ=J(i)+35; else JJ=J(i)-35; end
%     
%     flipped(i)=d_inc(II,JJ);
%     
% end
% 
% scatter(flipped,incs)


%%
% clear, clc, clf
% 
% etc.dataset     = 'BLSA50';
% etc.dir_data    = '~/Research/data/MRI/BLSA/';
% etc.dir_results = [etc.dir_data 'results/'];
% etc.dir_figs    = '~/Research/figs/MRI/BLSA/';
% etc.fname       = '_Lhats_loo_vs_m';
% 
% load([etc.dir_results etc.dataset '_data'])
% s=length(uniq_keeps);
% As=nan(70,70,s);
% targs=nan(70);
% for i=1:s
%     AA = uniq_keeps(i).Adj;
%     AA(isnan(AA))=0;
%     %     AA(AA<0.2)=0;
%     %     AA(AA>0.2)=1;
%     As(:,:,i)=AA;
%     targs(i) = uniq_keeps(i).sex;
% end
% 
% constants = get_constants(As,targs);     % get constants to ease classification code
% 


%%
% AAflip=fliplr(flipud(d_pval))';
% 
% subplot(131), imagesc(d_pval)
% subplot(132), imagesc(AAflip)
% subplot(133), scatter(AAflip(:),d_pval(:))


%%
%
% quad=AA(1:n/2,n/2+1:end);
% quadT=fliplr(flipud(quad))';
%
% subplot(323), image(quad); colorbar
% subplot(324), image(quadT); colorbar
%
%
% AATT=AAT;
% AATT(1:n/2,n/2+1:end)=quadT;
%
% subplot(325), image(AA); colorbar
% subplot(326), image(AATT); colorbar
%
%
%
% AAflip=fliplr(flipud(AA))';
% subplot(121), imagesc(AA)
% subplot(122), imagesc(AAflip)