clear, clc, close all
load('~/Research/necog/data/BLSA/MR_14')
figdir='/Users/joshyv/Research/necog/figs/MR/';

%% plot images
figure(111), clf, colormap('default')
s0=2; s1=1;
amax=max(As(:));
As1=60*As/amax;
for i=1:G.s
    if ys(i)==0
        subplot(ceil(G.s/2),2,s0)
        s0=s0+2;
    else
        subplot(ceil(G.s/2),2,s1)
        s1=s1+2;
    end
    image(As1(:,:,i))
    title(num2str(i))
    set(gca,'XTick',[],'YTick',[])
end

wh=[3 8];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname=[figdir 'MR14_images'];
print('-dpdf',figname)

%% plot images
figure(118), clf, colormap('default')
s0=2; s1=1;
amax=max(As(:));
As1=60*As/amax;
for i=1:G.s
    if ys(i)==0
        subplot(ceil(G.s/2),2,s0)
        s0=s0+2;
    else
        subplot(ceil(G.s/2),2,s1)
        s1=s1+2;
    end
    spy(As1(:,:,i))
%     title(num2str(i))
    set(gca,'XTick',[],'YTick',[])
end

wh=[3 8];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname=[figdir 'MR14_spy'];
print('-dpdf',figname)


%% plot histograms
figure(112)
k=50;
A=As(:,:,1);

h(1)=subplot(411);
hist(A(:),k)
title('hist of subject 1')

A0      = As(:,:,G.y0);
A1      = As(:,:,G.y1);

h(2)=subplot(412);
hist(A0(:),k)
title('hist class 0')
h(3)=subplot(413);
hist(A1(:),k)
title('hist class 1')

h(4)=subplot(414);
hist(As(:),k)
title('total hist')

linkaxes(h,'x')

wh=[6 8];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname=[figdir 'hist_FA'];
print('-dpdf',figname)

%% plot mean FA for subjects, which is just weird
var0=[]; mean0=[]; var1=[]; mean1=[];
xticklabel0=[]; xticklabel1=[];
for i=1:G.s
    A=As(:,:,i);
    if ys(i)==0
        var0=[var0 var(A(:))];
        mean0=[mean0 mean(A(:))];
        xticklabel0=[xticklabel0 {ids{i}}];
    else
        var1=[var1 var(A(:))];
        mean1=[mean1 mean(A(:))];
        xticklabel1=[xticklabel1 {ids{i}}];
    end
end


figure(4), clf
h(1)=subplot(311); hold all
plot(mean0,'k')
% plot(var0,'--k')
set(gca,'XTickLabel',xticklabel0)
h(2)=subplot(312); hold all
plot(mean1,'r')
% plot(var1,'--r')
set(gca,'XTickLabel',xticklabel1)
h(3)=subplot(313); hold all
plot(mean0,'k')
plot(mean1,'r')

% plot(var0,'--k')
% plot(var1,'--r')

linkaxes(h,'y')
axis([1 7 0.15 0.4])

wh=[5 7];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname=[figdir 'mean_FA'];
print('-dpdf',figname)
% print('-depsc',figname)
% saveas(fig,figname)

%% E0 and E1

figure(8), clf
pmax=max(max([P.E0(:) P.E1(:)]));
subplot(231), image(60*P.E0/pmax), title('E0')
subplot(232), image(60*P.E1/pmax), title('E1')
subplot(233), image(60*P.del/pmax), title('deltahat')

k=50; clear h
h(1)=subplot(234); hist(P.E0(:),k), title('E0')
h(2)=subplot(235); hist(P.E1(:),k), title('E1')
h(3)=subplot(236); hist(P.del(:),k), title('deltahat')
linkaxes(h,'x')

wh=[7 5];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname=[figdir 'E0_E1'];
print('-dpdf',figname)

%% block histograms

G.b = 2;
P   = MR_get_params(As,G.y0,G.y1,G);

B0{1,1} = P.E0(1:G.n/2      , 1:G.n/2);
B0{1,2} = P.E0(1:G.n/2      , G.n/2+1:G.n);
B0{2,1} = P.E0(G.n/2+1:G.n  , 1:G.n/2);
B0{2,2} = P.E0(G.n/2+1:G.n  , G.n/2+1:G.n);

B1{1,1} = P.E1(1:G.n/2      , 1:G.n/2);
B1{1,2} = P.E1(1:G.n/2      , G.n/2+1:G.n);
B1{2,1} = P.E1(G.n/2+1:G.n  , 1:G.n/2);
B1{2,2} = P.E1(G.n/2+1:G.n  , G.n/2+1:G.n);

figure(5), clf
k=20;
subplot(241)
hist(B0{1,1}(:),k)
subplot(242)
hist(B0{1,2}(:),k)
subplot(245)
hist(B0{2,1}(:),k)
subplot(246)
hist(B0{2,2}(:),k)

subplot(243)
hist(B1{1,1}(:),k)
subplot(244)
hist(B1{1,2}(:),k)
subplot(247)
hist(B1{2,1}(:),k)
subplot(248)
hist(B1{2,2}(:),k)

wh=[7 5];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname=[figdir 'block_histograms'];
print('-dpdf',figname)



%% plot blockmodels


figure(7), clf
BB0 = 60*P.B0/max(P.B0(:));
BB1 = 60*P.B1/max(P.B1(:));
subplot(221), image(BB0); colorbar; title('mean B0')
subplot(222), image(BB1); colorbar; title('mean B1')

% figure(8),
for i=1:2
    for j=1:2
        BB0=B0{i,j};
        BB1=B1{i,j};
        V0(i,j) = var(BB0(:));
        V1(i,j) = var(BB1(:));
    end
end

vmax = max(max(V0(:)),max(V1(:)));
V0 = 60*V0/vmax;
V1 = 60*V1/vmax;

subplot(223), image(V0); colorbar; title('var B0')
subplot(224), image(V1); colorbar; title('var B1')
        
wh=[5 7];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname=[figdir 'mean_var_block'];
print('-dpdf',figname)
% print('-depsc',figname)
% saveas(fig,figname)



































