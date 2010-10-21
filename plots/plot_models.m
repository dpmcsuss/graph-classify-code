


clear; clc

alg.datadir = '~/Research/data/graph_sims/unlabeled/';
alg.figdir  = '~/Research/figs/graph_sims/unlabeled/';
alg.fname   = 'models';
alg.save    = 1;

n   = 10;   % # of vertices
p   = 0.5;  % prob of connection for kidney
q0  = 0.25; % prob of connection for egg
q1  = 0.75; % prob of connection for egg
egg = 1:3;  % vertices in egg    

E0=p*ones(n);   % params in class 0
E0(egg,egg)=q0; % egg params in class 0

E1=p*ones(n);   % params in class 1
E1(egg,egg)=q1; % egg params in class 1

params1.E0=E0; params1.E1=E1;


n   = 10;   % # of vertices

E0=rand(n);   % params in class 0
E1=rand(n);   % params in class 1


params2.E0=E0; params2.E1=E1;


%% plot params

figure(1), clf
fs=6;
w=0.3;
h=0.3;
s=(1-3*w)/15;


emax=max(max([params1.E0(:) params1.E1(:) abs(params1.E0(:)-params1.E1(:))]));

% subplot('Position',[s .1 w h])
subplot(231)
image(60*params1.E0/emax)
colormap('gray')
title('class 0 mean','fontsize',fs)
% xlabel('vertices','fontsize',fs)
% ylabel('vertices','fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])

% subplot('Position',[2*s+w .1 w h])
subplot(232)
image(60*params1.E1/emax)
title('class 1 mean','fontsize',fs)
set(gca,'fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])

% subplot('Position',[4*s+2*w .1 w h])
subplot(233)
image(60*abs(params1.E0-params1.E1)/emax)
colormap('gray')
title('difference','fontsize',fs)
set(gca,'fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])

Position=get(gca,'Position');
pos(1)=Position(1)+Position(3)+0.01;
pos(2)=Position(2)+.05;
pos(3)=0.01;
pos(4)=Position(4)*.75;
colorbar('Position',pos,'YTick',[0 15 30 45 60],'YTickLabel',[0 25 50 75 100],'fontsize',fs)

% subplot('Position',[s .1 w h])
subplot(234)
image(60*params2.E0/emax)
colormap('gray')
title('class 0 mean','fontsize',fs)
xlabel('vertices','fontsize',fs)
ylabel('vertices','fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])

% subplot('Position',[2*s+w .1 w h])
subplot(235)
image(60*params2.E1/emax)
title('class 1 mean','fontsize',fs)
set(gca,'fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])

% subplot('Position',[4*s+2*w .1 w h])
subplot(236)
image(60*abs(params2.E0-params2.E1)/emax)
colormap('gray')
title('difference','fontsize',fs)
set(gca,'fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])

Position=get(gca,'Position');
pos(1)=Position(1)+Position(3)+0.01;
pos(2)=Position(2)*1.5;
pos(3)=0.01;
pos(4)=Position(4)*.75;
colorbar('Position',pos,'YTick',[0 15 30 45 60],'YTickLabel',[0 25 50 75 100],'fontsize',fs)

if alg.save
    wh=[3 2];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_params'];
    print('-dpdf',figname)
    print('-deps',figname)
    saveas(gcf,figname)
end

