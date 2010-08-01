function plot_params(est_params,alg,params)
%% plot params

figure(1), clf
if nargin==3, nrows=2; else nrows=1; end
ncols=3;


if nargin==3
    emax=max(max([est_params.E0(:) est_params.E1(:) abs(est_params.E0(:)-est_params.E1(:))]));
    subplot(nrows,ncols,1), cla
    image(60*params.E0/emax)
    colormap('gray')
    title('class 0 mean')
    xlabel('vertices')
    ylabel('vertices')

    subplot(nrows,ncols,2), cla
    image(60*params.E1/emax)
    title('class 1 mean')

    subplot(nrows,ncols,3), cla
    image(60*abs(params.E0-params.E1)/emax)
    colormap('gray')
    title('difference')
end

emax=max(max([est_params.E0(:) est_params.E1(:) abs(est_params.E0(:)-est_params.E1(:))]));
subplot(nrows,ncols,1+ncols*(nrows-1)), cla
image(60*est_params.E0/emax)
colormap('gray')
title('estimated class 0 mean')
xlabel('vertices')
ylabel('vertices')

subplot(nrows,ncols,2+ncols*(nrows-1)), cla
image(60*est_params.E1/emax)
title('estimated class 1 mean')

subplot(nrows,ncols,3+ncols*(nrows-1)), cla
image(60*est_params.d_pos/emax)
colormap('gray')
title('estimated difference')

if alg.save
    wh=[6 1.8];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_params'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end

