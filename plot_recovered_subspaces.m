function plot_recovered_subspaces(constants,est_params,alg)
%% plot example recovery image

figure(3), clf
nrows=3;
ncols=1;
i=0;

% true signal
i=i+1; subplot(nrows,ncols,i), cla
imagesc(est_params.d_opt)
colormap('gray')
title('naive Bayes signal subgraph')
% ylabel('truth')

% incoherent
i=i+1; subplot(nrows,ncols,i), cla
delsort=sort(est_params.d_pos(:),'descend');
delsparse=est_params.d_opt;
delsparse(delsparse<delsort(alg.num_inc_edges))=0;
imagesc(delsparse);
title('incoherent signal subgraph')

% max degree
i=i+1; subplot(nrows,ncols,i), cla
constants.Nmax=constants.n/10;
ind = get_max_edges(est_params.d_pos);
delmax=zeros(constants.n);
delmax(ind)=est_params.d_opt(ind);
imagesc(delmax);
title('coherent signal subgraph')


if alg.save
    wh=[2 6];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_subspaces'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end
