if alg.save==1
    figname=[alg.figdir alg.fname alg.ext];
    wh=[6 3];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    print('-dpdf',figname)
    print('-dpng',figname)
    saveas(gcf,figname)
    save([alg.postdir alg.fname alg.ext])
end