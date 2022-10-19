function Plot_DBIM(H,sub_H,image,iter_i,test_name,fctrs)
%H is the handle of figure, which decide in which figure to plot or create
%a new figure.
figure(H)
subplot(2,2,sub_H)
imagesc(squeeze(image'));
axis xy
axis square
%%%%%%%
Matn=sprintf('%s%d%s','iteration ',iter_i,' of EpsInf');
title(Matn);
caxis('auto')
colorbar;
axis image;

ax = axes('position',[0,0,1,1],'visible','off');
tx = text(0.4,0.95,[ test_name ' in the frequency of ' num2str(fctrs/1e9) ' GHz' ]);
set(tx,'fontweight','bold');
end