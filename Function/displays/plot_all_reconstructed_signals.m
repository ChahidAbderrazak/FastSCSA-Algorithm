%% plot figure

cnt_lgd=1;lgnd{cnt_lgd}= 'Input signal'; figure; plot(f, yf,'k'); hold on
cnt_lgd=cnt_lgd+1; lgnd{cnt_lgd}= 'fastSCSA'; plot(yf_fast,'r','LineWidth',2.5); hold on
cnt_lgd=cnt_lgd+1; lgnd{cnt_lgd}= 'Grid search';plot(yf_grid,'b','LineWidth',2); hold on

legend(lgnd,'Location','northeast');
xlabel('t')
ylabel('Intensity')
set(gca,'YTickLabel',[],'XTickLabel',[])
%     xlim([t(1) t(end)])
title('Signal reconstruction') 
set(gcf,'color','w') 
%     set(gca,'Xdir','reverse');
set(gca,'fontsize',16)
%     text(1.8,60,'NAA');text(1.40,40,'Lac(1)');text(1.2,30,'Lac(2)');text(4.7,30,'Water residue');
box 
