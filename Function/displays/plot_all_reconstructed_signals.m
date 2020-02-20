%% plot figure
close all ;figr=1; figure(figr);
cnt_lgd=1;lgnd{cnt_lgd}= 'Input signal';  plot(f, yf,'k','LineWidth',2.5); hold on
cnt_lgd=cnt_lgd+1; lgnd{cnt_lgd}= 'fastSCSA'; plot(f,yf_fast,'r','LineWidth',2.3); hold on
cnt_lgd=cnt_lgd+1; lgnd{cnt_lgd}= 'Grid search';plot(f,yf_grid,'b','LineWidth',2); hold on

legend(lgnd,'Location','northeast');
xlabel('Time (s)')
ylabel('Intensity')
set(gca,'YTickLabel',[],'XTickLabel',[])
%     xlim([t(1) t(end)])
% title('Signal reconstruction') 
set(gcf,'color','w') 
%     set(gca,'Xdir','reverse');
set(gca,'fontsize',16)
%     text(1.8,60,'NAA');text(1.40,40,'Lac(1)');text(1.2,30,'Lac(2)');text(4.7,30,'Water residue');
box 

%% save figure
name=strcat(tag_experimt,'_',num2str(noise_level),'_compare');
save_figure(Results_path,figr,name)