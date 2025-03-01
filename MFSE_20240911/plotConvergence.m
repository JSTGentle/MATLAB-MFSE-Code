function plotConvergence(values , volfrac)
iter = 1:length(values);
%% 绘制柔度
axes1 = gca;
yyaxis(axes1,'left');
plot(iter, values,'b-','LineWidth',1.5);
ylabel('Structural compliance','FontSize',14,'FontName','Times New Roman','Color',[0 0 1]);
set(axes1,'YColor','b','FontSize',14,'FontName','Times New Roman');

%% Y right 加体分比
yyaxis(axes1,'right');
plot(iter,volfrac,'r-.','LineWidth',1.5);hold on;
set(axes1,'ylim',[0 1]);set(axes1,'ytick',0:.1:1);
ylabel('Volume constraint','FontSize',14,'FontName','Times New Roman','Color',[1 0 0]);
set(axes1,'YColor','r','FontSize',14,'FontName','Times New Roman');

%% 坐标轴
xlabel('Number of iterations','FontSize',14,'FontName','Times New Roman');
drawnow;
end