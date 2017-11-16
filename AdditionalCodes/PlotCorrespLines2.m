function PlotCorrespLines2(I1,I2,Pt1,Pt2,L1,L2,corr)


[r c]=size(corr);
if size(Pt1,1)==6
    Pt1=Pt1([1 2 4 5],:);
    Pt2=Pt2([1 2 4 5],:);
end


% cc: is just for plotting with multi colors
cc=lines(size(corr,2));
% ..................

p1=Pt1(1:2,:);  p2=Pt1(3:4,:);
p3=Pt2(1:2,:);  p4=Pt2(3:4,:);

figure;

[figsize]= get(0,'ScreenSize');
figsize = figsize +[10 10 -100 -100];
set(gcf,'Position',figsize)

%winsize = get(gcf,'Position');

%winsize(1:2) = [0 0];
% Mov = moviein(gcf,winsize);
% set(fig1,'NextPlot','replacechildren')


for i=1:r
    subplot(1,2,1),imshow(I1);title('\bf{Image 1}','fontsize',10); hold on;
    subplot(1,2,2),imshow(I2);title('\bf{Image 2}','fontsize',10); hold on;

    subplot(1,2,1),plot([p1(1,corr(i,1)); p2(1,corr(i,1))],[p1(2,corr(i,1)); p2(2,corr(i,1))],'g','LineWidth',3),...
                   text(p1(1,corr(i,1)),p1(2,corr(i,1)),num2str(corr(i,1)),'color','r','FontSize',20);
                   hline(L1(:,corr(i,1)),'y');
    for j=2:size(corr,2)
        subplot(1,2,2),plot([p3(1,corr(i,j)); p4(1,corr(i,j))],[p3(2,corr(i,j)); p4(2,corr(i,j))],'color',cc(j,:),'LineWidth',3),...
                       text(p3(1,corr(i,j)),p3(2,corr(i,j)),num2str(j-1),'color',cc(j,:),'FontSize',16)
    end
       
  pause;
  hold off;
 
end

end

