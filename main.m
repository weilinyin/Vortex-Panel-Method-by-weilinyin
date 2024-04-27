%计算主程序
tic;
alpha=[-5:16];
s=1000;
c=5;
[cp,l,cl,x]=vortex_panel_method(2412,alpha,c,s,10);
plot(alpha,cl);

%[X,Y]=meshgrid(alpha,x);
%s=surf(X,Y,cp,'FaceAlpha',0.5);
%s.EdgeColor= 'none';


toc;