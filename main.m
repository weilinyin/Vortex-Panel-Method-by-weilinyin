%计算主程序
tic;
alpha=9;
s=1000;
c=5;
[cp,l,cl,x]=vortex_panel_method(0012,alpha,c,s,10);
%plot(alpha,cl);

plot(x(1:s),cp(1:s));
hold;
plot(x(s+1:end),cp(s+1:end));

toc;