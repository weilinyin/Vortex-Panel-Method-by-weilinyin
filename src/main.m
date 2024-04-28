%计算主程序
tic;

n=2412;%你想计算的翼型，4,5,6位NACA均可
alpha=[5,6];%攻角 °
s=1000;%上、下翼面的板块数量
c=5;%弦长 m
v=10;%速度 m/s
isgpu=1;%是否使用gpu求解
if isgpu == 1
    [cp,l,cl,x]=vortex_panel_method_gpu(n,alpha,c,s,v);
else
    [cp,l,cl,x]=vortex_panel_method(n,alpha,c,s,v);
end
f_2=figure(2);
plot(alpha,cl);%画cl-AoA图
title('cl-AoA');
xlabel('alpha (°)');
ylabel('cl');
grid on;

f_3=figure(3);%画cp-x图
f_3.Position(3)=f_3.Position(3)*size(alpha,2);
tiledlayout(1,size(alpha,2));
for i=1:size(alpha,2)
    nexttile;
    plot(x(1:s),cp(1:s,i));
    hold on;
    plot(x(s+1:end),cp(s+1:end,i));
    ax=gca;
    ax.YDir='reverse';%y轴翻转
    grid on;
    xlim([0,1]);
    ylim([-2.5,1.5]);
    title(['alpha=',num2str(alpha(i))]);
    xlabel('x/c');
    ylabel('cp');
end
legend('upper','lower','Location','northeast');

toc;