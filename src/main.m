%计算主程序
tic;

n=2412;%你想计算的翼型，4,5,6位NACA均可
alpha=[1:5];%输入你想要的攻角（°）
s=1000;%上、下翼面的板块数量
c=5;%弦长
v=10;%速度
[cp,l,cl,x]=vortex_panel_method(n,alpha,c,s,v);

f_1=figure(1);
plot(alpha,cl);%画cl-AoA二维图

f_2=figure(2);
f_2.Position(3)=f_2.Position(3)*size(alpha,2);
tiledlayout(1,size(alpha,2));
for i=1:size(alpha,2)
    nexttile;
    plot(x,cp(:,i));
    title(['alhpa=',num2str(alpha(i))]);
    xlabel('x/c');
    ylabel('cp');
end

toc;