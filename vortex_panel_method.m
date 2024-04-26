%% 涡板块法计算函数
%输出数据
%cp 压强系数
%l  升力
%cl 升力系数
%m  力矩
%cm 力矩系数
%输入数据
%n  翼型号（4 5 6位）
%alpha 攻角（°）,支持数组
%c  弦长(m)
%s  板块数的一半
%v  来流速度
function[cp,l,cl,x]=vortex_panel_method(n,alpha,c,s,v)

[x_u, x_l, y_u, y_l]=naca(n,0,c,s+1,1,1);  
%构建几何形状

x_l=(fliplr(x_l))';
y_l=(fliplr(y_l))';
%下翼面坐标顺序倒置

x=[x_u'; x_l(2:end-1)];
y=[y_u'; y_l(2:end-1)];
%几何点数据合并



data=zeros(2*s-1,4);
%   预分配内存
%   x 控制点x坐标 1
%   y 控制点y坐标 2
%   length of panel 3
%   theta (rad) 板块与x轴夹角 4
%   每行为一个板块的数据

data(:,1)=(x(1:end-1)+x(2:end))./2;
data(:,2)=(y(1:end-1)+y(2:end))./2;
%根据翼型几何数据求出涡板块控制点坐标

data(:,3)=sqrt((x(1:end-1)-x(2:end)).^2+(y(1:end-1)-y(2:end)).^2);
%求出板块长度

data(:,4)=atan((y(1:end-1)-y(2:end))./(x(1:end-1)-x(2:end)));
%求出板块与x轴夹角

J=zeros(2*s-1);
%J矩阵预分配显存

x_con=(data(:,1))';
y_con=(data(:,2))';
%控制点x y坐标向量化

J=(1+((data(:,2)-y_con)./(data(:,1)-x_con)).^2).^(-1).*(sin(data(:,4)).*(data(:,2)-y_con)./(data(:,1)-x_con).^2+cos(data(:,4)).*(data(:,1)-x_con).^(-1));
%求出Jij

b=v.*cos(pi/2-(data(:,4)-(alpha.*pi./180))).*(2*pi);
A=J.*(data(:,3))';

A(isnan(A)) = 0;
%方程组系数矩阵

b(s,:)=0;
A(s,:)=zeros(1,2*s-1);
A(s,s:s+1)=[1,1];
%库塔后缘条件，舍去第s个方程

gamma=A\b;
%解出gamma

Gamma=zeros(2*s-1,size(alpha,2));
Gamma(1:end-1,:)=(gamma(1:end-1,:)+gamma(2:end,:))./2;
Gamma(end,:)=(gamma(end,:)+gamma(1,:))./2;
%去除噪声

l=-1.23*v*sum(gamma.*data(:,3));
%升力

cl=2*l/(1.23*v^2*c);
%升力系数

cp=1-(Gamma./v).^2;
for i=2:(2*s-1)
    if abs(cp(i))>10
        cp(i)=cp(i-1);
    end
end
%压力系数，抑制尖峰

x=data(:,1)./c;
%输出x/c

end