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
function[cp,l,cl,x]=vortex_panel_method_gpu(n,alpha,c,s,v)

[x_u, x_l, y_u, y_l]=naca(n,0,c,s+1,1,1);  
%构建几何形状

x_l=(fliplr(x_l))';
y_l=(fliplr(y_l))';
%下翼面坐标顺序倒置

x=gpuArray([x_u'; x_l(2:end-1)]);
y=gpuArray([y_u'; y_l(2:end-1)]);
%几何点数据合并



data=gpuArray.zeros(2*s-1,4);
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

J=gpuArray.zeros(2*s-1);
%J矩阵预分配内存

J=(1+((data(:,2)-(data(:,2))')./(data(:,1)-(data(:,1))')).^2).^(-1).*(sin(data(:,4)).*(data(:,2)-(data(:,2))')./(data(:,1)-(data(:,1))').^2+cos(data(:,4)).*(data(:,1)-(data(:,1))').^(-1));
%求出方向导数

b=v.*cos(pi/2-(data(:,4)-(alpha.*pi./180))).*(2*pi);
A=J.*(data(:,3))';
A(isnan(A)) = 0;
%构建方程组系数矩阵

b(s,:)=0;
A(s,:)=gpuArray.zeros(1,2*s-1);
A(s,s:s+1)=[1,1];
%库塔后缘条件，替代第s个方程

gamma=A\b;
%解出gamma

Gamma=gpuArray.zeros(2*s-1,size(alpha,2));
Gamma(1:end-1,:)=(gamma(1:end-1,:)+gamma(2:end,:))./2;
Gamma(end,:)=(gamma(end,:)+gamma(1,:))./2;
%去除噪声

l=gather(-1.23*v*sum(gamma.*data(:,3)));
%升力

cl=2*l/(1.23*v^2*c);
%升力系数

cp=gather(1-(Gamma./v).^2);
%压力系数

x=gather(data(:,1)./c);
%输出x/c

end