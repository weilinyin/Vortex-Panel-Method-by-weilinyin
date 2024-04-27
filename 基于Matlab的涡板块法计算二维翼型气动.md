# 基于Matlab的涡板块法计算二维翼型气动

## 概述

### 涡板块法(Vortex Panel Method)

涡板块法是一种用于计算流体力学问题的数值方法，主要应用于机翼的气动力学分析。这种方法通过在固体物体的表面离散分布一系列的涡线，并假设每个涡线产生的涡旋在固体表面产生的诱导速度与实际流体的速度相匹配，从而近似地模拟流体的流动。

涡板块法的基本思想是将物体的表面分割成许多小板块，每个板块上分布有一个涡旋（涡线），这些涡旋在流场中产生诱导速度。根据库塔-儒科夫斯基定理，物体的升力可以通过计算物体表面涡旋所产生的诱导速度来求得。

## 核心算法

设 $P(x,y)$为流场中一点，则该点由第$j$个面板诱导的速度势为
$$
\Delta \phi_{j}=-\frac{1}{2 \pi} \int_{j} \theta_{p j} \gamma_{j} d s_{j}
\tag{1}
$$
其中，$\gamma_i$在第$j$个面板上，为仅取决于$j$的常数，而
$$
\theta_{pj}=\arctan\frac{y-y_j}{x-x_j} \tag{2}
$$
故，$P$点处由所有板块诱导的速度势为

$$
\phi(P)=-\sum_{j=1}^{n}\frac{\gamma_j}{2\pi}\int_j\theta_{pj}ds_j \tag{3}
$$
将$P$置于第$i$个板块的控制点上，则

$$
\phi(x_i,y_i)=-\sum_{j=1}^n\frac{\gamma_j}{2\pi}\int_j\theta_{ij}ds_j \tag{4}
$$
板块上速度对板块的法向分量为$0$，故
$$
V_\infty\cos\beta_i-\sum_{j=1}^n\frac{\gamma_j}{2\pi}\int_j \frac{\partial\theta_{ij}}{\partial n_i}ds_j=0\\ 
其中， \beta_i=\frac{\pi}{2}-(\alpha-\theta_i) \tag{5}
$$
式(5)中用控制点近似替代板块，则有

$$
V_\infty\cos\beta_i-\sum_{j=1}^n\frac{\gamma_j s_j}{2\pi}\frac{\partial\theta_{ij}}{\partial n_i}=0 \tag{6}
$$
由方向导数定义，有

$$
\frac{\partial \theta_{ij}}{\partial n_i}=\frac{1}{1+(\frac{y_i-y_j}{x_i-x_j})^2}\cdot[-\frac{y_i-y_j}{(x_i-x_j)^2}\frac{\partial x}{\partial n_i}+\frac{1}{x_i-x_j} \frac{\partial y}{\partial n_i}] \tag{7}
$$
其中，

$$
\frac{\partial x}{\partial n_i}=-\sin\theta_i \tag{8}
$$

$$
\frac{\partial x}{\partial n_i}=-\sin\theta_i \tag{9}
$$

当$i=j$ 时，

$$
\frac{\partial \theta_{ij}}{\partial n_i}=0 \tag{10}
$$
上述各式构成了一个关于$\gamma_i$的n元一次方程组

此外，翼型后缘应满足库塔后缘条件，在数值计算中用两个距离极近的板块替代后缘，这两个板块应满足

$$
\gamma_i=\gamma_{i-1} \tag{11}
$$
这使得方程数量增加为n+1，而未知量仅有n个，因此需要舍弃之前的一个方程。

构成方程组：
$$
\left[
\begin{matrix}
0 & s_2\frac{\partial \theta_{1\,2}}{\partial n_1} & \cdots\ & s_{n-1}\frac{\partial \theta_{1\,n-1}}{\partial n_1}&s_n\frac{\partial \theta_{1\,n}}{\partial n_1}\\
s_1\frac{\partial \theta_{2\,1}}{\partial n_2} & 0 & \cdots &  s_{n-1}\frac{\partial \theta_{2\,n-1}}{\partial n_2} &s_{n}\frac{\partial \theta_{2\,n}}{\partial n_2} \\
\vdots & \vdots & \ddots & \vdots & \vdots\\
s_{1}\frac{\partial \theta_{n-1\,1}}{\partial n_{n-1}} & s_{2}\frac{\partial \theta_{n-1\,2}}{\partial n_{n-1}} & \cdots & 0 & s_{n}\frac{\partial \theta_{n-1\,n}}{\partial n_{n-1}}\\
1 & 0 & \cdots & 0 & 1
\end{matrix}
\right]
\left[
\begin{matrix}
\gamma_1 \\
\gamma_2 \\
\vdots \\
\gamma_{n-1}\\
\gamma_n
\end{matrix}
\right]
=
\left[
\begin{matrix}
V_\infty2\pi\cos\beta_1\\
V_\infty2\pi\cos\beta_2\\
\vdots\\
V_\infty2\pi\cos\beta_{n-1}\\
V_\infty2\pi\cos\beta_n
\end{matrix}
\right]
\\ \tag{12}
$$
在求解出$\gamma_i$后，各板块处流速有
$$
V_{si}=-\gamma_i \tag{13}
$$
因此，升力由库塔-儒科夫斯基定理给出，

$$
L'=\rho_\infty V_\infty \Gamma = -\rho_\infty V_\infty\sum_{i=1}^n\gamma_i s_i
\tag{14}
$$

$$
c_l= \frac{L'}{\frac{1}{2}\rho_\infty V_\infty^2 c} \tag{15}
$$

各板块压力系数有， 
$$
C_{pi}=1-(\frac{V_{si}}{V_\infty})^2=1-(\frac{\gamma_i}{V_\infty})^2 \tag{16}
$$


## Matlab脚本实现



## NVIDIA GPU加速计算

## 总结
