%/*梯度投影法求7-DOF机械臂逆解
% *版本：Matlab 2019b
% *2020.09.17
% *Written by Rot_Tianers*/

clear;clc;set(0,'defaultfigurecolor','w');
%%设定直线初始点和末端点（曲线类型可自由修改,注意始末点的选取要在工作空间之内）
Initial = [-0.8025,0,0.2767];
End = [-0.6551,0.2,0.6037];
%%设置步长
deta = 0.01;
%%求直线长度
lline = sqrt(((Initial(1)-End(1))^2+(Initial(2)-End(2))^2+(Initial(3)-End(3))^2));
%%迭代步数
if rem(lline,deta) == 0
    STEP = lline/deta;
else
    STEP = floor(lline/deta)+1;
end
posx(1) = Initial(1);
posy(1) = Initial(2);
posz(1) = Initial(3);
for i=1:STEP
    posx(i+1) = Initial(1)+(End(1)-Initial(1))/STEP*i;
    posy(i+1) = Initial(2)+(End(2)-Initial(2))/STEP*i;
    posz(i+1) = Initial(3)+(End(3)-Initial(3))/STEP*i;
end
%%=========预期轨迹===================================================%%

%%iiwa机械臂DH参数
BS = 0.34;SE = 0.4;EW = 0.4;WT = 0.1266;
L(1) = Link('d', BS, 'a', 0, 'alpha', 0,'modified');
L(2) = Link('d', 0, 'a', 0, 'alpha', pi/2,'modified');
L(3) = Link('d', SE, 'a', 0, 'alpha', -pi/2,'modified');
L(4) = Link('d', 0, 'a', 0, 'alpha', pi/2,'modified');
L(5) = Link('d', EW, 'a', 0, 'alpha', -pi/2,'modified');
L(6) = Link('d', 0, 'a', 0, 'alpha', pi/2,'modified');
L(7) = Link('d', WT, 'a', 0, 'alpha', -pi/2,'modified');
bot = SerialLink([L(1) L(2) L(3) L(4) L(5) L(6) L(7)],'name','KUKA iiwa7');
%%初始关节角
Initial_q = [0 pi/3 0 pi/3 0 0 0];
dx = [0;0;0;0;0;0];

%% 设置关节极限范围
qmax = deg2rad([170 120 170 120 170 120 175]);
qmin = deg2rad([-170 -120 -170 -120 -170 -120 -175]);

%%姿态变化
Ori_ini = [0,0,0];
%%终止点位姿
Ori_end = [0,pi/3,0];

%%赋值准备迭代
q = Initial_q;
Theta(1,:) = rad2deg(q);
for i=1:STEP
    %%位置微分步
    dx(1) = posx(i+1)-posx(i);
    dx(2) = posy(i+1)-posy(i);
    dx(3) = posz(i+1)-posz(i);
    %%姿态微分步（姿态等量变化即可）
    dx(4) = (Ori_end(1)-Ori_ini(1))/STEP;
    dx(5) = (Ori_end(2)-Ori_ini(2))/STEP;
    dx(6) = (Ori_end(3)-Ori_ini(3))/STEP;
    %%求雅可比矩阵
    J = bot.jacob0(q);
    sigma = 0.02;  %%最小奇异值边界,取值范围要尽量接近0
    lambda0 = 5;   %%阻尼系数
    %%阻尼最小二乘法避奇异点
    sv = svd(J);
    sor = sort(sv); %%求雅可比矩阵的奇异值
    theta = sor(1); %%得到最小奇异值
    %%变阻尼系数的阻尼最小二乘法
    if theta < sigma
        lambda=lambda0^2*(1-(theta/sigma)^2);
    else
        lambda=0;
    end
    MPJ = J'/(J*J'+lambda*eye(6,6));   %%MPJ是采用阻尼最小二乘法得出的伪逆（也成为奇异鲁棒性逆），当lanmda2=0时，退化为普通的伪逆法
    %%避极限指标
    I = eye(7,7);
    K = -10;  %%放大系数，目标函数取最小时为负，取最大时为正
    for j = 1:7
       a(j) = (qmax(j)+qmin(j))/2;  %%第i个关节中值
       dH(j) = 2*(q(j)-a(j))/(7*(a(j)-qmax(j))^2);  %%避关节极限优化函数的梯度
    end
    %%梯度投影法：零空间优化指标
    NJ = K*(I-MPJ*J)*dH';  %%dH是7×1列向量
    dq = MPJ*dx+NJ;
    q = q+dq';
    %%记录每一次关节角
    Theta(i+1,:) = rad2deg(q);
end

%%画七个关节曲线
plot(1:STEP+1,Theta(:,1)','-d','Markersize',8,'linewidth',1.5);
hold on;
plot(1:STEP+1,Theta(:,2)','-o','Markersize',8,'linewidth',1.5);
hold on;
plot(1:STEP+1,Theta(:,3)','-*','Markersize',8,'linewidth',1.5);
hold on;
plot(1:STEP+1,Theta(:,4)','--h','Markersize',8,'linewidth',1.5);
hold on;
plot(1:STEP+1,Theta(:,5)',':x','Markersize',8,'linewidth',1.5);
hold on;
plot(1:STEP+1,Theta(:,6)',':s','Markersize',8,'linewidth',1.5);
hold on;
plot(1:STEP+1,Theta(:,7)',':p','Markersize',8,'linewidth',1.5);
hold on;
legend('\theta1','\theta2','\theta3','\theta4','\theta5','\theta6','\theta7','orientation','horizontal');

%%验证末端坐标值
bot.fkine(q);

%%可视化始末两个状态的机械臂
figure;
robot = importrobot('iiwa7.urdf');
robot.DataFormat='column';
%%初始关节角
q0 = [Initial_q(1);Initial_q(2);Initial_q(3);Initial_q(4);Initial_q(5);Initial_q(6);Initial_q(7)];
show(robot,q0,'PreservePlot',true);
hold on;
%%末端关节角
qn = [q(1);q(2);q(3);q(4);q(5);q(6);q(7)];
show(robot,qn,'PreservePlot',true);
%%图片显示的一些设置
axis([-0.6 0.6 -0.6 0.6 0 1.35]);%%坐标范围
box on;
camva('auto');  %%设置相机视角
daspect([1 1 1]); %%控制沿每个轴的数据单位长度,要在所有方向上采用相同的数据单位长度，使用 [1 1 1] 
set(gca,'FontSize',20,'Fontname', 'Times New Roman','linewidth',1.2);
xlabel('x(m)','FontSize',20);
ylabel('y(m)','FontSize',20);
zlabel('z(m)','FontSize',20);
set(gca,'XColor','k','YColor','k');
%%程序结束
