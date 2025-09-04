clear
clc
close all
%% Laplacian Matrix
L1=[
 1/2   0    0   -1/2	0    0;    
-1/4  3/4 -1/4	-1/4	0    0;    
-1/4   0   3/4	-1/4	0  -1/4;    
  0  -1/3 -1/3	 2/3	0    0;
  0  -1/2   0     0    1/2   0;
  0    0    0     0   -1/2  1/2];
L2=[
 1/2   0    0   -1/2	0    0;  
-1/4  3/4 -1/4	-1/4	0    0;
-1/3   0   2/3	-1/3	0    0;    
  0  -1/3 -1/3	 2/3	0    0;    
  0    0    0     0     0    0;
  0    0    0     0     0    0];
L3=[
 1/2   0    0   -1/2	0    0;    
-1/4  3/4 -1/4	-1/4	0    0;    
-1/4   0   3/4	-1/4  -1/4   0;    
  0  -1/3 -1/3	 2/3	0    0;    
  0  -1/2   0     0    1/2   0;
  0    0    0     0     0    0];
%% Surplus Matrix
S1=[ 
 1/3   0    0    1/4    0    0;    
 1/3  1/3  1/3	 1/4    0    0;    
 1/3   0   1/3	 1/4    0   1/2;    
  0   1/3  1/3	 1/4    0    0;    
  0   1/3   0     0    1/2   0;
  0    0    0     0    1/2  1/2];
S2=[ 
 1/3   0    0    1/4    0    0;    
 1/3  1/2  1/3	 1/4    0    0;    
 1/3   0   1/3	 1/4    0    0;    
  0   1/2  1/3	 1/4    0    0;    
  0    0    0     0     0    0;
  0    0    0     0     0    0];
S3=[ 
 1/3   0    0    1/4    0    0;    
 1/3  1/3  1/3	 1/4    0    0;    
 1/3   0   1/3	 1/4   1/2   0;    
  0   1/3  1/3	 1/4    0    0;    
  0   1/3   0     0    1/2   0;
  0    0    0     0     0    0];
E1=[
1	0	0	0	0   0;
0	1	0	0	0   0;
0	0	1	0	0   0;
0	0	0	1	0   0;
0	0	0	0	1   0;
0	0	0	0	0   1
    ];
E2=[
1	0	0	0	0   0;
0	1	0	0	0   0;
0	0	1	0	0   0;
0	0	0	1	0   0;
0	0	0	0	0   0;
0	0	0	0	0   0
    ];
E3=[
1	0	0	0	0   0;
0	1	0	0	0   0;
0	0	1	0	0   0;
0	0	0	1	0   0;
0	0	0	0	1   0;
0	0	0	0	0   0
    ];
% gains
alpha = 0.5;
beta = 1.0;
gama = 0.1;
psi = 0.01;
varphi = 0.2;
epsilon = 0.15;
 
%% Initial States
x_x(:,:) = 0.1*[5 3 -4 -3 -1 -2]';
x_y(:,:) = 0.1*[2 3 -3  2 -1 -2]';
v_x(:,:) = [0 0 0 0 0 0]';
v_y(:,:) = [0 0 0 0 0 0]';
z_x(:,:) = 0.1*[1 -1 2 -1 1 2]';
z_y(:,:) = 0.1*[2 -1 2 -1 2 1]';
p_x(:,:) = [ 0 0 0 0 0 0]';
p_y(:,:) = [ 0 0 0 0 0 0]';
q_x(:,:) = [ 0 0 0 0 0 0]';
q_y(:,:) = [ 0 0 0 0 0 0]'; 
s_x(:,:) = [ 0 0 0 0 0 0]';
s_y(:,:) = [ 0 0 0 0 0 0]';
%u1(:,:) = [1 2 3 4 0 0 ]';
z1_x=0.1*[0;0;0;0;-1;0];
z1_y=0.1*[0;0;0;0;-1;0];
z2_x=0.1*[0;0;0;0;1;0];
z2_y=0.1*[0;0;0;0;2;0];
%% Time Parameters
tBegin = 0;
tFinal1 = 100;
tFinal2 = 200;
tFinal3 = 300;
dT = 1;
 h = 0.1;
times1 = (tFinal1-tBegin)/dT;
times2 = (tFinal2-tBegin)/dT;
times3 = (tFinal3-tBegin)/dT;
t(1,1) = 0;

%% Iteration Calculate
for k=1:times1
    % record time
    t(:,k+1) = t(:,k) + dT;
    
    % calculate control inputs
    u1_x(:,k) =(gama+1) * [-alpha*L1 -beta*L1] * [x_x(:,k); z_x(:,k)]+ E1* gama * z_x(:,k) + gama * E1 * epsilon*s_x(:,k);
    u1_y(:,k) =(gama+1) * [-alpha*L1 -beta*L1] * [x_y(:,k); z_y(:,k)]+ E1* gama * z_y(:,k) + gama * E1 * epsilon*s_y(:,k);
    % update statues
    x_x(:,k+1)= E1 * x_x(:,k)+E1 * dT*v_x(:,k)+E1 * 1/2*dT*dT*u1_x(:,k)+E1 * psi*p_x(:,k);
    v_x(:,k+1)= E1 * v_x(:,k)+E1 * dT*u1_x(:,k)+E1 * varphi*q_x(:,k);
    z_x(:,k+1)= E1 * z_x(:,k)+E1 * dT*[-alpha*L1 -beta*L1] * [x_x(:,k); z_x(:,k)]+E1 * epsilon*s_x(:,k);
    p_x(:,k+1)= (S1-psi*E1)*p_x(:,k)-E1* 1/2 * (gama+1) * [-alpha*L1 -beta*L1] * [x_x(:,k); z_x(:,k)]+E1*q_x(:,k)+E1*1/2*gama*s_x(:,k);
    q_x(:,k+1)= (S1-varphi*E1)*q_x(:,k)-E1 * (gama+1) * [-alpha*L1 -beta*L1] * [x_x(:,k); z_x(:,k)]+E1*gama*s_x(:,k);
    s_x(:,k+1)= (S1-epsilon*E1)*s_x(:,k)-E1*[-alpha*L1 -beta*L1] * [x_x(:,k); z_x(:,k)];

    x_y(:,k+1)= E1 * x_y(:,k)+E1 * dT*v_y(:,k)+E1 * 1/2*dT*dT*u1_y(:,k)+E1 * psi*p_y(:,k);
    v_y(:,k+1)= E1 * v_y(:,k)+E1 * dT*u1_y(:,k)+E1 * varphi*q_y(:,k);
    z_y(:,k+1)= E1 * z_y(:,k)+E1 * dT*[-alpha*L1 -beta*L1] * [x_y(:,k); z_y(:,k)]+E1 * epsilon*s_y(:,k);
    p_y(:,k+1)= (S1-psi*E1)*p_y(:,k)-E1* 1/2 * (gama+1) * [-alpha*L1 -beta*L1] * [x_y(:,k); z_y(:,k)]+E1*q_y(:,k)+E1*1/2*gama*s_y(:,k);
    q_y(:,k+1)= (S1-varphi*E1)*q_y(:,k)-E1 * (gama+1) * [-alpha*L1 -beta*L1] * [x_y(:,k); z_y(:,k)]+E1*gama*s_y(:,k);
    s_y(:,k+1)= (S1-epsilon*E1)*s_y(:,k)-E1*[-alpha*L1 -beta*L1] * [x_y(:,k); z_y(:,k)];
    
end
t1=0:0.1:(tFinal1)*h;
t1_big=0:0.1:10*h;
t2=0:0.1:tFinal1*h-h;
z3_x=[0;0;0;0;x_x(5:5,end);x_x(6:6,end)];
z3_y=[0;0;0;0;x_y(5:5,end);x_y(6:6,end)];
z4_x=[0;0;0;0;z_x(5:5,end);z_x(6:6,end)];
z4_y=[0;0;0;0;z_y(5:5,end);z_y(6:6,end)];
x_x(:,end)=x_x(:,end)+(L2-L1)*z3_x;
x_y(:,end)=x_y(:,end)+(L2-L1)*z3_y;
z_x(:,end)=z_x(:,end)+(L2-L1)*z4_x;
z_y(:,end)=z_y(:,end)+(L2-L1)*z4_y;
for k=times1+1:times2
    % record time
    t(:,k+1) = t(:,k) + dT;
  
    % calculate control inputs
    u1_x(:,k) =(gama+1) * [-alpha*L2 -beta*L2] * [x_x(:,k); z_x(:,k)]+ E2* gama * z_x(:,k) + gama * E2 * epsilon*s_x(:,k);
    u1_y(:,k) =(gama+1) * [-alpha*L2 -beta*L2] * [x_y(:,k); z_y(:,k)]+ E2* gama * z_y(:,k) + gama * E2 * epsilon*s_y(:,k);
    % update statues
    x_x(:,k+1)= E2*x_x(:,k)+E2*dT*v_x(:,k)+E2*1/2*dT*dT*u1_x(:,k)+E2*psi*p_x(:,k);
    v_x(:,k+1)= E2*v_x(:,k)+E2*dT*u1_x(:,k)+E2*varphi*q_x(:,k);
    z_x(:,k+1)= E2*z_x(:,k)+E2*dT*[-alpha*L2 -beta*L2] * [x_x(:,k); z_x(:,k)]+E2*epsilon*s_x(:,k);
    p_x(:,k+1)= (S2-psi*E2)*p_x(:,k)-E2* 1/2 * (gama+1) * [-alpha*L2 -beta*L2] * [x_x(:,k); z_x(:,k)]+E2*q_x(:,k)+E2*1/2*gama*s_x(:,k);
    q_x(:,k+1)= (S2-varphi*E2)*q_x(:,k)-E2 * (gama+1) * [-alpha*L2 -beta*L2] * [x_x(:,k); z_x(:,k)]+E2*gama*s_x(:,k);
    s_x(:,k+1)= (S2-epsilon*E2)*s_x(:,k)-E2*[-alpha*L2 -beta*L2] * [x_x(:,k); z_x(:,k)];

    x_y(:,k+1)= E2*x_y(:,k)+E2*dT*v_y(:,k)+E2*1/2*dT*dT*u1_y(:,k)+E2*psi*p_y(:,k);
    v_y(:,k+1)= E2*v_y(:,k)+E2*dT*u1_y(:,k)+E2*varphi*q_y(:,k);
    z_y(:,k+1)= E2*z_y(:,k)+E2*dT*[-alpha*L2 -beta*L2] * [x_y(:,k); z_y(:,k)]+E2*epsilon*s_y(:,k);
    p_y(:,k+1)= (S2-psi*E2)*p_y(:,k)-E2* 1/2 * (gama+1) * [-alpha*L2 -beta*L2] * [x_y(:,k); z_y(:,k)]+E2*q_y(:,k)+E2*1/2*gama*s_y(:,k);
    q_y(:,k+1)= (S2-varphi*E2)*q_y(:,k)-E2 * (gama+1) * [-alpha*L2 -beta*L2] * [x_y(:,k); z_y(:,k)]+E2*gama*s_y(:,k);
    s_y(:,k+1)= (S2-epsilon*E2)*s_y(:,k)-E2*[-alpha*L2 -beta*L2] * [x_y(:,k); z_y(:,k)];
end
t3=(tFinal1)*h+h:h:tFinal2*h;
t_leave=tFinal1*h+h:h:110*h-h;
% t4=(tFinal1)*h+1:h:tFinal2*h-1;
x_x(:,end)=x_x(:,end)+z1_x;
x_y(:,end)=x_y(:,end)+z1_y;
z_x(:,end)=z_x(:,end)+z2_x;
z_y(:,end)=z_y(:,end)+z2_y;
for k=times2+1:times3
    % record time
    t(:,k+1) = t(:,k) + dT;
  
    % calculate control inputs
    u1_x(:,k) =(gama+1) * [-alpha*L3 -beta*L3] * [x_x(:,k); z_x(:,k)]+ E3* gama * z_x(:,k) + gama * E3 * epsilon*s_x(:,k);
    u1_y(:,k) =(gama+1) * [-alpha*L3 -beta*L3] * [x_y(:,k); z_y(:,k)]+ E3* gama * z_y(:,k) + gama * E3 * epsilon*s_y(:,k);
    
    % update statues
    x_x(:,k+1)= E3*x_x(:,k)+E3*dT*v_x(:,k)+E3*1/2*dT*dT*u1_x(:,k)+E3*psi*p_x(:,k);
    v_x(:,k+1)= E3*v_x(:,k)+E3*dT*u1_x(:,k)+E3*varphi*q_x(:,k);
    z_x(:,k+1)= E3*z_x(:,k)+E3*dT*[-alpha*L3 -beta*L3] * [x_x(:,k); z_x(:,k)]+E3*epsilon*s_x(:,k);
    p_x(:,k+1)= (S3-psi*E3)*p_x(:,k)-E3* 1/2 * (gama+1) * [-alpha*L3 -beta*L3] * [x_x(:,k); z_x(:,k)]+E3*q_x(:,k)+E3*1/2*gama*s_x(:,k);
    q_x(:,k+1)= (S3-varphi*E3)*q_x(:,k)-E3 * (gama+1) * [-alpha*L3 -beta*L3] * [x_x(:,k); z_x(:,k)]+E3*gama*s_x(:,k);
    s_x(:,k+1)= (S3-epsilon*E3)*s_x(:,k)-E3*[-alpha*L3 -beta*L3] * [x_x(:,k); z_x(:,k)];

    x_y(:,k+1)= E3*x_y(:,k)+E3*dT*v_y(:,k)+E3*1/2*dT*dT*u1_y(:,k)+E3*psi*p_y(:,k);
    v_y(:,k+1)= E3*v_y(:,k)+E3*dT*u1_y(:,k)+E3*varphi*q_y(:,k);
    z_y(:,k+1)= E3*z_y(:,k)+E3*dT*[-alpha*L3 -beta*L3] * [x_y(:,k); z_y(:,k)]+E3*epsilon*s_y(:,k);
    p_y(:,k+1)= (S3-psi*E3)*p_y(:,k)-E3* 1/2 * (gama+1) * [-alpha*L3 -beta*L3] * [x_y(:,k); z_y(:,k)]+E3*q_y(:,k)+E3*1/2*gama*s_y(:,k);
    q_y(:,k+1)= (S3-varphi*E3)*q_y(:,k)-E3 * (gama+1) * [-alpha*L3 -beta*L3] * [x_y(:,k); z_y(:,k)]+E3*gama*s_y(:,k);
    s_y(:,k+1)= (S3-epsilon*E3)*s_y(:,k)-E3*[-alpha*L3 -beta*L3] * [x_y(:,k); z_y(:,k)];
end
t4=(tFinal2)*h+h:h:tFinal3*h;
t5=(tFinal2)*h+h:h:tFinal3*h-h;
t11 = [t1,t3,t4];
t12 = [t1,t3];
t22 = t2;
t33 = [t1,t3,t5];
t1leave = [t1,t_leave];
% x_x(5:6,times1+1)=0; 
% v_x(5:6,times1+1)=0;
% z_x(5:6,times1+1)=0;
% p_x(5:6,times1+1)=0; 
% q_x(5:6,times1+1)=0;
% s_x(5:6,times1+1)=0;
% rowSums_x_x = sum(x_x)+sum(p_x);
% rowSums_v_x = sum(v_x)+sum(q_x);
% rowSums_z_x = sum(z_x)+sum(s_x);
% results_x_x = [];
% results_v_x = [];
% results_z_x = [];
% ave_x_x = [];
% ave_v_x = [];
% ave_z_x = [];
% for s=1:times1
%     total_x_x = x_x(:,1)+(s-1)*(v_x(:,1))+gama*((s-1)^2+1/2*(s-1))*(z_x(:,1));
%     total_v_x = v_x(:,1)+gama*(s-1)*(z_x(:,1));
%     total_z_x = z_x(:,1);
%     results_x_x = [results_x_x, sum(total_x_x)];
%     results_v_x = [results_v_x, sum(total_v_x)];
%     results_z_x = [results_z_x, sum(total_z_x)];
% %     ave_x_x = [ave_x_x, sum(total_x_x)/6];
% %     ave_v_x = [ave_v_x, sum(total_v_x)/6];
% end
% for s=times1+1:times2
%     total_x_x = x_x(:,times1+1)+(s-times1-1)*(v_x(:,times1+1))+gama*((s-times1-1)^2+1/2*(s-times1-1))*(z_x(:,times1+1));
%     total_v_x = v_x(:,times1+1)+gama*(s-times1-1)*(z_x(:,times1+1));
%     total_z_x = z_x(:,times1+1);
%     results_x_x = [results_x_x, sum(total_x_x)];
%     results_v_x = [results_v_x, sum(total_v_x)];
%     results_z_x = [results_z_x, sum(total_z_x)];
% %     ave_x_x = [ave_x_x, sum(total_x_x)/4];
% %     ave_v_x = [ave_v_x, sum(total_v_x)/4];
% end
% for s=times2+1:times3+1
%     total_x_x = x_x(:,times2+1)+(s-times2-1)*(v_x(:,times2+1))+gama*((s-times2-1)^2+1/2*(s-times2-1))*(z_x(:,times2+1));
%     total_v_x = v_x(:,times2+1)+gama*(s-times2-1)*(z_x(:,times2+1));
%     total_z_x = z_x(:,times2+1);
%     results_x_x = [results_x_x, sum(total_x_x)];
%     results_v_x = [results_v_x, sum(total_v_x)];
%     results_z_x = [results_z_x, sum(total_z_x)];
% %     ave_x_x = [ave_x_x, sum(total_x_x)/5];
% %     ave_v_x = [ave_v_x, sum(total_v_x)/5];
% end
% disp(results_x);
% disp(results_v);
x_x(5:5,length(t22)+1:length(t12)-1)=NaN; 
v_x(5:5,length(t22)+1:length(t12)-1)=NaN;
z_x(5:5,length(t22)+1:length(t12)-1)=NaN;
p_x(5:5,length(t22)+1:length(t12)-1)=NaN; 
q_x(5:5,length(t22)+1:length(t12)-1)=NaN;
s_x(5:5,length(t22)+1:length(t12)-1)=NaN;
x_y(5:5,length(t22)+1:length(t12)-1)=NaN; 
v_y(5:5,length(t22)+1:length(t12)-1)=NaN;
z_y(5:5,length(t22)+1:length(t12)-1)=NaN;
p_y(5:5,length(t22)+1:length(t12)-1)=NaN; 
q_y(5:5,length(t22)+1:length(t12)-1)=NaN;
s_y(5:5,length(t22)+1:length(t12)-1)=NaN;
%% Draw graphs
% figure(1)
% plot(t11,results_v_x(1,:),'linewidth',1.5)
% hold on;
% plot(t11,rowSums_v_x(1,:),'linewidth',1.5,'LineStyle','--')
% xlabel('$$k$$','Interpreter','latex')
% ylabel('Sum of Values','Interpreter','latex')
% h=legend('$\sum(\bar{x}_i+k\bar{v}_i)$','$\sum(x_{i,k}+p_{i,k})$');
% set(h,'interpreter','latex')
% grid on; set(gca,'GridLineStyle',':','GridColor',[0 0 0],'GridAlpha',1)
% figure(2)
% plot(t11,results_z_x(1,:),'linewidth',1.5)
% hold on;
% plot(t11,rowSums_z_x(1,:),'linewidth',1.5,'LineStyle','--')
% xlabel('$$k$$','Interpreter','latex')
% ylabel('Sum of Values','Interpreter','latex')
% h=legend('$\sum\bar{v}_i$','$\sum(v_{i,k}+q_{i,k})$');
% set(h,'interpreter','latex')
% grid on; set(gca,'GridLineStyle',':','GridColor',[0 0 0],'GridAlpha',1)

figure(1)
% z=tiledlayout(2,1,'TileSpacing','tight','Padding','compact');
% nexttile
plot(t11,x_x(1:5,:),'linewidth',1.5)
hold on;
plot(t22,x_x(6:6,1:length(t22)),'linewidth',1.5)
% plot(t11,ave_x(1,:),'linewidth',1.5,'Color',[0.5, 0.5, 0.5],'LineStyle','--')
xline(10,'--b','LineWidth',1.5,'Label','UAVs 5&6 leave','FontName','Times New Roman');
xline(20,'--r','LineWidth',1.5,'Label','UAV 5 rejoins','FontName','Times New Roman');
  xlabel('$k$','Interpreter','latex')
 ylabel('Position $\xi^k_{i_x}$','Interpreter','latex','FontName','Times New Roman')
 h=legend('$\xi^k_{1_x}$','$\xi^k_{2_x}$','$\xi^k_{3_x}$','$\xi^k_{4_x}$','$\xi^k_{5_x}$','$\xi^k_{6_x}$');
set(h,'interpreter','latex')
grid on; set(gca,'GridLineStyle',':','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
% 插入局部放大 (离开)
    ax1 = axes('Position',[0.44 0.56 0.2 0.2]);
    box on; hold on;
    plot(t_leave,x_x(:,length(t1):length(t1leave)-1),'LineWidth',1.5);
    xlim([10+0.1,11-0.1]); 
    ylim([30,50]);
    set(gca,'GridLineStyle','none','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
    grid on;
    title('Zoom near leave','FontName','Times New Roman');

    % 插入局部放大 (加入)
    ax2 = axes('Position',[0.16 0.32 0.2 0.2]);
    box on; hold on;
        plot(t1_big,x_x(:,1:length(t1_big)),'LineWidth',1.5);
    xlim([0,1]); 
    ylim([-1,1]);
    set(gca,'GridLineStyle','none','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
    grid on;
    title('Zoom near original','FontName','Times New Roman');

figure(2)
plot(t11,x_y(1:5,:),'linewidth',1.5)
hold on;
plot(t22,x_y(6:6,1:length(t22)),'linewidth',1.5)
xline(10,'--b','LineWidth',1.5,'Label','UAVs 5&6 leave','FontName','Times New Roman');
xline(20,'--r','LineWidth',1.5,'Label','UAV 5 rejoins','FontName','Times New Roman');
  xlabel('$k$','Interpreter','latex')
 ylabel('Position $\xi^k_{i_y}$','Interpreter','latex')
 h=legend('$\xi^k_{1_y}$','$\xi^k_{2_y}$','$\xi^k_{3_y}$','$\xi^k_{4_y}$','$\xi^k_{5_y}$','$\xi^k_{6_y}$');
set(h,'interpreter','latex')
grid on; set(gca,'GridLineStyle',':','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
% 插入局部放大 (离开)
    ax1 = axes('Position',[0.44 0.57 0.2 0.2]);
    box on; hold on;
    plot(t_leave,x_y(:,length(t1):length(t1leave)-1),'LineWidth',1.5);
    xlim([10+0.1,11-0.1]); 
    ylim([40,60]);
    set(gca,'GridLineStyle','none','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
    grid on;
    title('Zoom near leave','FontName','Times New Roman');

    % 插入局部放大 (加入)
    ax2 = axes('Position',[0.16 0.32 0.2 0.2]);
    box on; hold on;
        plot(t1_big,x_y(:,1:length(t1_big)),'LineWidth',1.5);
    xlim([0,1]); 
    ylim([-1,1]);
    set(gca,'GridLineStyle','none','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
    grid on;
    title('Zoom near original','FontName','Times New Roman');
% 
figure(3)
plot(t11,v_x(1:5,:),'linewidth',1.5)
hold on;
plot(t22,v_x(6:6,1:length(t22)),'linewidth',1.5)
xline(10,'--b','LineWidth',1.5,'Label','UAVs 5&6 leave','FontName','Times New Roman');
xline(20,'--r','LineWidth',1.5,'Label','UAV 5 rejoins','FontName','Times New Roman');
xlabel('$k$','Interpreter','latex')
ylabel('Velocity $\zeta^k_{i_x}$','Interpreter','latex')
h=legend('$\zeta^k_{1_x}$','$\zeta^k_{2_x}$','$\zeta^k_{3_x}$','$\zeta^k_{4_x}$','$\zeta^k_{5_x}$','$\zeta^k_{6_x}$');
set(h,'interpreter','latex')
grid on; set(gca,'GridLineStyle',':','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
% 插入局部放大 (离开)
    ax1 = axes('Position',[0.44 0.57 0.2 0.2]);
    box on; hold on;
    plot(t_leave,v_x(1:4,length(t1):length(t1leave)-1),'LineWidth',1.5);
    xlim([10+0.1,11-0.1]); 
    ylim([-3,3]);
    set(gca,'GridLineStyle','none','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
    grid on;
    title('Zoom near leave','FontName','Times New Roman');

    % 插入局部放大 (加入)
    ax2 = axes('Position',[0.17 0.4 0.2 0.2]);
    box on; hold on;
        plot(t1_big,v_x(:,1:length(t1_big)),'LineWidth',1.5);
    xlim([0,1]); 
    ylim([-0.5,0.5]);
    set(gca,'GridLineStyle','none','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
    grid on;
    title('Zoom near original','FontName','Times New Roman');

figure(4)
plot(t11,v_y(1:5,:),'linewidth',1.5)
hold on;
plot(t22,v_y(6:6,1:length(t22)),'linewidth',1.5)
% plot(t11,ave_v(1,:),'linewidth',1.5,'Color',[0.5, 0.5,0.5],'LineStyle','--')
xline(10,'--b','LineWidth',1.5,'Label','UAVs 5&6 leave','FontName','Times New Roman');
xline(20,'--r','LineWidth',1.5,'Label','UAV 5 rejoins','FontName','Times New Roman');
xlabel('$$k$$','Interpreter','latex')
ylabel('Velocity $\zeta^k_{i_y}$','Interpreter','latex')
h=legend('$\zeta^k_{1_y}$','$\zeta^k_{2_y}$','$\zeta^k_{3_y}$','$\zeta^k_{4_y}$','$\zeta^k_{5_y}$','$\zeta^k_{6_y}$');
set(h,'interpreter','latex')
grid on; set(gca,'GridLineStyle',':','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
% 插入局部放大 (离开)
    ax1 = axes('Position',[0.44 0.57 0.2 0.2]);
    box on; hold on;
    plot(t_leave,v_y(1:4,length(t1):length(t1leave)-1),'LineWidth',1.5);
    xlim([10+0.1,11-0.1]); 
    ylim([-3,3]);
    set(gca,'GridLineStyle','none','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
    grid on;
    title('Zoom near leave','FontName','Times New Roman');

    % 插入局部放大 (加入)
    ax2 = axes('Position',[0.17 0.46 0.2 0.2]);
    box on; hold on;
        plot(t1_big,v_y(:,1:length(t1_big)),'LineWidth',1.5);
    xlim([0,1]); 
    ylim([-0.2,0.2]);
    set(gca,'GridLineStyle','none','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
    grid on;
    title('Zoom near original','FontName','Times New Roman');
% 
figure(5)
plot(t11,z_x(1:5,:),'linewidth',1.5)
hold on;
plot(t22,z_x(6:6,1:length(t22)),'linewidth',1.5)
% plot(t11,ave_v(1,:),'linewidth',1.5,'Color',[0.5, 0.5, 0.5],'LineStyle','--')
xline(10,'--b','LineWidth',1.5,'Label','UAVs 5&6 leave','FontName','Times New Roman');
xline(20,'--r','LineWidth',1.5,'Label','UAV 5 rejoins','FontName','Times New Roman');
xlabel('$$k$$','Interpreter','latex')
ylabel('Filter $\varsigma^k_{i_x}$','Interpreter','latex')
h=legend('$\varsigma^k_{1_x}$','$\varsigma^k_{2_x}$','$\varsigma^k_{3_x}$','$\varsigma^k_{4_x}$','$\varsigma^k_{5_x}$','$\varsigma^k_{6_x}$');
set(h,'interpreter','latex')
grid on; set(gca,'GridLineStyle',':','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
% 插入局部放大 (离开)
    ax1 = axes('Position',[0.44 0.5 0.2 0.2]);
    box on; hold on;
    plot(t_leave,z_x(1:4,length(t1):length(t1leave)-1),'LineWidth',1.5);
    xlim([10+0.1,11-0.1]); 
    ylim([-4,4]);
    set(gca,'GridLineStyle','none','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
    grid on;
    title('Zoom near leave','FontName','Times New Roman');

    % 插入局部放大 (加入)
    ax2 = axes('Position',[0.18 0.4 0.2 0.2]);
    box on; hold on;
        plot(t1_big,z_x(:,1:length(t1_big)),'LineWidth',1.5);
    xlim([0,1]); 
    ylim([-0.5,0.5]);
    set(gca,'GridLineStyle','none','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
    grid on;
    title('Zoom near original','FontName','Times New Roman');
% 
figure(6)
plot(t11,z_y(1:5,:),'linewidth',1.5)
hold on;
plot(t22,z_y(6:6,1:length(t22)),'linewidth',1.5)
% plot(t11,ave_v(1,:),'linewidth',1.5,'Color',[0.5, 0.5, 0.5],'LineStyle','--')
xline(10,'--b','LineWidth',1.5,'Label','UAVs 5&6 leave','FontName','Times New Roman');
xline(20,'--r','LineWidth',1.5,'Label','UAV 5 rejoins','FontName','Times New Roman');
xlabel('$$k$$','Interpreter','latex')
ylabel('Filter $\varsigma^k_{i_y}$','Interpreter','latex')
h=legend('$\varsigma^k_{1_y}$','$\varsigma^k_{2_y}$','$\varsigma^k_{3_y}$','$\varsigma^k_{4_y}$','$\varsigma^k_{5_y}$','$\varsigma^k_{6_y}$');
set(h,'interpreter','latex')
grid on; set(gca,'GridLineStyle',':','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
% 插入局部放大 (离开)
    ax1 = axes('Position',[0.44 0.5 0.2 0.2]);
    box on; hold on;
    plot(t_leave,z_y(1:4,length(t1):length(t1leave)-1),'LineWidth',1.5);
    xlim([10+0.1,11-0.1]); 
    ylim([-4,4]);
    set(gca,'GridLineStyle','none','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
    grid on;
    title('Zoom near leave','FontName','Times New Roman');

    % 插入局部放大 (加入)
    ax2 = axes('Position',[0.18 0.45 0.2 0.2]);
    box on; hold on;
        plot(t1_big,z_y(:,1:length(t1_big)),'LineWidth',1.5);
    xlim([0,1]); 
    ylim([-0.5,0.5]);
    set(gca,'GridLineStyle','none','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
    grid on;
    title('Zoom near original','FontName','Times New Roman');
 
figure(7)
plot3(x_x(1:1,:),v_x(1:1,:),t11,'linewidth',1.5)
hold on;
plot3(x_x(2:2,:),v_x(2:2,:),t11,'linewidth',1.5)
hold on;
plot3(x_x(3:3,:),v_x(3:3,:),t11,'linewidth',1.5)
hold on;
plot3(x_x(4:4,:),v_x(4:4,:),t11,'linewidth',1.5)
hold on;
plot3(x_x(5:5,:),v_x(5:5,:),t11,'linewidth',1.5)
hold on;
plot3(x_x(6:6,1:length(t22)),v_x(6:6,1:length(t22)),t22,'linewidth',1.5)
%plot(x(2:4,1:length(t22)),v(2:4,1:length(t22)),'linewidth',1.5)
% hold on;
% plot(x(5:6,length(t1)+1:length(t22)),v(5:6,length(t1)+1:length(t22)),'linewidth',1.5)
 xlabel('$x(k)$','Interpreter','latex')
 ylabel('$y(k)$','Interpreter','latex')
 zlabel('$k$','Interpreter','latex')
 h=legend('UAV 1','UAV 2','UAV 3','UAV 4','UAV 5','UAV 6','FontName','Times New Roman');
set(h,'interpreter','latex')
grid on; set(gca,'GridLineStyle',':','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')

figure(8)
z=tiledlayout(3,1,'TileSpacing','tight','Padding','compact');
nexttile
plot(t11,x_x(1:5,:),'linewidth',1.5)
hold on;
plot(t22,x_x(6:6,1:length(t22)),'linewidth',1.5)
% plot(t11,ave_x(1,:),'linewidth',1.5,'Color',[0.5, 0.5, 0.5],'LineStyle','--')
xline(10,'--b','LineWidth',1.5,'Label','UAVs 5&6 leave','FontName','Times New Roman');
xline(20,'--r','LineWidth',1.5,'Label','UAV 5 rejoins','FontName','Times New Roman');
  %xlabel('$k$','Interpreter','latex')
 ylabel('Positions','Interpreter','latex','FontName','Times New Roman')
 h=legend('UAV 1','UAV 2','UAV 3','UAV 4','UAV 5','UAV 6');
% h=legend('$\xi^k_{1}$','$\xi^k_{2}$','$\xi^k_{3_x}$','$\xi^k_{4_x}$','$\xi^k_{5_x}$','$\xi^k_{6_x}$');
%set(h,'interpreter','latex')
grid on; set(gca,'GridLineStyle',':','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')

nexttile;
plot(t11,v_x(1:5,:),'linewidth',1.5)
hold on;
plot(t22,v_x(6:6,1:length(t22)),'linewidth',1.5)
xline(10,'--b','LineWidth',1.5,'Label','UAVs 5&6 leave','FontName','Times New Roman');
xline(20,'--r','LineWidth',1.5,'Label','UAV 5 rejoins','FontName','Times New Roman');
%xlabel('$k$','Interpreter','latex')
ylabel('Velocities','Interpreter','latex','FontName','Times New Roman')
%h=legend('$\zeta^k_{1_x}$','$\zeta^k_{2_x}$','$\zeta^k_{3_x}$','$\zeta^k_{4_x}$','$\zeta^k_{5_x}$','$\zeta^k_{6_x}$');
%set(h,'interpreter','latex')
grid on; set(gca,'GridLineStyle',':','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')

nexttile;
plot(t11,z_x(1:5,:),'linewidth',1.5)
hold on;
plot(t22,z_x(6:6,1:length(t22)),'linewidth',1.5)
% plot(t11,ave_v(1,:),'linewidth',1.5,'Color',[0.5, 0.5, 0.5],'LineStyle','--')
xline(10,'--b','LineWidth',1.5,'Label','UAVs 5&6 leave','FontName','Times New Roman');
xline(20,'--r','LineWidth',1.5,'Label','UAV 5 rejoins','FontName','Times New Roman');
%xlabel('$$k$$','Interpreter','latex')
ylabel('Filters','Interpreter','latex','FontName','Times New Roman')
%h=legend('$\varsigma^k_{1_x}$','$\varsigma^k_{2_x}$','$\varsigma^k_{3_x}$','$\varsigma^k_{4_x}$','$\varsigma^k_{5_x}$','$\varsigma^k_{6_x}$');
%set(h,'interpreter','latex')
grid on; set(gca,'GridLineStyle',':','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')

xlabel('$$k$$','Interpreter','latex')

