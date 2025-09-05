clear
clc
close all
%% Laplacian Matrix
L = [
 1/2   0    0   -1/2	0    0;    
-1/4  3/4 -1/4	-1/4	0    0;    
-1/4   0   3/4	-1/4	0  -1/4;    
  0  -1/3 -1/3	 2/3	0    0;
  0  -1/2   0     0    1/2   0;
  0    0    0     0   -1/2  1/2];
% gains
n = 6; 
alpha = 0.4; 
a = 3; 
b = 8;

%% Initial states
x0 = 0.1*[5 3 -4 -3 -1 -2]'; 
v0 = [0 0 0 0 0 0]';
w0 = 0.1*[1 -1 2 -1 1 2]';
theta0 = zeros(n,1);
z0 = [x0; v0; w0; theta0];

%% Dynamics
f = @(t,z) dynamics_double(t,z,n,alpha,a,b,L);
Tend = 20;
[t,z] = ode45(f,[0 Tend],z0);
x = z(:,1:n); 
v = z(:,n+1:2*n);
w = z(:,2*n+1:3*n);

%% Draw graphs
figure;
tiledlayout(3,1);
vars = {x, v, w};
ylabs = {'Positions','Velocities','Filters','FontName','Times New Roman'};
for k = 1:3
    nexttile;
    plot(t, vars{k}, 'LineWidth', 1.5);
    grid on;
    ylabel(ylabs{k});
    set(gca,'GridLineStyle',':','GridColor',[0 0 0],'GridAlpha',1,'FontName','Times New Roman')
    if k == 1
        legend(compose('UAV %d',1:n), 'Location','bestoutside','FontName','Times New Roman');
    end
    if k == 3
        xlabel('Time (s)','FontName','Times New Roman');
    end
end

%% Dynamics function
function dz = dynamics_double(~,z,n,alpha,a,b,L)
    x = z(1:n); 
    v = z(n+1:2*n); 
    w = z(2*n+1:3*n);
    theta = z(3*n+1:4*n);

    % Filter
    dtheta = -a*w;
    w_new = theta + b*L*x;

    % Control
    u = -alpha*(L*x + w);

    % Dynamics
    dx = v;
    dv = u;
    dw = -a*w + b*L*v;
    dz = [dx; dv; dw; dtheta];
end
