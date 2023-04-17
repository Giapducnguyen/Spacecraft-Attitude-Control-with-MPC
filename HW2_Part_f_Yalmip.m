% % HW 2 - Part f
clear all
close all
clc
%% Time data
Tsim = 200;
Ts = 2;
Nsim = floor(Tsim/Ts);
h = Ts/10;

%% Model data
x0 = zeros(6,1); u0 = zeros(3,1);
[Ac,Bc] = symLin(x0,u0);
sysc = ss(Ac,Bc,eye(6),zeros(6,3));
sysd = c2d(sysc,Ts);
Ad = sysd.A;
Bd = sysd.B;

%% Constraints              Gu * u <= gu
u_min = [-0.1;-0.1;-0.1];
u_max = -u_min;

%% MPC data
N1 = 2;     % CASE 1
N2 = 20;    % CASE 2
Q = blkdiag(eye(3),0.01*eye(3));
R = 0.01*eye(3);
[Qf,~,~] = idare(Ad,Bd,Q,R);

%% Yalmip formulation
yalmip clear; yalmip('clear');
nx = 6;
nu = 3;
% % State-input stacked vectors
% CASE 1
Xvec_1 = sdpvar(repmat(nx,1,N1+1),repmat(1,1,N1+1));
Uvec_1 = sdpvar(repmat(nu,1,N1),repmat(1,1,N1));
% CASE 2
Xvec_2 = sdpvar(repmat(nx,1,N2+1),repmat(1,1,N2+1));
Uvec_2 = sdpvar(repmat(nu,1,N2),repmat(1,1,N2));

% % 
constraints_1 = [];
objective_1 = 0;
for k = 1:N1
%     objective_1 = objective_1 + norm(Q*Xvec_1{k},2) + norm(R*Uvec_1{k},2);
    objective_1 = objective_1 + Xvec_1{k}'*Q*Xvec_1{k} + Uvec_1{k}'*R*Uvec_1{k};
    constraints_1 = [constraints_1, Xvec_1{k+1} == Ad*Xvec_1{k} + Bd*Uvec_1{k}];
    constraints_1 = [constraints_1, u_min <= Uvec_1{k} <= u_max];
end
objective_1 = objective_1 + Xvec_1{N1+1}'*Qf*Xvec_1{N1+1};

% % 
constraints_2 = [];
objective_2 = 0;
for k = 1:N2
%     objective_2 = objective_2 + norm(Q*Xvec_2{k},2) + norm(R*Uvec_2{k},2);
    objective_2 = objective_2 + Xvec_2{k}'*Q*Xvec_2{k} + Uvec_2{k}'*R*Uvec_2{k};
    constraints_2 = [constraints_2, Xvec_2{k+1} == Ad*Xvec_2{k} + Bd*Uvec_2{k}];
    constraints_2 = [constraints_2, u_min <= Uvec_2{k} <= u_max];
end
objective_2 = objective_2 + Xvec_2{N2+1}'*Qf*Xvec_2{N2+1};

% % Controllers
options = sdpsettings('solver', 'osqp');
% CASE 1
controller_1 = optimizer(constraints_1,objective_1,options,Xvec_1{1},Uvec_1{1});
% CASE 2
controller_2 = optimizer(constraints_2,objective_2,options,Xvec_2{1},Uvec_2{1});

%% Main simulation loop
% % initial condition
x0 = [-0.4; -0.8; 1.2; -0.02; -0.02; 0.02];
% % Data logging initialization
Time_log = 0:Ts:Tsim;       % time logging
X1_log = zeros(6,Nsim+1);   % state logging (case 1)
X1_log(:,1) = x0;
X2_log = zeros(6,Nsim+1);   % state logging (case 2)
X2_log(:,1) = x0;
U1_log = zeros(3,Nsim);     % control logging (case 1)
U2_log = zeros(3,Nsim);     % control logging (case 2)

% % Execution time loggings
execTime1 = zeros(1,Nsim);
execTime2 = zeros(1,Nsim);

% % Main Simulation Loop - case 1: N1 = 2;
for k = 1:Nsim
    tic;
    % Control calculation
    U1_log(:,k) = controller_1(X1_log(:,k));
    execTime1(k) = toc;
    % State update
    [~,X1] = ode45(@(t,x) scdynamics(t,x,U1_log(:,k)), (k-1)*Ts:h:k*Ts, X1_log(:,k));
    X1 = X1'; X1_log(:,k+1) = X1(:,end);
end

% % Main Simulation Loop - case 2: N2 = 20;
for k = 1:Nsim
    tic;
    % Control calculation
    U2_log(:,k) = controller_2(X2_log(:,k));
    execTime2(k) = toc;
    % State update
    [~,X2] = ode45(@(t,x) scdynamics(t,x,U2_log(:,k)), (k-1)*Ts:h:k*Ts, X2_log(:,k));
    X2 = X2'; X2_log(:,k+1) = X2(:,end);
end

%% Plots
figure(1)
fig1 = tiledlayout(3,1);

nexttile
hold on
plot(Time_log,X1_log(1,:),'LineStyle','-','Color','m','LineWidth',1);
plot(Time_log,X2_log(1,:),'LineStyle','--','Color','b','LineWidth',1);
ylabel('$\phi$ (rad)','interpreter','latex')
legend({'$N_1= 2$','$N_2= 20$'},'Interpreter','latex','location','best')
title('Orientation Angles','Interpreter','latex')
% ylim([-0.02 0.6])
box on
grid on
ax = gca; ax.FontSize = 14;

nexttile
hold on
plot(Time_log,X1_log(2,:),'LineStyle','-','Color','m','LineWidth',1);
plot(Time_log,X2_log(2,:),'LineStyle','--','Color','b','LineWidth',1);
ylabel('$\theta$ (rad)','interpreter','latex')
% legend({'$N_1= 2$','$N_2= 20$'},'Interpreter','latex','location','best')
% title('Pitch Angle','Interpreter','latex')
% ylim([-0.02 0.6])
box on
grid on
ax = gca; ax.FontSize = 14;

nexttile
hold on
plot(Time_log,X1_log(3,:),'LineStyle','-','Color','m','LineWidth',1);
plot(Time_log,X2_log(3,:),'LineStyle','--','Color','b','LineWidth',1);
ylabel('$\psi$ (rad)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
% legend({'$N_1= 2$','$N_2= 20$'},'Interpreter','latex','location','best')
% title('Yaw Angle','Interpreter','latex')
% ylim([-0.02 0.6])
box on
grid on
ax = gca; ax.FontSize = 14;

fig1.TileSpacing = 'compact';
fig1.Padding = 'compact';
set(gcf,'Units','points','position',[300, 50, 400, 400])

figure(2)
fig2 = tiledlayout(3,1);

nexttile
hold on
plot(Time_log,X1_log(4,:),'LineStyle','-','Color','m','LineWidth',1);
plot(Time_log,X2_log(4,:),'LineStyle','--','Color','b','LineWidth',1);
ylabel('$\omega_1$ (rad/s)','interpreter','latex')
legend({'$N_1= 2$','$N_2= 20$'},'Interpreter','latex','location','best')
title('Angular Velocity Components','Interpreter','latex')
% ylim([-0.02 0.6])
box on
grid on
ax = gca; ax.FontSize = 14;

nexttile
hold on
plot(Time_log,X1_log(5,:),'LineStyle','-','Color','m','LineWidth',1);
plot(Time_log,X2_log(5,:),'LineStyle','--','Color','b','LineWidth',1);
ylabel('$\omega_2$ (rad/s)','interpreter','latex')
% legend({'$N_1= 2$','$N_2= 20$'},'Interpreter','latex','location','best')
% title('Angular Velocity Components','Interpreter','latex')
% ylim([-0.02 0.6])
box on
grid on
ax = gca; ax.FontSize = 14;

nexttile
hold on
plot(Time_log,X1_log(6,:),'LineStyle','-','Color','m','LineWidth',1);
plot(Time_log,X2_log(6,:),'LineStyle','--','Color','b','LineWidth',1);
ylabel('$\omega_3$ (rad/s)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
% legend({'$N_1= 2$','$N_2= 20$'},'Interpreter','latex','location','best')
% title('Angular Velocity Components','Interpreter','latex')
% ylim([-0.02 0.6])
box on
grid on
ax = gca; ax.FontSize = 14;

fig2.TileSpacing = 'compact';
fig2.Padding = 'compact';
set(gcf,'Units','points','position',[300, 50, 400, 400])

figure(3)
fig3 = tiledlayout(3,1);

nexttile
hold on
stairs(Time_log(1,1:end-1),U1_log(1,:),'LineStyle','-','Color','m','LineWidth',1);
stairs(Time_log(1,1:end-1),U2_log(1,:),'LineStyle','--','Color','b','LineWidth',1);
stairs(Time_log(1,1:end-1),(-0.1)*ones(size(Time_log(1,1:end-1))),'LineStyle',':','Color','k','LineWidth',1);
stairs(Time_log(1,1:end-1),(0.1)*ones(size(Time_log(1,1:end-1))),'LineStyle',':','Color','k','LineWidth',1);
ylabel('$M_1$','interpreter','latex')
legend({'$N_1= 2$','$N_2= 20$','Constraints'},'Interpreter','latex','location','best')
title('Control Moments','Interpreter','latex')
ylim([-0.12 0.12])
box on
grid on
ax = gca; ax.FontSize = 14;

nexttile
hold on
stairs(Time_log(1,1:end-1),U1_log(2,:),'LineStyle','-','Color','m','LineWidth',1);
stairs(Time_log(1,1:end-1),U2_log(2,:),'LineStyle','--','Color','b','LineWidth',1);
stairs(Time_log(1,1:end-1),(-0.1)*ones(size(Time_log(1,1:end-1))),'LineStyle',':','Color','k','LineWidth',1);
stairs(Time_log(1,1:end-1),(0.1)*ones(size(Time_log(1,1:end-1))),'LineStyle',':','Color','k','LineWidth',1);
ylabel('$M_2$','interpreter','latex')
% legend({'$N_1= 2$','$N_2= 20$'},'Interpreter','latex','location','best')
% title('Control Moments','Interpreter','latex')
ylim([-0.12 0.12])
box on
grid on
ax = gca; ax.FontSize = 14;

nexttile
hold on
stairs(Time_log(1,1:end-1),U1_log(3,:),'LineStyle','-','Color','m','LineWidth',1);
stairs(Time_log(1,1:end-1),U2_log(3,:),'LineStyle','--','Color','b','LineWidth',1);
stairs(Time_log(1,1:end-1),(-0.1)*ones(size(Time_log(1,1:end-1))),'LineStyle',':','Color','k','LineWidth',1);
stairs(Time_log(1,1:end-1),(0.1)*ones(size(Time_log(1,1:end-1))),'LineStyle',':','Color','k','LineWidth',1);
ylabel('$M_3$','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
% legend({'$N_1= 2$','$N_2= 20$'},'Interpreter','latex','location','best')
% title('Control Moments','Interpreter','latex')
ylim([-0.12 0.12])
box on
grid on
ax = gca; ax.FontSize = 14;

fig3.TileSpacing = 'compact';
fig3.Padding = 'compact';
set(gcf,'Units','points','position',[300, 50, 400, 400])

figure(4)
hold on
plot(Time_log(1,1:end-1),execTime1,'LineStyle','-','Color','m','LineWidth',1);
plot(Time_log(1,1:end-1),execTime2,'LineStyle','--','Color','b','LineWidth',1);
xlabel('Time (sec)','interpreter','latex')
ylabel('t (sec)','interpreter','latex')
legend({'$N_1= 2$','$N_2= 20$'},'Interpreter','latex','location','best')
title('Computation Time','Interpreter','latex')
% ylim([-0.02 0.6])
box on
grid on
ax = gca; ax.FontSize = 14;

%% Function helpers

% % Spacecraft Attitude Dynamics
function xdot = scdynamics(t,x,u)
J1 = 120;
J2 = 100;
J3 = 80;
% % x = [\phi, \theta, \psi, omega_1, omega_2, omega_3]';
phi = x(1); theta = x(2); psi = x(3);
omega1 = x(4); omega2 = x(5); omega3 = x(6);
% % u = [M1, M2, M3]';
M1 = u(1); M2 = u(2); M3 = u(3);

xdot = [omega1 + (1/cos(theta))*sin(phi)*sin(theta)*omega2 + (1/cos(theta))*cos(phi)*sin(theta)*omega3;
        cos(phi)*omega2 + (-sin(phi))*omega3;
        (1/cos(theta))*sin(phi)*omega2 + (1/cos(theta))*cos(phi)*omega3;
        ((J2-J3)/J1)*omega2*omega3 + (1/J1)*M1;
        ((J3-J1)/J2)*omega3*omega1 + (1/J2)*M2;
        ((J1-J2)/J3)*omega1*omega2 + (1/J3)*M3];
end

% % Linearized model
function [Ac,Bc] = symLin(x_0,u_0)
x = sym('x',[6;1]); 
u = sym('u',[3;1]);
t = sym('t',1);
Ac = double(subs(jacobian(scdynamics(t,x,u),x),[x;u],[x_0;u_0]));
Bc = double(subs(jacobian(scdynamics(t,x,u),u),[x;u],[x_0;u_0]));
end