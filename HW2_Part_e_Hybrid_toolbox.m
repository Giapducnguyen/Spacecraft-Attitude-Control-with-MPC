% % HW 2 - Part e
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
sysd = c2d(sysc,Ts,'zoh');
Ad = sysd.A;
Bd = sysd.B;

%% MPC data
N1 = 2;     % CASE 1
N2 = 20;    % CASE 2
Q = blkdiag(eye(3),0.01*eye(3));
R = 0.01*eye(3);
[Qf,~,~] = idare(Ad,Bd,Q,R);

%% Constraints
u_min = [-0.1;-0.1;-0.1];
u_max = -u_min;

%% Formulation
% Define cost, constraints and horizon
clear cost1 constraints1 horizon1

cost1.Q = Q;
cost1.R = R;
cost1.P = Qf;
% cost.K = Bd'*0;
Smpc = ss(Ad,Bd,eye(6),zeros(6,3),Ts) ;
% cost1.rho = 0.1*Inf; % set to Inf for hard constraints (if present)
horizon1.N = N1 ; % output horizon \sum_{k=0}^{Ny-1}
% Horizon.Nu = 12; % input horizon u(0), …, u(Nu-1)
% horizon.Ncu = 11; % input constraints horizon k=0,…, Ncu
% horizon.Ncu = 11; % output constraints horizon k=0,…, Ncy
constraints1.umax = u_max ;
constraints1.umin = u_min ;
% constraints1.ymax = [1.1 0.2 0.1 1.1] ; % x, u, r
% constraints1.ymin = -constraints1.ymax ;
MPCctrl1 = lincon(Smpc, 'reg', cost1, horizon1, constraints1, 'qpact',0) ;

% Define cost, constraints and horizon
clear cost2 constraints2 horizon2

cost2.Q = Q;
cost2.R = R;
cost2.P = Qf;
% cost.K = Bd'*0;
Smpc = ss(Ad,Bd,eye(6),zeros(6,3),Ts) ;
% cost2.rho = 0.1*Inf; % set to Inf for hard constraints (if present)
horizon2.N = N2 ; % output horizon \sum_{k=0}^{Ny-1}
% Horizon.Nu = 12; % input horizon u(0), …, u(Nu-1)
% horizon.Ncu = 11; % input constraints horizon k=0,…, Ncu
% horizon.Ncu = 11; % output constraints horizon k=0,…, Ncy
constraints2.umax = u_max ;
constraints2.umin = u_min ;
% constraints1.ymax = [1.1 0.2 0.1 1.1] ; % x, u, r
% constraints1.ymin = -constraints1.ymax ;
MPCctrl2 = lincon(Smpc, 'reg', cost2, horizon2, constraints2, 'qpact',0) ;

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
    try
        U1_log(:,k) = eval(MPCctrl1,X1_log(:,k));
    catch
        U1_log(:,k) = zeros(3,1);
    end
    execTime1(k) = toc;
    % State update
    [~,X1] = ode45(@(t,x) scdynamics(t,x,U1_log(:,k)), (k-1)*Ts:h:k*Ts, X1_log(:,k));
    X1 = X1'; X1_log(:,k+1) = X1(:,end);
end

% % Main Simulation Loop - case 2: N2 = 20;
for k = 1:Nsim
    tic;
    % Control calculation
    try
        U2_log(:,k) = eval(MPCctrl2,X2_log(:,k));
    catch
        U2_log(:,k) = zeros(3,1);
    end
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