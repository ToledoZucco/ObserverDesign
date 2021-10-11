%In this code, a closed-loop simulation of the vibrating string is presented.
%The string is attached at one side and actuated at the other side with a
%force actuator.
%The control law is given by energy shaping and damping injection
%Refer to Machelli2017 to more details of the control law

close all
clear all
clc

N = 10;%Number of state variables to described the infinite-dimensional system
L = 1;  %Lengh of the string
rho = 1;    %Mass density
T = 1;      %Young's modulus
[A,B,C,D,Q,h] = AttachedActuatedString(N,L,rho,T);
Cw = h*[ones(1,N/2),zeros(1,N/2)];

%% Controller Design
qc = 200;
dc = 50;

Ac = 0;
Bc = 1;
Cc = qc;
Dc = dc;

%% Closed-Loop matrix
ACL = [A-B*Dc*C,-B*Cc;Bc*C,Ac];

%Time discretization
t0 = 0;
dt = 0.0001;
t = t0:dt:1000*dt;
Nt = length(t);
I = eye(length(ACL));
ACLd = inv(I-dt/2*ACL)*(I+dt/2*ACL);



x0 = [ones(N/2,1);zeros(N/2,1)];
xc0 = 0.8;
xc0 = Cw*x0;

z0 = [x0;xc0];
z = zeros(N+1,Nt);
z(:,1) = z0;

for k = 1:Nt
    
    %Closed loop Dynamic
    z(:,k+1) = ACLd*z(:,k);
end
z = z(:,1:end-1);
x = z(1:N,:);
xc = z(N+1:end,:);


%%

%Power preserving interconnection
y = C*x;
uc = y;
yc = Cc*xc+Dc*uc;
u = -yc;  

%Energies
H = 1/2*x'*Q*x;
Hc = 1/2*xc'*qc*xc;

%End-tip position
wb = Cw*x;


%% Figures
x0screen=100;y0screen=50;width=1000;height=600;font=35;lw=4;ms = 15;

%Output and input
figure
subplot(2,1,1)
hold on
plot(t,u,'LineWidth',lw)
legend({'$u(t)$'},'Interpreter','latex','FontSize',font)
grid on
set(gca,'FontSize',font);

subplot(2,1,2)
hold on
plot(t,y,'LineWidth',lw)
legend({'$y(t)$'},'Interpreter','latex','FontSize',font)
grid on
set(gca,'FontSize',font);


% Energies
figure
hold on
plot(t,H,'LineWidth',lw)
plot(t,Hc,'LineWidth',lw)
plot(t,H+Hc,'LineWidth',lw)
legend({'$H(t)$','$H_c(t)$','$V(t)$'},'Interpreter','latex','FontSize',font)
grid on
set(gca,'FontSize',font);

% Lyapunov function
figure
hold on
plot(t,H+Hc,'LineWidth',lw)
legend({'$V(t)$'},'Interpreter','latex','FontSize',font)
grid on
set(gca,'FontSize',font);


% End-tip position
figure
hold on
plot(t,wb,'LineWidth',lw)
plot(t,zeros(1,Nt),'--','LineWidth',lw)
legend({'$w(b,t)$','$w_{desired}(b,t)$'},'Interpreter','latex','FontSize',font)
grid on
set(gca,'FontSize',font);