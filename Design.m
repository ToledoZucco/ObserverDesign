close all
clear all
clc

%false: only show the simulation
%true: save the figures in png
saveFigures = false;

%String parameters
%spatial variable zeta in [a,b]
a = 0;
b = 1;
rho = 1;    %mass density
T = 1;      %modulous of elasticity

%Model
[A,B,C,D,Q,h] = AttachedActuatedString(N,long,rho_c,T_c);
np = N/2;
nq = N/2;

%State feeback design
QKlqr = 0.1*eye(nd);
RKlqr = eye(2);
NKlqr = 0;
[K,SK,PK] = lqr(Ad,Bd,QKlqr,RKlqr,NKlqr);

AK = Ad-Bd*K;
EK = eig(AK);


n = length(Ad);
CK = -(K'*Cd+Cd'*K);
alpha = 0.1;
Rc = alpha*eye(n);
% Rc(np+1:n,np+1:n) = diag(1*ones(np,1));

Rc(1:nqd,1:nqd) = diag(0.01*ones(nqd,1));
Rc(nqd+1:end,nqd+1:end) = diag(0.01*ones(npd,1));

% 
% Rc(1,1) = 10;
% Rc(2,2) = 10;
% Rc(3,3) = 10;
% % % Rc(10,10) = 100;
% Rc(npd-2,npd-2) = 10;
% Rc(nqd-1,nqd-1) = 10;
% Rc(nqd,nqd) = 10;
% 
gainRc = 1500;

Rc(nqd+1,nqd+1) = gainRc;
% Rc(nqd+2,nqd+2) = 10;
% Rc(nqd+3,nqd+3) = 25;
% Rc(nqd+4,nqd+4) = 12;
% 
% Rc(end-3,end-3) = 12;
% Rc(end-2,end-2) = 25;
% Rc(end-1,end-1) = 10;
Rc(end,end) = gainRc;


% Rc(1:np,1:np) = diag(linspace(0.0001,0.1,np));
% Rc(np+1:n,np+1:n) = diag(linspace(1,10,np)).^4;
% Rc(np+1:n,np+1:n) = diag(410*ones(np,1));
% Rc = diag([linspace(0.5,0.01,np),linspace(10,1000,np)]);
Hm = [AK,2*Rc;
      -CK,-AK'];
Em = eig(Hm);

figure
plot(real(Em),imag(Em),'x')
grid on




%%% Solving the ARE
Kc = K;
Gd = -CK;
Add = AK';
B1 = chol(1.1*(abs(min(eig(Gd))))*eye(n))';
B2 = chol(Gd+B1*B1')';
min(min(Gd - B2*B2' + B1*B1')) %this must be zero

PHI = SolveARE(Add,B1,B2,Rc);
Qci = PHI;
Qc = inv(Qci);
min(min(AK'*Qc + Qc*AK + 2*Qc*Rc*Qc + CK))
min(eig(Qc))



Jc = 1/2*(AK*Qci-Qci*AK'-Qci*(Kc'*Cd-Cd'*Kc)*Qci);
Bc = Qci*Kc';
Lc = Bc;


ME = (Jc-Rc)*Qc-(Ad-Bd*Kc-Lc*Cd);
mME = min(min(ME))

AL = Ad-Lc*Cd;
EL = eig(AL);
Ed = eig(Ad);



font=24; lw=2; ms = 10;
x0screen=100;y0screen=100;WidthScreen=1000;HeightScreen=950;


figure
hold on
% plot(real(Ed),imag(Ed),'xk','color',[0,0,0])
% plot(real(EK),imag(EK),'xk','color',[0,0,0]+0.40)
% plot(real(EL),imag(EL),'xk','color',[0,0,0]+0.60)
plot(real(Ed),imag(Ed),'xk','LineWidth',lw)
plot(real(EK),imag(EK),'xb','LineWidth',lw)
plot(real(EL),imag(EL),'xg','LineWidth',lw)
title({'Eigenvalues'},'Interpreter','latex','FontSize',font)
legend({'$\lambda{(A)}$','$\lambda{(A_{K})}$','$\lambda{(A_{L})}$'},'Location','northwest','Interpreter','latex','FontSize',font)
xlabel({'Real axis'},'Interpreter','latex','FontSize',font)
% xlim([-1.2,0.2])
ylabel({'Imag axis'},'Interpreter','latex','FontSize',font)
grid on



L = Lc;

save('OBSFController','K','L','Ad','Bd','Cd','npd','nqd','zpd','zqd','hd')

AL = Ad-L*Cd;
EL = eig(AL);





%Save a file with the PNG figures
figureDirectory = 'Eig';
filename = 'Eig';
N = 30;

SaveEigsFigures(Ad,Bd,Cd,K,L,npd,a,b,rho,T,N,figureDirectory,filename,saveFigures) 
