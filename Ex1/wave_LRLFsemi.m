% Semilinear wave equation with two nonlinear terms, 
% discetized with second-order finite differences and
% zero boundary conditions.
% 
clear;close all;
clc;

xmin = 0; xmax = 1;
ymin = xmin; ymax = xmax;
alf = 1; beta = 1e-1; gama = 1e-3; delta_numb = 1; 
fu = @(u,w) u.^2 + sin(w);
T = 0.1;  % end of time
M = 10; % time grids
N = 512; % space grids
h = (xmax - xmin)/N; % hx = hy
tau = T/M;
x = xmin + (0:N)'*h; y = x;
t = (0:M)'*tau;
[u0,v0] = initials(x(2:N),y(2:N));

tic;
ra = 13; % approx rank of A(t)
rb = ra; % approx rank of B(t)
omega = [98;1;1]/100; % [1;1;1]/3; % The weights of splitting 
% v --> v = u'(x,y,t)     
[Uu,Su,Vu] = svd(u0);                 
[Uv,Sv,Vv] = svd(v0); 
Uu = Uu(:,1:ra); Su = Su(1:ra,1:ra); Vu = Vu(:,1:ra);
Uv = Uv(:,1:rb); Sv = Sv(1:rb,1:rb); Vv = Vv(:,1:rb);
% u0 = Uu*Su*Vu'; 
clearvars u0 v0;
% the matrices in space
% x-direction
e1 = -ones(N - 1,1);
Lx = 1/h^2*spdiags([e1 -2*e1 e1],[-1,0,1],N - 1,N - 1); % approx (-\Delta)
% y-direction
Ly = Lx; % L2
clearvars e1;
% initial set
Is = speye(N - 1);
% half step for linear left subproblem
expleft = expm(tau/2*[sparse(N - 1,N - 1),omega(1)*Is;...
                      -alf*Lx - 1/2*delta_numb*Is, -beta*Lx - 1/2*gama*Is]); 
% half step for linear right subproblem
expright = expm(tau/2*[sparse(N - 1,N - 1),-alf*Ly - 1/2*delta_numb*Is;...
                       omega(2)*Is, -beta*Ly - 1/2*gama*Is]);
clearvars Is Lx Ly;

for j = 1:M
% Linear left: tau/2
% delta = exp(tau/2*[0,omega1^2*I; Lx, 0])*[A0;B0]
    delta = expleft*[Uu*Su*Vu';Uv*Sv*Vv'];  
    [Uu,Su,Vu,Uv,Sv,Vv] = lo_linear_left(Uu,Su,Vu,Uv,Sv,Vv,delta,N);
% Linear right: tau/2  
% delta = [A0,B0]*exp(tau*[0,Ly; omega2^2*I, 0])
    delta = [Uu*Su*Vu',Uv*Sv*Vv']*expright; 
    [Uu,Su,Vu,Uv,Sv,Vv] = lo_linear_right(Uu,Su,Vu,Uv,Sv,Vv,delta,N);
% Nonlinear: tau
    [Uu,Su,Vu,Uv,Sv,Vv] = stlo_var_step1(Uu,Su,Vu,Uv,Sv,Vv,fu,tau,omega(3),[t(j);t(j) + tau/2]);
% Linear right: tau/2  
% delta = [A0,B0]*exp(tau*[0,Ly; omega2^2*I, 0])
    delta = [Uu*Su*Vu',Uv*Sv*Vv']*expright; 
    [Uu,Su,Vu,Uv,Sv,Vv] = lo_linear_right(Uu,Su,Vu,Uv,Sv,Vv,delta,N);
% Linear left: tau/2
% delta = exp(tau/2*[0,omega1^2*I; Lx, 0])*[A0;B0]
    delta = expleft*[Uu*Su*Vu';Uv*Sv*Vv'];  
    [Uu,Su,Vu,Uv,Sv,Vv] = lo_linear_left(Uu,Su,Vu,Uv,Sv,Vv,delta,N);
% low-rank approx
%     u(2:N,2:N,j + 1) = Uu*Su*Vu';    
end
toc
u = zeros(N + 1,N + 1); % low-rank approx solution 
ut = u;
u(2:N,2:N) = Uu*Su*Vu';
ut(2:N,2:N) = Uv*Sv*Vv';
clear Is Uu Su Vu Uv Sv Vv delta Lx Ly expleft expright;
load('u_0.1_80000_512.mat'); 
ratio = 512/N;
err = norm(u - u_ref(1:ratio:end,1:ratio:end),'fro')/norm(u_ref(1:ratio:end,1:ratio:end),'fro');
fprintf('The err: %.4e\n',err);
% [X,Y] = meshgrid(x,y);
% subplot(1,2,1);
% mesh(X,Y,u);
% subplot(1,2,2);
% mesh(X,Y,ut);
% for j = 1:M + 1
%     [du(2:N + 1,2:N + 1,j),~,~,~] = ...
%           initials_homwave(x,y,t(j),'travel2');
% end

