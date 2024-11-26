function [U1,S1,V1,L1] = Lie_split_1st(fu,U,S,V,tspan) % for Runge-Kutta 4
% Lie-Trotter splitting
% First order projector splitting integrator (KSL)
% for dynamical low-rank approximation 
% of the solution of a matrix differential equation 
% Each subproblem is solved by explicit Runge-Kutta with order 4 or 2
%
% INPUT:
%      U,S,V: initial values for the factors 
%      tspan: time interval
% OUTPUT:
%      U,S,V: factors of the solution after one step

%% Runge-Kutta order 4/2
ns_int = 5; % numer of intermediate time steps

% K step: Update U
K = U*S;
K_fun = @(t,X) fu(X*V')*V;
K = rk_2(K_fun,K,tspan,ns_int);
[U1,~] = qr(K,0);
M = U1'*U;
% L step: Update S and V
L = V*S';
L_fun = @(t,X) fu(U*X')'*U;
L = rk_2(L_fun,L,tspan,ns_int);
[V1,~] = qr(L,0);
N = V1'*V;
% S step
S_fun = @(t,X) U1'*fu(U1*X*V1')*V1;
S1 = rk_2(S_fun,M*S*N',tspan,ns_int);
L1 = V1*S1';
end

function [Y,YY] = rk_2(f,Y,tspan,ns) 
% explicit Runge-Kutta method of order 2: Heun method
% INPUT
% f: function (rhs of the differential equation)
% tspan: interval of integration
% ns: number of time steps
% Y: initial condition
%
% OUTPUT
% Y: approximate solution at T = tspan(2)
% YY: tensor structure which collects the solution at each time step

k = (tspan(2) - tspan(1))/ns; % time step

for i = 1:ns
%     Y_1 = f(tspan(1) + (i - 1)*k,Y);
%     Y_2 = f(tspan(1) + (i - 1)*k + 1/3*k,Y + 1/3*k*Y_1);
%     Y_3 = f(tspan(1) + (i - 1)*k + 2/3*k,Y - 1/3*k*Y_1 + k*Y_2);
%     Y_4 = f(tspan(1) + (i - 1)*k + k,Y + k*Y_1 - k*Y_2+ k*Y_3);
%     Y = Y + 1/8*k*(Y_1 + 3*Y_2 + 3*Y_3 + Y_4);
    Y_1 = f(tspan(1) + (i - 1)*k,Y);
    Y_2 = f(tspan(1) + (i - 1)*k + 2/3*k,Y + 2/3*k*Y_1);
%     Y_3 = f(tspan(1) + (i - 1)*k + 1/2*k,Y + 1/2*k*Y_2);
%     Y_4 = f(tspan(1) + (i - 1)*k + k,Y + k*Y_3);
    Y = Y + 1/4*k*(Y_1 + 3*Y_2);
    YY(:,:,i) = Y;
end
end

