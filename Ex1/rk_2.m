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