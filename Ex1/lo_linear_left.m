function [U1,S1,V1,T1,R1,W1] = lo_linear_left(U0,S0,V0,T0,R0,W0,delta,N)
% in delta, tau already included
% project-splitting integrator with explicit Euler method
% adjusted for application to second-order problems
% A = U*S*V^H; B = T*R*W^H
deltaA = delta(1:N - 1,:);
deltaB = delta(N:end,:);

% A-step
US1 = deltaA*V0; % U*S
% K1 = U0 * S0 + Kt;
[U1,~] = qr(US1,0);
% tilde_S0 = hat_S1 - U1'* Kt;
% L1 = V0 * S0' + delta'*U1;
VSH1 = deltaA'*U1; % V*S^H
[V1,S1] = qr(VSH1,0);
S1 = S1';

% B-step
TR1 = deltaB*W0; % T*R
% K1 = U0 * S0 + Kt;
[T1,~] = qr(TR1,0);
% tilde_S0 = hat_S1 - U1'* Kt;
% L1 = V0 * S0' + delta'*U1;
WRH1 = deltaB'*T1; % W*R^H
[W1,R1] = qr(WRH1,0);
R1 = R1';

