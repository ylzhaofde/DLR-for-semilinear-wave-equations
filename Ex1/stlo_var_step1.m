function [U1,S1,V1,T1,R1,W1] = stlo_var_step1(U0,S0,V0,T0,R0,W0,fu,tau,omega,tspan)
% fu --> the nonlinear term, functional handle

%% B-step: tau/2, delta_B = tau/2*f(U*S*V')
% delta = tau/2*fu(U0*S0*V0', T0*R0*W0');
% [T,R,W,L] = unconv_step(T0,R0,W0,delta);
[T,R,W,L] = Lie_split_1st(@(w) fu(U0*S0*V0',w),T0,R0,W0,tspan);

%% A-step:tau, delta_A = omega*tau*T*L'
[U1,S1,V1,La] = unconv_step(U0,S0,V0,T,L',omega*tau);

%% B-step: tau/2, delta_B = tau/2*f(U*S*V')
% delta = tau/2*fu(U1*La',T*R*W');
% [T1,R1,W1,~] = unconv_step(T,R,W,delta);
[T1,R1,W1,~] = Lie_split_1st(@(w) fu(U1*La',w),T,R,W,tspan);

end

function [U1,S1,V1,L1] = unconv_step(U,S,V,T,LL,delta)
% in delta, tau already included
% adjusted for application to second-order problems
Kt = delta*T*(LL*V);
K = U * S + Kt;
[U1,~]= qr(K,0);
M = U1' * U;
L = V * S' + delta*LL'*(T'*U);
[V1,~] = qr(L,0);
N = V1'* V;
S0 = M * S * N';
Kt = delta*T*(LL*V1);
S1 = S0 + U1' * Kt;
L1 = V1 * S1';
end

