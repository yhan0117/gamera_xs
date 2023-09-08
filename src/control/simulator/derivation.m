% Dynamics of Quadrotor with slung load
clc;clear;close all

% Generalized Coordinates
syms x(t) y(t) z(t) phi(t) theta(t) psi(t)
syms dq [1 8]
syms ddq [1 8]

q = [x y z phi theta psi];
q = q(t);


%% System Parameters
% Physical Parameters
% I1 I2 I3 mg mp d l g Cd Cl
% 
% Coordinates
% x y z phi theta psi alpha beta
% 
% Inputs
% u1~4 = w^2

syms I_1 I_2 I_3
syms m d g C_d C_l
param = [I_1 I_2 I_3 m d g C_d C_l];
I_g = [I_1   0   0;
       0    I_2  0;
       0    0   I_3];


%% Kinematics
r_go = [x,y,z].';
v_go = diff(r_go, t);  % wrt to I

% 3 Rotations all wrt to body frame
% Roll->Pitch->Yaw
R_IB = R(psi,3)*R(theta,2)*R(phi,1);

w_IB = (R(theta,2)*R(phi,1)).'*[0 0 diff(psi,t)].' +  R(phi,1).'*[0 diff(theta,t) 0].' + [diff(phi, t) 0 0].';


%% Total Energy
T_quad = m_g*sum(v_go.*v_go)/2 + w_IB.'*I_g*w_IB/2 ;
U_quad = m_g*g*z;

%% Lagrangian Formulation
% Lagrangian
L = T_quad - U_quad;

% LHS of Euler Lagrangian
LHS = t*zeros(1,8);
for i = 1:8
    LHS(i) = simplify(diff(diff(L,diff(q(i),t)), t) - diff(L,q(i)), "Steps", 30);
end

% Generalized Forces
syms u [1 4]
syms T
syms F tau1 tau2 tau3

B = [C_l  C_l  C_l  C_l;
     0   C_l  0   -C_l;
     -C_l 0   C_l  0;
     C_d  -C_d C_d  -C_d];

r = r_go + d*R_IB*[[1 0 0]', [0 1 0]', -[1 0 0]', -[0 1 0]'];
r = r(t);

RHS = t*zeros(1,6);
for i = 1:6
    RHS(i) = simplify(Q(q(i), r, u, C_l, R_IB, T, r_go), "Steps", 5);
end

% Euler-Lagrange Eqns
eqns = t*zeros(1,8);
for i = 1:6
    eqns(i) = LHS(i) == RHS(i);
end
eqns = subs(eqns, diff(q,t,2), ddq);
eqns = subs(eqns, diff(q,t), dq);

% intermediate EoM prior to separating 2nd degree terms
int_ode = t*zeros(1,8);
int_ode(1) = simplify(solve(eqns(1), ddq1), "Steps",100);
int_ode(2) = simplify(solve(eqns(2), ddq2), "Steps",100);
int_ode(3) = simplify(solve(eqns(3), ddq3), "Steps",100);
int_ode(4) = simplify(solve(eqns(4), ddq4), "Steps",100);
int_ode(5) = simplify(solve(eqns(5), ddq5), "Steps",100);
int_ode(6) = simplify(solve(eqns(6), ddq6), "Steps",100);
int_ode(7) = simplify(solve(eqns(7), ddq7), "Steps",100);
int_ode(8) = simplify(solve(eqns(8), ddq8), "Steps",100);

int_ode = subs(int_ode, dq, diff(q,t));
int_ode = subs(int_ode, ddq, diff(q,t,2));
int_ode = subs(int_ode, C_l*sum(u), F);
int_ode = subs(int_ode, d*C_l*(u1-u3), tau2);
int_ode = subs(int_ode, d*C_l*(u2-u4), tau1);
int_ode = subs(int_ode, C_l*sum(u), F);

int_ode = diff(q,t,2) == int_ode;
int_ode.'
%%
% Equations of Motion
sol = solve(eqns, ddq);
ode(1) = subs(sol.ddq1, dq, diff(q,t));
ode(2) = subs(sol.ddq2, dq, diff(q,t));
ode(3) = subs(sol.ddq3, dq, diff(q,t));
ode(4) = subs(sol.ddq4, dq, diff(q,t));
ode(5) = subs(sol.ddq5, dq, diff(q,t));
ode(6) = subs(sol.ddq6, dq, diff(q,t));
ode(7) = subs(sol.ddq7, dq, diff(q,t));
ode(8) = subs(sol.ddq8, dq, diff(q,t));

ode = simplify(ode, "Steps", 150);
ode = diff(q,t,2) == ode;

latex()
% % clearvars -except int_ode ode r_po R_IB q u param t % clean up workspace
% % save("Mats\EOM.mat")


%%
clear;
sys_parameters = [0.0140 0.0140 0.0300 1.6450 0.3 0.2500 1 9.8000 0.0005 0.0012];
save("Mats\sys_param.mat")


%--------------------------------------Functions--------------------------------------%
% Generalized Forces
function Gen_F = Q(var, r, u, Cl, R_IB, T, r_go)
    Gen_F = 0;
    for i = 1:4
        Gen_F = Gen_F + sum( (u(i)*Cl*R_IB*[0 0 1]') .* (diff(r(:,i), var)) );
    end
    Gen_F = Gen_F + sum( (T*R_IC*[0 0 -1].') .* (diff(r_go, var))) + sum( (T*R_IC*[0 0 1].') .* (diff(r_po, var)));
end


% Rotation Matrices
function Rot_Mat = R(x, y)
    if y == 3
        Rot_Mat = [cos(x)   -sin(x) 0;
                   sin(x)   cos(x)  0;
                   0        0       1];
    elseif y == 2      
        Rot_Mat = [cos(x)   0       sin(x);
                   0        1       0;
                   -sin(x)  0       cos(x)];
    elseif y == 1
        Rot_Mat = [1        0       0;
                   0        cos(x)  -sin(x)
                   0        sin(x)  cos(x)];
    end
end

