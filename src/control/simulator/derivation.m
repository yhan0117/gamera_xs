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
w_IB = [[1 0 0].' R(phi,1).'*[0 1 0].' (R(theta,2)*R(phi,1)).'*[0 0 1].'];
w_IB = (R(theta,2)*R(phi,1)).'*[0 0 diff(psi,t)].' +  R(phi,1).'*[0 diff(theta,t) 0].' + [diff(phi, t) 0 0].';


%% Total Energy
T_quad = m_g*sum(v_go.*v_go)/2 + w_IB.'*I_g*w_IB/2 ;
U_quad = m_g*g*z;

%% Lagrangian 
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

int_ode = subs(int_ode, dq, diff(q,t));
int_ode = subs(int_ode, ddq, diff(q,t,2));
int_ode = subs(int_ode, C_l*sum(u), F);
int_ode = subs(int_ode, d*C_l*(u1-u3), tau2);
int_ode = subs(int_ode, d*C_l*(u2-u4), tau1);
int_ode = subs(int_ode, C_l*sum(u), F);

int_ode = diff(q,t,2) == int_ode;

%% Equations of Motion
sol = solve(eqns, ddq);
ode(1) = subs(sol.ddq1, dq, diff(q,t));
ode(2) = subs(sol.ddq2, dq, diff(q,t));
ode(3) = subs(sol.ddq3, dq, diff(q,t));
ode(4) = subs(sol.ddq4, dq, diff(q,t));
ode(5) = subs(sol.ddq5, dq, diff(q,t));
ode(6) = subs(sol.ddq6, dq, diff(q,t));

ode = simplify(ode, "Steps", 150);
ode = diff(q,t,2) == ode;

%% Linearization
clear;clc

% states (independent vars)
syms z(t) [1 9]
syms u I [1 3]
syms t
z = z.';
dz = diff(z,t);


% euler body rate conversion
W = [[1 0 0].' R(z4,1).'*[0 1 0].' (R(z5,2)*R(z4,1)).'*[0 0 1].'];

% integral and 0th order terms
f1 = diff(z1,t) - z4;
f2 = diff(z2,t) - z5;
f3 = diff(z3,t) - z6;
f4 = diff(z4,t) - z7;
f5 = diff(z5,t) - z8;
f6 = diff(z6,t) - z9;

% rewrite equations of motion
syms w_ib
w_ib = W*[z7 z8 z9].';
w_ib = w_ib(t);

LHS = simplify(diff(w_ib,t), "Steps", 10);
RHS = simplify(([w_ib(2)*w_ib(3)*(I2 - I3) ; w_ib(1)*w_ib(3)*(I3 - I1) ; w_ib(1)*w_ib(2)*(I1 - I2)] + u.')./I.',"Steps", 10);
f7 = LHS(1)-RHS(1);
f8 = LHS(2)-RHS(2);
f9 = LHS(3)-RHS(3);

f = [f1;f2;f3;f4;f5;f6;f7;f8;f9];

E = jacobian(f,dz);
E = subs(E(t),dz(t),zeros(9,1))
E = subs(E,z(t),zeros(9,1))
F = jacobian(f,z);
F = subs(F(t),z(t),zeros(9,1))
G = jacobian(f,u.');
G = subs(G(t),z(t),zeros(9,1))

A = double(-inv(E)*F)
B = double(subs(-inv(E)*G,I, [2611.52 2611.52 5160.31]));
Q_ = diag([0.01 0.01 0.01 15 15 5 0.5 0.5 0.1]);
R_ = diag([0.01,0.01,0.1]);
lqr(A,B,Q_,R_)
%% functions
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

