% Quadcopter simulation
clc; clear; clear global; close all;

record = 0; 
trial = 1;


%% Parameters
d = 1;      % arm length
p.m = 4;    % quadcopter mass
p.g = 1;    % gravity 
p.I = [1 1 1];    % moment of inertia, 3 vector

% actuator dynamics matrix
Cl = 1;     % lift coefficients of propellers
Cd = 1;     % drag coefficients of propellers
p.K = [ Cl Cl Cl Cl;
        0 d*Cl 0 -d*Cl;
        -d*Cl 0 d*Cl 0;
        Cd -Cd Cd -Cd];

% rotation matrices
syms a b c
R1 = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
R2 = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
R3 = [cos(c) -sin(c) 0; sin(c) cos(c) 0; 0 0 1];

p.R = matlabFunction(simplify(R3*R2*R1));   % rotation matrix B_R_I
L = [eye(3)*[1;0;0] R1^-1*[0;1;0] R1^-1*R2^-1*[0;0;1]];
p.L = matlabFunction(simplify(L));          % euler rate to body rate mapping matrix B_L_I

% simulation parameters
t_span = 10;                    % duration of simulation
fps = 30;                       % frames per second
t_s = linspace(0,t_span,t_span*fps);    % time stamps
z0 = zeros(12,1);               % initial state


%% Controls 
% desired states
k.zd = [1 0 0 0];

% control gains
k.kz = 6;
k.kdz = 0.3;
k.kp = 1;
k.kdp = 1;
k.kq = 1;
k.kdq = 1;
k.kr = 1;
k.kdr = 1;

% dummy input for testing
% u = @(t,z)[1.01,1,0.99,1]';

% feedback controller u = f(t,z) by calling the control function
u = @(t,z)control(t,z,k,p);


%% Simultation
% solves for the trajectory numerically through an RK integrater
disp("Producing simulation")
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~, z] = ode45(@(t,z)dynamic(t,z,u,p), t_s, z0, options);  % <-- produce simulation


%% Save data
if record
    if ~exist("data", 'dir')
       mkdir("data")    % make sure directory exist
    end
    save("data\trial" + trial + ".mat", 'z')
end


%% Plot 
% TODO:
%   Plot the states and input trajectory (also reference signal u(t))
% figure
% zdirect = [z(:,9), z(:,3), k.zd(1)*ones(size(z(:,3)))];
% plot(t_s,zdirect);
% xlabel('Time', 'FontSize', 14);
% ylabel('Height', 'FontSize', 14);
% % title('Nice-Looking MATLAB Plot', 'FontSize', 16);
% grid on;
% %   Try to make this its own file


%% Animation
% TODO:
%   Han: change the scaling and timing of the animation 
filename = "animations\trial" + trial + ".gif";
if record 
    if ~exist("animations", 'dir')
        mkdir("animations") % make sure directory exist
    elseif exist(filename, "file") 
        delete(filename)    % clear existing file
    end
end

% create and configurate a figure
figure(2); clf;
set(gcf,"position", [0,0,900,600])  % set window size
movegui(gcf, 'center');             % center animation
view(90,0);                     % initial view angle, adjustable

animate(t_s,z,p,filename,record,fps)    % <-- produce animation

disp("Done!!")
%close all


