% This function computes the state derivates of the quadcopter through the state equations
function dz = dynamic(t,z,u,p)

    % INPUTS:
    %   t = simulation time 
    %   z = [x;y,;z;psi;theta;phi;dx;dy;dz;p;q;r] 
    %     = state of the system, n = 12
    %   u = system input (pwm to each motor), m = 4
    %   p = parameter struct
    %       m = quadcopter mass
    %       g = gravity 
    %       I = moment of inertia, 3 vector
    %       K = actuator dynamics matrix
    %       L = body rate to euler rate matrix B_L_I
    %       R = body frame to inertial frame rotation matrix B_R_I
    %
    % OUTPUTS:
    %   dz = dz/dt = time derivative of states
    
    % control law
    u = u(t,z);
    
    % TODO: 
    %   add saturation, aka max torque the motors can realistically achieve

    % unpack physical parameters
    m = p.m;    % quadcopter mass
    g = p.g;    % gravity 
    I = p.I;    % moment of inertia, 3 vector
    K = p.K;    % actuator dynamics matrix
    L = p.L;    % matrix B_L_I
    R = p.R;    % rotation matrix B_R_I

    % actuator dynamics -> transformed inputs
    % U = K*(u.^2);
    U = K*u;

    % evaluate rotation matrices
    R = R(z(4),z(5),z(6));
    L = L(z(4),z(5));

    %% equations of motion
    % translational state equations
    dz = zeros(12,1);
    dz(1) = z(7);
    dz(2) = z(8);
    dz(3) = z(9);
    dz(7:9) = U(1)*R(:,3)/m;
    dz(9) = dz(9) - g;

    % rotational state equations
    % euler rates obtained from body rates
    dz(4:6) = L\z(10:12);
    dz(10) = (z(11)*z(12)*(I(2)-I(3)) + U(2))/I(1);
    dz(11) = (z(10)*z(12)*(I(3)-I(1)) + U(3))/I(2);
    dz(12) = (z(10)*z(11)*(I(2)-I(3)) + U(4))/I(3);
end