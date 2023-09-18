% This function represents the control law that computes the input
% based on full state feedback
function u = control(t,z,c,p)

    %   document here the input outputs, see the "dynamic.m" file for
    %   reference

    % inner loop states
    zi = [z(3) z(4) z(5) z(6) z(9) z(10) z(11) z(12)];
        
    % unpack control gains
    kz = c.kz;
    kdz = c.kdz;
    kp = c.kp;
    kdp = c.kdp;
    kq = c.kq;
    kdq = c.kdq;
    kr = c.kr;
    kdr = c.kdr;
    
    % error term
    ez = zi(1) - c.zd(1);
    edz = zi(5) - 0;
    ep = zi(2) - c.zd(2);
    edp = zi(6) - 0;
    eq = zi(3) - c.zd(3);
    edq = zi(7) - 0;
    er = zi(4) - c.zd(4);
    edr = zi(8) - 0;
    
    % PD control 
    u = zeros(4,1);
    u(1) = kz*ez + kdz*edz;
    u(2) = kp*ep + kdp*edp;
    u(3) = kq*eq + kdq*edq;
    u(4) = kr*er + kdr*edr;

    % saturation
    u = p.K\(-u);

    % example dummy control law 
    % u = repmap(c.kp*z(1) + c.kd*z(2),4);
end

