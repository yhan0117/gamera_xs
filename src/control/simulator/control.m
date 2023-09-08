% This function represents the control law that computes the input
% based on full state feedback
function u = control(t,z,c)

    % TODO:
    %   figure out control law
    %   document here the input outputs, see the "dynamic.m" file for
    %   reference

    % example dummy control law 
    u = repmap(c.kp*z(1) + c.kd*z(2),4);
end

