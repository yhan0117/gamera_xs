function animate(t_s,z_,p,filename,record,fps)
    hold on % <-- very important
    
    % extract states and rotation matrix
    R = @(i) p.R(z_(i,4),z_(i,5),z_(i,6));
    x = z_(:,1);
    y = z_(:,2);
    z = z_(:,3);
    
    % axis and labels
    title(sprintf('Trajectory\nTime: %0.2f sec', t_s(1)), 'Interpreter', 'Latex');
    xlabel('x', 'Interpreter', 'Latex')
    ylabel('y', 'Interpreter', 'Latex')
    zlabel('z', 'Interpreter', 'Latex')
    
    grid minor; axis equal; rotate3d on;
    xlim([-1 1]); ylim([-1 1]); zlim([-1 1]); 

    % plotting with no color to set axis limits
    plot3(x,y,z,'Color','none');
    
    % plotting the first iteration
    p = plot3(x(1),y(1),z(1),'b');
    m = scatter3(x(1),y(1),z(1),'filled','b','square');
    
    % draw the quad
    l = 0.125;
    r = 0.04;
    h = 0.02;
    r_mg = [-l 0 0;0 -l 0;l 0 0;0 l 0];
    mtr = R(1)*r_mg.'+ [x(1) y(1) z(1)].';
    v = [r -r -h; r r -h; -r r -h; -r -r -h;r -r h; r r h; -r r h; -r -r h];
    V1 = R(1)*v.'+ mtr(:,1);
    V2 = R(1)*v.'+ mtr(:,2);
    V3 = R(1)*v.'+ mtr(:,3);
    V4 = R(1)*v.'+ mtr(:,4);
    
    face = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
    box1 = patch('Vertices',V1.','Faces',face,'FaceColor',[0.1 0.1 0.8],'FaceAlpha',0.2);
    box2 = patch('Vertices',V2.','Faces',face,'FaceColor',[0.1 0.1 0.8],'FaceAlpha',0.2);
    box3 = patch('Vertices',V3.','Faces',face,'FaceColor',[0.8 0.1 0.1],'FaceAlpha',0.2);
    box4 = patch('Vertices',V4.','Faces',face,'FaceColor',[0.1 0.8 0.1],'FaceAlpha',0.2);
    
    arm1 = plot3([mtr(1,1), mtr(1,3)], [mtr(2,1), mtr(2,3)], [mtr(3,1), mtr(3,3)],'k', "LineWidth", 0.5);
    arm2 = plot3([mtr(1,2), mtr(1,4)], [mtr(2,2), mtr(2,4)], [mtr(3,2), mtr(3,4)],'k', "LineWidth", 0.5);
    
    
    % iterating through the length of the time array
    for k = 1:numel(t_s)
        % updating the trajectory trace
        p.XData = x(1:k);
        p.YData = y(1:k);
        p.ZData = z(1:k);
    
        % updating quad
        mtr = R(k)*r_mg.'+ [x(k) y(k) z(k)].';
        V1 = R(k)*v.'+ mtr(:,1);
        V2 = R(k)*v.'+ mtr(:,2);
        V3 = R(k)*v.'+ mtr(:,3);
        V4 = R(k)*v.'+ mtr(:,4);
    
        set(box1,'Faces',face,'Vertices',V1.','FaceColor',[0.1 0.1 0.8]);
        set(box2,'Faces',face,'Vertices',V2.','FaceColor',[0.1 0.1 0.8]);
        set(box3,'Faces',face,'Vertices',V3.','FaceColor',[0.8 0.1 0.1]);
        set(box4,'Faces',face,'Vertices',V4.','FaceColor',[0.1 0.8 0.1]);
        
        arm1.XData = [mtr(1,1), mtr(1,3)];
        arm2.XData = [mtr(1,2), mtr(1,4)];
        arm1.YData = [mtr(2,1), mtr(2,3)];
        arm2.YData = [mtr(2,2), mtr(2,4)];
        arm1.ZData = [mtr(3,1), mtr(3,3)];
        arm2.ZData = [mtr(3,2), mtr(3,4)];
    
        % updating the title
        title(sprintf('Trajectory\nTime: %0.2f sec', t_s(k)),'Interpreter','Latex');
    
        % saving the figure
        if record
            frame = getframe(gcf);
        
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if k == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime', 1/fps);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', 1/fps);
            end
        else
            pause(0.01)
        end
    end
end

