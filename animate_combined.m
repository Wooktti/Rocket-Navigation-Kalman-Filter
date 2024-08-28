function animate_combined(t, ENU_pos_history, quat_history, viewpoint)
% animation for trajectory, attitude.
N = length(t);
N_pos = length(ENU_pos_history);

nanArray = nan * ones(3, N - N_pos);
ENU_pos_history = [ENU_pos_history, nanArray];


xE = ENU_pos_history(1, :);
xN = ENU_pos_history(2, :);
xU = ENU_pos_history(3, :);

hfg = figure; % use single window
set(hfg, 'Position',  [0, 0, 1920, 1080])

%% Setup for trajectory animation plot
hax1 = axes('OuterPosition', [0, 0, 0.5, 0.9]);

ENU_axis_length = 50; % m
body_axis_length = 20; % m

r = 3;
h = 20;
r = 0.104;
h = 2;
[Z1, Y1, X1] = cylinder(r);
X1 = X1*h;
X1 = X1-h/2;

C1 = ones(2, 21, 3);
hp_vehicle_1 = surf(hax1, X1, Y1, Z1, C1, 'DisplayName', 'TU-1.f Rocket');
hold(hax1, 'on')
grid(hax1, 'on')

hp_trajectory = plot(hax1, 0, 0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Trajectory');

plot3(hax1, [0 ENU_axis_length], [0 0], [0 0], 'r--', 'LineWidth', 1.5, 'DisplayName', 'E axis');
plot3(hax1, [0 0], [0 ENU_axis_length], [0 0], 'g--', 'LineWidth', 1.5, 'DisplayName', 'N axis');
plot3(hax1, [0 0], [0 0], [0 ENU_axis_length], 'b--', 'LineWidth', 1.5, 'DisplayName', 'U axis');

hpx1 = plot3(hax1, [0 body_axis_length], [0 0], [0 0], 'r-', 'LineWidth', 1.5, 'DisplayName', 'Body X axis');
hpy1 = plot3(hax1, [0 0], [0 body_axis_length], [0 0], 'g-', 'LineWidth', 1.5, 'DisplayName', 'Body Y axis');
hpz1 = plot3(hax1, [0 0], [0 0], [0 body_axis_length], 'b-', 'LineWidth', 1.5, 'DisplayName', 'Body Z axis');

legend(hax1, 'off')

xlim(hax1, [min(xE)-30, max(xE)+30]);
ylim(hax1, [min(xN)-30, max(xN)+30]);
zlim(hax1, [min(xU)-30, max(xU)+30]);

xlabel(hax1, "pos E [m]")
ylabel(hax1, "pos N [m]")
zlabel(hax1, "pos U [m]")

fontsize(hax1, 14, "points")
daspect(hax1, [1 1 1])
title(hax1, "Trajectory")
view(hax1, viewpoint)

%% Setup for attitude animation plot
hax2 = axes('OuterPosition', [0.5, 0.45, 0.25, 0.45]);

r = 0.104;
h = 1.5;

[Z2, Y2, X2] = cylinder(r);
X2 = X2*h;
X2 = X2-h/2;

C2 = ones(2, 21, 3);
hp_vehicle_2 = surf(hax2, X2, Y2, Z2, C2, 'DisplayName', 'TU-1.f Rocket');
hold(hax2, 'on')
grid(hax2, 'on')

plot3(hax2, [0 2], [0 0], [0 0], 'r--', 'LineWidth', 1.5, 'DisplayName', 'E axis');
plot3(hax2, [0 0], [0 2], [0 0], 'g--', 'LineWidth', 1.5, 'DisplayName', 'N axis');
plot3(hax2, [0 0], [0 0], [0 2], 'b--', 'LineWidth', 1.5, 'DisplayName', 'U axis');

hpx2 = plot3(hax2, [0 1.2], [0 0], [0 0], 'r-', 'LineWidth', 1.5, 'DisplayName', 'Body X axis');
hpy2 = plot3(hax2, [0 0], [0 0.5], [0 0], 'g-', 'LineWidth', 1.5, 'DisplayName', 'Body Y axis');
hpz2 = plot3(hax2, [0 0], [0 0], [0 0.5], 'b-', 'LineWidth', 1.5, 'DisplayName', 'Body Z axis');

legend(hax2, 'off')

xlim(hax2, [-1 1])
ylim(hax2, [-1 1])
zlim(hax2, [-1 1])

fontsize(hax2, 14, "points")
daspect(hax2, [1 1 1])
title(hax2, "Attitude")

%% Setup for 2D trajectory animation plot
% EN Plane
hax_EN = axes('OuterPosition', [0.75, 0.45, 0.25, 0.45]);
hold(hax_EN, 'on')
grid(hax_EN, 'on')
hp_trajectory_EN = plot(hax_EN, 0, 0, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Trajectory');

legend(hax_EN, 'off')

xlim(hax_EN, [min(xE), max(xE)]);
ylim(hax_EN, [min(xN), max(xN)]);

xlabel(hax_EN, "pos E [m]")
ylabel(hax_EN, "pos N [m]")

fontsize(hax_EN, 14, "points")
daspect(hax_EN, [1 1 1])
title(hax_EN, "EN (XY) Plane")

% NU Plane
hax_NU = axes('OuterPosition', [0.5, 0, 0.25, 0.45]);
hold(hax_NU, 'on')
grid(hax_NU, 'on')
hp_trajectory_NU = plot(hax_NU, 0, 0, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Trajectory');

legend(hax_NU, 'off')

xlim(hax_NU, [min(xN), max(xN)]);
ylim(hax_NU, [min(xU), max(xU)]);

xlabel(hax_NU, "pos N [m]")
ylabel(hax_NU, "pos U [m]")

fontsize(hax_NU, 14, "points")
daspect(hax_NU, [1 1 1])
title(hax_NU, "NU (YZ) Plane")


% EU Plane
hax_EU = axes('OuterPosition', [0.75, 0, 0.25, 0.45]);
hold(hax_EU, 'on')
grid(hax_EU, 'on')
hp_trajectory_EU = plot(hax_EU, 0, 0, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Trajectory');

legend(hax_EU, 'off')

xlim(hax_EU, [min(xE), max(xE)]);
ylim(hax_EU, [min(xU), max(xU)]);

xlabel(hax_EU, "pos E [m]")
ylabel(hax_EU, "pos U [m]")

fontsize(hax_EU, 14, "points")
daspect(hax_EU, [1 1 1])
title(hax_EU, "EU (XZ) Plane")


%% Animation
for i=1:1:N
    if mod(i, 5) ~= 0
        continue;
    end

    q = quat_history(:, i); % attitude at time t(i)
    dcm = quat_2_dcm(q);
    dcm = transpose(dcm);

    X_rot_trajectory = zeros(2, 21);
    Y_rot_trajectory = zeros(2, 21);
    Z_rot_trajectory = zeros(2, 21);

    X_rot_attitude = zeros(2, 21);
    Y_rot_attitude = zeros(2, 21);
    Z_rot_attitude = zeros(2, 21);

    for j=1:2
        for k=1:21
            % calculate transformed coordinates of cylinder for trajectory
            % plot
            rotated_pt_trajectory = dcm * [X1(j, k); Y1(j, k); Z1(j, k)];
            X_rot_trajectory(j, k) = rotated_pt_trajectory(1) + xE(i);
            Y_rot_trajectory(j, k) = rotated_pt_trajectory(2) + xN(i);
            Z_rot_trajectory(j, k) = rotated_pt_trajectory(3) + xU(i);

            % calculate transformed coordinates of cylinder for attitude
            % plot
            rotated_pt_attitude = dcm * [X2(j, k); Y2(j, k); Z2(j, k)];
            X_rot_attitude(j, k) = rotated_pt_attitude(1);
            Y_rot_attitude(j, k) = rotated_pt_attitude(2);
            Z_rot_attitude(j, k) = rotated_pt_attitude(3);
        end
    end

    rotated_xb_trajectory = dcm * [body_axis_length; 0; 0] + [xE(i); xN(i); xU(i)];
    rotated_yb_trajectory = dcm * [0; body_axis_length; 0] + [xE(i); xN(i); xU(i)];
    rotated_zb_trajectory = dcm * [0; 0; body_axis_length] + [xE(i); xN(i); xU(i)];

    rotated_xb_attitude = dcm * [1.2; 0; 0];
    rotated_yb_attitude = dcm * [0; 0.5; 0];
    rotated_zb_attitude = dcm * [0; 0; 0.5];
    
    % update 3d trajectory plot
    set(hp_trajectory, 'XData', xE(1:i), 'YData', xN(1:i), 'ZData', xU(1:i))
    set(hp_vehicle_1, "XData", X_rot_trajectory, "YData", Y_rot_trajectory, "ZData", Z_rot_trajectory)
    set(hpx1, "XData", [xE(i) rotated_xb_trajectory(1)], "YData", [xN(i) rotated_xb_trajectory(2)], "ZData", [xU(i) rotated_xb_trajectory(3)])
    set(hpy1, "XData", [xE(i) rotated_yb_trajectory(1)], "YData", [xN(i) rotated_yb_trajectory(2)], "ZData", [xU(i) rotated_yb_trajectory(3)])
    set(hpz1, "XData", [xE(i) rotated_zb_trajectory(1)], "YData", [xN(i) rotated_zb_trajectory(2)], "ZData", [xU(i) rotated_zb_trajectory(3)])
    
    if mod(i, 30) == 0
        plot3(hax1, [xE(i) rotated_xb_trajectory(1)], [xN(i) rotated_xb_trajectory(2)], [xU(i) rotated_xb_trajectory(3)], 'r-', 'HandleVisibility','off')
        plot3(hax1, [xE(i) rotated_yb_trajectory(1)], [xN(i) rotated_yb_trajectory(2)], [xU(i) rotated_yb_trajectory(3)], 'g-', 'HandleVisibility','off')
        plot3(hax1, [xE(i) rotated_zb_trajectory(1)], [xN(i) rotated_zb_trajectory(2)], [xU(i) rotated_zb_trajectory(3)], 'b-', 'HandleVisibility','off')
    end

    % update attitude plot
    set(hp_vehicle_2, "XData", X_rot_attitude, "YData", Y_rot_attitude, "ZData", Z_rot_attitude)
    set(hpx2, "XData", [0 rotated_xb_attitude(1)], "YData", [0 rotated_xb_attitude(2)], "ZData", [0 rotated_xb_attitude(3)])
    set(hpy2, "XData", [0 rotated_yb_attitude(1)], "YData", [0 rotated_yb_attitude(2)], "ZData", [0 rotated_yb_attitude(3)])
    set(hpz2, "XData", [0 rotated_zb_attitude(1)], "YData", [0 rotated_zb_attitude(2)], "ZData", [0 rotated_zb_attitude(3)])
    
    % update 2d trajectory plot
    set(hp_trajectory_EN, 'XData', xE(1:i), 'YData', xN(1:i))
    set(hp_trajectory_NU, 'XData', xN(1:i), 'YData', xU(1:i))
    set(hp_trajectory_EU, 'XData', xE(1:i), 'YData', xU(1:i))

    % update title
    if t(i) < 58.230
        STATE = "BOOT";
    elseif t(i) < 528.779
        STATE = "STANDBY";
    elseif t(i) < 531.279
        STATE = "BURN";
    elseif t(i) < 538.780
        STATE = "COAST";
    elseif t(i) < 650
        STATE = "DESCENT (PARACHUTE DEPLOYED) ";
    else
        STATE = "LANDED (RECOVERING VEHICLE) ";
    end
        
    sgtitle("[State: " + STATE + "], Time after launch detection: " + num2str(t(i) - 528.779, '%.2f') + " sec", 'FontSize', 20)
    pause(0.1)
end

end
