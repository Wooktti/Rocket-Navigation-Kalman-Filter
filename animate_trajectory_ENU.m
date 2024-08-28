function animate_trajectory_ENU(t, ENU_pos_history, quat_history)
% animation for trajectory and attitude visualization given pos and
% attitude history

% Settings
ENU_axis_length = 50; % m
body_axis_length = 20; % m

xE = ENU_pos_history(1, :);
xN = ENU_pos_history(2, :);
xU = ENU_pos_history(3, :);

N = length(t);
r = 3;
h = 20;

[Z, Y, X] = cylinder(r);
X = X*h;
X = X-h/2;
hfg = figure;

set(hfg, 'Position',  [0, 0, 1920, 1080])

hax = axes(hfg);
C = ones(2, 21, 3);
hp = surf(X, Y, Z, C, 'DisplayName', 'TU-1.f Rocket');
hold on
grid on

hp_trajectory = plot(0, 0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Trajectory');

E_axis = plot3([0 ENU_axis_length], [0 0], [0 0], 'r--', 'LineWidth', 1.5, 'DisplayName', 'E axis');
N_axis = plot3([0 0], [0 ENU_axis_length], [0 0], 'g--', 'LineWidth', 1.5, 'DisplayName', 'N axis');
U_axis = plot3([0 0], [0 0], [0 ENU_axis_length], 'b--', 'LineWidth', 1.5, 'DisplayName', 'U axis');

hpx = plot3([0 body_axis_length], [0 0], [0 0], 'r-', 'LineWidth', 1.5, 'DisplayName', 'Body X axis');
hpy = plot3([0 0], [0 body_axis_length], [0 0], 'g-', 'LineWidth', 1.5, 'DisplayName', 'Body Y axis');
hpz = plot3([0 0], [0 0], [0 body_axis_length], 'b-', 'LineWidth', 1.5, 'DisplayName', 'Body Z axis');

legend show

xlim([min(xE)-30, max(xE)+30]);
ylim([min(xN)-30, max(xN)+30]);
zlim([min(xU)-30, max(xU)+30]);

xlabel("pos E [m]")
ylabel("pos N [m]")
zlabel("pos U [m]")

fontsize(14,"points")

daspect([1 1 1])

for i=1:1:N
    % if mod(i, 5) ~= 0
    %     continue;
    % end

    q = quat_history(:, i); % attitude at time t(i)
    dcm = quat_2_dcm(q);
    dcm = transpose(dcm);
    X_rot = zeros(2, 21);
    Y_rot = zeros(2, 21);
    Z_rot = zeros(2, 21);
    for j=1:2
        for k=1:21
            rotated_pt = dcm * [X(j, k); Y(j, k); Z(j, k)];
            X_rot(j, k) = rotated_pt(1) + xE(i);
            Y_rot(j, k) = rotated_pt(2) + xN(i);
            Z_rot(j, k) = rotated_pt(3) + xU(i);
        end
    end
    rotated_xb = dcm * [body_axis_length; 0; 0] + [xE(i); xN(i); xU(i)];
    rotated_yb = dcm * [0; body_axis_length; 0] + [xE(i); xN(i); xU(i)];
    rotated_zb = dcm * [0; 0; body_axis_length] + [xE(i); xN(i); xU(i)];
    
    set(hp_trajectory, 'XData', xE(1:i), 'YData', xN(1:i), 'ZData', xU(1:i))
    set(hp, "XData", X_rot, "YData", Y_rot, "ZData", Z_rot)
    set(hpx, "XData", [xE(i) rotated_xb(1)], "YData", [xN(i) rotated_xb(2)], "ZData", [xU(i) rotated_xb(3)])
    set(hpy, "XData", [xE(i) rotated_yb(1)], "YData", [xN(i) rotated_yb(2)], "ZData", [xU(i) rotated_yb(3)])
    set(hpz, "XData", [xE(i) rotated_zb(1)], "YData", [xN(i) rotated_zb(2)], "ZData", [xU(i) rotated_zb(3)])
    
    if mod(i, 30) == 0
        plot3([xE(i) rotated_xb(1)], [xN(i) rotated_xb(2)], [xU(i) rotated_xb(3)], 'r-', 'HandleVisibility','off')
        plot3([xE(i) rotated_yb(1)], [xN(i) rotated_yb(2)], [xU(i) rotated_yb(3)], 'g-', 'HandleVisibility','off')
        plot3([xE(i) rotated_zb(1)], [xN(i) rotated_zb(2)], [xU(i) rotated_zb(3)], 'b-', 'HandleVisibility','off')
    end

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
        
    title(hax, "[STATE: " + STATE + "], Time after launch detection: " + num2str(t(i) - 528.779, '%.2f') + " sec")
    pause(0.1)
end

end

