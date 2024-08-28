function animate_with_cylinder(t, quat_history)
% animation for attitude visualization given quaternion history
N = length(quat_history);
r = 0.104;
h = 1.5;

[Z, Y, X] = cylinder(r);
X = X*h;
X = X-h/2;
hfg = figure;
hax = axes(hfg);
C = ones(2, 21, 3);
hp = surf(X, Y, Z, C);
hold on

E_axis = plot3([0 2], [0 0], [0 0], 'k--', 'LineWidth', 1.5);
N_axis = plot3([0 0], [0 2], [0 0], 'm--', 'LineWidth', 1.5);
U_axis = plot3([0 0], [0 0], [0 2], 'y--', 'LineWidth', 1.5);

hpx = plot3([0 1.2], [0 0], [0 0], 'r-', 'LineWidth', 1.5);
hpy = plot3([0 0], [0 0.5], [0 0], 'g-', 'LineWidth', 1.5);
hpz = plot3([0 0], [0 0], [0 0.5], 'b-', 'LineWidth', 1.5);
daspect([1 1 1])
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
legend("TU-1.f", "East", "North", "Up", "xb", "yb", "zb")


for i=1:5:N
    q = quat_history(:, i); % attitude at time t(i)
    dcm = quat_2_dcm(q);
    dcm = transpose(dcm);
    X_rot = zeros(2, 21);
    Y_rot = zeros(2, 21);
    Z_rot = zeros(2, 21);
    for j=1:2
        for k=1:21
            rotated_pt = dcm * [X(j, k); Y(j, k); Z(j, k)];
            X_rot(j, k) = rotated_pt(1);
            Y_rot(j, k) = rotated_pt(2);
            Z_rot(j, k) = rotated_pt(3);
        end
    end
    rotated_xb = dcm * [1.2; 0; 0];
    rotated_yb = dcm * [0; 0.5; 0];
    rotated_zb = dcm * [0; 0; 0.5];

    set(hp, "XData", X_rot, "YData", Y_rot, "ZData", Z_rot)
    set(hpx, "XData", [0 rotated_xb(1)], "YData", [0 rotated_xb(2)], "ZData", [0 rotated_xb(3)])
    set(hpy, "XData", [0 rotated_yb(1)], "YData", [0 rotated_yb(2)], "ZData", [0 rotated_yb(3)])
    set(hpz, "XData", [0 rotated_zb(1)], "YData", [0 rotated_zb(2)], "ZData", [0 rotated_zb(3)])
    
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
        
    title(hax, "[STATE: " + STATE + "], Time after boot up: " + num2str(t(i), '%.2f') + " sec")
    pause(0.01)
end


end