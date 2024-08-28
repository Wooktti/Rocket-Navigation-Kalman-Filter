function [x_curr, P_curr] = KalmanFilter(dt, x_prev, u_prev, P_prev, z_curr)
   
    % x: [xE; xN; xU; vE; vN; vU] in m and m/s
    % z: measurement of xU (i.e. AGL altitude) in m
    % u: [accE; accN, accU] in m/s2
    % y: measurement residual

    % P: 6 x 6 Covariance of the state
    % A: 6 x 6 state transition matrix
    % B: 6 x 3 control input matrix
    % Q: 6 x 6 Estimated Process Error Covariance
    % H: 1 x 6 observation matrix
    % R: 1 x 1 Estimated Measurement Error Covariance

    A = [1, 0, 0, dt, 0, 0;
         0, 1, 0, 0, dt, 0;
         0, 0, 1, 0, 0, dt;
         0, 0, 0, 1, 0, 0;
         0, 0, 0, 0, 1, 0;
         0, 0, 0, 0, 0, 1];

    B = [dt^2 / 2, 0, 0;
         0, dt^2 / 2, 0;
         0, 0, dt^2 / 2;
         dt, 0, 0;
         0, dt, 0;
         0, 0, dt];

    H = [0, 0, 1, 0, 0, 0]; 

    Q = 0.0001 * eye(6);
    
    R = 0.092 * 1.3; % standard deviation of barometer altitude
        
    % Prediction
    x_predict = A * x_prev + B * u_prev; % predict x
    P_predict = A * P_prev * A' + Q; % predict P matrix
    
    % Correction
    y_curr = z_curr - H * x_predict;
    S_k = H * P_predict * H' + R; 

    K = P_predict * H' * inv(S_k); % calculate Kamlan gain matrix
    x_curr = x_predict + K * y_curr; % correct x
    P_curr = (eye(6) - K * H) * P_predict; % update P matrix
    P_curr = P_predict - K * S_k * K';
end