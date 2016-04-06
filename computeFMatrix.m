function [F_Q, F_Omg] = computeFMatrix(Xv, dt)
%%-----------------------------------------------------------------------------
% Compute F
%std::cout << "[OBS_COMPUTOR]  Get F" << std::endl;
omegaOld = Xv(10:12, :);
qOld = Xv(3:6, :);
v = omegaOld * dt;
theta = norm(v, 2);
if(theta < 1e-6)
    q = [1,0,0,0];
else
    v_n = v / theta;
    v_n_reshape =  sin(theta / 2.0) * (v_n / norm(v_n, 2));
    q = [cos(theta / 2.0) v_n_reshape(0,0) v_n_reshape(0,1) v_n_reshape(0,2)];
end

% F matrix subblock:  F_Q
R = q(0,0); X = q(0,1); Y = q(0,2); Z = q(0,3);
F_Q = [
    R -X -Y -Z
    X R Z -Y
    Y -Z R X
    Z Y -X R
    ];

% F matrix subblock:   F_Omg
R = qOld(0,0); X = qOld(1,0); Y = qOld(2,0); Z = qOld(3,0);
dq3_by_dq1 = [
    R -X -Y -Z
    X R -Z Y
    Y Z R -X
    Z -Y X R
    ];

F_Omg = dq3_by_dq1 * dqomegadt_by_domega(omegaOld,dt);