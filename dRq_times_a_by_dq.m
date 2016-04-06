function [RES] = dRq_times_a_by_dq(q,  aMat)
%%
assert(aMat.n_rows == 3 && aMat.n_cols == 1);

q0 = q(1,1);
qx = q(2,1);
qy = q(3,1);
qz = q(4,1);

dR_by_dq0 = [
    2.0*q0 -2.0*qz 2.0*qy
    2.0*qz 2.0*q0 -2.0*qx
    -2.0*qy 2.0*qx 2.0*q0
    ];

dR_by_dqx = [
    2.0*qx 2.0*qy 2.0*qz
    2.0*qy -2.0*qx -2.0*q0
    2.0*qz 2.0*q0 -2.0*qx
    ];

dR_by_dqy = [
    -2.0*qy 2.0*qx 2.0*q0
    2.0*qx 2.0*qy 2.0*qz
    -2.0*q0 2.0*qz -2.0*qy
    ];

dR_by_dqz = [
    -2.0*qz -2.0*q0 2.0*qx
    2.0*q0 -2.0*qz 2.0*qy
    2.0*qx 2.0*qy 2.0*qz
    ];

RES = zeros(3,4);
RES(:, 1) = dR_by_dq0 * aMat;
RES(:, 2) = dR_by_dqx * aMat;
RES(:, 3) = dR_by_dqy * aMat;
RES(:, 4) = dR_by_dqz * aMat;
