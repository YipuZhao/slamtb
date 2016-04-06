function [H13, H47] = compute_H_subblock (Xv, yi, zi)
%% 
rw = Xv(1:3, :);
qwr = Xv(4:7, :);
Rrw = inv(q2r(qwr));
RelRw = yi - rw;

% dhd_dhu
ud = zi(1,1);
vd = zi(1,2);
xd = (zi(1,1) - Cx) * dx;
yd = (zi(1,2) - Cy) * dy;
rd2 = (xd * xd) + (yd * yd);
rd4 = rd2 * rd2;
uu_ud = (1+k1*rd2+k2*rd4)+(ud-Cx)*(k1+2*k2*rd2)*(2*(ud-Cx)*dx*dx);
vu_vd = (1+k1*rd2+k2*rd4)+(vd-Cy)*(k1+2*k2*rd2)*(2*(vd-Cy)*dy*dy);
uu_vd = (ud-Cx)*(k1+2*k2*rd2)*(2*(vd-Cy)*dy*dy);
vu_ud = (vd-Cy)*(k1+2*k2*rd2)*(2*(ud-Cx)*dx*dx);

J_undistor = [
    uu_ud uu_vd
    vu_ud vu_vd
    ];
dhd_dhu = inv(J_undistor);

% dhu_dhrl
hrl = Rrw * RelRw;
if ( abs(hrl(3,1)) < 1e-6 )
    dhu_dhrl  = zeros(2,3);
else
    dhu_dhrl = [
        f*ku/(hrl(3,1)) 0 -hrl(1,1)*f*ku/( pow(hrl(3,1), 2.0))
        0 f*kv/(hrl(3,1)) -hrl(2,1)*f*kv/( pow(hrl(3,1), 2.0))
        ];
end

dh_dhrl = dhd_dhu *dhu_dhrl;
qwr_conj = Qconj(qwr);
dhrl_dqwr = dRq_times_a_by_dq( qwr_conj ,  RelRw) * dqbar_by_dq;

% H matrix subblock (cols 1~3): H13
H13 = -1.0 * (dh_dhrl *  Rrw);

% H matrix subblock (cols 4~7): H47
H47 = dh_dhrl * dhrl_dqwr;

