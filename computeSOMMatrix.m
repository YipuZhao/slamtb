function [SOM] = computeSOMMatrix(Obs, Sen, Lmk, triplet, Xv, F_Q, F_Omg, dt  )
%%
% triplet are the indices of selected landmarks in Lmk
global Map

tic;
%    %-----------------------------------------------------------------------------
%    % Compute H and update Obs information for each feature

% Feature position
Y = zeros(3,1);
lr = Lmk(triplet(1)).state.r;        % lmk range in Map
featPos  = Map.x(lr) ;               % lmk mean
Y(1) = featPos(1);
Y(2) = featPos(2);
Y(3) = featPos(3);

% Measurement
Z = zeros(1,2);
kpUn = Obs(Sen.sen, triplet(1)).meas;
Z(1,1) = kpUn.y(1);
Z(1,2) = kpUn.y(2);

% Compute the H matrix based on a triplet (different from the GF-SLAM)
[H13, H47] = compute_H_subblock(Xv, Y, Z);

%std::cout << "[OBS_COMPUTOR]  LinObsMat" << std::endl;
LinObsMat = zeros(26, 13);

% Copy first 3 columns
LinObsMat(:, 1:3) = repmat(H13, 13, 1);
% Columns 7~9
LinObsMat(:, 8:10) = LinObsMat(:, 1:3); % H13_mul_fac;
if(abs(dt - 1.0) > 1e-6)
    LinObsMat(:, 8:10) = LinObsMat(:, 8:10) * dt;
end

% First segment in linObsMat: just H:
LinObsMat(1:2, 4:7)= H47;

% Segment 2 to 13
rollingQAdd= eye(4, 4);  % (Q^(n-1) + Q^(n-2) +... + I)
rollingQFac = F_Q;  % Q^n
for j = 1:13
    % [0,1,2,   3,4,5,6,   7,8,9,   10,11,12 ]
    
    % 2nd block:  (H47 * Q^n)
    LinObsMat(j*2+1:(j+1)*2, 4:7) = H47 * rollingQFac;
    
    % 4th block: Q^(n-1) + Q^(n-2) +... + I) * Omega
    LinObsMat(j*2+1:(j+1)*2, 11:13) = H47 * (rollingQAdd * F_Omg);
    
    % Update
    rollingQAdd = rollingQAdd + rollingQFac;
    rollingQFac = rollingQFac * F_Q;
end

SOM = LinObsMat;
