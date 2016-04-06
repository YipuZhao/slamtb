function [score] = computeSOMMatrix(Obs, Sen, Lmk, triplet, Xv, F_Q, F_Omg, dt  )
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

% Compute the H matrix

compute_H_subblock(Xv, Y, Z, H13, H47);

            %std::cout << "[OBS_COMPUTOR]  LinObsMat" << std::endl;
            arma::mat LinObsMat = arma::zeros<arma::mat>(26, 13);
            % Copy first 3 columns
            LinObsMat.cols(0, 2) = arma::repmat(H13, 13, 1);
            % Columns 7~9
            LinObsMat.cols(7,9) = LinObsMat.cols(0, 2) % H13_mul_fac;
            if(fabs(dt - 1.0) > 1e-6) {
                LinObsMat.cols(7,9) = LinObsMat.cols(7,9) * dt;
            }
            % First segment in linObsMat: just H:
            LinObsMat(arma::span(0,1), arma::span(3,6))= H47;

            % Segment 2 to 13
            arma::mat rollingQAdd= arma::eye<arma::mat>(4, 4);  % (Q^(n-1) + Q^(n-2) +... + I)
            arma::mat rollingQFac = F_Q;  % Q^n
            for (int j = 1; j < 13; j++) {
                    % [0,1,2,   3,4,5,6,   7,8,9,   10,11,12 ]

                    % 2nd block:  (H47 * Q^n)
                    LinObsMat( arma::span(j*2, ((j+1)*2-1)), arma::span(3,6)) = H47 * rollingQFac;

                    % 4th block: Q^(n-1) + Q^(n-2) +... + I) * Omega
                    LinObsMat( arma::span(j*2, ((j+1)*2-1)), arma::span(10,12)) = H47 * (rollingQAdd * F_Omg);

                    % Update
                    rollingQAdd = rollingQAdd + rollingQFac;
                    rollingQFac = rollingQFac * F_Q;
            }

            % Update the SOM
            arma::mat SOM_old = pMP->ObsMat;
            int Len_SOM_old = static_cast<int>(SOM_old.n_rows);

            assert(Len_SOM_old% 26 == 0);
            if(Len_SOM_old == 0) {
                %/ Case 1: first time observed, the SOM_old is an empty matrix
                pMP->ObsMat  = LinObsMat;
                pMP->ObsScore= 0.0;
                pMP->ObsRank = arma::rank(LinObsMat);
            }
            else if(Len_SOM_old == 26*tmprLen) {
                %/ Case 2: phase 4, "fully loaded"
                % SOM block order: row 0 -- oldest, row End -- latest.

                arma::mat SOM = arma::zeros<arma::mat>(SOM_old.n_rows, SOM_old.n_cols);
                SOM.rows(0, 26*(tmprLen-1)-1) = SOM_old.rows(26, 26*tmprLen-1);
                SOM.rows(26*(tmprLen-1), 26*tmprLen-1) = LinObsMat;

                % Save
                arma::vec s = arma::svd(SOM);
                pMP->ObsMat  = SOM;
                pMP->ObsScore= s(12);
                pMP->ObsRank = arma::rank(SOM);
            }
            else if(Len_SOM_old< 26*tmprLen) {
                %/ Case 3: phase 2/3, "NOT fully loaded, but has been initialized"
                arma::mat SOM = arma::join_vert(SOM_old, LinObsMat);
                pMP->ObsMat  = SOM;
                pMP->ObsRank = arma::rank(SOM);


                if (SOM.n_rows >= 32) {
                    % Full rank
                    arma::vec s = arma::svd(SOM);
                     pMP->ObsScore= s(12);
                } else {
                     pMP->ObsScore = 0.0;
                }
            }

        } %Perform regular temporal obs update!
     } % For: all pMP
     
toc