% ENPM 673 Project 2 -  Visual Odometry
% Function to Estimate Rotation and translation
% Author : Pranav Inani
function [Rnew,tnew] = EstimateRotTrans(sPts1,sPts2,E,K)
    
    [U,~,V] = svd(E);

    W = [0 -1 0;
         1  0 0;
         0  0 1];
    % Compute 4 possible H matrices 
    H_choices(:,:,1) = [ U*W*V'   U(:,3) ; 0 0 0 1];
    H_choices(:,:,2) = [ U*W*V'  -U(:,3) ; 0 0 0 1];
    H_choices(:,:,3) = [ U*W'*V'  U(:,3) ; 0 0 0 1];
    H_choices(:,:,4) = [ U*W'*V' -U(:,3) ; 0 0 0 1];

    % Ensure the rotaions are legal by making the
    % deteminants of the rotation matrices positive
    for k=1:4
        if det(H_choices(1:3,1:3,k)) < 0
            H_choices(1:3,1:3,k) = -H_choices(1:3,1:3,k);
        end
    end

    H_posZ = [];

    % Throw away H matrices with negative z componenet
    for i = 1: 4
        if (H_choices(3,4,i) > 0)
            H_posZ = cat( 3, H_posZ, H_choices(:,:,i));
        end
    end

    % 3-d reconstruction: 
    % take two point correspondences and check that
    % the z componenets of both their reprojections is positive


    M1 = [ 1 0 0 0;
           0 1 0 0;
           0 0 1 0];
       
    H_est = [];
 
    % Continue to reproject if 3-d reconstruction fails for
    % one set of point correspondances
    for t=1:length(sPts1)
        if isempty(H_est) == false
            break;
        end
    
    pt1 = inv(K) * [sPts1(t,:) 1]';
    pt2 = inv(K) * [sPts2(t,:) 1]';

    p1_skew = [ 0 -pt1(3,1) pt1(2,1);
               -pt1(3,1)  0 -pt1(1,1);
             -pt1(2,1) pt1(1,1) 0];

    p2_skew = [ 0 -pt2(3,1) pt2(2,1);
               -pt2(3,1) 0 -pt2(1,1);
             -pt2(2,1) pt2(1,1) 0];
    
    for j=1:2
        H_posZ_1 = inv(H_posZ(:,:,j));
        M2 = H_posZ_1(1:3,1:4);
        A = [ p1_skew * M1; p2_skew * M2 ];
        [~,~,V] = svd(A);
        P = V(:,4);                     % get last column of V
        P1est = P/P(4)     ;            % normalize

        P2est = H_posZ_1 * P1est;

        if P1est(3) >= 0 && P2est(3) >= 0
            % Solution Found. Break
             H_est = H_posZ(:,:,j);
            break;      
        end
    end
    end
    
    % If the 3-d reconstruction fails, 
    % pass the identity as rotation and 
    % 0 vector as translation
    if isempty(H_est)
        H_est = [1 0 0 0;
                 0 1 0 0;
                 0 0 1 0;
                 0 0 0 1];
    end
    % Return best choices of Rotation and translation
    Rnew = H_est(1:3,1:3);
    tnew = H_est(1:3,4);
end