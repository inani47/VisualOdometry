% ENPM 673 Project 2 -  Visual Odometry
% Function to Estimate Fundamental Matrix
% Author : Pranav Inani
function E = EstimateEssentialMatrix(sPts1,sPts2,K)

% Preprocessing of Points

sz = size(sPts1,1);
sPts1_x = sPts1(:,1);
sPts2_x = sPts2(:,1);
sPts1_y = sPts1(:,2);
sPts2_y = sPts2(:,2);

% Scale and translate image points so that the centroid of
% the points is at the origin, and the average distance of the points to the
% origin is equal to sqrt(2).

centroid_sPts1_x = mean(sPts1_x);
centroid_sPts1_y = mean(sPts1_y);
sPts1_x = sPts1_x - centroid_sPts1_x * ones(sz,1);
sPts1_y = sPts1_y - centroid_sPts1_y * ones(sz,1);
avg_dist = sqrt(sum(sPts1_x.^2  + sPts1_y.^2)) / sz;
scale_factor = sqrt(2) / avg_dist;
sPts1(:,1) = scale_factor * sPts1_x;
sPts1(:,2) = scale_factor * sPts1_y;
T1 = [scale_factor 0 -scale_factor*centroid_sPts1_x;
       0 scale_factor -scale_factor*centroid_sPts1_y;
       0 0 1];  

cent_x2_x = mean(sPts2_x);
cent_x2_y = mean(sPts2_y);
sPts2_x = sPts2_x - cent_x2_x * ones(sz,1);
sPts2_y = sPts2_y - cent_x2_y * ones(sz,1);
avg_dist = sqrt(sum(sPts2_x.^2  + sPts2_y.^2)) / sz;
scale_factor = sqrt(2) / avg_dist;
sPts2(:,1) = scale_factor * sPts2_x;
sPts2(:,2) = scale_factor * sPts2_y;
T2 = [scale_factor 0 -scale_factor*cent_x2_x;
       0 scale_factor -scale_factor*cent_x2_y;
       0 0 1];
   
% Compute Fundamental Matrix F from point correspondences.
A = [sPts1(:,1).*sPts2(:,1) sPts1(:,1).*sPts2(:,2) sPts1(:,1) sPts1(:,2).* sPts2(:,1) sPts1(:,2) .* sPts2(:,2) sPts1(:,2) sPts2(:,1) sPts2(:,2) ones(8,1)];
% A = [sPts2(:,1).*sPts1(:,1) sPts2(:,1).*sPts1(:,2) sPts2(:,1) sPts2(:,2).*sPts1(:,1) sPts2(:,2).*sPts1(:,2) sPts2(:,2) sPts1(:,1) sPts1(:,2) ones(8,1)];
[~,~,V] = svd(A);
x = V(:,end);
% Transpose because reshape arranges elements in column-first order
F = reshape(x,3,3)';

% Enforce Rank 2 Condition on Fundamental Matrix
[U,D,V] = svd(F);
D(3,3) = 0;
F_r2 = U*D*V';

% Undo Preprocessing
F_est = T1' * F_r2 * T2;

F_est = F_est/norm(F_est) ;

% Estimate E
E_est = K'*F_est*K;
[u,~,v] = svd(E_est);
% Enforce Rank 2 on E Condition and make eigen values equal to 1
E = u*diag([1,1,0])*v';





end