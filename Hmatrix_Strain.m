function h = Hmatrix_Strain(Bf, Br, detJ)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the element h matrix (Connesson et al. 2015)
%
% Inputs: 1) B matrix (strain matrix)
%         2) detJ - element volume (determinant of Jacobian matrix)
%
% Output: h matrix - assume constant variance in strain noise
%
% Renee Miller
% 8 March 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate H matrix for current element
% V(Shear) ~ E(fg^2(e,eVF)) = (1/2)*stdU^2*UeVF'*H*UeVF

% t = 1/3 * Tr(e*)
t = 1/3 * sum(Br(1:3,:),1);

% Calculate each component of h
h1 = (Bf(1,:) - t)'*(Bf(1,:) - t);
h2 = (Bf(2,:) - t)'*(Bf(2,:) - t);
h3 = (Bf(3,:) - t)'*(Bf(3,:) - t);
h12 = 4*(0.5*Bf(4,:))'*(0.5*Bf(4,:));
h13 = 4*(0.5*Bf(5,:))'*(0.5*Bf(5,:));
h23 = 4*(0.5*Bf(6,:))'*(0.5*Bf(6,:));

h = (detJ^2)*(h1 + h2 + h3 + h12 + h13 + h23); %Hg