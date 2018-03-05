function [C1, C2, C3, C4, C5, H] = C3D8F_5p(U, nodesSubZone, elemSubZone, globalNodeNums, GaussPoints, paramInit, orientation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function returns Ak, Ag and H matrices for constructing constraints on 
% virtual displacement field as well as strain (for validating methods with 
% Abaqus output)
%
% Inputs: 1) U - displacement vector 
%         2) nodesSubZone - nodal coordinates (4 x # of nodes)
%         3) elemSubZone - element connectivity (9 x # of elements)
%         4) globalNodeNums - list of global node numbers
%         5) GaussPoints - number of Gauss points (per direction)
%         6) paramInit - initial parameter estimate (1 x 3)
%         7) orientation - element material orientations (three vectors - 
% cross product of first two vectors gives fibre orientation in element)
%
% Outputs: 1) Constraints: C1, C2, C3, C4 and C5
%          2) H - matrix ((# nodes x nodal DOFs) x (# nodes x nodal DOFs))
%
% Note: Full integration, number of Gauss points depends on user input.
%
% Renee Miller
% 28 June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodes per element
nodesPerElem = size(elemSubZone,2) - 1;

% Nodal DOF
DOF = size(nodesSubZone,2) - 1;

% Gauss point coordinates and weights
switch GaussPoints
    case 1
        zeta = 0;
        w = 2;  
        GP = 1;
    case 2
        zeta = [-sqrt(1/3); sqrt(1/3)];
        w = [1; 1];
        GP = 8;
    case 3
        zeta = [-sqrt(3/5); 0; sqrt(3/5)];
        w = [5/9; 8/9; 5/9];
        GP = 27;
end

% Initialise vectors for conditions on the virtual field
C1 = zeros(1,length(nodesSubZone)*DOF);
C2 = zeros(1,length(nodesSubZone)*DOF);
C3 = zeros(1,length(nodesSubZone)*DOF);
C4 = zeros(1,length(nodesSubZone)*DOF);
C5 = zeros(1,length(nodesSubZone)*DOF);

% Variable to save strain at Gauss points (to comopare to Abaqus output)
%strain = zeros(size(elemSubZone,1)*length(zeta)^3,8); % 7 = 1 element label + gauss point + 6 strain components (e11, e22, e33, e12, e23, e13)

% Initialise row and col vectors
row_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);
col_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);
h_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);

% Create waitbar
WH = waitbar(0, 'Assembling global matrices...');

% Loop through each element
for i = 1:size(elemSubZone,1)

% 	count = 0; % make counter to assign strain values
    
    % Update waitbar to give user an indication of time
    percentComplete = i/size(elemSubZone,1);
    waitbar(percentComplete, WH, sprintf('%.2f%% complete...', 100*percentComplete))
    
    % Initialise variables
    C1_1 = 0; C2_1 = 0; C3_1 = 0; C4_1 = 0; C5_1 = 0; h_1 = 0; 
    C1_2 = 0; C2_2 = 0; C3_2 = 0; C4_2 = 0; C5_2 = 0; h_2 = 0; 
    C1_3 = 0; C2_3 = 0; C3_3 = 0; C4_3 = 0; C5_3 = 0; h_3 = 0;
    
    % Get vectors of nodal coordinates
    [X, Y, Z] = getElemCoords(nodesSubZone, elemSubZone, i);
    
    % Get rotation matrix
    L = getRotMat(orientation, elemSubZone(i));
    
    %% Indexing - there are two types of node indexing - local and global. 
    % Local refers to the node index in the subzone and global refers to 
    % the node index in the whole model.
    
    % Get indics of nodal DOF in global list
    globalNodeIdcs = getNodeIdcs(globalNodeNums, elemSubZone, i, DOF);
    
    % Get indics of nodal DOF in global list
    localNodeIdcs = getNodeIdcs(nodesSubZone(:,1), elemSubZone, i, DOF);
       
    % Get displacements for particular element using indices
    Ue = U(globalNodeIdcs);
  
    for j = 1:length(zeta)
        o = zeta(j);
        
        for k = 1:length(zeta)
            n = zeta(k);
            
            for l = 1:length(zeta)
                m = zeta(l);
                
                % Gauss integration coordinate = (m,n,o) %
                
                % Evaluate the derivative of the shape functions at m, n, o
                dNdm = 0.125 * [-1*(1-n)*(1-o) (1-n)*(1-o) (1+n)*(1-o) -1*(1+n)*(1-o) -1*(1-n)*(1+o) (1-n)*(1+o) (1+n)*(1+o) -1*(1+n)*(1+o)];
                dNdn = 0.125 * [-1*(1-m)*(1-o) -1*(1+m)*(1-o) (1+m)*(1-o) (1-m)*(1-o) -1*(1-m)*(1+o) -1*(1+m)*(1+o) (1+m)*(1+o) (1-m)*(1+o)];
                dNdo = 0.125 * [-1*(1-m)*(1-n) -1*(1+m)*(1-n) -1*(1+m)*(1+n) -1*(1-m)*(1+n) (1-m)*(1-n) (1+m)*(1-n) (1+m)*(1+n) (1-m)*(1+n)];
                
                dN = [dNdm; dNdn; dNdo];
                
                % Calculate Jacobian for current element
                jac = dN*[X Y Z];
                
                % Multiply inverse of jacobian times the derivative of shape functions
                dNdXYZ = jac\dN;
                                
                % Determinant of Jacobian
                detJ = det(jac);
                                
                % Calculate B matrix (strain matrix)
                B = [];
                for c = 1:nodesPerElem %Loop through number of nodes per element
                    Bi = [dNdXYZ(1,c)     0              0; ...
                          0            dNdXYZ(2,c)       0; ...
                          0               0         dNdXYZ(3,c); ...
                          dNdXYZ(2,c)  dNdXYZ(1,c)       0; ...
                          dNdXYZ(3,c)     0         dNdXYZ(1,c); ...
                          0            dNdXYZ(3,c)  dNdXYZ(2,c)];
                    B = [B Bi];
                end
                                
                % Calculate strain: B*U
                eV = B*Ue; % Strain of measured displacements using full integration
                                
%                 % Save strain from measured displacement field - compare to Abaqus (CHECK)
%                 count = count + 1;
%                 idx = (i-1)*GP + count;
%                 strain(idx,:) = [i count eV.'];
                
                % Convert strains to square strain tensor
                e = [eV(1) 0.5*eV(4) 0.5*eV(5); ...
                    0.5*eV(4) eV(2) 0.5*eV(6);...
                    0.5*eV(5) 0.5*eV(6) eV(3)];
                
                % Rotate strain tensor
                eRot = L * e * L';
                
                % Put strain tensor back into vector form
                eV = [eRot(1,1); eRot(2,2); eRot(3,3); 2*eRot(1,2); 2*eRot(1,3); 2*eRot(2,3)];
                
                % Calculate rotated B matrix for calculation of virtual displacement fields
                eVF11r = L(1,1)*L(1,1)*B(1,:) + L(1,2)*L(1,1)*0.5*B(4,:) + L(1,1)*L(1,3)*0.5*B(5,:) + L(1,1)*L(1,2)*0.5*B(4,:) + L(1,2)*L(1,2)*B(2,:) + ...
                    L(1,3)*L(1,2)*0.5*B(6,:) + L(1,1)*L(1,3)*0.5*B(5,:) + L(1,2)*L(1,3)*0.5*B(6,:) + L(1,3)*L(1,3)*B(3,:);
                eVF12r = L(2,1)*L(1,1)*B(1,:) + L(2,2)*L(1,1)*0.5*B(4,:) + L(2,3)*L(1,1)*0.5*B(5,:) + L(2,1)*L(1,2)*0.5*B(4,:) + L(2,2)*L(1,2)*B(2,:) + ...
                    L(2,3)*L(1,2)*0.5*B(6,:) + L(2,1)*L(1,3)*0.5*B(5,:) + L(2,2)*L(1,3)*0.5*B(6,:) + L(2,3)*L(1,3)*B(3,:);
                eVF13r = L(3,1)*L(1,1)*B(1,:) + L(3,2)*L(1,1)*0.5*B(4,:) + L(3,3)*L(1,1)*0.5*B(5,:) + L(3,1)*L(1,2)*0.5*B(4,:) + L(3,2)*L(1,2)*B(2,:) + ...
                    L(3,3)*L(1,2)*0.5*B(6,:) + L(3,1)*L(1,3)*0.5*B(5,:) + L(3,2)*L(1,3)*0.5*B(6,:) + L(3,3)*L(1,3)*B(3,:);
                eVF22r = L(2,1)*L(2,1)*B(1,:) + L(2,2)*L(2,1)*0.5*B(4,:) + L(2,3)*L(2,1)*0.5*B(5,:) + L(2,1)*L(2,2)*0.5*B(4,:) + L(2,2)*L(2,2)*B(2,:) + ...
                    L(2,3)*L(2,2)*0.5*B(6,:) + L(2,1)*L(2,3)*0.5*B(5,:) + L(2,2)*L(2,3)*0.5*B(6,:) + L(2,3)*L(2,3)*B(3,:);
                eVF23r = L(3,1)*L(2,1)*B(1,:) + L(3,2)*L(2,1)*0.5*B(4,:) + L(3,3)*L(2,1)*0.5*B(5,:) + L(3,1)*L(2,2)*0.5*B(4,:) + L(3,2)*L(2,2)*B(2,:) + ...
                    L(3,3)*L(2,2)*0.5*B(6,:) + L(3,1)*L(2,3)*0.5*B(5,:) + L(3,2)*L(2,3)*0.5*B(6,:) + L(3,3)*L(2,3)*B(3,:);
                eVF33r = L(3,1)*L(3,1)*B(1,:) + L(3,2)*L(3,1)*0.5*B(4,:) + L(3,3)*L(3,1)*0.5*B(5,:) + L(3,1)*L(3,2)*0.5*B(4,:) + L(3,2)*L(3,2)*B(2,:) + ...
                    L(3,3)*L(3,2)*0.5*B(6,:) + L(3,1)*L(3,3)*0.5*B(5,:) + L(3,2)*L(3,3)*0.5*B(6,:) + L(3,3)*L(3,3)*B(3,:);
                
                % Rotated B matrix - used to derive constraints
                B_rot = [eVF11r; eVF22r; eVF33r; 2*eVF12r; 2*eVF13r; 2*eVF23r];
                              
                % Calculate constraints: c1, c2, c3, c4 and c5
                f_c11 = detJ*( eV(1)*B_rot(1,:) + eV(2)*B_rot(2,:) +  eV(1)*B_rot(2,:) + eV(2)*B_rot(1,:) );
                f_c33 = detJ*( eV(3)*B_rot(3,:) );
                f_c44 = detJ*( eV(4)*B_rot(4,:) - 2*eV(1)*B_rot(2,:) - 2*eV(2)*B_rot(1,:) ); 
                f_c66 = detJ*( eV(5)*B_rot(5,:) + eV(6)*B_rot(6,:) );
                f_c13 = detJ*( eV(1)*B_rot(3,:) + eV(3)*B_rot(1,:) + eV(2)*B_rot(3,:) + eV(3)*B_rot(2,:) );
                                
                                %% Optimisation of virtual fields %%%
                % Calculate the H matrix              
                N11 = (detJ^2)*( 2*B_rot(1,:).'*B_rot(1,:) + 4*B_rot(1,:).'*B_rot(2,:) + 2*B_rot(2,:).'*B_rot(2,:) );
                N22 = (detJ^2)*( B_rot(3,:).'*B_rot(3,:) );
                N33 = (detJ^2)*( B_rot(4,:).'*B_rot(4,:) + 4*B_rot(1,:).'*B_rot(1,:) + 4*B_rot(2,:).'*B_rot(2,:)); 
                N44 = (detJ^2)*( B_rot(5,:).'*B_rot(5,:) + B_rot(6,:).'*B_rot(6,:)); 
                N55 = (detJ^2)*( 2*B_rot(3,:).'*B_rot(3,:) + 2*B_rot(1,:).'*B_rot(2,:) + B_rot(1,:).'*B_rot(1,:) + B_rot(2,:).'*B_rot(2,:) );
                N13 = (detJ^2)*( -2*B_rot(1,:).'*B_rot(1,:) - 2*B_rot(2,:).'*B_rot(2,:) - 4*B_rot(1,:).'*B_rot(2,:) );
                N15 = (detJ^2)*( 2*B_rot(1,:).'*B_rot(3,:) + 2*B_rot(2,:).'*B_rot(3,:) );
                N25 = (detJ^2)*( B_rot(1,:).'*B_rot(3,:) + B_rot(2,:).'*B_rot(3,:) );
                N35 = (detJ^2)*( -2*B_rot(1,:).'*B_rot(3,:) - 2*B_rot(2,:).'*B_rot(3,:) );
                
                % Parameter estimates
                C11_app = paramInit(1);
                C33_app = paramInit(2);
                C44_app = paramInit(3);
                C66_app = paramInit(4);
                C13_app = paramInit(5);
                
                % Compile H matrix
                h = C11_app^2 * N11 + C33_app^2 * N22 + C44_app^2 * N33 + C66_app^2 * N44 + C13_app^2 * N55 + ...
                    2*C11_app*C44_app * N13 + 2*C11_app*C13_app * N15 + 2*C33_app*C13_app * N25 + 2*C44_app*C13_app * N35;
%                 h = conj(C11_app)*C11_app * N11 + conj(C33_app)*C33_app * N22 + conj(C44_app)*C44_app * N33 + conj(C66_app)*C66_app * N44 + conj(C13_app)*C13_app * N55 + ...
%                     conj(C11_app)*C44_app * N13 + C11_app*conj(C44_app) * N13 + conj(C11_app)*C13_app * N15 + C11_app*conj(C13_app) * N15 + conj(C33_app)*C13_app * N25 + ...
%                     C33_app*conj(C13_app) * N25 + conj(C44_app)*C13_app * N35 + C44_app*conj(C13_app) * N35;
                
                % Sum the weighted functions
                C1_1 = C1_1 + w(l).*f_c11;
                C2_1 = C2_1 + w(l).*f_c33;
                C3_1 = C3_1 + w(l).*f_c44;
                C4_1 = C4_1 + w(l).*f_c66;
                C5_1 = C5_1 + w(l).*f_c13;
                h_1 = h_1 + w(l).*h;

            end
            
            % Sum the weighted functions
            C1_2 = C1_2 + w(k).*C1_1;
            C2_2 = C2_2 + w(k).*C2_1;
            C3_2 = C3_2 + w(k).*C3_1;
            C4_2 = C4_2 + w(k).*C4_1;
            C5_2 = C5_2 + w(k).*C5_1;
            h_2 = h_2 + w(k).*h_1;

            % Clear the inner sums
            C1_1 = C1_1.*0;
            C2_1 = C2_1.*0;
            C3_1 = C3_1.*0;
            C4_1 = C4_1.*0;
            C5_1 = C5_1.*0;
            h_1 = h_1.*0;

        end
        
        % Sum the weighted functions
        C1_3 = C1_3 + w(j).*C1_2;
        C2_3 = C2_3 + w(j).*C2_2;
        C3_3 = C3_3 + w(j).*C3_2;
        C4_3 = C4_3 + w(j).*C4_2;
        C5_3 = C5_3 + w(j).*C5_2;
        h_3 = h_3 + w(j).*h_2;
   
        % Clear the inner sums
        C1_2 = C1_2.*0;
        C2_2 = C2_2.*0;
        C3_2 = C3_2.*0;
        C4_2 = C4_2.*0;
        C5_2 = C5_2.*0;
        h_2 = h_2.*0;

    end
    
    % Assemble constraint vectors
    C1(localNodeIdcs) = C1(localNodeIdcs) + C1_3;
    C2(localNodeIdcs) = C2(localNodeIdcs) + C2_3;
    C3(localNodeIdcs) = C3(localNodeIdcs) + C3_3;
    C4(localNodeIdcs) = C4(localNodeIdcs) + C4_3;
    C5(localNodeIdcs) = C5(localNodeIdcs) + C5_3;

    % Assemble H matrix
    c = 0; % counter
    for a = 1:length(h_3)
        for b = 1:length(h_3)
            c = c + 1;
            %H(localNodeIdcs(a), localNodeIdcs(b)) = H(localNodeIdcs(a), localNodeIdcs(b)) + h_3(a,b); %%%%%%% SLOW %%%%%%
            row_ind(i,c) = localNodeIdcs(a);
            col_ind(i,c) = localNodeIdcs(b);
            h_ind(i,c) = h_3(a,b);
        end
    end
    
end

% Put together sparse matrix for H
H = sparse(row_ind, col_ind, h_ind);

% Close wait bar
close(WH);
