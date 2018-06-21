function [A, C, K, strain] = solveVFM_3p_C3D8(u, uVF1, uVF2, uVF3, rho, omega, nodesSubZone, elemSubZone, globalNodeNums, orientation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the RHS (B) and LHS (A and K) for a transversely
% isotropic material described by three parameters: G12, G13 and T (G12*E3/E1)
%
% Inputs: 1) u - MRE displacement field
%         2) uVF1, uVF2 and uVF3 - virtual displacement fields
%		  3) rho - density
%         4) omega - angular frequency
%         5) DOF - number of degrees of freedom
%         6) nodesSubZone - list of nodal coordinates in subzone (numNodes x 3)
%         7) elemSubZone - (numElems x 8) - 8 node numbers that make up the
%         elements in the subzone
%         8) globalNodeNums - node numbers in the entire model (not just subzone)
%         9) orientation - vectors: a, b and f for each element (a x b = f)
%         10) GaussPoints - number of Gauss points per direction (1, 2 or 3)
%
% Outputs: 1) A (LHS)
%          2) B (RHS)
%          3) K (matrix of bulk stress values)
%          4) strain - to validate with Abaqus
%
% Written by: Renee Miller
% Updated: 30 June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss point coordinates
zeta = [-sqrt(1/3); sqrt(1/3)];
% Integration point weights
w = [1; 1];
% Number of Gauss points per element
GP = 8;

% Initialise components of matrices
A = zeros(3,3);
C = zeros(3,1);
K = zeros(3,1);

% Create waitbar
WH = waitbar(0, 'Assembling global matrices...');

% Initialise variable to create list of elements with negative jacobians
elemNegJac = [];

% Nodes per elem
nodesPerElem = size(elemSubZone,2) - 1;

% Nodal DOF
DOF = size(nodesSubZone,2) - 1;

% Variable to save strain at Gauss points (to comopare to Abaqus output)
strain = zeros(size(elemSubZone,1)*size(zeta,1)^3,8); % 7 = 1 element label + gauss point + 6 strain components (e11, e22, e33, e12, e23, e13)

% Loop through each element
for i = 1:size(elemSubZone,1)
    
    count = 0; % make counter to assign strain values
    
    % Update waitbar to give user an iclar allndication of time
    percentComplete = i/size(elemSubZone,1);
    waitbar(percentComplete, WH, sprintf('%.2f%% complete...', 100*percentComplete))
    
    % Initialise variables and matrices
    A_1 = zeros(3,3); K_1 = zeros(3,1); C_1 = zeros(3,1);
    A_2 = zeros(3,3); K_2 = zeros(3,1); C_2 = zeros(3,1);
    A_3 = zeros(3,3); K_3 = zeros(3,1); C_3 = zeros(3,1);
    
    % Get vectors of nodal coordinates
    [X, Y, Z] = getElemCoords(nodesSubZone, elemSubZone, i);
    
    %% Indexing - there are two types of node indexing - local and global. 
    % Local refers to the node index in the subzone and global refers to 
    % the node index in the whole model.
    
    % Get indics of nodal DOF in global list
    globalNodeIdcs = getNodeIdcs(globalNodeNums, elemSubZone, i, DOF);
    
    % Get indics of nodal DOF in global list
    localNodeIdcs = getNodeIdcs(nodesSubZone(:,1), elemSubZone, i, DOF);
    
    % Get displacements for particular element using indices
    Ue = u(globalNodeIdcs);
    UeVF1 = uVF1(localNodeIdcs);
    UeVF2 = uVF2(localNodeIdcs);
    UeVF3 = uVF3(localNodeIdcs);
    
    % Get rotation matrix
    L = getRotMat(orientation, elemSubZone(i));
      
    % Rotate displacement vectors
    % Initialise variables
    UeRot = zeros(size(Ue));
    UeVF1Rot = zeros(size(UeVF1));
    UeVF2Rot = zeros(size(UeVF2));
    UeVF3Rot = zeros(size(UeVF3));
    % Loop through vector of all displacements and rotate each individual nodal displacement vector
    for a = 1:DOF:length(Ue)
        UeRot(a:(a+2)) = L*Ue(a:(a+2));
        UeVF1Rot(a:(a+2)) = L*UeVF1(a:(a+2));
        UeVF2Rot(a:(a+2)) = L*UeVF2(a:(a+2));
        UeVF3Rot(a:(a+2)) = L*UeVF3(a:(a+2));
    end
    
    %% Compute reduced integration B matrix - Uniform strain method - Flanagan 1981
    
    % Make B matrix - Flanagan 1981 - Appendix I
    B = makeBmatrix(X, Y, Z);
    
    % Calculate element volume using uniform strain formulation - Flanagan 1981
    detJ = calcVolumeUniformStrain(X, Y, Z);
    
    % Divide B matrix by element volume to get B matrix for element
    Br = B/detJ; % Reduced integration B matrix
    
    % Compute the dilatational part of Br
    tmp = sum(Br(1:3,:),1);
    Br_dil = 1/3 * [tmp; tmp; tmp; zeros(1, length(tmp)); zeros(1, length(tmp)); zeros(1, length(tmp))];

    
    % Loop through Gauss points
    for j = 1:length(zeta)
        o = zeta(j);
        
        for k = 1:length(zeta)
            n = zeta(k);
            
            for l = 1:length(zeta)
                m = zeta(l);
                
                % Gauss integration coordinate = (m,n,o) %
                
                % Shape functions
                N1 = 1/8*(1-m)*(1-n)*(1-o);
                N2 = 1/8*(1+m)*(1-n)*(1-o);
                N3 = 1/8*(1+m)*(1+n)*(1-o);
                N4 = 1/8*(1-m)*(1+n)*(1-o);
                N5 = 1/8*(1-m)*(1-n)*(1+o);
                N6 = 1/8*(1+m)*(1-n)*(1+o);
                N7 = 1/8*(1+m)*(1+n)*(1+o);
                N8 = 1/8*(1-m)*(1+n)*(1+o);
                
                % Compile N - shape functions
                ShapeFuns = [N1 N2 N3 N4 N5 N6 N7 N8];
                
                % Evaluate the derivative of the shape functions at m, n
                % and o
                DN = 0.125*[-1*(1-n)*(1-o) (1-n)*(1-o) (1+n)*(1-o) -1*(1+n)*(1-o) -1*(1-n)*(1+o) (1-n)*(1+o) (1+n)*(1+o) -1*(1+n)*(1+o);...
                    -1*(1-m)*(1-o) -1*(1+m)*(1-o) (1+m)*(1-o) (1-m)*(1-o) -1*(1-m)*(1+o) -1*(1+m)*(1+o) (1+m)*(1+o) (1-m)*(1+o);
                    -1*(1-m)*(1-n) -1*(1+m)*(1-n) -1*(1+m)*(1+n) -1*(1-m)*(1+n) (1-m)*(1-n) (1+m)*(1-n) (1+m)*(1+n) (1-m)*(1+n)];
                
                % Calculate components of stiffness matrix
                jac = DN*[X Y Z];
                
                % Determinant of Jacobian  matrix
                detJ = det(jac);
                
                % Check for negative jacobians (distorted elements)
                if inv(jac) <= 0
                    %disp('Error: negative Jacobian in element');
                    elemNegJac = [elemNegJac i];
                end
                
                %Multiply inverse of jacobian times the derivative of shape functions
                dNdXYZ = jac\DN;
                
                % Calculate fully integrated B matrix
                Bf = [];
                for x = 1:nodesPerElem %Loop through number of nodes per element
                    Bi = [dNdXYZ(1,x)       0           0; ...
                            0           dNdXYZ(2,x)     0; ...
                            0               0       dNdXYZ(3,x); ...
                        dNdXYZ(2,x)     dNdXYZ(1,x)     0; ...
                        dNdXYZ(3,x)         0       dNdXYZ(1,x); ...
                        0               dNdXYZ(3,x) dNdXYZ(2,x)];
                    Bf = [Bf Bi];
                end
           
                
                %%%%%%%%% Selectively Reduced Stiffness Matrix %%%%%%%%%%%%
                % Full integration used for deviatoric part and reduced
                % integration used for dilational part (Hughes 1980)
                
                % Compute the dilatational part of Bf and Br
                tmp = sum(Bf(1:3,:),1);
                Bf_dil = 1/3 * [tmp; tmp; tmp; zeros(1, length(tmp)); zeros(1, length(tmp)); zeros(1, length(tmp))];
                
                % Compute the deviatoric part of Bf
                Bf_dev = Bf - Bf_dil;
                
                % Compute B_bar (the final B matrix to use)
                B_bar = Bf_dev + Br_dil;
                
                
                % Calculate strain: B*U
                eV = B_bar*Ue; % Strain of measured displacements
                eVF1_V = B_bar*UeVF1; % Strain of virtual displacement field
                eVF2_V = B_bar*UeVF2; % Strain of virtual displacement field
                eVF3_V = B_bar*UeVF3; % Strain of virtual displacement field
                
                % Convert strains to square strain tensors
                e = [eV(1) 0.5*eV(4) 0.5*eV(5); ...
                    0.5*eV(4) eV(2) 0.5*eV(6);...
                    0.5*eV(5) 0.5*eV(6) eV(3)];
                
                eVF1 = [eVF1_V(1) 0.5*eVF1_V(4) 0.5*eVF1_V(5); ...
                    0.5*eVF1_V(4) eVF1_V(2) 0.5*eVF1_V(6);...
                    0.5*eVF1_V(5) 0.5*eVF1_V(6) eVF1_V(3)];
                
                eVF2 = [eVF2_V(1) 0.5*eVF2_V(4) 0.5*eVF2_V(5); ...
                    0.5*eVF2_V(4) eVF2_V(2) 0.5*eVF2_V(6);...
                    0.5*eVF2_V(5) 0.5*eVF2_V(6) eVF2_V(3)];
                
                eVF3 = [eVF3_V(1) 0.5*eVF3_V(4) 0.5*eVF3_V(5); ...
                    0.5*eVF3_V(4) eVF3_V(2) 0.5*eVF3_V(6);...
                    0.5*eVF3_V(5) 0.5*eVF3_V(6) eVF3_V(3)];
                
                % Rotate: rotMat * strain tensor * rotMat'
                eRot = L * e * L';
                eVF1_Rot = L * eVF1 * L';
                eVF2_Rot = L * eVF2 * L';
                eVF3_Rot = L * eVF3 * L';
                
                % Put tensor back into vector form
                eV = [eRot(1,1); eRot(2,2); eRot(3,3); 2*eRot(1,2); 2*eRot(1,3); 2*eRot(2,3)];
                eVF1_V = [eVF1_Rot(1,1); eVF1_Rot(2,2); eVF1_Rot(3,3); 2*eVF1_Rot(1,2); 2*eVF1_Rot(1,3); 2*eVF1_Rot(2,3)];
                eVF2_V = [eVF2_Rot(1,1); eVF2_Rot(2,2); eVF2_Rot(3,3); 2*eVF2_Rot(1,2); 2*eVF2_Rot(1,3); 2*eVF2_Rot(2,3)];
                eVF3_V = [eVF3_Rot(1,1); eVF3_Rot(2,2); eVF3_Rot(3,3); 2*eVF3_Rot(1,2); 2*eVF3_Rot(1,3); 2*eVF3_Rot(2,3)];
             
                % Save strain from measured displacement field -
                % compare to Abaqus (CHECK)
                count = count + 1;
                idx = (i-1)*GP + count;
                strain(idx,:) = [i count eV.'];
                
                % Calculate each component of the K matrix
                c_k1 = detJ*trace(eRot)*trace(eVF1_Rot);
                c_k2 = detJ*trace(eRot)*trace(eVF2_Rot);
                c_k3 = detJ*trace(eRot)*trace(eVF3_Rot);
                c_k = [c_k1; c_k2; c_k3];
                
                % Calculate each component of the A matrix
  
                % mu12 
                c_mu12_1 = detJ*(8/9.*(eV(1)*eVF1_V(1) + eV(2)*eVF1_V(2)) - 10/9.*(eV(1)*eVF1_V(2) + eV(2)*eVF1_V(1)) - ...
                    4/9.*eV(3)*eVF1_V(3) + 2/9.*(eV(1)*eVF1_V(3) + eV(3)*eVF1_V(1) + eV(2)*eVF1_V(3) + eV(3)*eVF1_V(2)) + 2*eV(4)*0.5*eVF1_V(4)); 
                c_mu12_2 = detJ*(8/9.*(eV(1)*eVF2_V(1) + eV(2)*eVF2_V(2)) - 10/9.*(eV(1)*eVF2_V(2) + eV(2)*eVF2_V(1)) - ...
                    4/9.*eV(3)*eVF2_V(3) + 2/9.*(eV(1)*eVF2_V(3) + eV(3)*eVF2_V(1) + eV(2)*eVF2_V(3) + eV(3)*eVF2_V(2)) + 2*eV(4)*0.5*eVF2_V(4));
                c_mu12_3 = detJ*(8/9.*(eV(1)*eVF3_V(1) + eV(2)*eVF3_V(2)) - 10/9.*(eV(1)*eVF3_V(2) + eV(2)*eVF3_V(1)) - ...
                    4/9.*eV(3)*eVF3_V(3) + 2/9.*(eV(1)*eVF3_V(3) + eV(3)*eVF3_V(1) + eV(2)*eVF3_V(3) + eV(3)*eVF3_V(2)) + 2*eV(4)*0.5*eVF3_V(4));
                                
                % mu13
                c_mu13_1 = detJ*(2*eV(6)*0.5*eVF1_V(6) + 2*eV(5)*0.5*eVF1_V(5));
                c_mu13_2 = detJ*(2*eV(6)*0.5*eVF2_V(6) + 2*eV(5)*0.5*eVF2_V(5));
                c_mu13_3 = detJ*(2*eV(6)*0.5*eVF3_V(6) + 2*eV(5)*0.5*eVF3_V(5));
                                            
                % T
                c_T_1 = detJ*( 4/9.*(eV(1)*eVF1_V(1) + eV(2)*eVF1_V(2) + eV(1)*eVF1_V(2) + eV(2)*eVF1_V(1)) + ...
                    16/9.*eV(3)*eVF1_V(3) - 8/9.*(eV(1)*eVF1_V(3) + eV(3)*eVF1_V(1) + eV(2)*eVF1_V(3) + eV(3)*eVF1_V(2)));
                c_T_2 = detJ*( 4/9.*(eV(1)*eVF2_V(1) + eV(2)*eVF2_V(2) + eV(1)*eVF2_V(2) + eV(2)*eVF2_V(1)) + ...
                    16/9.*eV(3)*eVF2_V(3) - 8/9.*(eV(1)*eVF2_V(3) + eV(3)*eVF2_V(1) + eV(2)*eVF2_V(3) + eV(3)*eVF2_V(2)));
                c_T_3 = detJ*( 4/9.*(eV(1)*eVF3_V(1) + eV(2)*eVF3_V(2) + eV(1)*eVF3_V(2) + eV(2)*eVF3_V(1)) + ...
                    16/9.*eV(3)*eVF3_V(3) - 8/9.*(eV(1)*eVF3_V(3) + eV(3)*eVF3_V(1) + eV(2)*eVF3_V(3) + eV(3)*eVF3_V(2)));
                
                % A matrix
                a = [c_mu12_1 c_mu13_1 c_T_1; ...
                    c_mu12_2 c_mu13_2 c_T_2; ...
                    c_mu12_3 c_mu13_3 c_T_3];
                                              
                %%%%%%% RHS %%%%%%%
                
                N = []; % Initialise the list of shape functions 
                for sf = 1:length(ShapeFuns)
                    N = [N eye(DOF).*ShapeFuns(sf)]; % Append shape functions
                end
                
                % RHS = Integral(rho*omega^2*u*uVF*dV)
                sh = N'*N;
                sh = diag(sum(sh,1)); % Lumped mass matrix
                b1 = UeVF1Rot.'*rho*(omega^2)*sh*UeRot*detJ;
                b2 = UeVF2Rot.'*rho*(omega^2)*sh*UeRot*detJ;
                b3 = UeVF3Rot.'*rho*(omega^2)*sh*UeRot*detJ;
                b = [b1; b2; b3];
                
                % Sum the weighted functions
                A_1 = A_1 + w(l).*a;
                K_1 = K_1 + w(l).*c_k;
                C_1 = C_1 + w(l).*b;
                
            end
            
            % Sum the weighted functions
            A_2 = A_2 + w(k).*A_1;
            K_2 = K_2 + w(k).*K_1;
            C_2 = C_2 + w(k).*C_1;

            % Clear the inner sums
            A_1 = A_1.*0;
            K_1 = K_1.*0;
            C_1 = C_1.*0;

        end
        
        % Sum the weighted functions
        A_3 = A_3 + w(j).*A_2;
        K_3 = K_3 + w(j).*K_2;
        C_3 = C_3 + w(j).*C_2;

        % Clear the inner sums
        A_2 = A_2.*0;
        K_2 = K_2.*0;
        C_2 = C_2.*0;

    end
    
    % Final element LHS and RHS
    A = A + A_3;
    C = C + C_3;
    K = K + K_3;

end

% Get number of elements with negative jacobians
countNegJac = length(unique(elemNegJac));
fprintf('There were %d elements with negative Jacobians.', countNegJac);

% Close wait bar
close(WH);
