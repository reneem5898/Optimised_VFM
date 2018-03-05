function [C1, C2, C3, C4, H] = C3D8R_3p(U, nodesSubZone, elemSubZone, globalNodeNums, paramInit, orientation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function returns Ak, Ag and H matrices for constructing constraints on
% virtual displacement field as well as strain (for validating methods with
% Abaqus output)
%
% Inputs: 1) U - displacement vector
%         2) nodesSubZone - nodal coordinates (4 x # of nodes)
%         3) elemSubZone - element connectivity (9 x # of elements)
%         4) globalNodeNums - list of global node numbers
%         5) paramInit - initial parameter estimate (1 x 3)
%         6) orientation - element material orientations (three vectors -
% cross product of first two vectors gives fibre orientation in element)
%
% Outputs: 1) Constraints: C1, C2, C3 and C4
%          2) H - matrix ((# nodes x nodal DOFs) x (# nodes x nodal DOFs))
%
% Note: This reduced integration method utilises the uniform strain
% formulation from Flanagan 1981. This is the same method used by Abaqus.
%
% Renee Miller
% 28 June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodes per element
nodesPerElem = size(elemSubZone,2) - 1;

% Nodal DOF
DOF = size(nodesSubZone,2) - 1;

% Weighting for "single gauss point" - just used for constraints
w = 2;

% Initialise vectors for conditions on the virtual field
C1 = zeros(1,length(nodesSubZone)*DOF);
C2 = zeros(1,length(nodesSubZone)*DOF);
C3 = zeros(1,length(nodesSubZone)*DOF);
C4 = zeros(1,length(nodesSubZone)*DOF);

% Variable to save strain at Gauss points (to comopare to Abaqus output)
%strain = zeros(size(elemSubZone,1),7); % 7 = 1 element label + 6 strain components (e11, e22, e33, e12, e23, e13)

% Initialise row and col vectors
row_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);
col_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);
h_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);

% Create waitbar
WH = waitbar(0, 'Assembling global matrices...');

%count = 0; % make counter to assign strain values

% Loop through each element
for i = 1:size(elemSubZone,1)
    
    % Update waitbar to give user an indication of time
    percentComplete = i/size(elemSubZone,1);
    waitbar(percentComplete, WH, sprintf('%.2f%% complete...', 100*percentComplete))
    
    % Get rotation matrix
    L = getRotMat(orientation, elemSubZone(i));
        
    %% Indexing - there are two types of node indexing - local and global.
    % Local refers to the node index in the subzone and global refers to
    % the node index in the whole model.
    
    % Get indics of nodal DOF in global list
    globalNodeIdcs = getNodeIdcs(globalNodeNums, elemSubZone, i, DOF);
    
    % Get indics of nodal DOF in global list
    localNodeIdcs = getNodeIdcs(nodesSubZone(:,1), elemSubZone, i, DOF);
    
    % Update waitbar to give user an indication of time
    percentComplete = i/size(elemSubZone,1);
    waitbar(percentComplete, WH, sprintf('%.2f%% complete...', 100*percentComplete))
    
    % Get vectors of nodal coordinates
    [X, Y, Z] = getElemCoords(nodesSubZone, elemSubZone, i);
    
    % Get displacements for particular element using indices
    Ue = U(globalNodeIdcs);
    
    %% Uniform strain method - Flanagan 1981
    % This element type was implemented since it is the same as C3D8R element
    % type in Abaqus
    
    % Make B matrix - Flanagan 1981 - Appendix I
    B = makeBmatrix(X, Y, Z);
    
    % Calculate element volume using uniform strain formulation - Flanagan 1981
    detJ = calcVolumeUniformStrain(X, Y, Z);
    
    % Divide B matrix by element volume to get B matrix for element
    B = B/detJ; 
    
    % Calculate strain: B*U
    eV = B*Ue; % Strain of measured displacements using full integration
    
    % Convert strains to square strain tensor
    e = [eV(1) 0.5*eV(4) 0.5*eV(5); ...
        0.5*eV(4) eV(2) 0.5*eV(6);...
        0.5*eV(5) 0.5*eV(6) eV(3)];
    
    % Rotate strain tensor
    eRot = L * e * L';
    
    % Put strain tensor back into vector form
    eV = [eRot(1,1); eRot(2,2); eRot(3,3); 2*eRot(1,2); 2*eRot(1,3); 2*eRot(2,3)];
    
%     % Save strain from measured displacement field -
%     % compare to Abaqus (CHECK)
%     count = count + 1;
%     strain(count,:) = [i eV.'];
    
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
    
    % Calculate constraints: c1, c2, c3, and c4
    c_k = trace(e)*sum(B_rot(1:3,:),1)*detJ;
    c_G12 = detJ*(8/9.*(eV(1)*B_rot(1,:) + eV(2)*B_rot(2,:)) - 10/9.*(eV(1)*B_rot(2,:) + eV(2)*B_rot(1,:)) - ...
        4/9.*eV(3)*B_rot(3,:) + 2/9.*(eV(1)*B_rot(3,:) + eV(3)*B_rot(1,:) + eV(2)*B_rot(3,:) + eV(3)*B_rot(2,:)) + 2*eV(4)*0.5*B_rot(4,:));
    c_G13 = detJ*(2*eV(6)*0.5*B_rot(6,:) + 2*eV(5)*0.5*B_rot(5,:));
    c_T = detJ*( 4/9.*(eV(1)*B_rot(1,:) + eV(2)*B_rot(2,:) + eV(1)*B_rot(2,:) + eV(2)*B_rot(1,:)) + ...
        16/9.*eV(3)*B_rot(3,:) - 8/9.*(eV(1)*B_rot(3,:) + eV(3)*B_rot(1,:) + eV(2)*B_rot(3,:) + eV(3)*B_rot(2,:)));
    
    % Calculate the H matrix
    N11 = (detJ^2)*(56/27*(B_rot(1,:).'*B_rot(1,:) + B_rot(2,:).'*B_rot(2,:)) - 104/27*(B_rot(1,:).'*B_rot(2,:)) - 8/27*(B_rot(1,:).'*B_rot(3,:) + B_rot(2,:).'*B_rot(3,:)) ...
        + 8/27*(B_rot(3,:).'*B_rot(3,:)) + 4*0.5*B_rot(4,:).'*0.5*B_rot(4,:));
    N22 = (detJ^2)*(4*0.5*B_rot(5,:).'*0.5*B_rot(5,:) + 4*0.5*B_rot(6,:).'*0.5*B_rot(6,:));
    N33 = (detJ^2)*(128/27*(B_rot(1,:).'*B_rot(3,:) + B_rot(2,:).'*B_rot(3,:) + B_rot(3,:).'*B_rot(3,:)) + 32/27*(B_rot(1,:).'*B_rot(1,:) + B_rot(2,:).'*B_rot(2,:) + 2*B_rot(1,:).'*B_rot(2,:)));
    N13 = (detJ^2)*(-8/27*(B_rot(1,:).'*B_rot(1,:) + 2*B_rot(1,:).'*B_rot(2,:) + B_rot(2,:).'*B_rot(2,:)) + 64/81*(B_rot(1,:).'*B_rot(3,:) + B_rot(2,:).'*B_rot(3,:)) - 32/27*(B_rot(3,:).'*B_rot(3,:)));
       
    G12_app = paramInit(1);
    G13_app = paramInit(2);
    T_app = paramInit(3);
      
    h = G12_app^2 * N11  +  G12_app*T_app * 2 * N13  +  G13_app^2 * N22  +  T_app^2 * N33;
%     h = conj(G12_app)*G12_app * N11  +  G12_app*conj(T_app) * N13  + conj(G12_app)*T_app * N13 +  conj(G13_app)*G13_app * N22  +  conj(T_app)*T_app * N33;
    
    
    % Assemble constraint vectors
    C1(localNodeIdcs) = C1(localNodeIdcs) + c_k*w;
    C2(localNodeIdcs) = C2(localNodeIdcs) + c_G12*w;
    C3(localNodeIdcs) = C3(localNodeIdcs) + c_G13*w;
    C4(localNodeIdcs) = C4(localNodeIdcs) + c_T*w;
    
    % Assemble H matrix
    c = 0; % counter
    for a = 1:length(h)
        for b = 1:length(h)
            c = c + 1;
            %H(localNodeIdcs(a), localNodeIdcs(b)) = H(localNodeIdcs(a), localNodeIdcs(b)) + h_3(a,b); %%%%%%% SLOW %%%%%%
            row_ind(i,c) = localNodeIdcs(a);
            col_ind(i,c) = localNodeIdcs(b);
            h_ind(i,c) = h(a,b);
        end
    end
    
end

% Put together sparse matrix for H
H = sparse(row_ind, col_ind, h_ind);

% Close wait bar
close(WH);
