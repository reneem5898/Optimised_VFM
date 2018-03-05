function [A, C, strain] = solveVFM_5p_C3D8R(u, uVF1, uVF2, uVF3, uVF4, uVF5, rho, omega, nodesSubZone, elemSubZone, globalNodeNums, orientation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the RHS (B) and LHS (A and K) for a transversely
% isotropic material described by three parameters: G12, G13 and T (G12*E3/E1)
%
% Inputs: 1) u - MRE displacement field
%         2) uVF1, uVF2, uVF3, uVF4 and uVF5 - virtual displacement fields
%		  3) rho - density
%         4) omega - angular frequency
%         5) DOF - number of degrees of freedom
%         6) nodesSubZone - list of nodal coordinates in subzone (numNodes x 3)
%         7) elemSubZone - (numElems x 8) - 8 node numbers that make up the
%         elements in the subzone
%         8) globalNodeNums - node numbers in the entire model (not just subzone)
%         9) orientation - vectors: a, b and f for each element (a x b = f)
%
% Outputs: 1) A (LHS)
%          2) B (RHS)
%          3) K (matrix of bulk stress values)
%          4) strain - to validate with Abaqus
%
% Written by: Renee Miller
% Updated: 30 June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodal DOF
DOF = size(nodesSubZone,2) - 1;

% Weighting for "single gauss point" - just used for constraints
w = 2;

% Flanagan 1981 - Uniform Strain Paper
% Components of shape functions
Sigma = [1, 1, 1, 1, 1, 1, 1, 1];
Lambda1 = [-1, 1, 1, -1, -1, 1, 1, -1];
Lambda2 = [-1, -1, 1, 1, -1, -1, 1, 1];
Lambda3 = [-1, -1, -1, -1, 1, 1, 1, 1];
Gamma1 = [1, 1, -1, -1, -1, -1, 1, 1];
Gamma2 = [1, -1 -1, 1, -1, 1, 1, -1];
Gamma3 = [1, -1, 1, -1 1, -1, 1, -1];
Gamma4 = [-1, 1, -1, 1, 1, -1, 1, -1];

% Initialise components of matrices
A = zeros(5,5);
C = zeros(5,1);

% Create waitbar
WH = waitbar(0, 'Assembling global matrices...');

% Initialise variable to create list of elements with negative jacobians
elemNegJac = [];

% Variable to save strain at Gauss points (to comopare to Abaqus output)
strain = zeros(size(elemSubZone,1),7); % 7 = 1 element label + gauss point + 6 strain components (e11, e22, e33, e12, e23, e13)

count = 0; % make counter to assign strain values

% Loop through each element
for i = 1:size(elemSubZone,1)
    
    % Update waitbar to give user an iclar allndication of time
    percentComplete = i/size(elemSubZone,1);
    waitbar(percentComplete, WH, sprintf('%.2f%% complete...', 100*percentComplete))
    
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
	UeVF4 = uVF4(localNodeIdcs);
	UeVF5 = uVF5(localNodeIdcs);
    
    % Get rotation matrix
    L = getRotMat(orientation, elemSubZone(i));
    
    % Rotate displacement vectors
    % Initialise variables
    UeRot = zeros(size(Ue));
    UeVF1Rot = zeros(size(UeVF1));
    UeVF2Rot = zeros(size(UeVF2));
    UeVF3Rot = zeros(size(UeVF3));
	UeVF4Rot = zeros(size(UeVF4));
	UeVF5Rot = zeros(size(UeVF5));
	
    % Loop through vector of all displacements and rotate each individual nodal displacement vector
    for a = 1:DOF:length(Ue)
        UeRot(a:(a+2)) = L*Ue(a:(a+2));
        UeVF1Rot(a:(a+2)) = L*UeVF1(a:(a+2));
        UeVF2Rot(a:(a+2)) = L*UeVF2(a:(a+2));
        UeVF3Rot(a:(a+2)) = L*UeVF3(a:(a+2));
		UeVF4Rot(a:(a+2)) = L*UeVF4(a:(a+2));
        UeVF5Rot(a:(a+2)) = L*UeVF5(a:(a+2));
    end
    
    %% Compute reduced integration B matrix - Uniform strain method - Flanagan 1981
    
    % Make B matrix - Flanagan 1981 - Appendix I
    B = makeBmatrix(X, Y, Z);
    
    % Calculate element volume using uniform strain formulation - Flanagan 1981
    detJ = calcVolumeUniformStrain(X, Y, Z);
    
    % Divide B matrix by element volume to get B matrix for element
    B = B/detJ; % Reduced integration B matrix
    
    % Calculate strain: B*U - reduced integration
    eV = B*Ue; % Strain of measured displacements using reduced integration
    eVF1_V = B*UeVF1; % Strain of virtual displacement field #1
    eVF2_V = B*UeVF2; % Strain of virtual displacement field #2
    eVF3_V = B*UeVF3; % Strain of virtual displacement field #3
	eVF4_V = B*UeVF4; % Strain of virtual displacement field #4
	eVF5_V = B*UeVF5; % Strain of virtual displacement field #5
    
    % Convert strain to square tensor
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
		
	eVF4 = [eVF4_V(1) 0.5*eVF4_V(4) 0.5*eVF4_V(5); ...
        0.5*eVF4_V(4) eVF4_V(2) 0.5*eVF4_V(6);...
        0.5*eVF4_V(5) 0.5*eVF4_V(6) eVF4_V(3)];
                
    eVF5 = [eVF5_V(1) 0.5*eVF5_V(4) 0.5*eVF5_V(5); ...
        0.5*eVF5_V(4) eVF5_V(2) 0.5*eVF5_V(6);...
        0.5*eVF5_V(5) 0.5*eVF5_V(6) eVF5_V(3)];
    
    % Rotate strain tensor
    eRot = L * e * L';
    eVF1Rot = L * eVF1 * L';
    eVF2Rot = L * eVF2 * L';
    eVF3Rot = L * eVF3 * L';
	eVF4Rot = L * eVF4 * L';
	eVF5Rot = L * eVF5 * L';
    
    % Put strain tensor back into vector
    eV = [eRot(1,1); eRot(2,2); eRot(3,3); 2*eRot(1,2); 2*eRot(1,3); 2*eRot(2,3)];
    eVF1_V = [eVF1Rot(1,1); eVF1Rot(2,2); eVF1Rot(3,3); 2*eVF1Rot(1,2); 2*eVF1Rot(1,3); 2*eVF1Rot(2,3)];
    eVF2_V = [eVF2Rot(1,1); eVF2Rot(2,2); eVF2Rot(3,3); 2*eVF2Rot(1,2); 2*eVF2Rot(1,3); 2*eVF2Rot(2,3)];
    eVF3_V = [eVF3Rot(1,1); eVF3Rot(2,2); eVF3Rot(3,3); 2*eVF3Rot(1,2); 2*eVF3Rot(1,3); 2*eVF3Rot(2,3)];
	eVF4_V = [eVF4Rot(1,1); eVF4Rot(2,2); eVF4Rot(3,3); 2*eVF4Rot(1,2); 2*eVF4Rot(1,3); 2*eVF4Rot(2,3)];
	eVF5_V = [eVF5Rot(1,1); eVF5Rot(2,2); eVF5Rot(3,3); 2*eVF5Rot(1,2); 2*eVF5Rot(1,3); 2*eVF5Rot(2,3)];          
    
    
    % Save strain from measured displacement field -
    % compare to Abaqus (CHECK)
    count = count + 1;
    strain(count,:) = [i eV.'];
    
    % Calculate each component of the A matrix
    
    % Calculate each component of the A matrix
	% f_c11
	f_c11_vf1 = detJ*( eV(1)*eVF1_V(1) + eV(2)*eVF1_V(2) +  eV(1)*eVF1_V(2) + eV(2)*eVF1_V(1) );
	f_c11_vf2 = detJ*( eV(1)*eVF2_V(1) + eV(2)*eVF2_V(2) +  eV(1)*eVF2_V(2) + eV(2)*eVF2_V(1) );
	f_c11_vf3 = detJ*( eV(1)*eVF3_V(1) + eV(2)*eVF3_V(2) +  eV(1)*eVF3_V(2) + eV(2)*eVF3_V(1) );
	f_c11_vf4 = detJ*( eV(1)*eVF4_V(1) + eV(2)*eVF4_V(2) +  eV(1)*eVF4_V(2) + eV(2)*eVF4_V(1) );
	f_c11_vf5 = detJ*( eV(1)*eVF5_V(1) + eV(2)*eVF5_V(2) +  eV(1)*eVF5_V(2) + eV(2)*eVF5_V(1) );
					
	% f_c33
	f_c33_vf1 = detJ*(eV(3)*eVF1_V(3));
	f_c33_vf2 = detJ*(eV(3)*eVF2_V(3));
	f_c33_vf3 = detJ*(eV(3)*eVF3_V(3));
	f_c33_vf4 = detJ*(eV(3)*eVF4_V(3));
	f_c33_vf5 = detJ*(eV(3)*eVF5_V(3));              
	
	% f_c44
	f_c44_vf1 = detJ*(eV(4)*eVF1_V(4) - 2*eV(2)*eVF1_V(1) - 2*eV(1)*eVF1_V(2)); 
	f_c44_vf2 = detJ*(eV(4)*eVF2_V(4) - 2*eV(2)*eVF2_V(1) - 2*eV(1)*eVF2_V(2)); 
	f_c44_vf3 = detJ*(eV(4)*eVF3_V(4) - 2*eV(2)*eVF3_V(1) - 2*eV(1)*eVF3_V(2)); 
	f_c44_vf4 = detJ*(eV(4)*eVF4_V(4) - 2*eV(2)*eVF4_V(1) - 2*eV(1)*eVF4_V(2)); 
	f_c44_vf5 = detJ*(eV(4)*eVF5_V(4) - 2*eV(2)*eVF5_V(1) - 2*eV(1)*eVF5_V(2)); 
	
	% f_c66
	f_c66_vf1 = detJ*(eV(6)*eVF1_V(6) + eV(5)*eVF1_V(5)); 
	f_c66_vf2 = detJ*(eV(6)*eVF2_V(6) + eV(5)*eVF2_V(5)); 
	f_c66_vf3 = detJ*(eV(6)*eVF3_V(6) + eV(5)*eVF3_V(5)); 
	f_c66_vf4 = detJ*(eV(6)*eVF4_V(6) + eV(5)*eVF4_V(5)); 
	f_c66_vf5 = detJ*(eV(6)*eVF5_V(6) + eV(5)*eVF5_V(5));                 
	
	% f_c13
	f_c13_vf1 = detJ*( eV(1)*eVF1_V(3) + eV(2)*eVF1_V(3) + eV(3)*eVF1_V(1) + eV(3)*eVF1_V(2));
	f_c13_vf2 = detJ*( eV(1)*eVF2_V(3) + eV(2)*eVF2_V(3) + eV(3)*eVF2_V(1) + eV(3)*eVF2_V(2));
	f_c13_vf3 = detJ*( eV(1)*eVF3_V(3) + eV(2)*eVF3_V(3) + eV(3)*eVF3_V(1) + eV(3)*eVF3_V(2));
	f_c13_vf4 = detJ*( eV(1)*eVF4_V(3) + eV(2)*eVF4_V(3) + eV(3)*eVF4_V(1) + eV(3)*eVF4_V(2));
	f_c13_vf5 = detJ*( eV(1)*eVF5_V(3) + eV(2)*eVF5_V(3) + eV(3)*eVF5_V(1) + eV(3)*eVF5_V(2));
				   
	% A matrix
	a = [f_c11_vf1 f_c33_vf1 f_c44_vf1 f_c66_vf1 f_c13_vf1; ...
		f_c11_vf2 f_c33_vf2 f_c44_vf2 f_c66_vf2 f_c13_vf2; ...
		f_c11_vf3 f_c33_vf3 f_c44_vf3 f_c66_vf3 f_c13_vf3; ...
		f_c11_vf4 f_c33_vf4 f_c44_vf4 f_c66_vf4 f_c13_vf4; ...
		f_c11_vf5 f_c33_vf5 f_c44_vf5 f_c66_vf5 f_c13_vf5];
    
    %%%%%%% RHS %%%%%%%
    
    % Single integration point - only for determining the volume of the element
    m = 0; n = 0; o = 0;
    
    % Shape functions - Flanagan 1981
    ShapeFuns = 1/8 * Sigma + 1/4 * m * Lambda1 + 1/4 * n * Lambda2 + 1/4 * o * Lambda3 + ...
        1/2 * n * o * Gamma1 + 1/2 * m * o * Gamma2 + 1/2 * m * n * Gamma3 + ...
        1/2 * m * n * o * Gamma4;
    
    % Compile N
    %  _                                                       _
    % | N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8 0 0 |
    % | 0 N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8 0 |
    % | 0 0 N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8 |
    %  -                                                       -
    N = []; % Initialise the list of shape functions with a column of zeros
    for sf = 1:length(ShapeFuns) % Loop through each individual shape function
        N = [N eye(DOF).*ShapeFuns(sf)]; % Append identity shape functions (3 x 3)
    end    
    
    % RHS = Integral(rho*omega^2*u*uVF*dV)
    sh = N'*N;
    sh = diag(sum(sh,1)); % Lumped mass matrix
    b1 = UeVF1Rot.'*rho*(omega^2)*sh*UeRot*detJ;
    b2 = UeVF2Rot.'*rho*(omega^2)*sh*UeRot*detJ;
    b3 = UeVF3Rot.'*rho*(omega^2)*sh*UeRot*detJ;
	b4 = UeVF4Rot.'*rho*(omega^2)*sh*UeRot*detJ;
    b5 = UeVF5Rot.'*rho*(omega^2)*sh*UeRot*detJ;
    c = [b1; b2; b3; b4; b5];
    
    
    % Final element LHS and RHS
    A = A + w .* a;
    C = C + w .* c;
    
end

% Get number of elements with negative jacobians
countNegJac = length(unique(elemNegJac));
disp(sprintf('There were %d elements with negative Jacobians.', countNegJac));

% Close wait bar
close(WH);
