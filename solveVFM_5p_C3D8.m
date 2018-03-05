function [A, C, strain] = solveVFM_5p_C3D8(u, uVF1, uVF2, uVF3, uVF4, uVF5, rho, omega, nodesSubZone, elemSubZone, globalNodeNums, orientation)

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
%          2) C (RHS)
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
A = zeros(5,5);
C = zeros(5,1);

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
    A_1 = zeros(5,5); C_1 = zeros(5,1);
    A_2 = zeros(5,5); C_2 = zeros(5,1);
    A_3 = zeros(5,5); C_3 = zeros(5,1);
    
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
    Br = B/detJ; % Reduced integration B matrix
    
    % Calculate strain: B*U - reduced integration
    eV_r = Br*Ue; % Strain of measured displacements using reduced integration
    eVF1_Vr = Br*UeVF1; % Strain of virtual displacement field #1
    eVF2_Vr = Br*UeVF2; % Strain of virtual displacement field #2
    eVF3_Vr = Br*UeVF3; % Strain of virtual displacement field #3
	eVF4_Vr = Br*UeVF4; % Strain of virtual displacement field #4
	eVF5_Vr = Br*UeVF5; % Strain of virtual displacement field #5
    
    % Convert strain to square tensor
    eR = [eV_r(1) 0.5*eV_r(4) 0.5*eV_r(5); ...
        0.5*eV_r(4) eV_r(2) 0.5*eV_r(6);...
        0.5*eV_r(5) 0.5*eV_r(6) eV_r(3)];
    
    eVF1r = [eVF1_Vr(1) 0.5*eVF1_Vr(4) 0.5*eVF1_Vr(5); ...
        0.5*eVF1_Vr(4) eVF1_Vr(2) 0.5*eVF1_Vr(6);...
        0.5*eVF1_Vr(5) 0.5*eVF1_Vr(6) eVF1_Vr(3)];
    
    eVF2r = [eVF2_Vr(1) 0.5*eVF2_Vr(4) 0.5*eVF2_Vr(5); ...
        0.5*eVF2_Vr(4) eVF2_Vr(2) 0.5*eVF2_Vr(6);...
        0.5*eVF2_Vr(5) 0.5*eVF2_Vr(6) eVF2_Vr(3)];
    
    eVF3r = [eVF3_Vr(1) 0.5*eVF3_Vr(4) 0.5*eVF3_Vr(5); ...
        0.5*eVF3_Vr(4) eVF3_Vr(2) 0.5*eVF3_Vr(6);...
        0.5*eVF3_Vr(5) 0.5*eVF3_Vr(6) eVF3_Vr(3)];
		
	eVF4r = [eVF4_Vr(1) 0.5*eVF4_Vr(4) 0.5*eVF4_Vr(5); ...
        0.5*eVF4_Vr(4) eVF4_Vr(2) 0.5*eVF4_Vr(6);...
        0.5*eVF4_Vr(5) 0.5*eVF4_Vr(6) eVF4_Vr(3)];
		
	eVF5r = [eVF5_Vr(1) 0.5*eVF5_Vr(4) 0.5*eVF5_Vr(5); ...
        0.5*eVF5_Vr(4) eVF5_Vr(2) 0.5*eVF5_Vr(6);...
        0.5*eVF5_Vr(5) 0.5*eVF5_Vr(6) eVF5_Vr(3)];
        
    % Rotate strain tensor
    eRotR = L * eR * L';
    eVF1_RotR = L * eVF1r * L';
    eVF2_RotR = L * eVF2r * L';
    eVF3_RotR = L * eVF3r * L';
	eVF4_RotR = L * eVF4r * L';
	eVF5_RotR = L * eVF5r * L';
    
    % Put strain tensor back into vector
    eVr = [eRotR(1,1); eRotR(2,2); eRotR(3,3); 2*eRotR(1,2); 2*eRotR(1,3); 2*eRotR(2,3)];
    eVF1_Vr = [eVF1_RotR(1,1); eVF1_RotR(2,2); eVF1_RotR(3,3); 2*eVF1_RotR(1,2); 2*eVF1_RotR(1,3); 2*eVF1_RotR(2,3)];
    eVF2_Vr = [eVF2_RotR(1,1); eVF2_RotR(2,2); eVF2_RotR(3,3); 2*eVF2_RotR(1,2); 2*eVF2_RotR(1,3); 2*eVF2_RotR(2,3)];
    eVF3_Vr = [eVF3_RotR(1,1); eVF3_RotR(2,2); eVF3_RotR(3,3); 2*eVF3_RotR(1,2); 2*eVF3_RotR(1,3); 2*eVF3_RotR(2,3)];
	eVF4_Vr = [eVF4_RotR(1,1); eVF4_RotR(2,2); eVF4_RotR(3,3); 2*eVF4_RotR(1,2); 2*eVF4_RotR(1,3); 2*eVF4_RotR(2,3)];
	eVF5_Vr = [eVF5_RotR(1,1); eVF5_RotR(2,2); eVF5_RotR(3,3); 2*eVF5_RotR(1,2); 2*eVF5_RotR(1,3); 2*eVF5_RotR(2,3)];
    
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
                
                %%%%%%%% LHS matrix %%%%%%%%
                
                %Multiply inverse of jacobian times the derivative of shape functions
                dNdXYZ = jac\DN;
                
                % Calculate B matrix
                B = [];
                for x = 1:nodesPerElem %Loop through number of nodes per element
                    Bi = [dNdXYZ(1,x) 0 0; 0 dNdXYZ(2,x) 0; 0 0 dNdXYZ(3,x); ...
                        dNdXYZ(2,x) dNdXYZ(1,x) 0; dNdXYZ(3,x) 0 dNdXYZ(1,x); 0 dNdXYZ(3,x) dNdXYZ(2,x)];
                    B = [B Bi];
                end
                
                % Full integration used for deviatoric part of strain
                Bf = B;
                
                % Calculate strain: B*U - full integrated strain
                eVf = Bf*Ue; % Strain of measured displacements
                eVF1_Vf = Bf*UeVF1; % Strain of virtual displacement field #1
                eVF2_Vf = Bf*UeVF2; % Strain of virtual displacement field #2
                eVF3_Vf = Bf*UeVF3; % Strain of virtual displacement field #3
				eVF4_Vf = Bf*UeVF4; % Strain of virtual displacement field #3
				eVF5_Vf = Bf*UeVF5; % Strain of virtual displacement field #3
                
                % Rotate strain matrix
                % Convert strains to square strain tensors
                eF = [eVf(1) 0.5*eVf(4) 0.5*eVf(5); ...
                    0.5*eVf(4) eVf(2) 0.5*eVf(6);...
                    0.5*eVf(5) 0.5*eVf(6) eVf(3)];
                
                eVF1f = [eVF1_Vf(1) 0.5*eVF1_Vf(4) 0.5*eVF1_Vf(5); ...
                    0.5*eVF1_Vf(4) eVF1_Vf(2) 0.5*eVF1_Vf(6);...
                    0.5*eVF1_Vf(5) 0.5*eVF1_Vf(6) eVF1_Vf(3)];
                
                eVF2f = [eVF2_Vf(1) 0.5*eVF2_Vf(4) 0.5*eVF2_Vf(5); ...
                    0.5*eVF2_Vf(4) eVF2_Vf(2) 0.5*eVF2_Vf(6);...
                    0.5*eVF2_Vf(5) 0.5*eVF2_Vf(6) eVF2_Vf(3)];
                
                eVF3f = [eVF3_Vf(1) 0.5*eVF3_Vf(4) 0.5*eVF3_Vf(5); ...
                    0.5*eVF3_Vf(4) eVF3_Vf(2) 0.5*eVF3_Vf(6);...
                    0.5*eVF3_Vf(5) 0.5*eVF3_Vf(6) eVF3_Vf(3)];
					
				eVF4f = [eVF4_Vf(1) 0.5*eVF4_Vf(4) 0.5*eVF4_Vf(5); ...
                    0.5*eVF4_Vf(4) eVF4_Vf(2) 0.5*eVF4_Vf(6);...
                    0.5*eVF4_Vf(5) 0.5*eVF4_Vf(6) eVF4_Vf(3)];
					
				eVF5f = [eVF5_Vf(1) 0.5*eVF5_Vf(4) 0.5*eVF5_Vf(5); ...
                    0.5*eVF5_Vf(4) eVF5_Vf(2) 0.5*eVF5_Vf(6);...
                    0.5*eVF5_Vf(5) 0.5*eVF5_Vf(6) eVF5_Vf(3)];
                
                % Rotate: rotMat * strain tensor * rotMat'
                eRotF = L * eF * L';
                eVF1_RotF = L * eVF1f * L';
                eVF2_RotF = L * eVF2f * L';
                eVF3_RotF = L * eVF3f * L';
				eVF4_RotF = L * eVF4f * L';
				eVF5_RotF = L * eVF5f * L';
                
                % Put tensor back into vector form
                eVf = [eRotF(1,1); eRotF(2,2); eRotF(3,3); 2*eRotF(1,2); 2*eRotF(1,3); 2*eRotF(2,3)];
                eVF1_Vf = [eVF1_RotF(1,1); eVF1_RotF(2,2); eVF1_RotF(3,3); 2*eVF1_RotF(1,2); 2*eVF1_RotF(1,3); 2*eVF1_RotF(2,3)];
                eVF2_Vf = [eVF2_RotF(1,1); eVF2_RotF(2,2); eVF2_RotF(3,3); 2*eVF2_RotF(1,2); 2*eVF2_RotF(1,3); 2*eVF2_RotF(2,3)];
                eVF3_Vf = [eVF3_RotF(1,1); eVF3_RotF(2,2); eVF3_RotF(3,3); 2*eVF3_RotF(1,2); 2*eVF3_RotF(1,3); 2*eVF3_RotF(2,3)];
				eVF4_Vf = [eVF4_RotF(1,1); eVF4_RotF(2,2); eVF4_RotF(3,3); 2*eVF4_RotF(1,2); 2*eVF4_RotF(1,3); 2*eVF4_RotF(2,3)];
				eVF5_Vf = [eVF5_RotF(1,1); eVF5_RotF(2,2); eVF5_RotF(3,3); 2*eVF5_RotF(1,2); 2*eVF5_RotF(1,3); 2*eVF5_RotF(2,3)];
             
                % Save strain from measured displacement field -
                % compare to Abaqus (CHECK)
                count = count + 1;
                idx = (i-1)*GP + count;
                strain(idx,:) = [i count eVr(1:3).' eVf(4:6).'];
                
                % Calculate each component of the A matrix
                % f_c11 - REDUCED INTEGRATION
                f_c11_vf1 = detJ*( eVr(1)*eVF1_Vr(1) + eVr(2)*eVF1_Vr(2) +  eVr(1)*eVF1_Vr(2) + eVr(2)*eVF1_Vr(1) );
                f_c11_vf2 = detJ*( eVr(1)*eVF2_Vr(1) + eVr(2)*eVF2_Vr(2) +  eVr(1)*eVF2_Vr(2) + eVr(2)*eVF2_Vr(1) );
                f_c11_vf3 = detJ*( eVr(1)*eVF3_Vr(1) + eVr(2)*eVF3_Vr(2) +  eVr(1)*eVF3_Vr(2) + eVr(2)*eVF3_Vr(1) );
                f_c11_vf4 = detJ*( eVr(1)*eVF4_Vr(1) + eVr(2)*eVF4_Vr(2) +  eVr(1)*eVF4_Vr(2) + eVr(2)*eVF4_Vr(1) );
                f_c11_vf5 = detJ*( eVr(1)*eVF5_Vr(1) + eVr(2)*eVF5_Vr(2) +  eVr(1)*eVF5_Vr(2) + eVr(2)*eVF5_Vr(1) );
                                
                % f_c33 - REDUCED INTEGRATION
                f_c33_vf1 = detJ*(eVr(3)*eVF1_Vr(3));
                f_c33_vf2 = detJ*(eVr(3)*eVF2_Vr(3));
                f_c33_vf3 = detJ*(eVr(3)*eVF3_Vr(3));
                f_c33_vf4 = detJ*(eVr(3)*eVF4_Vr(3));
                f_c33_vf5 = detJ*(eVr(3)*eVF5_Vr(3));              
                
                % f_c44 - FULL INTEGRATION
                f_c44_vf1 = detJ*(eVf(4)*eVF1_Vf(4) - 2*eVf(2)*eVF1_Vf(1) - 2*eVf(1)*eVF1_Vf(2)); 
                f_c44_vf2 = detJ*(eVf(4)*eVF2_Vf(4) - 2*eVf(2)*eVF2_Vf(1) - 2*eVf(1)*eVF2_Vf(2)); 
                f_c44_vf3 = detJ*(eVf(4)*eVF3_Vf(4) - 2*eVf(2)*eVF3_Vf(1) - 2*eVf(1)*eVF3_Vf(2)); 
                f_c44_vf4 = detJ*(eVf(4)*eVF4_Vf(4) - 2*eVf(2)*eVF4_Vf(1) - 2*eVf(1)*eVF4_Vf(2)); 
                f_c44_vf5 = detJ*(eVf(4)*eVF5_Vf(4) - 2*eVf(2)*eVF5_Vf(1) - 2*eVf(1)*eVF5_Vf(2)); 
                
                % f_c66 - FULL INTEGRATION
                f_c66_vf1 = detJ*(eVf(6)*eVF1_Vf(6) + eVf(5)*eVF1_Vf(5)); 
                f_c66_vf2 = detJ*(eVf(6)*eVF2_Vf(6) + eVf(5)*eVF2_Vf(5)); 
                f_c66_vf3 = detJ*(eVf(6)*eVF3_Vf(6) + eVf(5)*eVF3_Vf(5)); 
                f_c66_vf4 = detJ*(eVf(6)*eVF4_Vf(6) + eVf(5)*eVF4_Vf(5)); 
                f_c66_vf5 = detJ*(eVf(6)*eVF5_Vf(6) + eVf(5)*eVF5_Vf(5));                 
                
                % f_c13 - REDUCED INTEGRATION
                f_c13_vf1 = detJ*( eVr(1)*eVF1_Vr(3) + eVr(2)*eVF1_Vr(3) + eVr(3)*eVF1_Vr(1) + eVr(3)*eVF1_Vr(2));
                f_c13_vf2 = detJ*( eVr(1)*eVF2_Vr(3) + eVr(2)*eVF2_Vr(3) + eVr(3)*eVF2_Vr(1) + eVr(3)*eVF2_Vr(2));
                f_c13_vf3 = detJ*( eVr(1)*eVF3_Vr(3) + eVr(2)*eVF3_Vr(3) + eVr(3)*eVF3_Vr(1) + eVr(3)*eVF3_Vr(2));
                f_c13_vf4 = detJ*( eVr(1)*eVF4_Vr(3) + eVr(2)*eVF4_Vr(3) + eVr(3)*eVF4_Vr(1) + eVr(3)*eVF4_Vr(2));
                f_c13_vf5 = detJ*( eVr(1)*eVF5_Vr(3) + eVr(2)*eVF5_Vr(3) + eVr(3)*eVF5_Vr(1) + eVr(3)*eVF5_Vr(2));
                               
                % A matrix
                a = [f_c11_vf1 f_c33_vf1 f_c44_vf1 f_c66_vf1 f_c13_vf1; ...
                    f_c11_vf2 f_c33_vf2 f_c44_vf2 f_c66_vf2 f_c13_vf2; ...
                    f_c11_vf3 f_c33_vf3 f_c44_vf3 f_c66_vf3 f_c13_vf3; ...
                    f_c11_vf4 f_c33_vf4 f_c44_vf4 f_c66_vf4 f_c13_vf4; ...
                    f_c11_vf5 f_c33_vf5 f_c44_vf5 f_c66_vf5 f_c13_vf5];
                                              
                %%%%%%% RHS Vector %%%%%%%
                
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
                b4 = UeVF4Rot.'*rho*(omega^2)*sh*UeRot*detJ;
                b5 = UeVF5Rot.'*rho*(omega^2)*sh*UeRot*detJ;
                b = [b1; b2; b3; b4; b5];
                
                % Sum the weighted functions
                A_1 = A_1 + w(l).*a;
                C_1 = C_1 + w(l).*b;
                
            end
            
            % Sum the weighted functions
            A_2 = A_2 + w(k).*A_1;
            C_2 = C_2 + w(k).*C_1;

            % Clear the inner sums
            A_1 = A_1.*0;
            C_1 = C_1.*0;

        end
        
        % Sum the weighted functions
        A_3 = A_3 + w(j).*A_2;
        C_3 = C_3 + w(j).*C_2;

        % Clear the inner sums
        A_2 = A_2.*0;
        C_2 = C_2.*0;

    end
    
    % Final element LHS and RHS
    A = A + A_3;
    C = C + C_3;
    
end

% Get number of elements with negative jacobians
countNegJac = length(unique(elemNegJac));
disp(sprintf('There were %d elements with negative Jacobians.', countNegJac));

% Close wait bar
close(WH);
