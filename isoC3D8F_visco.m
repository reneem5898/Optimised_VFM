function [AkR, AgR, AkI, AgI, H] = isoC3D8F_visco(Ur, Ui, nodesSubZone, elemSubZone, globalNodeNums)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function returns Ak, Ag and H matrices for constructing constraints on 
% virtual displacement field as well as strain (for validating methods with 
% Abaqus output)
%
% Inputs: 1) Ur - real part of displacement field
%         2) Ui - imaginary part of displacement field
%         3) nodesSubZone - nodal coordinates (4 x # of nodes)
%         4) elemSubZone - element connectivity (9 x # of elements)
%         5) globalNodeNums - node numbers in entire model
%
% Outputs: 1) Ak - row vector (1 x (# nodes x nodal DOFs))
%          2) Ag - row vectors (1 x (# nodes x nodal DOFs))
%          3) H - matrix ((# nodes x nodal DOFs) x (# nodes x nodal DOFs))
%          4) strain - matrix (# elements x 7)
%
% Note: This reduced integration method utilises the uniform strain
% formulation from Flanagan 1981. This is the same method used by Abaqus. 
% This function is used only with isotropic materials.
%
% Renee Miller
% 28 June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodes per element
nodesPerElem = size(elemSubZone,2) - 1;

% Nodal DOF
DOF = size(nodesSubZone,2) - 1;

% Gauss point coordinates
zeta = [-sqrt(1/3); sqrt(1/3)];
% Integration point weights
w = [1; 1];
% Number of Gauss points per element
GP = 8;

% Initialise vectors for Ag and Ak - conditions on the virtual field
AkR = zeros(1,length(nodesSubZone)*DOF);
AgR = zeros(1,length(nodesSubZone)*DOF);
AkI = zeros(1,length(nodesSubZone)*DOF);
AgI = zeros(1,length(nodesSubZone)*DOF);

% Variable to save strain at Gauss points (to comopare to Abaqus output)
strain = zeros(size(elemSubZone,1)*length(zeta)^3,8); % 7 = 1 element label + 6 strain components (e11, e22, e33, e12, e23, e13)

% Initialise row and col vectors
row_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);
col_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);
h_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);

% Create waitbar
WH = waitbar(0, 'Assembling global matrices...');

% Loop through each element
for i = 1:size(elemSubZone,1)
    
    count = 0; % make counter to assign strain values
    
    % Update waitbar to give user an indication of time
    percentComplete = i/size(elemSubZone,1);
    waitbar(percentComplete, WH, sprintf('%.2f%% complete...', 100*percentComplete))
    
    % Initialise variables and matrices
    A1_1r = 0; A2_1r = 0; A1_1i = 0; A2_1i = 0; h1_1 = 0;
    A1_2r = 0; A2_2r = 0; A1_2i = 0; A2_2i = 0; h1_2 = 0;
    A1_3r = 0; A2_3r = 0; A1_3i = 0; A2_3i = 0; h1_3 = 0;
    
    % Get vectors of nodal coordinates
    [X, Y, Z] = getElemCoords(nodesSubZone, elemSubZone, i);
    
    % Calculate delX, delY and delZ - length of sides of element
    [delX, delY, delZ] = calcHexSides(X, Y, Z);
    
    %% Indexing - there are two types of node indexing - local and global. 
    % Local refers to the node index in the subzone and global refers to 
    % the node index in the whole model.
    
    % Get indics of nodal DOF in global list
    globalNodeIdcs = getNodeIdcs(globalNodeNums, elemSubZone, i, DOF);
    
    % Get indics of nodal DOF in global list
    localNodeIdcs = getNodeIdcs(nodesSubZone(:,1), elemSubZone, i, DOF);
       
    % Get displacements for particular element using indices
    UeR = Ur(globalNodeIdcs);
    UeI = Ui(globalNodeIdcs);
    
    % Loop through integration points 
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
                eV_R = B*UeR; % Strain of measured real displacements 
                eV_I = B*UeI; % Strain of measured imaginary displacements 
                                
                % Convert strains to square strain tensor
                eR = [eV_R(1) 0.5*eV_R(4) 0.5*eV_R(5); ... %% Fully integrated - real component
                    0.5*eV_R(4) eV_R(2) 0.5*eV_R(6);...
                    0.5*eV_R(5) 0.5*eV_R(6) eV_R(3)];
                
                eI = [eV_I(1) 0.5*eV_I(4) 0.5*eV_I(5); ... %% Fully integrated - imaginary component
                    0.5*eV_I(4) eV_I(2) 0.5*eV_I(6);...
                    0.5*eV_I(5) 0.5*eV_I(6) eV_I(3)];
                
                % Ak term for current element (row vector)
                % fk = ak * UeVF = 0
                akR = trace(eR)*sum(B(1:3,:),1)*detJ;
                akI = trace(eI)*sum(B(1:3,:),1)*detJ;
                
                % Ag term for current element (row vector)
                % fg = ag * UeVF = 1
                tmpR = eV_R(1)*B(1,:) + eV_R(2)*B(2,:) + eV_R(3)*B(3,:) + 2*0.5*eV_R(4)*0.5*B(4,:) ...
                    + 2*0.5*eV_R(5)*0.5*B(5,:) + 2*0.5*eV_R(6)*0.5*B(6,:); % e:eVF
                tmpI = eV_I(1)*B(1,:) + eV_I(2)*B(2,:) + eV_I(3)*B(3,:) + 2*0.5*eV_I(4)*0.5*B(4,:) ...
                    + 2*0.5*eV_I(5)*0.5*B(5,:) + 2*0.5*eV_I(6)*0.5*B(6,:); % e:eVF
                % ag = 2 * (e:e* - 1/3 * Tr(e) * Tr(e*)) * dV
                agR = 2*(tmpR - (1/3)* trace(eR)* sum(B(1:3,:),1)) * detJ;
                agI = 2*(tmpI - (1/3)* trace(eI)* sum(B(1:3,:),1)) * detJ;
                                
                %% H Matrix - Optimisation matrix for isotropic case (Connesson et al. 2015)
    
                % Construct H matrix
                h = Hmatrix(B, detJ, delX, delY, delZ);  
                h = h + h*1i; % ADDED 30 JUNE 2017 - RENEE MILLER
                
                % Sum the weighted functions
                A1_1r = A1_1r + w(l).*akR;
                A1_1i = A1_1i + w(l).*akI;
                A2_1r = A2_1r + w(l).*agR;
                A2_1i = A2_1i + w(l).*agI;
                h1_1 = h1_1 + w(l).*h;
                
            end
            
            % Sum the weighted functions
            A1_2r = A1_2r + w(k).*A1_1r;
            A1_2i = A1_2i + w(k).*A1_1i;
            A2_2r = A2_2r + w(k).*A2_1r;
            A2_2i = A2_2i + w(k).*A2_1i;
            h1_2 = h1_2 + w(k).*h1_1;
            
            % Clear the inner sums
            A1_1r = A1_1r.*0;
            A1_1i = A1_1i.*0;
            A2_1r = A2_1r.*0;
            A2_1i = A2_1i.*0;
            h1_1 = h1_1.*0;
            
        end
        
        % Sum the weighted functions
        A1_3r = A1_3r + w(j).*A1_2r;
        A1_3i = A1_3i + w(j).*A1_2i;
        A2_3r = A2_3r + w(j).*A2_2r;
        A2_3i = A2_3i + w(j).*A2_2i;
        h1_3 = h1_3 + w(j).*h1_2;
        
        % Clear the inner sums
        A1_2r = A1_2r.*0;
        A1_2i = A1_2i.*0;
        A2_2r = A2_2r.*0;
        A2_2i = A2_2i.*0;
        h1_2 = h1_2.*0;
        
    end
    
    % Assemble Ak and Ag vectors
    AkR(nodeIdcs) = AkR(nodeIdcs) + A1_3r;
    AkI(nodeIdcs) = AkI(nodeIdcs) + A1_3i;
    AgR(nodeIdcs) = AgR(nodeIdcs) + A2_3r;
    AgI(nodeIdcs) = AgI(nodeIdcs) + A2_3i;
    
    % Assemble H matrix
    c = 0;
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
