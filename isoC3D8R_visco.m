function [AkR, AgR, AkI, AgI, H] = isoC3D8R_visco(Ur, Ui, nodesSubZone, elemSubZone, globalNodeNums)

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

% Weighting for "single gauss point" - just used for constraints
w = 2;

% Initialise vectors for Ag and Ak - conditions on the virtual field
AkR = zeros(1,length(nodesSubZone)*DOF);
AgR = zeros(1,length(nodesSubZone)*DOF);
AkI = zeros(1,length(nodesSubZone)*DOF);
AgI = zeros(1,length(nodesSubZone)*DOF);

% Initialise row and col vectors
row_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);
col_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);
h_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);

% Create waitbar
WH = waitbar(0, 'Assembling global matrices...');

% Loop through each element
for i = 1:size(elemSubZone,1)
    
    % Update waitbar to give user an indication of time
    percentComplete = i/size(elemSubZone,1);
    waitbar(percentComplete, WH, sprintf('%.2f%% complete...', 100*percentComplete))
    
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
    eVr = B*UeR; % Strain of measured real displacements
    eVi = B*UeI; % Strain of measured imaginary displacements
    
    % Convert strains to square strain tensor
    eR = [eVr(1) 0.5*eVr(4) 0.5*eVr(5); ...
        0.5*eVr(4) eVr(2) 0.5*eVr(6);...
        0.5*eVr(5) 0.5*eVr(6) eVr(3)];
    % Convert strains to square strain tensor
    eI = [eVi(1) 0.5*eVi(4) 0.5*eVi(5); ...
        0.5*eVi(4) eVi(2) 0.5*eVi(6);...
        0.5*eVi(5) 0.5*eVi(6) eVi(3)];
    
    %% FK and FG constraints on the virtual field
    
    % Ak term for current element (row vector)
    % fk = ak * UeVF = 0
    akR = trace(eR)*sum(B(1:3,:),1)*detJ;
    akI = trace(eI)*sum(B(1:3,:),1)*detJ;
    
    % Ag term for current element (row vector)
    % fg = ag * UeVF = 1
    tmpR = eVr(1)*B(1,:) + eVr(2)*B(2,:) + eVr(3)*B(3,:) + 2*0.5*eVr(4)*0.5*B(4,:) + 2*0.5*eVr(5)*0.5*B(5,:) + 2*0.5*eVr(6)*0.5*B(6,:); % e:eVF
    tmpI = eVi(1)*B(1,:) + eVi(2)*B(2,:) + eVi(3)*B(3,:) + 2*0.5*eVi(4)*0.5*B(4,:) + 2*0.5*eVi(5)*0.5*B(5,:) + 2*0.5*eVi(6)*0.5*B(6,:); % e:eVF
    agR = 2*(tmpR - (1/3)*trace(eR)*sum(B(1:3,:),1))*detJ;
    agI = 2*(tmpI - (1/3)*trace(eI)*sum(B(1:3,:),1))*detJ;
    
    % Assemble Ak and Ag vectors
    AkR(localNodeIdcs) = AkR(localNodeIdcs) + akR*w;
    AgR(localNodeIdcs) = AgR(localNodeIdcs) + agR*w;
    AkI(localNodeIdcs) = AkI(localNodeIdcs) + akI*w;
    AgI(localNodeIdcs) = AgI(localNodeIdcs) + agI*w;
    
    %% H Matrix - Optimisation matrix for isotropic case (Connesson et al. 2015)
    
    % Construct H matrix
    h = Hmatrix(B, detJ, delX, delY, delZ);
    h = h + h*1i;
    
    % Assemble H matrix
    c = 0;
    for a = 1:length(h)
        for b = 1:length(h)
            c = c + 1;
            %H(localNodeIdcs(a), localNodeIdcs(b)) = H(localNodeIdcs(a), localNodeIdcs(b)) + h_3(a,b); %%%%%%% SLOW %%%%%%
            row_ind(i,c) = localNodeIdcs(a);
            col_ind(i,c) = localNodeIdcs(b);
            h_ind(i,c) = h(a,b)*w;
        end
    end
    
end

% Put together sparse matrix for H
H = sparse(row_ind, col_ind, h_ind);

% Close wait bar
close(WH);
