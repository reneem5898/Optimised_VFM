function [Ak, Ag, H, strain] = isoC2D4(U, nodesSubZone, elemSubZone, globalNodeNums, GP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function returns Ak, Ag and H matrices for constructing constraints on
% virtual displacement field as well as strain (for validating methods with
% Abaqus output)
%
% Inputs: 1) U - complex displacement vector 
%         2) nodesSubZone - nodal coordinates 
%         3) elemSubZone - element connectivity 
%         4) globalNodeNums - global node numbers
%         5) GaussPoints - number of Gauss points to use within an element
%
% Outputs: 1) Ak - row vector (1 x (# nodes x nodal DOFs))
%          2) Ag - row vectors (1 x (# nodes x nodal DOFs))
%          3) H - matrix ((# nodes x nodal DOFs) x (# nodes x nodal DOFs))
%          4) strain - matrix (# elements x 7)
%
% Renee Miller (renee.miller@auckland.ac.nz)
% 2 March 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodes per element
nodesPerElem = size(elemSubZone,2) - 1;

% Nodal DOF
DOF = size(nodesSubZone,2) - 1;

% Gauss point coordinates and weights
switch GP
    case 1
        zeta = 0;
        w = 2;
    case 2
        zeta = [-sqrt(1/3); sqrt(1/3)];
        w = [1; 1];
    case 3
        zeta = [-sqrt(3/5); 0; sqrt(3/5)];
        w = [5/9; 8/9; 5/9];
end

% Initialise vectors for Ag and Ak - conditions on the virtual field
Ak = zeros(1,length(nodesSubZone)*DOF);
Ag = zeros(1,length(nodesSubZone)*DOF);

% Variable to save strain at Gauss points (to comopare to Abaqus output)
strain = zeros(size(elemSubZone,1)*length(zeta)^3,5); % 4 = 1 element label + 3 strain components (e11, e22, e12)

% Initialise row and col vectors
row_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);
col_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);
h_ind = zeros(size(elemSubZone,1), (nodesPerElem*DOF)^2);

% Create waitbar
WH = waitbar(0, 'Assembling global matrices...');

count = 0; % make counter to assign strain values

% Loop through each element
for i = 1:size(elemSubZone,1)
    
    % Update waitbar to give user an indication of time
    percentComplete = i/size(elemSubZone,1);
    waitbar(percentComplete, WH, sprintf('%.2f%% complete...', 100*percentComplete))
    
    % Initialise variables and matrices
    A1_1 = 0; A2_1 = 0; h_1 = 0;
    A1_2 = 0; A2_2 = 0; h_2 = 0;
    
    % Get vectors of nodal coordinates
    [X, Y] = getElemCoords2D(nodesSubZone, elemSubZone, i);
    
    % Calculate delX, delY and delZ - length of sides of element
    [delX, delY] = calcQuadSides(X, Y);
    
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
        n = zeta(j);
        
        for k = 1:length(zeta)
            m = zeta(k);
            
            % Gauss integration coordinate = (m,n) %
            
            % Evaluate the derivative of the shape functions at m, n
            dNdm = 0.25 * [-1*(1-n) (1-n) (1+n) -1*(1+n)];
            dNdn = 0.25 * [-1*(1-m) -1*(1+m) (1+m) (1-m)];
                        
            DN = [dNdm; dNdn];
            
            % Calculate Jacobian for current element
            jac = DN*[X Y];
            
            % Multiply inverse of jacobian times the derivative of shape functions
            dNdXY = jac\DN;
            
            % Determinant of Jacobian
            detJ = det(jac);
            
            % Calculate B matrix (strain matrix)
            B = [];
            for c = 1:nodesPerElem %Loop through number of nodes per element
                Bi = [dNdXY(1,c)     0        ; ...
                      0            dNdXY(2,c) ; ...
                      dNdXY(2,c)  dNdXY(1,c)];
                B = [B Bi];
            end
            
            % Calculate strain: B*U
            eV = B*Ue; % Strain of measured displacements using full integration
            
            % Save strain from measured displacement field - compare to Abaqus (CHECK)
            count = count + 1;
            idx = (i-1)*GP + count;
            strain(idx,:) = [i count eV'];
            
            % Convert strains to square strain tensor
            e = [eV(1) 0.5*eV(3); ...
                0.5*eV(3) eV(2) ];
            
            % Ak term for current element (row vector)
            % fk = ak * UeVF = 0
            ak = trace(e)*sum(B(1:2,:),1)*detJ;
            
            % Ag term for current element (row vector)
            % fg = ag * UeVF = 1
            tmp = eV(1)*B(1,:) + eV(2)*B(2,:) + 2*0.5*eV(3)*0.5*B(3,:); % e:eVF
            % ag = 2 * (e:e* - 1/3 * Tr(e) * Tr(e*)) * dV
            ag = 2*(tmp - (1/3)* trace(e)* sum(B(1:2,:),1)) * detJ;
            
            %% H Matrix - Optimisation matrix for isotropic case (Connesson et al. 2015)
            
            % Construct H matrix
            h = Hmatrix2D(B, detJ, delX, delY);
            
            % Sum the weighted functions
            A1_1 = A1_1 + w(k).*ak;
            A2_1 = A2_1 + w(k).*ag;
            h_1 = h_1 + w(k).*h;
            
        end
        
        % Sum the weighted functions
        A1_2 = A1_2 + w(j).*A1_1;
        A2_2 = A2_2 + w(j).*A2_1;
        h_2 = h_2 + w(j).*h_1;
        
        % Clear the inner sums
        A1_1 = A1_1.*0;
        A2_1 = A2_1.*0;
        h_1 = h_1.*0;
        
    end
    
    
    % Assemble Ak and Ag vectors
    Ak(localNodeIdcs) = Ak(localNodeIdcs) + A1_2;
    Ag(localNodeIdcs) = Ag(localNodeIdcs) + A2_2;
    
    % Assemble H matrix
    c = 0;
    for a = 1:length(h_2)
        for b = 1:length(h_2)
            c = c + 1;
            %H(localNodeIdcs(a), localNodeIdcs(b)) = H(localNodeIdcs(a), localNodeIdcs(b)) + h_3(a,b); %%%%%%% SLOW %%%%%%
            row_ind(i,c) = localNodeIdcs(a);
            col_ind(i,c) = localNodeIdcs(b);
            h_ind(i,c) = h_2(a,b);
        end
    end
    
end

% Put together sparse matrix for H
H = sparse(row_ind, col_ind, h_ind);

% Close wait bar
close(WH);
