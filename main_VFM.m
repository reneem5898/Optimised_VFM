function [Results] = main_VFM(modelDir, outDir, rho, freq, refParams, maxIter, subzones, elemType)

%% Optimised Virtual Fields Method %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function utilises the optimised virutal fields method to solve for
% the global linear elastic stiffness properties of a 3D structure from
% a complex harmonic displacement field given:
% 1) node coordinates, 2) element connectivitiy and 3) nodal displacements
%
% Allows for estimating 1, 3 or 5 parameters for an isotropic or transversely
% isotropic material. Number of parameters to estimate is obtained by the
% length of "refParams". If refParams is empty, the program will estimate
% 1 parameter (assume isotropy).
%
% Inputs: 1) modelDir - model directory where to find node, element, etc. data
%         2) outDir - output directory
%         3) rho - density
%         4) freq - frequency (Hz)
%         5) refParams - reference or initial parameters to use in estimation
%         6) maxIter - maximum number of iterations allowable (default = 30)
%         7) subzones = struct with three fields: x, y and z where subzones.x is a N x 2 matrix (N = number of subzones, 2 = range of x coords in each suzbone)
%         8) elemType - string controlling integration type: 'C3D8' (selectively reduced integration elements),
% 'C3D8R' (fully reduced integration) and 'C3D8F' (fully integrated - giving the user the choice of the number of Gauss points)
%
% Output: Results - struct containing results from VFM analysis
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Input file formats and locations:
%
% <modelDir>/nodeCoords.txt
% 2,	25.2733307,	9.80368423,	-0.286143214
% 3,	26.25741,	-0.0725992396,	-12.8894196
% 4,	26.3661003,	-11.5723,	-0.198858097
% 5,	26.1123409,	-0.00780306477,	9.79632473
% 6,	14.2937803,	12.3374205,	-0.0866783634
% 7,	15.3667002,	0.146519497,	-17.0545006
% ...
%
% <modelDir>/elems.txt
%    1,    2,  131, 2584, 2582,  389,  390, 2585, 2583
%    2,    3,  133, 2588, 2586,  391,  392, 2589, 2587
%    3,    4,  135, 2592, 2590,  393,  394, 2593, 2591
%    4,    5,  137, 2596, 2594,  395,  396, 2597, 2595
%    5,    6,  139, 2600, 2598,  397,  398, 2601, 2599
% ...
%
% <outDir>/uComplex.txt
% -0.46153978	0.53962666
% 0.39781991	0.81265861
% -2.3007324	-0.80385208
% 0.40850705	-1.3562853
% -0.0043678368	-0.89747274
% -1.6591556	-0.24158493
% ...
%
% <modelDir>/Dir1.txt   (vectors defining the FIRST material orientation at the centroid of each element)
% 0.6245934367,	 0.0000000000,	 -0.7809501886
% 0.5084568262,	 -0.0000000298,	 -0.8610874414
% 0.7487915158,	 0.0000000000,	 0.6628056169
% 0.0843197405,	 -0.0000000298,	 0.9964387417
% 0.5443310142,	 0.0000000000,	 -0.8388704658
% ...
%
% <modelDir>/Dir2.txt   (vectors defining the SECOND material orientation at the centroid of each element)
% 0.2296404690,	 0.9557892680,	 0.1836633533
% -0.5473996401,	 0.7719300389,	 -0.3232297897
% -0.1096123829,	 0.9862305522,	 0.1238324270
% 0.8808015585,	 0.4675824940,	 -0.0745343864
% 0.1499492675,	 0.9838942289,	 0.0972999260
% ...
%
% <modelDir>/Dir3.txt   (vectors defining the THIRD material orientation - FIBRE ORIENTATION - at the centroid of each element)
% 0.7464237213,	 -0.2940526605,	 0.5969796777
% 0.6646993160,	 0.6357074380,	 0.3924930990
% -0.6536791921,	 -0.1653763652,	 0.7384810448
% -0.4659173191,	 0.8839495182,	 0.0394264758
% 0.8253598809,	 -0.1787513793,	 0.5355641842
% ...
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%
%
% Written: Renee Miller (reneem5898@gmail.com)
% Updated: 6 March 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants

% Material properties
omega = freq*2*pi; % angular frequency in rad/s

% Number of parameters to solve for: 1, 3, or 5
if ~isempty(refParams)
    if length(refParams) ~= 3 && length(refParams) ~= 5
        disp('You have not entered an appropriate number of reference parameters. To estimate transversely isotropic material properties, either 3 or 5 reference parameters are needed. The material will be assumed isotropic.');
        numParam = 1;
    else
        numParam = length(refParams);
    end
else
    numParam = 1;
end

% Number of Gauss points in each direction - used for fully integrated elements only
GaussPoints = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create diary file
tmp = clock;
diaryFile = sprintf('%s/logFile_optimisedVF_%dparam_%d%d%d_%d%d.txt', outDir, numParam, tmp(1), tmp(2), tmp(3), tmp(4), tmp(5));
diary(diaryFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load model data

disp('Loading model node and element data...');

% Node coordinates
n = load(strcat(modelDir,'/nodeCoords.txt'));

% Each element = 1 row with eight node numbers
e = load(strcat(modelDir,'/elems.txt'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove any nodes from node list which are not in elements

% Get list of nodes in elements
elemNodes = e(:,2:end);
elemNodes = unique(elemNodes(:));

for i = 1:length(n)
    if ~any(elemNodes == n(i,1))
        n(i,:) = [0, 0, 0, 0];
    end
end
n( ~any(n,2), : ) = [];  %Remove zero rows

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Renumber nodes from 1

[nnodes, elems] = renumberNodes(n(:,1), e);
nodes = [nnodes' n(:,2:end)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load boundary nodes (if available)

% Load boundary condition nodes - if available
bcNodesFile = strcat(modelDir,'/surfNodes.txt');
if exist(bcNodesFile, 'file')
    bcNodes = load(bcNodesFile);
else
    bcNodes = [];
end

% Get bc node numbers in new renumbered nodes
bcNodesNew = bcNodes;
for i = 1:length(bcNodes)
    idx = find(n(:,1) == bcNodes(i));
    bcNodesNew(i) = idx;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Element Material Orientation

% If estimating anisotropic material properties, we need to know the material orientation
if numParam > 1
    
    % Load material orientation files from model directory
    if exist(strcat(modelDir,'/Dir1.txt'), 'f')
        
        dir1 = load(strcat(modelDir,'/Dir1.txt')); % num elements x 3: <u1, u2, u3>
        dir2 = load(strcat(modelDir,'/Dir2.txt')); % num elements x 3: <v1, v2, v3>
        dir3 = load(strcat(modelDir,'/Dir3.txt')); % num elements x 3: <f1, f2, f3> %%%% fibre direction %%%%
        
        
    % Load material orientation files from output directory - when varying material orientations was tested
    elseif exist(strcat(outDir,'/Dir1.txt'), 'f')
        
        dir1 = load(strcat(outDir,'/Dir1.txt')); % num elements x 3: <u1, u2, u3>
        dir2 = load(strcat(outDir,'/Dir2.txt')); % num elements x 3: <v1, v2, v3>
        dir3 = load(strcat(outDir,'/Dir3.txt')); % num elements x 3: <f1, f2, f3> %%%% fibre direction %%%%
                
    % If these don't exist, set element material orientation
    else
        
        disp('Could not find files with element material orientation. A default orientation: f = <0, 0, 1> will be applied.');
        
        numElems = size(elems,1);
        dir1 = [ones(numElems,1) zeros(numElems,1) zeros(numElems,1)];
        dir2 = [zeros(numElems,1) ones(numElems,1) zeros(numElems,1)];
        dir3 = [zeros(numElems,1) zeros(numElems,1) ones(numElems,1)];
        
    end
    
    % Combine orientations into one variable
    orientation = [elems(:,1) dir1 dir2 dir3];
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subzones

if isstruct(subzones)
    
    % Load x, y and z ranges for subzones from struct: subzones
    xRange = subzones.x;
    yRange = subzones.y;
    zRange = subzones.z;
else
    
    % Node range == entire model
    xRange = [min(nodes(:,2)) max(nodes(:,2))];
    yRange = [min(nodes(:,3)) max(nodes(:,3))];
    zRange = [min(nodes(:,4)) max(nodes(:,4))];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load displacements
dispFile = strcat(outDir, '/uComplex.txt');
U = load(dispFile);
[~, tmp] = size(U); % Get number of columns
if tmp == 2 % If data is complex (e.g. two columns)
    U = U(:,1) + U(:,2)*1i;
else % Else if data is real only
    U = U(:,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through subzones

% Count subzone
countZone = 0;

for m = 1:size(xRange,1)
    for n = 1:size(yRange,1)
        for o = 1:size(zRange,1)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Assemble subzone node and element lists
            disp('Assembling element subzone node and element lists...');
            
            % Counter for zone number
            countZone = countZone + 1;
            
            % Creating subzone node and element lists
            [nodesSubZone, elemSubZone] = getSubZone(nodes, elems, xRange, yRange, zRange, m, n, o);
            surfNodes = intersect(nodesSubZone(:,1), bcNodesNew);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% ISOTROPIC MATERIAL MODEL
            if numParam == 1
                
                disp(sprintf('Running Region #: %d\n', countZone));
                
                % Calculate the numeric virtual field
                disp('Calculating the numeric virtual field...')
                [uVF, eta, strain1] = numericVF_Iso(U, nodesSubZone, elemSubZone, nodes(:,1), surfNodes, elemType, GaussPoints);
                
                % Calculate shear modulus
                disp('Calculating the shear modulus...')
                if strcmp(elemType, 'C3D8R')
                    % Uniform strain elements - C3D8R
                    [fk, fg, b, strain2] = solveIsoVFM_C3D8R(U, uVF, rho, omega, nodesSubZone, elemSubZone, nodes(:,1));
                    
                elseif strcmp(elemType, 'C3D8')
                    % Selectively reduced integration type element - C3D8
                    [fk, fg, b, strain2] = solveIsoVFM_C3D8(U, uVF, rho, omega, nodesSubZone, elemSubZone, nodes(:,1));
                    
                else
                    % Fully integrated element
                    [fk, fg, b, strain2] = solveIsoVFM_C3D8F(U, uVF, rho, omega, nodesSubZone, elemSubZone, nodes(:,1), GaussPoints);
                    
                end
                
                % Calculate complex shear modulus
                G = b/fg
                
                % Save results as struct
                result = struct('G', G, 'FK', fk, 'FG', fg, 'B', b, 'uVF', uVF, 'eta', eta, 'strain1', strain1, 'strain2', strain2);
                Results(countZone) = result;
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% TRANSVERSELY ISOTROPIC MATERIAL MODEL - 3 PARAMETERS
            elseif numParam == 3
                
                % Initialise parameters
                diffPerc = 1;
                iter = 0;
                
                % Iteratively solve for parameters until change in parameters is less than 0.1% - Avril 2004
                while diffPerc > 0.1 && iter < maxIter
                    
                    % Count number of iterations
                    iter = iter + 1;
                    
                    disp(sprintf('Running Region #: %d\nIteration #: %d\n', countZone, iter));
                    
                    % Qapp - approximate parameters to use to calculate virtual displacement field
                    if ~exist('moduli','var')
                        paramEst = refParams;
                    else
                        paramEst = moduli; % Parameters calculated from previous subzone
                    end
                    
                    % Create numeric virtual displacement field
                    disp('Calculating the numeric virtual displacement fields...')
                    [uVF1, uVF2, uVF3, eta] = numericVF_3p(U, nodesSubZone, elemSubZone, nodes(:,1), paramEst, orientation, surfNodes, elemType, GaussPoints);
                    
                    % Calculate material properties
                    disp('Calculating the estimated parameters...')
                    if strcmp(elemType, 'C3D8R')
                        % Uniform strain elements - C3D8R
                        [A, B, K, strain] = solveVFM_3p_C3D8R(U, uVF1, uVF2, uVF3, rho, omega, nodesSubZone, elemSubZone, nodes(:,1), orientation);
                        
                    elseif strcmp(elemType, 'C3D8')
                        % Selectively reduced integration type element - C3D8
                        [A, B, K, strain] = solveVFM_3p_C3D8(U, uVF1, uVF2, uVF3, rho, omega, nodesSubZone, elemSubZone, nodes(:,1), orientation);
                        
                    else
                        % Fully integrated element
                        [A, B, K, strain] = solveVFM_3p_C3D8F(U, uVF1, uVF2, uVF3, rho, omega, nodesSubZone, elemSubZone, nodes(:,1), orientation, GaussPoints);
                        
                    end
                    
                    moduli = A\B
                    
                    % Calculate the percent change in estimated material parameters
                    diffPerc = max(abs((moduli - paramEst)./paramEst))*100;
                    disp(sprintf('Maximum change in parameters: %.2f%', diffPerc));
                    
                end
                
                % Get parameters from complex values of elasticity matrix (moduli)
                [G12, G13, damp, T] = getParams_3p(moduli, iter, maxIter);
                
                % Save results in struct
                result = struct('G12', G12, 'G13', G13, 'T', T, 'damp', damp, 'uVF1', uVF1, 'uVF2', uVF2, 'uVF3', uVF3, 'iter', iter, 'eta', eta, 'strain', strain, 'K', K, 'orientation', orientation);
                Results(countZone) = result;
                
                clear moduli
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% TRANSVERSELY ISOTROPIC MATERIAL MODEL - 5 PARAMETERS
            elseif numParam == 5
                
                % Initialise parameters
                diffPerc = 1;
                iter = 0;
                
                % Iteratively solve for parameters until change in parameters is less than 0.1% - Avril 2004
                while diffPerc > 0.1 && iter < maxIter
                    
                    % Count number of iterations
                    iter = iter + 1;
                    
                    disp(sprintf('Running Region #: %d\nIteration #: %d\n', countZone, iter));
                    
                    % Qapp - approximate parameters to use to calculate virtual displacement field
                    if ~exist('moduli','var')
                        paramEst = refParams; % Initial parameter guesses
                    else
                        paramEst = moduli; % Parameters calculated from previous subzone
                    end
                    
                    % Calculate engineering constants to use for calculating convergence
                    [Ec, G12c, G13c, E1c, E3c, v12c, v13c, v31c, dampc] = getParams_5p(paramEst, 1, maxIter);
                    paramComp = [G12c; G13c; E1c; E3c];
                    
                    % Create numeric virtual displacement field
                    disp('Calculating the numeric virtual displacement fields...')
                    [uVF1, uVF2, uVF3, uVF4, uVF5, eta] = numericVF_5p(U, nodesSubZone, elemSubZone, nodes(:,1), paramEst, orientation, surfNodes, elemType, GaussPoints);
                    
                    % Calculate material properties
                    disp('Calculating the estimated parameters...')
                    if strcmp(elemType, 'C3D8R')
                        % Uniform strain elements - C3D8R
                        [A, B, strain] = solveVFM_5p_C3D8R(U, uVF1, uVF2, uVF3, uVF4, uVF5, rho, omega, nodesSubZone, elemSubZone, nodes(:,1), orientation);
                        
                    elseif strcmp(elemType, 'C3D8')
                        % Selectively reduced integration type element - C3D8
                        [A, B, strain] = solveVFM_5p_C3D8(U, uVF1, uVF2, uVF3, uVF4, uVF5, rho, omega, nodesSubZone, elemSubZone, nodes(:,1), orientation);
                        
                    else
                        % Fully integrated element
                        [A, B, strain] = solveVFM_5p_C3D8F(U, uVF1, uVF2, uVF3, uVF4, uVF5, rho, omega, nodesSubZone, elemSubZone, nodes(:,1), orientation, GaussPoints);
                        
                    end
                    
                    moduli = A\B
                    
                    % Calculate parameters from elasticity matrix to calculate convergence
                    [Ei, G12i, G13i, E1i, E3i, v12i, v13i, v31i, dampi] = getParams_5p(moduli, 1, maxIter);
                    paramI = [G12i; G13i; E1i; E3i];
                    
                    % Calculate the percent change in estimated material parameters
                    %diffPerc = max(abs((moduli(3:4) - paramEst(3:4))./paramEst(3:4)))*100; % Only use parameters C44 and C66 to test for convergence - in incompressible media, C11, C33 and C13 are not reliably estimated
                    %diffPerc = max(abs((moduli - paramEst)./paramEst))*100;
                    diffPerc = max(abs((paramI - paramComp)./paramComp))*100;
                    disp(sprintf('Maximum change in parameters: %.2f%', diffPerc));
                    
                end
                
                % Get parameters from complex values of elasticity matrix (moduli)
                [E, G12, G13, E1, E3, v12, v13, v31, damp] = getParams_5p(moduli, iter, maxIter);
                
                % Save results as struct
                result = struct('E', E, 'G12', G12, 'G13', G13, 'E1', E1, 'E3', E3, 'v12', v12, 'v13', v13, 'v31', v31, ...
                    'damp', damp, 'uVF1', uVF1, 'uVF2', uVF2, 'uVF3', uVF3, 'uVF4', uVF4, 'uVF5', uVF5, 'iter', iter, 'eta', eta, 'strain', strain, 'orientation', orientation);
                Results(countZone) = result;
                
                clear moduli
                
                close all
                
            end
        end
    end
end

% Save results to mat file
matFile = sprintf('%s/optimisedVF_%dparam_%dsubzones_%s_BHmatrix.mat', outDir, numParam, countZone, elemType);
save(matFile, 'Results', 'maxIter');
disp(matFile)


diary off
% The end