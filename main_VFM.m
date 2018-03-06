function [subZoneResults] = main_VFM(modelDir, refParams, noise, num_MonteCarlo, subzones, rho, f, outDir, vs, elemType, GaussPoints, fibreNoise)

%% Virtual Fields Method %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function utilises the virutal fields method to solve for the material
% properties of a structure given: 1) node coordinates, 2) element connectivitiy
% and 3) nodal displacements
%
% Allows for estimating 1, 3 or 5 parameters for an isotropic or transversely
% isotropic material. Number of parameters to estimated is obtained by the length of "refParams".
% If refParams is empty, program will estimate 1 parameter (assume isotropy).
%
% Inputs: 1) modelDir - model directory where to find node, element, etc. data
%         2) refParams - reference or initial parameters to use in estimation
%         3) noise - list of amounts of noise to add to displacements (e.g. noise = linspace(0.15, 0.5, 11))
%         4) num_MonteCarlo - number of Monte Carlo simulations
% NB: analysis will run a Monte Carlo simulation (n = 30) for each amount of noise
%         4) subzones = 0/1 - 0, solve using whole model; 1, separate model into subzones
%         5) rho - density
%         6) f - frequency (Hz)
%         7) outDir - output directory
%         8) vs - varying start = 0/1 - 0, use reference parameters as initial estimate;
% 1, randomly choose initial estimate from Gaussian distribution centred at reference parameter
% with stdev = 20%*reference parameter
%         9) elemType - string: 'C3D8R', 'C3D8' or 'C3D8F'
%         10) GaussPoints - number of Gauss points - not necessary unless using C3D8F elements (GaussPoints must = 1 or 2)
%         11) fibreNoise - 0/1 - whether or not to add Gaussian noise to fibre orientations
%
% Output: subZoneResults - struct containing results from VFM analysis
%
% Written: Renee Miller
% Updated: 27 September 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants

% Type of noise added
noiseType = 'phaseNoise';
%noiseType = 'dispNoise';

% Material properties
omega = f*2*pi; % angular frequency in rad/s

% Maximum number of iterations per subzone - only relevant to anisotropic material parameter methods
maxIter = 30;

% Number of parameters to solve for: 1, 3, or 5
if ~isempty(refParams)
    numParam = length(refParams);
    
    if numParam ~= 3 && numParam ~= 5
        disp('You have not entered an appropriate number of parameters. Please enter 3 or 5 initial guesses.');
    end
else
    numParam = 1;
end

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
    
    orientation = [elems(:,1) dir1 dir2 dir3];
    
    % Add noise to fibre orientations
    if fibreNoise
        
        % Percentage of noise in fibre angles
        perc = 0.12; % 12% noise -- results in a variation in fibre angles with 2*std = 5degree
        
        orientation = addNoiseFibres(modelDir, orientation, perc);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gaussian noise

if ~isempty(noise)
    
    % Number of Monte-Carlo noise simulations
    %num_MCsim = 30;
    num_MCsim = 10;
    
    % Noise values - percentage of displacement stdev
    noiseVals = noise;
    
    disp('Creating noise copies...');
    
    % Add noise to phase images, then take fft
    if strcmp(noiseType,'phaseNoise')
        
        % Create directory to save results
        noiseDir = sprintf('%s/phaseNoise', outDir);
        
        if ~exist(noiseDir, 'dir')
            mkdir(noiseDir)
        end
        
        % Create noise copies
        addNoisePhase(outDir, noiseDir, num_MCsim, noiseVals);
        
        % Add noise directly to complex displacement data
    else
        % Create directory to save results
        noiseDir = sprintf('%s/dispNoise', outDir);
        
        if ~exist(noiseDir, 'dir')
            mkdir(noiseDir)
        end
        
        % Create noise copies
        addNoise(outDir, noiseDir, num_MCsim, noiseVals);
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subzones

if subzones
    
    %         xRange = [0 25; 25 50];
    %         yRange = [0 25; 25 50];
    %         zRange = [0 25; 25 50; 50 75; 75 100; 100 125; 125 150; 150 175; 175 200];
    
    % Ranges for LV
    xRange = [-14 16; 16 43];
    yRange = [-50 50];
    zRange = [-50 50];
    
    % Number of subzones
    numZones = size(xRange,1)*size(yRange,1)*size(zRange,1);
else
    % Node range == entire model
    xRange = [min(nodes(:,2)) max(nodes(:,2))];
    yRange = [min(nodes(:,3)) max(nodes(:,3))];
    zRange = [min(nodes(:,4)) max(nodes(:,4))];
    
    % Number of subzones
    numZones = size(xRange,1)*size(yRange,1)*size(zRange,1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve for parameters

%% Solving for parameters with Gausssian noise added to the reference displacements
if ~isempty(noise)
    
    % Loop through simulations with noise
    for p = 1:length(noiseVals)
        for q = 1:num_MCsim
            
            % Simulation number
            ns = (p-1)*num_MCsim + q;
            
            % Open displacement field with added Gaussian noise
            if strcmp(noiseType, 'phaseNoise')
                dispFile = strcat(noiseDir, '\uComplex_Noise', num2str(ns), '.txt');
            else
                dispFile = strcat(noiseDir, '\uComplex_noise_', num2str(ns), '.txt');
            end
            
            % Load displacements
            U = load(dispFile);
            U = U(:,1) + U(:,2)*1i;
                        
            % Count subzone
            countZone = 0;
            
            % Create initial parameters - if vs == 1, vary starting parameters, elseif vs == 0, use parameters given in refParams
            if vs
                paramInit = varyParamStart(refParams);
            else
                paramInit = refParams;
            end
            
            % Loop through subzones
            for m = 1:size(xRange,1)
                for n = 1:size(yRange,1)
                    for o = 1:size(zRange,1)
                        
                        
                        disp('Assembling element subzone node and element lists...');
                        % Counter for zone number
                        countZone = countZone + 1;
                        
                        % Creating subzone node and element lists
                        [nodesSubZone, elemSubZone] = getSubZone(nodes, elems, xRange, yRange, zRange, m, n, o);
                        surfNodes = intersect(nodesSubZone(:,1), bcNodesNew);
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % If solving for an isotropic shear modulus
                        if numParam == 1
                            
                            disp(sprintf('Running MC simulation #: %d\nSubzone #: %d\n', ns, countZone));
                            
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
                            
                            subZoneResult = struct('G', G, 'FK', fk, 'FG', fg, 'B', b, 'uVF', uVF, 'eta', eta, 'strain1', strain1, 'strain2', strain2);
                            subZoneResults(q,countZone) = subZoneResult;
                            
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % If solving for transversely isotropic material properties - 3 parameters
                        elseif numParam == 3
                            
                            % Iteratively solve for parameters until change in parameters is less than 0.1%
                            diffPerc = 1;
                            iter = 0;
                            while diffPerc > 0.1 && iter < maxIter
                                
                                % Count number of iterations
                                iter = iter + 1;
                                
                                disp(sprintf('Running MC simulation #: %d\nSubzone #: %d\nIteration #: %d\n', ns, countZone, iter));
                                
                                if ~exist('moduli','var')
                                    paramEst = paramInit;
                                else
                                    paramEst = moduli; % Parameters calculated from previous subzone
                                end
                                
                                % 3 parameter inversion
                                % Create numeric virtual displacement field
                                disp('Calculating the numeric virtual displacement fields...')
                                [uVF1, uVF2, uVF3, eta] = numericVF_3p(U, nodesSubZone, elemSubZone, nodes(:,1), paramEst, orientation, surfNodes, elemType, GaussPoints);
                                
                                % Calculate shear modulus
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
                                %diffPerc = max(abs((moduli(1:2) - paramEst(1:2))./paramEst(1:2)))*100;
                                disp(sprintf('Maximum change in parameters: %.2f%', diffPerc));
                                
                            end
                            
                            % Get parameters from complex values of elasticity matrix (moduli)
                            [G12, G13, damp, T] = getParams_3p(moduli, iter, maxIter);
                            
                            % Save results
                            subZoneResult = struct('G12', G12, 'G13', G13, 'T', T, 'damp', damp, 'uVF1', uVF1, 'uVF2', uVF2, 'uVF3', uVF3, ...
                                'iter', iter, 'eta', eta, 'K', K, 'strain', strain);
                            subZoneResults(q,countZone) = subZoneResult;
                            
                            clear moduli
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % If solving for transversely isotropic material properties - 5 parameters
                        elseif numParam == 5
                            
                            % Iteratively solve for parameters until change in parameters is less than 0.1%
                            diffPerc = 1;
                            iter = 0;
                            
                            while diffPerc > 0.1 && iter < maxIter
                                
                                % Count number of iterations
                                iter = iter + 1;
                                
                                disp(sprintf('Running MC simulation #: %d\nSubzone #: %d\nIteration #: %d\n', ns, countZone, iter));
                                
                                if ~exist('moduli','var')
                                    paramEst = paramInit; % Initial parameter guesses
                                else
                                    paramEst = moduli; % Parameters calculated from previous subzone
                                end
                                
                                % Calculate engineering constants to use for calculating convergence
                                [Ec, G12c, G13c, E1c, E3c, v12c, v13c, v31c, dampc] = getParams_5p(paramEst, 1, maxIter);
                                paramComp = [G12c; G13c; E1c; E3c];
                                
                                % 5 parameter inversion
                                % Create numeric virtual displacement field
                                disp('Calculating the numeric virtual displacement fields...')
                                [uVF1, uVF2, uVF3, uVF4, uVF5, eta] = numericVF_5p(U, nodesSubZone, elemSubZone, nodes(:,1), paramEst, orientation, surfNodes, elemType, GaussPoints);
                                
                                % Calculate shear modulus
                                disp('Calculating the estimated parameters...')
                                
                                if strcmp(elemType, 'C3D8R')
                                    % Uniform strain elements - C3D8R
                                    [A, B, strain2] = solveVFM_5p_C3D8R(U, uVF1, uVF2, uVF3, uVF4, uVF5, rho, omega, nodesSubZone, elemSubZone, nodes(:,1), orientation);
                                    
                                elseif strcmp(elemType, 'C3D8')
                                    % Selectively reduced integration type element - C3D8
                                    [A, B, strain2] = solveVFM_5p_C3D8(U, uVF1, uVF2, uVF3, uVF4, uVF5, rho, omega, nodesSubZone, elemSubZone, nodes(:,1), orientation);
                                    
                                else
                                    % Fully integrated element
                                    [A, B, strain2] = solveVFM_5p_C3D8F(U, uVF1, uVF2, uVF3, uVF4, uVF5, rho, omega, nodesSubZone, elemSubZone, nodes(:,1), orientation, GaussPoints);
                                    
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
                            
                            subZoneResult = struct('E', E, 'G12', G12, 'G13', G13, 'E1', E1, 'E3', E3, 'v12', v12, 'v13', v13, 'v31', v31, ...
                                'damp', damp, 'uVF1', uVF1, 'uVF2', uVF2, 'uVF3', uVF3, 'uVF4', uVF4, 'uVF5', uVF5, 'iter', iter, 'eta', eta, 'strain2', strain2);
                            subZoneResults(q,countZone) = subZoneResult;
                            
                            clear moduli
                            
                        end
                    end
                end
            end
        end
        
        % Save results
        if bcNodes
            save(sprintf('%s/numericVF_%dparam_%.3fnoise_%dsubzones_%s_vs%d_bcNodes.mat', noiseDir, numParam, noiseVals(p), countZone, elemType, vs), 'subZoneResults');
        else
            save(sprintf('%s/numericVF_%dparam_%.3fnoise_%dsubzones_%s_vs%d.mat', noiseDir, numParam, noiseVals(p), countZone, elemType, vs), 'subZoneResults');
        end
        
    end
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Solving for parameters WITHOUT noise
    
    % Open measured displacement field
    dispFile = strcat(outDir, '/uComplex.txt');
    U = load(dispFile);
    U = U(:,1) + U(:,2)*1i;
    
    % Count subzone
    countZone = 0;
    
    % Create initial parameters - if vs == 1, vary starting parameters, elseif vs == 0, use parameters given in refParams
    if vs
        paramInit = varyParamStart(refParams);
    else
        paramInit = refParams;
    end
    
    % Loop through subzones
    for m = 1:size(xRange,1)
        for n = 1:size(yRange,1)
            for o = 1:size(zRange,1)
                
                disp('Assembling element subzone node and element lists...');
                
                % Counter for zone number
                countZone = countZone + 1;
                
                % Creating subzone node and element lists
                [nodesSubZone, elemSubZone] = getSubZone(nodes, elems, xRange, yRange, zRange, m, n, o);
                surfNodes = intersect(nodesSubZone(:,1), bcNodesNew);
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % If solving for an isotropic shear modulus
                if numParam == 1
                    
                    disp(sprintf('Running Subzone #: %d\n', countZone));
                    
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
                    
                    subZoneResult = struct('G', G, 'FK', fk, 'FG', fg, 'B', b, 'uVF', uVF, 'eta', eta, 'strain1', strain1, 'strain2', strain2);
                    subZoneResults(countZone) = subZoneResult;
                    
                    
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % If solving for transversely isotropic material properties - 3 parameters
                elseif numParam == 3
                    
                    % Iteratively solve for parameters until change in parameters is less than 0.1%
                    diffPerc = 1;
                    iter = 0;
                    while diffPerc > 0.1 && iter < maxIter
                        
                        % Count number of iterations
                        iter = iter + 1;
                        
                        disp(sprintf('Running Subzone #: %d\nIteration #: %d\n', countZone, iter));
                        
                        if ~exist('moduli','var')
                            paramEst = paramInit;
                        else
                            paramEst = moduli; % Parameters calculated from previous subzone
                        end
                        
                        % 3 parameter inversion
                        % Create numeric virtual displacement field
                        disp('Calculating the numeric virtual displacement fields...')
                        [uVF1, uVF2, uVF3, eta] = numericVF_3p(U, nodesSubZone, elemSubZone, nodes(:,1), paramEst, orientation, surfNodes, elemType, GaussPoints);
                        
                        % Calculate shear modulus
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
                    
                    % Save results
                    subZoneResult = struct('G12', G12, 'G13', G13, 'T', T, 'damp', damp, 'uVF1', uVF1, 'uVF2', uVF2, 'uVF3', uVF3, 'iter', iter, 'eta', eta, 'strain', strain, 'K', K, 'orientation', orientation);
                    subZoneResults(countZone) = subZoneResult;
                    
                    clear moduli
                    
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % If solving for transversely isotropic material properties - 5 parameters
                elseif numParam == 5
                    
                    % Iteratively solve for parameters until change in parameters is less than 0.1%
                    diffPerc = 1;
                    iter = 0;
                    
                    while diffPerc > 0.1 && iter < maxIter
                        
                        % Count number of iterations
                        iter = iter + 1;
                        
                        disp(sprintf('Running Subzone #: %d\nIteration #: %d\n', countZone, iter));
                        
                        if ~exist('moduli','var')
                            paramEst = paramInit; % Initial parameter guesses
                        else
                            paramEst = moduli; % Parameters calculated from previous subzone
                        end
                        
                        % Calculate engineering constants to use for calculating convergence
                        [Ec, G12c, G13c, E1c, E3c, v12c, v13c, v31c, dampc] = getParams_5p(paramEst, 1, maxIter);
                        paramComp = [G12c; G13c; E1c; E3c];
                        
                        % 5 parameter inversion
                        % Create numeric virtual displacement field
                        disp('Calculating the numeric virtual displacement fields...')
                        [uVF1, uVF2, uVF3, uVF4, uVF5, eta] = numericVF_5p(U, nodesSubZone, elemSubZone, nodes(:,1), paramEst, orientation, surfNodes, elemType, GaussPoints);
                        
                        % Calculate shear modulus
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
                    
                    subZoneResult = struct('E', E, 'G12', G12, 'G13', G13, 'E1', E1, 'E3', E3, 'v12', v12, 'v13', v13, 'v31', v31, ...
                        'damp', damp, 'uVF1', uVF1, 'uVF2', uVF2, 'uVF3', uVF3, 'uVF4', uVF4, 'uVF5', uVF5, 'iter', iter, 'eta', eta, 'strain', strain, 'orientation', orientation);
                    subZoneResults(countZone) = subZoneResult;
                    
                    clear moduli
                    
                    close all
                    
                end
            end
        end
    end
    
    % Save results
    if bcNodes
        save(sprintf('%s/numericVF_%dparam_%dsubzones_%s_vs%d_bcNodes.mat', outDir, numParam, countZone, elemType, vs), 'subZoneResults');
    else
        save(sprintf('%s/numericVF_%dparam_%dsubzones_%s_vs%d.mat', outDir, numParam, countZone, elemType, vs), 'subZoneResults');
    end
    
end

% The end