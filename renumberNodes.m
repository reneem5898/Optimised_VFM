function [newNodes, newElems] = renumberNodes(nodes, elems)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function renumbers nodes from 1 and replaces old node numbers with
% new node numbers in elements
%
% Inputs: 1) nodes - # nodes x 4 matrix -- node num, x coord, y coord, z coord
%         2) elems - # elems x 9 matrix -- elem num, node1, node2, node3, etc.
%
% Outputs: 1) renumbered nodes
%          2) renumbered elements
%
% Written by: Renee Miller (reneem5898@gmail.com)
% Date Updated: 5 December 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create list of new node numbers from 1:length(nodeNums)
newNodeNums = linspace(1, size(nodes,1), size(nodes,1))';

% Initialise new elements list
newElems = elems;

% number of nodes per element
nodesPerElem = size(elems,2) - 1;

% Loop through elements and replace old node number with new node number (idx)
for i = 1:length(elems)
    newElems(i,1) = i;
    for j = 1:nodesPerElem
        idx = find(nodes(:,1) == elems(i,j+1));
        newElems(i,j+1) = idx;
    end
end

% Cat new node numbers with coordinates
newNodes = [newNodeNums nodes(:,2:end)];