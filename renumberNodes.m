function [newNodeNums, newElems] = renumberNodes(nodeNums, elems)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function renumbers nodes from 1 and replaces old node numbers with
% new node numbers in elements
%
% Inputs: 1) nodes - # nodes x 4 matrix -- node num, x coord, y coord, z coord
%         2) elems - # elems x 9 matrix -- elem num, node1, node2, node3, etc.
%
% Outputs: 1) renumbered node list
%          2) new element list
%
% Written by: Renee Miller
% Date: 22 June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create list of new node numbers from 1:length(nodeNums)
newNodeNums = linspace(1, length(nodeNums), length(nodeNums));

% Initialise new elements list
newElems = elems;

% number of nodes per element
nodesPerElem = size(elems,2) - 1;

% Loop through elements and replace old node number with new node number (idx)
for i = 1:length(elems)
    for j = 1:nodesPerElem
        idx = find(nodeNums == elems(i,j+1));
        newElems(i,j+1) = idx;
    end
end
    