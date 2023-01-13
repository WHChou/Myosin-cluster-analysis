function proj = stackZproj(stack, znum, thickness)
%%%%%%%%%%%%%%%%%%%
% Given a z stack and z number, compute average projection around the z num
% Input: 
%     stack: a 3-D stack
%     znum: The stack index to be averaged around
%     thickness (optional): how many slices around znum to average. Default
%                                         is 1
%%%%%%%%%%%%%%%%%%%
if nargin <3
    thickness = 1;
end
if size(stack, 3)==1
    error('Only z stacks can be processed!')
end

ub = znum+thickness;
lb = znum-thickness;
if ub>size(stack, 3)
    ub = size(stack, 3);
    warning('Upper bound exceed stack size');
end
if lb<1
    lb = 1;
    warning('Lower bound smaller than 1');
end
proj = uint16(mean(stack(:,:,lb:ub), 3));
