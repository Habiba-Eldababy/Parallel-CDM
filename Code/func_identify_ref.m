function [ref, ref_node_ID, ref_node_pos] = func_identify_ref(direction_load, fixnodes, ndof)
% This function identfies the reference dof to measure the displacement applied at each increment in the case of NR and FAL routines 

% Identify the first node with a displacement in the +ve Y direction
% (assuming displacement is applied along the Y axis)
for i=1:1:size(fixnodes,2)
    if fixnodes(3,i)>0 && fixnodes(2,i) == direction_load 
        ref_node_ID  = fixnodes(1,i);
        ref_node_pos = i; 
        break
    end
end

ref = (ref_node_ID - 1) * ndof + direction_load;

end