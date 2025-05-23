function [Beta] = func_calc_Beta(Constraint_type,eq_strain_type,dofs,Delastic,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,strain_tolerance,strain_var_mat_previousinc,IsProj,RoutineID)
% This functon calculates Beta

if Constraint_type == 1      % Cylindrical constraint
    Beta = 0;

elseif Constraint_type == 2  % Spherical constraint
    % Calculate stiffness
    [~,J, ~, ~, ~, ~, ~, ~,~,~] = parallel_func_globalstiffness(eq_strain_type,dofs,Delastic,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,strain_tolerance,strain_var_mat_previousinc,IsProj,RoutineID);

    % Calculate Beta
    X   = size(J,1);
    sum = 0;
    for i = 1:1:X
        sum = sum+J(i,i);
    end

    Beta = inv(sum/X);
end
end
