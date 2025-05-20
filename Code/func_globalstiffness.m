function [K, J, DRdu, Res_F, F_int, F_ext, history_var_mat, gausspoints_prop_mat, nodes_prop_mat, strain_var_mat] = func_globalstiffness(Damage_type,k_damage_parameter,eq_strain_type,dofs,Delastic,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,strain_tolerance,strain_mat_previousinc,IsProj,RoutineID)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====================== ASSEMBLE THE TANGENT MATRIX ======================
% ===================== ASSEMBLE THE RESIDUAL VECTOR ======================
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Include global variables
func_include_flags;

% Extract Poissons ratio
nu = materialprops(2);


% Variables strain_var_mat and strain_mat_previousinc are only active in
% the UAL local damage case
if RoutineID == 2 || SolverID == 2
    strain_var_mat = [];
    strain_mat_previousinc = zeros(size(history_var_mat_previousinc,1),size(history_var_mat_previousinc,2));
end

if IsProj == 0
    gausspoints_prop_mat = [];
    nodes_prop_mat = [];

    DRdu = zeros(ndof * nnodes, ndof * nnodes);

    % Create local copies of global variables to use inside the spmd block
    local_alpha_val = alpha_val; local_TangentID = TangentID; local_SolverID = SolverID; local_beta_val = beta_val; local_e_delta = e_delta; local_dmax = dmax; local_nelnodes = nelnodes; local_ndof = ndof;
    local_nnodes = nnodes; local_nelem = nelem; local_ncoord = ncoord; local_maxnodes = maxnodes; local_ndof2 = ndof2; local_connect_nds = connect_nds; local_coords  = coords; local_dofs = dofs;
    local_strain_mat_previousinc = strain_mat_previousinc; local_g = g;

    % Start spmd block
    spmd
        num_workers = spmdSize; % Get number of workers in the parallel pool
        local_elems = floor(local_nelem / num_workers); % Calculate number of elements per worker

        start_elem = (spmdIndex - 1) * local_elems + 1; % Set element range
        if spmdIndex == num_workers
            end_elem = local_nelem;
        else
            end_elem = spmdIndex * local_elems;
        end

        % Preallocate local matrices
        local_J = sparse(local_ndof*local_nnodes, local_ndof*local_nnodes);
        local_K = sparse(local_ndof2*local_nnodes, local_ndof2*local_nnodes);
        local_Res_F = sparse(local_ndof*local_nnodes, 1);
        local_F_int = sparse(local_ndof*local_nnodes, 1);
        local_F_ext = sparse(local_ndof*local_nnodes, 1);
        local_history_var_mat = sparse(local_nelem, 4);
        local_strain_var_mat = [];

        % Extract coords of nodes, DOF for the local elements
        for lmn = start_elem:end_elem

            lmncoord = zeros(local_ncoord, local_maxnodes);
            lmndof = zeros(local_ndof, local_maxnodes);

            for a = 1:local_nelnodes(lmn)
                for i = 1:local_ncoord
                    lmncoord(i, a) = local_coords(i, local_connect_nds(a, lmn));
                end
                for i = 1:local_ndof
                    lmndof(i, a) = local_dofs(local_ndof * (local_connect_nds(a, lmn) - 1) + i);
                end
            end

            % ANALYTICAL CALCULATION - TANGENT MATRIX (for each element):
            if local_TangentID == 1
                if local_SolverID == 1
                    [k_el, j_el, local_history_var_mat(lmn, :), local_strain_var_mat(lmn, :), ~, ~, Res_F_el, f_internal_el, f_external_el] = ...
                        func_elstif_Local(Damage_type, nu, k_damage_parameter,eq_strain_type, lmncoord, lmndof, Delastic, history_var_mat_previousinc(lmn, :), local_alpha_val, local_beta_val, local_e_delta, local_dmax, n_hood, weights, strain_tolerance, local_strain_mat_previousinc(lmn, :), 1, IsProj, RoutineID);
                elseif local_SolverID == 2
                    [k_el, j_el, local_history_var_mat(lmn, :), ~, ~, Res_F_el, f_internal_el, f_external_el] = ...
                        func_elstif_Nonlocgradient(Damage_type, nu, k_damage_parameter,eq_strain_type, lmncoord, lmndof, Delastic, history_var_mat_previousinc(lmn, :), local_g, local_alpha_val, local_beta_val, local_e_delta, local_dmax, strain_tolerance, n_hood, weights, 1, IsProj);
                else
                    disp("Check your SolverID - globalstiffness")
                end
            end

            % Assemble J and K matrices for this group of elements
            for a = 1:local_nelnodes(lmn)
                for i = 1:local_ndof
                    for b = 1:local_nelnodes(lmn)
                        for k = 1:local_ndof
                            rw = local_ndof * (local_connect_nds(a, lmn) - 1) + i;
                            cl = local_ndof * (local_connect_nds(b, lmn) - 1) + k;
                            local_J(rw, cl) = local_J(rw, cl) + j_el(local_ndof * (a - 1) + i, local_ndof * (b - 1) + k); % J matrix assembly
                        end
                    end
                end

                for i = 1:local_ndof2
                    for b = 1:local_nelnodes(lmn)
                        for k = 1:local_ndof2
                            rw = local_ndof2 * (local_connect_nds(a, lmn) - 1) + i;
                            cl = local_ndof2 * (local_connect_nds(b, lmn) - 1) + k;
                            local_K(rw, cl) = local_K(rw, cl) + k_el(local_ndof2 * (a - 1) + i, local_ndof2 * (b - 1) + k); % K matrix assembly
                        end
                    end
                end

            end

            % Calculate residual vector, internal force vector, external force vector
            for a = 1:local_nelnodes(lmn)
                for i = 1:local_ndof
                    rw = local_ndof * (local_connect_nds(a, lmn) - 1) + i;
                    local_Res_F(rw, 1) = local_Res_F(rw, 1) + Res_F_el(local_ndof * (a - 1) + i, 1); % residual vector
                    local_F_int(rw, 1) = local_F_int(rw, 1) + f_internal_el(local_ndof * (a - 1) + i, 1); % internal force vector
                    local_F_ext(rw, 1) = local_F_ext(rw, 1) + f_external_el(local_ndof * (a - 1) + i, 1); % external force vector
                end
            end
        end

        spmdBarrier; % Wait for all workers to complete their tasks

        % Combine results across workers
        J = spmdPlus(local_J, 1);
        K = spmdPlus(local_K, 1);
        Res_F = spmdPlus(local_Res_F, 1);
        F_int = spmdPlus(local_F_int, 1);
        F_ext = spmdPlus(local_F_ext, 1);
        history_var_mat = spmdCat(local_history_var_mat, 1, 1);
        strain_var_mat = spmdCat(local_strain_var_mat,1,1);
    end

    % After spmd block, extract the data from composite variables
    J = J{1};
    K = K{1};
    Res_F = Res_F{1};
    F_int = F_int{1};
    F_ext = F_ext{1};
    history_var_mat = history_var_mat{1};
    strain_var_mat = strain_var_mat{1};

elseif IsProj == 1
    % -------------------------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%% PLOTTING PROPERTIES ROUTINE %%%%%%%%%%%%%%%%%%%%%%%
    % -------------------------------------------------------------------------

    % Create empty entries for solver variables
    J               = [];
    K               = [];
    DRdu            = [];
    Res_F           = [];
    history_var_mat = [];
    F_int           = [];
    F_ext           = [];
    strain_var_mat  = [];

    % Create local copies of global variables for parfor loop
    local_alpha_val = alpha_val; local_SolverID = SolverID; local_beta_val = beta_val; local_e_delta = e_delta; local_dmax = dmax; local_nelnodes = nelnodes; local_connect_nds = connect_nds;
    local_strain_mat_previousinc = strain_mat_previousinc; local_nelem = nelem; local_g = g;
    
    % Initializing element/nodal properties and gather matrices at zero
    gausspoints_prop_mat = zeros(nelem * func_numberofintegrationpoints,9);
    nodes_prop_mat       = zeros(nnodes,9);

    % Preallocate temporary matrices for parallel computation
    temp_lmncoord             = zeros(nelem,ncoord,maxnodes);
    temp_lmndof               = zeros(nelem,ndof,maxnodes);
    temp_gausspoints_prop_mat_cell = cell(nelem,1); % 2500x1 matrix
    temp_nodes_prop_mat  = zeros(size(nodes_prop_mat));  % 2601x9 temp matrix

    % Extract coords of nodes, DOF for all elements
    for lmn = 1:nelem
        for a = 1:nelnodes(lmn)
            for i = 1:ncoord
                temp_lmncoord(lmn,i,a) = coords(i,connect_nds(a,lmn));
            end
            for i = 1:ndof
                temp_lmndof(lmn,i,a) = dofs(ndof*(connect_nds(a,lmn)-1)+i);
            end
        end
    end

    % parallel for loop to populate properties matrix before plotting
    parfor lmn = 1:local_nelem

        %Initialize temporary element matrices for the current element
        temp_gausspoints_prop_mat_elem = zeros(4, size(gausspoints_prop_mat, 2)); % 4x9 temporary matrix for each elem
        temp_nodes_prop_mat_elem = zeros(size(nodes_prop_mat)); % 2601x9 temp matrix

        lmncoord = squeeze(temp_lmncoord(lmn,:,:));
        lmndof = squeeze(temp_lmndof(lmn,:,:));

        % Compute element/nodal properties
        % -----------------------------------------------------------------
        if local_SolverID == 1
            [~, ~, ~, ~, temp_gausspoints_prop_mat_elem, nodes_prop_mat_elem, ~, ~, ~] = func_elstif_Local(Damage_type, nu, k_damage_parameter,eq_strain_type,lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),local_alpha_val,local_beta_val,local_e_delta,local_dmax,n_hood,weights,strain_tolerance,local_strain_mat_previousinc(lmn,:),1,IsProj,RoutineID);
            %------------------------------------------------------------------
        elseif local_SolverID == 2
            [~, ~, ~, temp_gausspoints_prop_mat_elem, nodes_prop_mat_elem, ~, ~, ~] = func_elstif_Nonlocgradient(Damage_type, nu, k_damage_parameter,eq_strain_type,lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),local_g,local_alpha_val,local_beta_val,local_e_delta,local_dmax,strain_tolerance,n_hood,weights,1,IsProj);
            % -----------------------------------------------------------------
        else
            disp("Check your SolverID - globalstiffness")
        end

        % Populate the properties matrix
        for a = 1:local_nelnodes(lmn)
            rw = local_connect_nds(a,lmn);
            temp_nodes_prop_mat_elem(rw,:) = nodes_prop_mat_elem(a,:);
        end

        % Aggregate results from this iteration into temporary arrays
        temp_nodes_prop_mat = temp_nodes_prop_mat + temp_nodes_prop_mat_elem;

        % Store results in temporary variables
        temp_gausspoints_prop_mat_cell{lmn} = temp_gausspoints_prop_mat_elem;
    end

    % Concatenate cell array into matrix
    temp_gausspoints_prop_mat = cat(1, temp_gausspoints_prop_mat_cell{:});

    % Calculate the final matrix of nodal properties
    nodes_prop_mat = temp_nodes_prop_mat ./ num_elem_at_node;

    % Assign temporary gausspoints_prop_mat to the original variable
    gausspoints_prop_mat = temp_gausspoints_prop_mat;
else

    disp("Check your IsProj variable - globalstiffness")

end
end
