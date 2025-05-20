function [kappa, omega, domega_dkappa] = func_geersmodel_Nonlocgradient(nonlocal_strain,kappa_previousinc,alpha_val,beta_val,e_delta,dmax,strain_tolerance,IsM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ========================== GEERS DAMAGE MODEL =========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computes the damage variable omega and its gradient wrt to the maximum
% nonlocal strain history variable (kappa) at this Gauss point based on:
% nonlocal_strain   : nonlocal equivalent strain
% e_delta           : threshold strain
% alpha             : model hyperparameter
% beta              : model hyperparameter

% IsM = 0: This is called only when we calculate the residual at plus/minus positions for the numerical tangent
if IsM == 0
    
    % -------------------------------------------------------------------------
    % Compute the nonlocal strain history variable
    kappa = nonlocal_strain; 
    
    % -------------------------------------------------------------------------
    % Compute the damage variable
     omega = max(0,dmax*(1 - (e_delta/kappa * ( (1 - alpha_val) + (alpha_val * exp(beta_val * (e_delta - kappa)) ) ) ) ) );
    
    % -------------------------------------------------------------------------
    % Compute the derivative of damage w.r.t. the nonlocal strain history variable
    domega_dkappa = 0;

    
% IsM = 1: This is called when we calculate the analytical tangent and the residual 
elseif IsM == 1

    % -------------------------------------------------------------------------
    % Compute the nonlocal strain history variable
    if (kappa_previousinc - nonlocal_strain > strain_tolerance)
        kappa = max(nonlocal_strain,kappa_previousinc);
    else 
        kappa = nonlocal_strain;
    end
    
    % -------------------------------------------------------------------------
    % Compute the damage variable
    omega = max(0,dmax*(1 - (e_delta/kappa * ( (1 - alpha_val) + (alpha_val * exp(beta_val * (e_delta - kappa)) ) ) ) ) );
    
    % -------------------------------------------------------------------------
    % Compute the derivative of damage w.r.t. the nonlocal strain history variable
    if (omega == 0) || (omega >= dmax) || (kappa_previousinc - nonlocal_strain > strain_tolerance)
        domega_dkappa = 0;
    else
        domega_dkappa = ((e_delta/(kappa^2)) * (1 - alpha_val + (alpha_val * exp(beta_val * (e_delta - kappa))))) + ((e_delta/kappa) * (alpha_val * beta_val * exp(beta_val * (e_delta - kappa))) ) ;        
    end

else
    
    disp("Check your IsM condition!")

end
end

