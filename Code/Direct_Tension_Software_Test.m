clc; clear; close all;

% This is a testing code which runs the direct tension test and compares
% the F-D results with reference data in the same folder. For more info, see the documentation
% This code is compatible with MATLAB version R2022b and later.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%============================ USER INPUTS ================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open the 'FEM_Main_Script_2D.m' and set model_name="Direct_tension_test".

% Enter location of current file path
main_file_path=''; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=================== NO CHANGES NEEDED AFTER THIS LINE ===================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the main script
run('FEM_Main_Script_2D.m');  

% Load the reference F–D data
ref_path = load(strcat(main_file_path,'\ref_fd_data.mat'),'plot_storage'); % Load ref data
ref_disp  = ref_path.plot_storage(:,1);                     % Reference displacements
ref_force = ref_path.plot_storage(:,2);                     % Reference forces

% Open the .mat of the latest increment from the current run of the code
files_path= dir(strcat(main_file_path,'\Saved Files\*.mat'));  % Saved Files of current run
[~, idx] = max([files_path.datenum]);
curr = load(fullfile(files_path(idx).folder, files_path(idx).name), 'plot_storage');
cur_disp  = curr.plot_storage(:,1);                    % Reference displacements
cur_force = curr.plot_storage(:,2);                    % Reference forces

% Error calculations
tol         = 1e-5;                                    % Tolerance to pass the test
max_err     = (abs(max(cur_force) - max(ref_force)));  % Find error between max forces
l2_err      = norm(cur_force - ref_force);             % L2 norm between the two curves

disp( repmat('-', 1, 50) );
fprintf("F–D Error Metrics: \n");
fprintf("  Max force error:      %.3e\n", max_err);
fprintf("  L2 norm error:       %.3e\n", l2_err);

if max_err > tol
    fprintf("  Test Failed (%.3e > tol %.1e)\n\n", max_err, tol);
else
    fprintf("  Test Passed (%.3e ≤ tol %.1e)\n\n", max_err, tol);
end

% Plot F-D curves
figure('Color','w');
plot(ref_disp, ref_force,  'k-','LineWidth',1.5); hold on;
plot(ref_disp, cur_force,  'r--','LineWidth',1.5);
xlabel('Displacement'); ylabel('Force');
legend('Reference','Current','Location','best');
title("Direct Tension Test F-D Curves", 'Interpreter','latex','FontSize',14,'Color','k');
grid on;
