function onepatch_func_plotmesh2(coords,connect_nds,nelem,damage,elocal,enonlocal,nelnodes,color,model_name,inc_success_counter,last_iteration,increment,loadfactor,solver_string,tangent_string,SolverID)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =============================== PLOTTING ================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure("Position",[400 200 1500 500]);
t = tiledlayout(1,3,'TileSpacing','compact');
% title(t, strcat(model_name + ": Increment #" + int2str(increment) + " - Iteration #" + int2str(last_iteration) + " - Loadfactor " + int2str(loadfactor)), 'interpreter','latex','fontsize',12,'Color','k');

% -------------------------------------------------------------------------
% ------------------------ PLOTTING DAMAGE CONTOUR ------------------------
% ------------------------------------------------------------------------- 
nexttile; % First figure - tile 1 Damage contour
hold on; box on;
colormap(jet);  %https://www.mathworks.com/help/matlab/ref/colormap.html
title("Damage contours", 'interpreter','latex','fontsize',14,'Color','k');
axis tight
axis off
colorbar

% Preallocating matrices
all_x = zeros(nelem,4,2);               % to store all element coordinates
all_x_reshaped = zeros(nelem*4,2);      % will contain the same as all_x but correctly ordered for patch function

all_col_damage= zeros(nelem, 4, 1);     % to store all damage 
all_col_elocal = zeros(nelem, 4, 1);
all_col_enonlocal = zeros(nelem, 4,1);

% Populating x, damage, elocal, enonlocal matrices for all elements
for lmn = 1:nelem
   for i = 1:nelnodes(lmn)
       for j = 1:2 
        all_x(lmn,i,j) = coords(j,connect_nds(i,lmn));
       end
       all_col_damage(lmn,i,1) = damage(connect_nds(i,lmn));
       all_col_elocal(lmn,i,1) = elocal(connect_nds(i,lmn));
       all_col_enonlocal(lmn,i,1) = enonlocal(connect_nds(i,lmn));
   end
end

% Organization of nodal points for all elements in patch accepted formats
% from size (num elements x 4 x 2) to ((num elements * 4) x 2) in chronological order

idx1 = 1; % matrix indeces for the reshaped coordinates matrix
idx2 = 1;
for lmn = 1:nelem
   if (nelnodes(lmn)~=4)
    disp("Check your nelnodes. Element number: ", lmn)
   end   
   for i = 1:nelnodes(lmn)
      for j = 1:2
        all_x_reshaped(idx1,idx2) = all_x(lmn,i,j);
        if idx2 ==1
            idx2=2;
        elseif idx2 ==2
            idx2=1;
        end
      end
      idx1 = idx1+1;
    end
end 

% Reshaping matrices to prepare for patch function 
all_faces = reshape(1:size(all_x, 1)*4, 4, []).';           % nelem x 4 where each row is [1 2 3 4] 

all_col_damage = all_col_damage';                           % transpose so that reshaping happens row-wise
all_col_damage_reshaped = reshape(all_col_damage, [], 1);   % reshape into (nelem * 4) x 1

all_col_elocal = all_col_elocal';                           
all_col_elocal_reshaped = reshape(all_col_elocal, [], 1);

all_col_enonlocal = all_col_enonlocal';
all_col_enonlocal_reshaped = reshape(all_col_enonlocal, [], 1);

% Plotting patch for all elements at once - damage contour
patch('Vertices', all_x_reshaped, 'Faces', all_faces, 'FaceVertexCData', all_col_damage_reshaped, 'FaceColor', 'interp', 'EdgeColor', 'k');

hold off

% -------------------------------------------------------------------------
% ----------------------- PLOTTING ELOCAL CONTOUR -------------------------
% -------------------------------------------------------------------------
nexttile
hold on; box on;
colormap(jet);  %https://www.mathworks.com/help/matlab/ref/colormap.html
colorbar('AxisLocation','in')

% Plotting patch for all elements at once - elocal contour
patch('Vertices', all_x_reshaped, 'Faces', all_faces, 'FaceVertexCData', all_col_elocal_reshaped, 'FaceColor', 'interp', 'EdgeColor', 'none');

colorbar
title("Local strain contours", 'fontweight','bold','interpreter','latex','fontsize',14,'Color','k');

axis tight
axis off
hold off

% -------------------------------------------------------------------------
% --------------------- PLOTTING ENONLOCAL CONTOUR ------------------------
% -------------------------------------------------------------------------
if SolverID == 2
    nexttile
    hold on; box on;
    colormap(jet);  %https://www.mathworks.com/help/matlab/ref/colormap.html
    colorbar('AxisLocation','in')

    % Plotting patch for all elements at once - enonlocal contour 
    patch('Vertices',all_x_reshaped,'Faces',all_faces,'FaceVertexCData',all_col_enonlocal_reshaped,'FaceColor','interp','EdgeColor','none');
 
    title("Nonlocal strain contours", 'interpreter','latex','fontsize',14,'Color','k');
    axis tight
    axis off
    hold off
end

saveas(f1, strcat(model_name + "_contours_" + solver_string + "_" + tangent_string + "_inc_" + int2str(increment) + ".png"))


% %--------------------------------------------------------------------------
% % --------------------- PLOTTING FORCE-DISPLACEMENT -----------------------
% % -------------------------------------------------------------------------
% f2 = figure("Position",[400 200 1500 800]);
% hold on;
% plot(plot_storage(:,1),plot_storage(:,2))
% title("Force-Displacement using UAL w/ nonlocal gradient damage")
% 
% saveas(f2, strcat(model_name + "_force_disp_" + solver_string + "_" + tangent_string + "_inc_" + int2str(increment) + ".png"))

% close all 
end


