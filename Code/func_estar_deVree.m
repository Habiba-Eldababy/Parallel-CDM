function [estar, s] = func_estar_deVree(exx, eyy, gxy, k_damage_parameter, nu)

% Compute the principal strains
e_princ_i = (exx + eyy)/2 + sqrt(((exx-eyy)/2)^2 + (gxy/2)^2);
e_princ_ii = (exx + eyy)/2 - sqrt(((exx-eyy)/2)^2 + (gxy/2)^2);

% Compute the strain invariants 
I1 = e_princ_i + e_princ_ii;                                                       % First strain invariant 
J2 = (2 * (e_princ_i^2)) + (2 * (e_princ_ii^2)) - (2 * (e_princ_i * e_princ_ii));  % Second deviatoric invariant of strain

% Compute the deVree equivalent strain
A = (k_damage_parameter - 1)/(1 - (2 * nu));
B = (12 * k_damage_parameter)/((1 + nu)^2);

estar = ((A/(2 * k_damage_parameter)) * I1) + ((1/(2 * k_damage_parameter)) * sqrt((A*I1)^2 + (B * J2)));


% Compute partial derivative components 

% a) destar_dI1
destar_dI1 = (A/(2 * k_damage_parameter))    +  ( (I1/k_damage_parameter) * (A^2) )/( 2 * sqrt( (A * I1)^2 + (B * J2) ) ) ; 

% b) destar_dJ2
destar_dJ2 = B / (4 * k_damage_parameter * sqrt( (A*I1)^2 + (B*J2) ));

% c) dI1_dep1
dI1_dep1 = 1;

% d) dI1_dep2
dI1_dep2 = 1;

% e) dJ2_dep1
dJ2_dep1 = (4 * e_princ_i) - (2 * e_princ_ii);

% f) dJ2_dep2
dJ2_dep2 = (4 * e_princ_ii) - (2 * e_princ_i);

% g) dep1_dexx 
dep1_dexx = (1/2) * ( 1 + (exx - eyy)/(sqrt( (exx - eyy)^2 + (gxy)^2 )));

% h) dep1_deyy 
dep1_deyy = (1/2) * ( 1 - (exx - eyy)/(sqrt( (exx - eyy)^2 + (gxy)^2 )));

% i) dep1_dgxy 
dep1_dgxy = (1/2) * ( (gxy) / (sqrt( (exx - eyy)^2 + (gxy)^2 )) );

% j) dep2_dexx 
dep2_dexx = (1/2) * ( 1 - (exx - eyy)/(sqrt( (exx - eyy)^2 + (gxy)^2 )));

% k) dep2_deyy 
dep2_deyy = (1/2) * ( 1 + (exx - eyy)/(sqrt( (exx - eyy)^2 + (gxy)^2 )));

% l) dep2_dgxy 
dep2_dgxy = -(1/2) * ( (gxy) / (sqrt( (exx - eyy)^2 + (gxy)^2 )) );




% Compute the partial derivatives of estar
destar_dexx = (destar_dI1 * ((dI1_dep1 * dep1_dexx) + (dI1_dep2 * dep2_dexx))) + (destar_dJ2 * ((dJ2_dep1 * dep1_dexx) + (dJ2_dep2 * dep2_dexx)));

destar_deyy = (destar_dI1 * ((dI1_dep1 * dep1_deyy) + (dI1_dep2 * dep2_deyy))) + (destar_dJ2 * ((dJ2_dep1 * dep1_deyy) + (dJ2_dep2 * dep2_deyy)));

destar_dgxy = (destar_dI1 * ((dI1_dep1 * dep1_dgxy) + (dI1_dep2 * dep2_dgxy))) + (destar_dJ2 * ((dJ2_dep1 * dep1_dgxy) + (dJ2_dep2 * dep2_dgxy)));

% Gather the derivatives in vector s
s = [destar_dexx; destar_deyy; destar_dgxy];

end