L1_range = linspace(10, 25, 10);
L4_range = linspace(10, 19, 10);
[LL1, LL4] = meshgrid(L1_range, L4_range);
L0 = [12,20,45,15]; 

global_manipulability_surface = zeros(size(LL1));
CPI_surface = zeros(size(LL1));


q2_max =0.4;
    q3_max = 0.3;

alpha_range = linspace(-20*pi/180, 20*pi/180,5);
theta_range = linspace(-30*pi/180, 90*pi/180,5);
[Alpha1, Theta] = meshgrid(alpha_range, theta_range);
W_v = zeros(size(Alpha1));
 
W = diag([1/q2_max, 1/q3_max]);
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

objective = @(L) -computeGlobalManipulability(L,Alpha1,Theta,W);

constraints = [];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [10,15,35,10]; 
ub = [25,25,55,19]; 
for i = 1:length(L1_range)
    for j = 1:length(L4_range)
        L = [L1_range(i), L0(2), L0(3), L4_range(j)]; 
        
        global_manipulability_surfacea(j, i) = computeGlobalManipulability(L, Alpha1, Theta, W);
    end
end


m = mean(global_manipulability_surfacea);
m = mean(m);


figure;
surf(LL1, LL4, global_manipulability_surfacea);
xlabel('L1');
ylabel('L4');
zlabel('Global Manipulability_a');
title('Surface Plot of Global Manipulability');
grid on;

[L_opt, fval] = fmincon(objective, L0, A, b, Aeq, beq, lb, ub, constraints, options);

disp(['Optimized Link Lengths: L1 = ', num2str(L_opt(1)), ', L2 = ', num2str(L_opt(2)),  ', L3 = ',num2str(L_opt(3)), ', L4 = ', num2str(L_opt(4))]);

function [global_manipulabilitya] = computeGlobalManipulability(L,Alpha1, Theta,W)
    global_manipulability_a = 0;global_manipulability_t = 0;
    num_samples = 0;
    n=0;
    for i = 1:length(Alpha1)
        for j = 1:length(Theta)
            alpha11 = Alpha1(j, i);
            theta = Theta(j, i);
            
           [q2, q3] = inverse_T(L,alpha11, theta);  % Note: x and q1 are not used in inverse kinematics here

            J = jacobian_T(L,q2, q3);
  
            [M,Me]=inertia_matrix_exo(L,q2,q3);
           B1=M*inv(J)+J'*Me;
              B2=M*inv(J)*inv(Me)+J';
            
            manipulability1 = sqrt(det(inv(B1'*W'*W*B1)));
            manipulability2 = sqrt(det(inv(B2'*W'*W*B2)));
            
            W_a(j, i) = manipulability1;
             W_t(j, i) = manipulability2;
            
           
            global_manipulability_a = global_manipulability_a + manipulability1;
             global_manipulability_t = global_manipulability_t + manipulability2;
            num_samples = num_samples + 1;
        end
    end
    
     
      global_manipulabilitya = global_manipulability_a / num_samples;
     global_manipulabilityt = global_manipulability_t / num_samples;
end

