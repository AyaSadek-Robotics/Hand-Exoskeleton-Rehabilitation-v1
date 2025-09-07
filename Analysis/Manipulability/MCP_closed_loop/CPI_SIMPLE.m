
L0 = [12,20,45,15];
q2_max = 20*pi/180;
    q3_max = 60*pi/180;
tq2_max = 0.4;
    tq3_max = 0.4;
% we define the ranges for joint angles
alpha_range = linspace(-20*pi/180, 20*pi/180, 5);
theta_range = linspace(-30*pi/180, 90*pi/180,5);
[Alpha1, Theta] = meshgrid(alpha_range, theta_range);
W_v = zeros(size(Alpha1));
 
   
W = diag([1/q2_max, 1/q3_max]);
Wd = diag([1/tq2_max, 1/tq3_max]);


options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

objective = @(L) -computeGlobalManipulability(L,Alpha1,Theta,W,Wd);

constraints = [];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [10,15,35,10];
ub = [25,25,55,19]; 
% Perform the optimization
[L_opt, fval] = fmincon(objective, L0, A, b, Aeq, beq, lb, ub, constraints, options);


disp(['Optimized Link Lengths: L1 = ', num2str(L_opt(1)), ', L2 = ', num2str(L_opt(2)), num2str(L_opt(3)), num2str(L_opt(4))]);

% Function to compute global manipulability
function CPI = computeGlobalManipulability(L,Alpha1, Theta,W,Wd)
    global_manipulability = 0;  global_manipulability_a = 0;  global_manipulability_t = 0;
    num_samples = 0;
    n=0;
    for i = 1:length(Alpha1)
        for j = 1:length(Theta)
            alpha11 = Alpha1(j, i);
            theta = Theta(j, i);
            
           [q2, q3] = inverse_T(L,alpha11, theta);  
            J = jacobian_T(L,q2, q3);
  
            JW_inv = inv(J)' *W';
            JW_inv1 = W*inv(J);
            
            [M,Me]=inertia_matrix_exo(L,q2,q3);
           B1=M*inv(J)+J'*Me;
              B2=M*inv(J)*inv(Me)+J';
            
            manipulability1 = sqrt(det(inv(B1'*W'*W*B1)));
            manipulability2 = sqrt(det(inv(B2'*W'*W*B2)));
            
            W_a(j, i) = manipulability1;
             W_t(j, i) = manipulability2;
            
            local_manipulability =  sqrt(det(inv(JW_inv * JW_inv1)));
             W_v(j, i) = local_manipulability;
            global_manipulability = global_manipulability + local_manipulability;
global_manipulability_a = global_manipulability_a + manipulability1;
             global_manipulability_t = global_manipulability_t + manipulability2;
            num_samples = num_samples + 1;
        end
    end
    
     
    global_manipulability = global_manipulability / num_samples;
    global_manipulability_a = global_manipulability / num_samples;
    CPI=global_manipulability+global_manipulability_t;
    
end

