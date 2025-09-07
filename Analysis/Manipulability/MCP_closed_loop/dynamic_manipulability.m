  L=[12,20,45,15];
Tq2_max = 2;
Tq3_max = 3;

alpha_range = linspace(-10*pi/180, 10*pi/180,5);
theta_range = linspace(-50*pi/180, 100*pi/180, 5);
[Alpha1, Theta] = meshgrid(alpha_range, theta_range);

W_a = zeros(size(Alpha1));
W_t = zeros(size(Alpha1));

W = diag([1/Tq2_max, 1/Tq3_max]);

for i = 1:length(alpha_range)
    for j = 1:length(theta_range)
        alpha = Alpha1(j, i);
        theta = Theta(j, i);

        [q2, q3] = inverse_T(L,alpha, theta) ; 
        J = jacobian_T(L,q2, q3);
        [M,Me]=inertia_matrix_exo(L,q2,q3);
         B1=M*inv(J)+J'*Me;
          B2=M*inv(J)*inv(Me)+J';

        manipulability1 = sqrt(det(inv(B1'*W'*W*B1)));
        manipulability2 = sqrt(det(inv(B2'*W'*W*B2)));

        W_a(j, i) = manipulability1;
         W_t(j, i) = manipulability2;
    end
end

figure;
surf( Theta,Alpha1, W_a);
xlabel('\theta');
ylabel('\alpha');
zlabel('Dynamic Manipulability Measure W_a');
title('Dynamic Manipulability Measure');
colorbar;
figure;
surf( Theta,Alpha1, W_t);
xlabel('\theta');
ylabel('\alpha');
zlabel('Dynamic Manipulability Measure W_t');
title('Dynamic Manipulability Measure');
colorbar;
