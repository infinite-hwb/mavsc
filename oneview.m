function [C1, C2, Y1, Y2, A] = centroid_MGSSC_1view(n, K, C1, C2, C_centroid, Y1, Y2, alpha, lambda1, mu)

% updating A

inv_A = inv(2*K+(mu(1)+mu(2))*eye(n));
A = inv_A * (2*K+mu(1)*(C1-Y1/mu(1))+mu(2)*(C2-Y2/mu(2)));


% updating C2 and C3
C1_new = soft_thresh(A+Y1/mu(1),lambda1/mu(1));
C1_new = C1_new - diag(diag(C1_new));
C2_new = 1/(2*alpha+mu(2))*(2*alpha*C_centroid+mu(2)*A+Y2);

%update variables
C1 = C1_new;
C2 = C2_new;

% updating Lagrange multipliers

Y1 = Y1 + mu(1) * (A - C1 + diag(diag(C1)));
Y2 = Y2 + mu(2) * (A - C2);
