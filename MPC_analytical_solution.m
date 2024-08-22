%% This function finds the matrices required for computing the optimal control action.
function [G_u, Matrix_F, Matrix_G] = MPC_analytical_solution(A, B, Np, Nc)

% extracting model order
% order of the system with respect to x (noise-free output)
% A(z-1) = 1 + a_1 z^-1 + a_2 z^-2 + ... + a_na z^-na
na  = numel(A) - 1;
% order of the system with respect to u (input) 
% B(z-1) = b_0 + b_1 z^-1 + b_2 z^-2 + ... + b_nb z^-nb
nb  = numel(B) - 1;

% making Ahat: Ahat(z-1) = A(z-1).Delta ,where Delta = 1 - z^-1
Delta = [1 -1];
Ahat  = conv(A, Delta);

% Making E(z-1), F(z-1), and G
E = zeros(Np+1);
F = zeros(Np+na);               % zeros(Np+1,na+1);
G = zeros(Np+1,Np+nb);

% if j=1 --> E(z-1) = e_{1,0} = a_0 = 1
E(1,1) = 1;
% if j=1 --> F(z-1) = f_{1,0} + f_{1,1} z^-1 + ... + f_{1,na} z^-na = -ahat_1 - ahat_2 z^-1 - ... - ahat_na z^-na
% notice: ahat_0 = 1
F(1,1:numel(Ahat)-1) = -Ahat(2:end);
% G_j(z-1) = E_j(z-1) B(z-1) --> G_1(z-1) = B(z-1) = [b_0, b_1, b_2, ..., b_nb]
G(1,1:numel(B)) = B;

% Recursion of the Diophantine equation: DOI: https://doi.org/10.1007/978-0-85729-398-5_2
for j = 1:Np
    % Recursion of F: f_{j+1,i} = f_{j,i+1} - f_{j,0} ahat_i+1; between i:0 to na-1
    % notice: numel(Ahat)-1 = numel(A)
    for i = 1:numel(Ahat)-1
        F(j+1,i) = F(j,i+1) - F(j,1)*Ahat(i+1);
    end

    % Recursion of E: E_j+1(z-1) = E_j(z-1) + e_{j+1,j} Z^-1; where e_{j+1,j} = f_{j,0}
    E(j+1,:)   = E(j,:);
    E(j+1,j+1) = F(j,1);
    
    % Recursion of G: G_j(z-1) = E_j(z-1) B(z-1) 
    G(j+1,1:numel(conv(B,E(j+1,:)))) = conv(B,E(j+1,:));    % G(j+1,1:j+nb) = conv(B,E(j+1,j+nb-1));
end

% Making matrices F and G from polynomial F and G
Matrix_F = F(1:Np,1:numel(Ahat)-1);  %F = F(1:Np,1:numel(Ahat)-1);

Matrix_G       = zeros(Np);
G_u            = zeros(Np,1);
for j = 2:Np+1
    % Matrix_G consists of polynomial coefficients G in a reverse order
    % G = [g0 0 ... 0;g1 g0 ...0;...;gNp-1 gNp1-2 ... g0]
    Matrix_G(j-1,:) = [G(j-1,j-1:-1:1),zeros(1,Np-(j-1))];
    % extracting G_u for free responce --> extracting first j-1 column of Matrix G
    G_u(j-1,:)      = G(j-1,j:j+nb-1);
end

% only send out Nu column of Matrix G for the sake of lower computational burden.
Matrix_G = Matrix_G(:,1:Nc);
end