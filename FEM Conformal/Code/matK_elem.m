function Kel = matK_elem(S1, S2, S3, Reftri)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matK_elem :
% calcul la matrices de raideur elementaire en P1 lagrange.
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%       * Reftri : reference du triangle.
%
% OUTPUT * Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul utilise une formule de quadrature de Gauss-Legendre.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture.
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% Points et poids de quadrature.
S_hat(:,1) = [1/6; 1/6];
S_hat(:,2) = [2/3; 1/6];
S_hat(:,3) = [1/6, 2/3];
w0 = 1/6;

% Gradients des fonctions de base sur le triangle de reference.
Delta_w(:,1)=[-1;-1];
Delta_w(:,2)=[1;0];
Delta_w(:,3)=[0;1];

% Transformation géométrique associée au triangle courant.
B_l = [x2-x1 x3-x1;y2-y1 y3-y1];
S_l = [x1;y1];
det_B = abs(det(B_l));
S=zeros(2,3);
for i=1:3
    S(:,i)= B_l*S_hat(:,i) +S_l;
end
B_inv = B_l^(-1);


Kel = zeros(3,3);
for i=1:3
    for j=1:3
        for k=1:3
            if Reftri==1
                Kel(i,j)=Kel(i,j) + sigma_1(S(:,k))*((B_inv')*Delta_w(:,i))'*((B_inv')*Delta_w(:,j))*det_B*w0;
            else 
                Kel(i,j)=Kel(i,j) + sigma_2(S(:,k))*((B_inv')*Delta_w(:,i))'*((B_inv')*Delta_w(:,j))*det_B*w0;
            end
        end
    end
end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021
