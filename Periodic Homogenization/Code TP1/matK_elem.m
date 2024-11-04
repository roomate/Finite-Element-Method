function [Kel] = matK_elem(S1, S2, S3, tenseur)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul la matrices de raideur elementaire en P1 lagrange
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exacte (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par 
%          les coordonnees barycentriques 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = zeros(3, 2);
norm(1, :) = [y2-y3, x3-x2];
norm(2, :) = [y3-y1, x1-x3];
norm(3, :) = [y1-y2, x2-x1];



%Point de quadrature
M = [1/3  1/5, 1/5, 3/5 ; 1/3, 1/5, 3/5, 1/5];
W = [-9/32; 25/96; 25/96; 25/96];

%Matrice de passage entre triangle et triangle de référence
C = [x2-x1, x3-x1; y2-y1, y3-y1]; 
det1 = abs(det(C));
O = transpose(inv(C));
B = [x1; y1];

%Fonctions de basses locales
lambda = [-1, 1, 0; -1, 0, 1];  

% D est, au signe pres, deux fois l'aire du triangle
D = abs(((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)));
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end

% calcul de la matrice de raideur
% -------------------------------
Kel = zeros(3,3);
%Cas où  c'est l'identité est identique au cas scalaire
if strcmp(tenseur, 'scalaire')
for i=1:3
  for j=i:3
    Kel(i,j) = (norm(i,1)*norm(j,1) + norm(i,2)*norm(j,2))/(2*D);
    Kel(j,i) = Kel(i,j);
  end % j
end % i
end

if strcmp(tenseur, 'tenseur')
    for i=1:3
        for j=1:3
            for l=1:4
                x = C*M(:,l) + B;
                Kel1 = A(x(1),x(2))*O*lambda(:,i);
                Kel2 = lambda(:,j)*det1;
                Kel(i,j) = Kel(i,j) + W(l)*dot(Kel1,O*Kel2);
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
