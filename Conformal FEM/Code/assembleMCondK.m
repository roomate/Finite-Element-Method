function [MCond, K] = assembleMCondK(Coorneu, Refneu, Numtri, Reftri)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembleMCondK :
% assemble les matrices de masse condensée et de raideur globales en P1 lagrange.
%
% SYNOPSIS [MCond, K] = assembleMCondK(Coorneu, Numtri, Reftri)
%          
% INPUT  * Coorneu : coordonnees (x, y) des sommets (matrice reelle Nbpt x 2)
%        * Refneu : reference des sommets (vecteur entier Nbpt x 1)
%        * Numtri : liste de triangles 
%                   (3 numeros de sommets - matrice entiere Nbtri x 3)
%        * Reftri : reference des triangles (matrice entiere Nbtri x 1)
%
% OUTPUT * M matrice de masse globale (vecteur Nbpt)
%        * K matrice de raideur globale (matrice NbptxNbpt)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nbpt = size(Coorneu, 1);
Nbtri = size(Numtri, 1);

% Declarations des matrices EF.
K = sparse(Nbpt,Nbpt);
MCondDiag = zeros(Nbpt, 1);

% Boucle d'assemblage sur les triangles.
for l=1:Nbtri
    
    S1 = Coorneu(Numtri(l,1), :);
    S2 = Coorneu(Numtri(l,2), :);
    S3 = Coorneu(Numtri(l,3), :);
    
    % Assemblage de la matrice de rigidité.
   Kel=matK_elem(S1,S2,S3,Reftri(l));
   Aire_Tl = 1/2 * abs((S2(1) - S3(1))*(S3(2) - S1(2)) - (S3(1) - S1(1))*(S2(2) - S3(2)));
     for i=1:3
      I=Numtri(l,i);
      for j=1:3
          J=Numtri(l,j);
          K(I,J)=K(I,J)+Kel(i,j);
      end
        
    
    % Assemblage de la diagonale de la matrice de masse.
      MCondDiag(I) = MCondDiag(I) + 1/3*Aire_Tl; %Ajout du triangle TL dont S1 est le sommet
     end
end % for l

% Application de la pseudo élimination.
for l=1:Nbpt
    if Refneu(l) == 1 || Refneu(l) == 2
        K(:,l) = zeros(Nbpt,1);
        K(l,:) = zeros(1,Nbpt);
        K(l,l) = 1 ;
    end
end

% Transformation de la diagonale en une matrice sparse.
MCond = spdiags(MCondDiag, 0, Nbpt, Nbpt);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021