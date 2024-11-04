function [M, K] = assembleMK(Coorneu, Refneu, Numtri, Reftri)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembleMK :
% assemble les matrices de masse et de raideur globales en P1 lagrange.
%
% SYNOPSIS [M, K] = assembleMK(Coorneu, Numtri, Reftri)
%
% INPUT  * Coorneu : coordonnees (x, y) des sommets (matrice reelle Nbpt x 2)
%        * Refneu : reference des sommets (vecteur entier Nbpt x 1)
%        * Numtri : liste de triangles
%                   (3 numeros de sommets - matrice entiere Nbtri x 3)
%        * Reftri : reference des triangles (matrice entiere Nbtri x 1)
%
% OUTPUT * M matrice de masse globale (matrice NbptxNbpt)
%        * K matrice de raideur globale (matrice NbptxNbpt)
%
% NOTE (1) les matrices de masse et de raideur sont définies au format
% sparse.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nbpt = size(Coorneu, 1);
Nbtri = size(Numtri, 1);

% Declarations des matrices EF.
KK = sparse(Nbpt,Nbpt);
MM = sparse(Nbpt,Nbpt);

%Déclaration des matrices EF après pseudo-élimination
K = sparse(Nbpt,Nbpt);
M = sparse(Nbpt,Nbpt);

% Boucle d'assemblage sur les triangles.
for l=1:Nbtri
    
    S1 = Coorneu(Numtri(l,1), :);
    S2 = Coorneu(Numtri(l,2), :);
    S3 = Coorneu(Numtri(l,3), :);
    
    % Calcul des matrices elementaires du triangle l.
   Kel=matK_elem(S1,S2,S3,Reftri(l));
   Mel=matM_elem(S1,S2,S3);
   
    % Assemblage des matrices globales.
     for i=1:3
      I=Numtri(l,i);
      for j=1:3
          J=Numtri(l,j);
          MM(I,J)=MM(I,J)+Mel(i,j);
          KK(I,J)=KK(I,J)+Kel(i,j);
      end
     end   
end

% Application de la pseudo élimination.
for l=1:Nbpt
    if Refneu(l)==1 || Refneu(l)==2
        for i=1:Npt
            if i==l
                K(i,i) = 1;
                M(i,i) = 1;
            else
                M(i,l) = 0;
                M(l,i) = 0;
                K(l,i) = 0;
                K(i,l) = 0;
            end
        end
    end     
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021