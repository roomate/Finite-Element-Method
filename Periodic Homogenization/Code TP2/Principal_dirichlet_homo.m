% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Dirichlet sur le maillage nom_maillage.msh
%
% | -\Delta u + u= f,   dans \Omega
% |         u = 0,   sur le bord
%
% =====================================================
%Code pour calculer le problème homogénéisé

% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomCarre_per.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);


% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre


% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  % Coordonnees des sommets du triangles
  S1=[Coorneu(Numtri(l,1),:)];
  S2=[Coorneu(Numtri(l,2),:)];
  S3=[Coorneu(Numtri(l,3),:)];
  % calcul des matrices elementaires du triangle l 
             
   Mel=matM_elem_TP2(S1, S2, S3);
   
   Kel=matK_elem_homo(S1, S2, S3,A_eff);
           
   
   for i=1:3
      for j=1:3
          KK(Numtri(l,i),Numtri(l,j)) = KK(Numtri(l,i),Numtri(l,j)) + Kel(i,j);
          MM(Numtri(l,i),Numtri(l,j)) = MM(Numtri(l,i),Numtri(l,j)) + Mel(i,j);
      end
   end
end % for l


% Calcul du second membre L
% -------------------------
% utiliser la routine f.m
FF = f_TP2(Coorneu(:,1),Coorneu(:,2));
LL = MM*FF;

%Nombre de noeud à l'intérieur
Noeud_I = find(Refneu==0);
N_noeud_I = size(Noeud_I,1);

% Projection sur l espace V_0
% ———————————————————
% matrice de projection 
PP = sparse(N_noeud_I,size(Refneu,1));
k=1;
for i=1:N_noeud_I
    PP(k,Noeud_I(i)) = 1;
    k = k + 1;
end 

KK0 = PP*KK*transpose(PP);
LL0 = PP*LL;

% inversion
% ----------
UU0 = KK0\LL0;

% Expression de la solution dans toute la base
% ———————
UU_exact = transpose(PP)*UU0; %Solution du poblème homogénéisé

% visualisation
% -------------
%affiche_3D(UU_exact, Numtri, Coorneu, sprintf('Périodique - %s', nom_maillage));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

