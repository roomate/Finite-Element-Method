% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions periodiques
% sur le maillage nom_maillage.msh
%
% | -div(A grad u ) u + u= f,   dans \Omega
% |         u periodique,   sur le bord
%
% =====================================================
%Code pour calculer les correcteurs

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
KK_2 = zeros(Nbpt,1); % Matrice de rigidité pour le second membre
MM = sparse(Nbpt,Nbpt); % matrice de rigidite

%Nombre de noeud à l'intérieur
N_noeud_I = size(find(Refneu==0),1);
%Il y a le même nombre de noeud pour tout les autres côtés
N_noeud_E = (size(Refneu,1) - N_noeud_I)/4 - 1;

x=Coorneu(:,1);
y=Coorneu(:,2);

%Paramètre eta
eta = 0.01; %Arbitraire

%Each edge has its own list
Noeud_bord1 = find(Refneu==1);
Noeud_bord2 = find(Refneu==2);
Noeud_bord3 = find(Refneu==3);
Noeud_bord4 = find(Refneu==4);
Noeud_coin = find(Refneu==5);
Noeud_I = find(Refneu==0);
k=2; % La première ligne est réservée pour les coins

% matrice de projection 
PP = sparse(1+N_noeud_I+2*N_noeud_E,size(Refneu,1));

for i=1:size(Refneu,1)
   if Refneu(i)==0
       PP(k,i) = 1;
       k = k + 1;
   elseif Refneu(i)==1
       for j=1:N_noeud_E
           if abs(Coorneu(i,1) - Coorneu(Noeud_bord3(j),1)) < 1e-10
               PP(k,i) = 1;
               PP(k,Noeud_bord3(j)) = 1;
               k = k + 1;
           end
       end
   elseif Refneu(i)==2
       for j=1:N_noeud_E
           if abs(Coorneu(i,2) - Coorneu(Noeud_bord4(j),2)) < 1e-10
               PP(k,i)=1;
               PP(k,Noeud_bord4(j))=1;
               k = k + 1;
           end
       end 
   elseif Refneu(i)==5
       for j=1:size(Noeud_coin,1)
           PP(1,Noeud_coin(j))=1;
       end
   end
end


% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  % Coordonnees des sommets du triangles
  S1=[Coorneu(Numtri(l,1),:)];
  S2=[Coorneu(Numtri(l,2),:)];
  S3=[Coorneu(Numtri(l,3),:)];
  % calcul des matrices elementaires du triangle l 
             
   Mel=matM_elem_TP2(S1, S2, S3);
    
   Kel=matK_elem_TP2(S1, S2, S3, 1);
   
   Kel_2 = matK_elem_TP2(S1,S2,S3, 1); 
  
   for i=1:3
      for j=1:3
          KK(Numtri(l,i),Numtri(l,j)) = KK(Numtri(l,i),Numtri(l,j)) + Kel(i,j);
          MM(Numtri(l,i),Numtri(l,j)) = MM(Numtri(l,i),Numtri(l,j)) + Mel(i,j);
      end
   end
end % for l

% Calcul du second membre L
% -------------------------


LL1 = -transpose(Coorneu(:,1))*KK;
LL2 = -transpose(Coorneu(:,2))*KK;

LL1 = transpose(LL1); %Pour avoir un vecteur colonne
LL2 = transpose(LL2); %Pour avoir un vecteur colonne

% Projection sur l espace V_p
% ———————————————————

AA = eta*MM + KK; %On n'oublie pas le paramètre eta

AAp = PP*AA*transpose(PP);
LLp1 = PP*LL1;
LLp2 = PP*LL2;
% inversion
% ----------
W_xp = AAp\LLp1;
W_yp = AAp\LLp2;

% Expression de la solution dans toute la base
% ———————
W_x = transpose(PP)*W_xp;
W_y = transpose(PP)*W_yp;
% visualisation
% -------------
%affiche(W_y, Numtri, Coorneu, sprintf('Periodique - %s', nom_maillage));
validation = 'oui';
% ----------

% %Calcul du tenseur généralisé A_eff
A_eff = A_effi(W_x,W_y,x,y,KK);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

