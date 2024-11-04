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

%Paramètre eta
eta = 0.01; %Arbitraire

%Each edge has its own list
Noeud_bord1 = find(Refneu==1);
Noeud_bord2 = find(Refneu==2);
Noeud_bord3 = find(Refneu==3);
Noeud_bord4 = find(Refneu==4);
Noeud_coin = find(Refneu==5);
Noeud_I = find(Refneu==0);

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  % Coordonnees des sommets du triangles
  S1=[Coorneu(Numtri(l,1),:)];
  S2=[Coorneu(Numtri(l,2),:)];
  S3=[Coorneu(Numtri(l,3),:)];
  % calcul des matrices elementaires du triangle l 
             
   Mel=matM_elem_TP2(S1, S2, S3);
    
   Kel=matK_elem_TP2(S1, S2, S3,0);
   
  
   for i=1:3
      for j=1:3
          KK(Numtri(l,i),Numtri(l,j)) = KK(Numtri(l,i),Numtri(l,j)) + Kel(i,j);
          MM(Numtri(l,i),Numtri(l,j)) = MM(Numtri(l,i),Numtri(l,j)) + Mel(i,j);
      end
      KK_2(Numtri(l,i)) = KK_2(Numtri(l,i)) + Kel_2(i,1);
   end
end % for l

% Calcul du second membre L
% -------------------------
LL = -KK_2;

% Projection sur l espace V_p
% ———————————————————
% matrice de projection 
PP = sparse(1+N_noeud_I+2*N_noeud_E,size(Refneu,1));
AA = eta*MM + KK; %On n'oublie pas le paramètre eta
k=2; % La première ligne est réservée pour les coins

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

AAp = PP*AA*transpose(PP);
LLp = PP*LL;

% inversion
% ----------
UUp = AAp\LLp;

% Expression de la solution dans toute la base
% ———————
UU = transpose(PP)*UUp;

% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Periodique - %s', nom_maillage));
validation = 'oui';
% ----------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

