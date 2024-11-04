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


% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomCarre.msh';
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
             
   Mel=matM_elem(S1, S2, S3);
    
   Kel=matK_elem(S1, S2, S3,'tenseur');
           
   
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
FF = f(Coorneu(:,1),Coorneu(:,2));
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

AA = MM+KK;
AA0 = PP*AA*transpose(PP);
LL0 = PP*LL;

% inversion
% ----------
UU0 = AA0\LL0;

% Expression de la solution dans toute la base
% ———————
UU = transpose(PP)*UU0;

% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));
UU_exact = sin(pi*Coorneu(:,1)).*sin(2*pi*Coorneu(:,2));
%affiche(UU_exact, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage))
validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
UU_exact = sin(pi*Coorneu(:,1)).*sin(2*pi*Coorneu(:,2));
err_L2 = 0;
for i=1:size(MM,1)
    for j=1:size(MM,1)
        err_L2 = err_L2 + (UU_exact(i) - UU(i))*MM(i,j)*(UU_exact(j) - UU(j));
    end
end

% Calcul de l erreur H1
err_H1 = 0;
for i=1:size(KK,1)
    for j=1:size(KK,1)
        err_H1 = err_H1 + (UU_exact(i) - UU(i))*KK(i,j)*(UU_exact(j) - UU(j));
    end
end
err_H1 = sqrt(err_H1);
% attention de bien changer le terme source (dans FF)
err_L2 = sqrt(err_L2); 
norm_u_L2 = pi/sqrt(2);
display(err_L2);
display(err_H1);
H = [1 log(1);1 log(2);1 log(4);1 log(8); 1 log(16)];
semi_norm_u = sqrt(3/2)*pi;
h = [1,2,4,8,16];
L2_err_L2 = [0.3414, 0.1054, 0.0295, 0.0076, 0.0019]/norm_u_L2;
L2_err_H1 = [2.7311, 1.0934, 0.3139, 0.0818, 0.0208]/semi_norm_u;
B_L2 = H\transpose(log(L2_err_L2));
%plot(-log(h),log(L2_err_L2));
hold on;
%plot(-log(h), H*B_L2);
legend({'Solution','Régression linéaire : pente = 1.8772'},'Location','northwest')
xlabel('-log(h)');
ylabel('erreur relative de la norme L2');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

