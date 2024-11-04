% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Neumann sur le maillage nom_maillage.msh
%
% | -Delta  u + u= f,   dans \Omega
% |         du/dn = 0,   sur le bord
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
MM = sparse(Nbpt,Nbpt); % matrice de masse
LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  % Coordonnees des sommets du triangles
  % A COMPLETER
  S1=[Coorneu(Numtri(l,1),1), Coorneu(Numtri(l,1),2)];
  S2=[Coorneu(Numtri(l,2),1), Coorneu(Numtri(l,2),2)];
  S3=[Coorneu(Numtri(l,3),1), Coorneu(Numtri(l,3),2)];
  % calcul des matrices elementaires du triangle l 
  
   Kel=matK_elem(S1, S2, S3,'tenseur');
           
   Mel=matM_elem(S1, S2, S3);
   
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

% inversion
% ----------
UU = (MM+KK)\LL;

% visualisation
% -------------

affiche(UU, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));

validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
UU_exact = cos(pi*Coorneu(:,1)).*cos(2*pi*Coorneu(:,2));
xlabel("x");
ylabel("y");
% Calcul de l erreur L2
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
h = [1,2,4,8, 16];
H = [1, log(1);1 log(2);1 log(4);1 log(8);1 log(16)];
norm_u = pi/sqrt(2);
semi_norm_u = sqrt(3/2)*pi;
L1_err_L2 = [0.3998, 0.1077, 0.0297, 0.0076, 0.0019]/norm_u;
L1_err_H1 = [1.7123, 0.6899, 0.1934, 0.0501, 0.0127]/semi_norm_u;
B_L2 = H\transpose(log(L1_err_H1));
L1_err_L2_t = [0.7830, 0.1589, 0.0414, 0.0105]/norm_u;
L1_err_H1_t = [2.7873, 1.1025, 0.3342, 0.0885]/semi_norm_u;
display(err_L2);
display(err_H1);
plot(-log(h),log(L1_err_H1));
hold on;
plot(-log(h), H*B_L2);
legend({'Solution','Régression linéaire : pente = 1.7933'},'Location','northwest')
xlabel('-log(h)');
ylabel('erreur relative de la semi-norme H1');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%