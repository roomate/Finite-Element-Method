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
%Code pour Q1

% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomCarre_per.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% Pas de la période
eps = 0.0001;

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
    
   Kel=matK_elem_TP2(S1, S2, S3, eps);
           
   
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
N_noeud = 0;
for i=1:size(Refneu,1)
    if Refneu(i)==0
        N_noeud = N_noeud + 1;
    end
end

%Projection sur l espace V_0
%———————————————————
%matrice de projection 
PP = sparse(N_noeud,size(Refneu,1));
k=1;
for i=1:size(Refneu,1)
   if Refneu(i)==0
       PP(k,i) = 1;
       k = k + 1;
   end 
end

KK0 = PP*KK*transpose(PP);
LL0 = PP*LL;

% inversion
% ----------
UU0 = KK0\LL0;

% Expression de la solution dans toute la base
% ———————
UU = transpose(PP)*UU0;
% visualisation
% -------------
%affiche_3D(UU, Numtri, Coorneu, sprintf('cas (iv) %s', nom_maillage));

validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
dderiv_U0 = deriv_U0(UU_exact,Coorneu,Numtri,Nbtri);
err_L2 = 0;
for i=1:size(MM,1)
    for j=1:size(MM,1)
        err_L2 = err_L2 + (UU_exact(i) + eps*(dderiv_U0(i,1)*W_x(i) + dderiv_U0(i,2)*W_y(i)) - UU(i))*MM(i,j)*(UU_exact(j) + eps*(dderiv_U0(j,1)*W_x(j) + dderiv_U0(j,2)*W_y(j)) - UU(j));
    end
end                                                         
% Calcul de l erreur H1
err_H1 = 0;
for i=1:size(KK,1)
    for j=1:size(KK,1)
        err_H1 = err_H1 + (UU_exact(i) + eps*(dderiv_U0(i,1)*W_x(i) + dderiv_U0(i,2)*W_y(i)) - UU(i))*KK(i,j)*(UU_exact(j) + eps*(dderiv_U0(j,1)*W_x(j) + dderiv_U0(j,2)*W_y(j)) - UU(j));
    end
end
% attention de bien changer le terme source (dans FF)
err_L2 = sqrt(err_L2); 
norm_u_L2 = pi; %Constante de normalisation de la solution pour la norme L2
display(err_L2);
display(sqrt(err_L2*err_L2 + err_H1));
norm_u_H1 = pi*sqrt(1+2*pi*pi); %Constante de normalisation pour la norme H1
eps = [5e-1,1e-1,5*1e-2,1e-2, 1e-3,1e-4];
L2_err_L2 = [0.3771, 0.0807, 0.0459, 0.0288, 0.0389, 0.0418]/norm_u_L2;
L2_err_H1 = [12.9678, 3.4427, 2.3122, 1.5256, 0.9693, 0.8438]/norm_u_H1;
plot(-log(eps),log(L2_err_H1));
xlabel('-log(eps)');
ylabel('erreur relative de la norme H1');
end

%affiche_3D(dderiv_U0(:,1),Numtri,Coorneu,sprintf('cas (iv) %s', nom_maillage))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

