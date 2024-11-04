% =====================================================
% principal_stationnaire;
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour :
%
% 1) l'equation suivante stationnaire, avec conditions de
% Dirichlet homogene
% | u - div(\sigma \grad u)= f,   dans \Omega=\Omega_1 U \Omega_2
% |         u = 0,   sur le bord
%
% avec
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
%
% =====================================================

% Donnees du probleme.
nom_maillage = 'geomRect_004.msh';
affichage = true; % false; %

% Lecture du maillage et affichage.
[Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri] = lecture_msh(nom_maillage);

% Declarations des matrices EF et vecteur second membre.
KK = sparse(Nbpt,Nbpt);
MM = sparse(Nbpt,Nbpt);
FF = zeros(Nbpt,1);

% Boucle d'assemblage sur chacun des triangles.
for l=1:Nbtri
    
    % Calcul des matrices elementaires du triangle l.
    S1 = Coorneu(Numtri(l,1), :);
    S2 = Coorneu(Numtri(l,2), :);
    S3 = Coorneu(Numtri(l,3), :);

    Kel=matK_elem(S1,S2,S3,Reftri(l));
    Mel=matM_elem(S1,S2,S3);    

    % Assemblage des matrices globales.
    for i=1:3
	I = Numtri(l,i);
	for j=1:3
	    J = Numtri(l,j);
	    MM(I,J) += Mel(i,j);
    	    KK(I,J) += Kel(i,j);
	end
    end
end

% Calcul du second membre par vectorisation.
FF = f(Coorneu);

% Matrice EF complète.
AA = MM + KK;

% Résolution du problème par inversion et sortie des résultats
%(validation et résolution numérique.
UU = AA\FF

% Visualisation de sigma et de la solution.
if affichage
    afficheSigma(Numtri, Reftri, Coorneu);
    affiche(UU, Numtri, Coorneu, sprintf('Stationnaire - %s', nom_maillage));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021

