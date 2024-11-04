function [Us, Kinetic, Potential, Times] = propage(M, K, interpU0, interpU1, dt, niter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propage :
% Propage les conditions initiales à partir du propagateur discret saute-mouton.
%
% SYNOPSIS [Us, Kinetic, Potential, Times] = propage(M, K, interpU0, interpU1, dt, niter)
%
% INPUT  * M, K : matrice de masse et de rigidité
%        * interpU0, interpU1 : interpolée en espace des conditions initiales
%        * dt : pas de temps du schéma.
%        * niter : nombre d'itérations
%
% OUTPUT * Us : saolution discrète à tous les pas de temps
%               avec les conditions initiales (matrice Nbpt x niter + 2)
%        * Kinetic, Potential : énergies cinétique et potentielle (vecteur niter)
%        * Times : vecteur temps (vecteur niter)
%
% NOTE (1) la matrice M^{-1} K est calculée et stockée avant la boucle en
% temps (peu efficace).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calcul et stockage de la matrice M^{-1} K.
% A COMPLETER

% Allocation.
Nbpts = length(interpU0);
Us = zeros(Nbpts, niter + 2);
Kinetic = zeros(niter, 1);
Potential = zeros(niter, 1);
Times = zeros(niter, 1);

% Conditions initiales.
Us(:, 1) = interpU0;
Us(:, 2) = 

% Iteriations.
for i = 1:niter
    
    % Calcul des énergies.
    % A COMPLETER
    
    % Calcul de la solution à l'itération suivante.
    % A COMPLETER
    
    Times(i) = i * dt;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021