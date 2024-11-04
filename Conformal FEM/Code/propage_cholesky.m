function [Us, Kinetic, Potential, Times] = propage_cholesky(M, K, interpU0, interpU1, dt, niter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propage_cholesky :
% Propage les conditions initiales à partir du propagateur discret saute-mouton.
%
% SYNOPSIS [Us, Kinetic, Potential, Times] = propage_cholesky(M, K, interpU0, interpU1, dt, niter)
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
% NOTE (1) On utilise la décomposition de Cholesky pour résoudre les
% sytèmes du type M X = b dans les itérations en temps.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calcul de la décomposition de Cholesky de la matrice de masse.
% A COMPLETER

% Allocation.
Nbpts = length(interpU0);
Us = zeros(Nbpts, niter + 2);
Kinetic = zeros(niter, 1);
Potential = zeros(niter, 1);
Times = zeros(niter, 1);

% Conditions initiales.
Us(:, 1) = % A COMPLETER
Us(:, 2) = % A COMPLETER

% Iteriations.
for i = 1:niter
    
    % Calcul des énergies.
    % A COMPLETER
    
    % Calcul de la solution à l'itération suivante par descente remontée.
    % A COMPLETER
    
    Times(i) = i * dt;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021