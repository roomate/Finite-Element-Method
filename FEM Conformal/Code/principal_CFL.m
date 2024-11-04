% =====================================================
% principal_CFL;
%
% une routine permettant d'évaluer la dépendance de la condition CFL
% en fonction du pas de maillage
%
% =====================================================

%% Definition of mesh files and associated mesh steps.

meshFilePath = [ 
    "C:\cygwin64\home\hugon\ANN201\Routines\GeomRect\0.04\geomRect_0.04.msh",
    "C:\cygwin64\home\hugon\ANN201\Routines\GeomRect\0.06\geomRect_0.06.msh",
    "C:\cygwin64\home\hugon\ANN201\Routines\GeomRect\0.08\geomRect_0.08.msh", 
    "C:\cygwin64\home\hugon\ANN201\Routines\GeomRect\0.10\geomRect_0.10.msh", 
    "C:\cygwin64\home\hugon\ANN201\Routines\GeomRect\0.12\geomRect_0.12.msh",
    "C:\cygwin64\home\hugon\ANN201\Routines\GeomRect\0.14\geomRect_0.14.msh",
    "C:\cygwin64\home\hugon\ANN201\Routines\GeomRect\0.16\geomRect_0.16.msh",
    "C:\cygwin64\home\hugon\ANN201\Routines\GeomRect\0.18\geomRect_0.18.msh",
    "C:\cygwin64\home\hugon\ANN201\Routines\GeomRect\0.20\geomRect_0.20.msh", 
    ];

steps = 0.04:0.02:0.2;


%% Computing CFL for every meshes.

nbMesh = length(steps);
if nbMesh ~= length(meshFilePath)
    error('Number of mesh files and number of steps differ.');
end

% Allocation.
cfl = zeros(nbMesh, 1);
X = ones(nbMesh, 2);

% Boucle sur les maillages.
for i=1:nbMesh
    X(i,2) = steps(i);
    [Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri] = lecture_msh(meshFilePath(i));
    
    % Assemblage de M et K et calcul de la CFL.
    [MCond, K] = assembleMCondK(Coorneu, Refneu, Numtri, Reftri);
    lambda = eigs(K,MCond,1);
    cfl(i) = 2/sqrt(abs(lambda));
end


%% Affichage de la cfl en fonction de h.

plot(steps,cfl);
hold on;
B = (X'*X)\X'*cfl;
disp(B);
plot(steps,X*B);
legend('Courbe simulé','Régression linéaire');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021
