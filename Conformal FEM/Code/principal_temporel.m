% =====================================================
% principal_temporel;
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour :
%
% 1) l'equation suivante des ondes en régime temporel, avec conditions de
% Dirichlet homogene
% | d^2_{tt} u - div(\sigma \grad u)= f,   dans \Omega=\Omega_1 U \Omega_2
% |         u = 0,   sur le bord
%
% avec
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
%
% =====================================================

%% Reading mesh file and assembling matrices.

meshFilePath = "C:\cygwin64\home\hugon\ANN201\Routines\GeomRect\0.04\geomRect_0.04.msh" ;
massType = 'Cond'; % 'Cond'; 'Exacte' %  

[Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri]=lecture_msh(meshFilePath);

if strcmp(massType, 'Exacte')
    [M, K] = assembleMK(Coorneu, Refneu, Numtri, Reftri) ;
elseif strcmp(massType, 'Cond')
    [M, K] = assembleMCondK(Coorneu, Refneu, Numtri, Reftri) ;
end
%% Computing CLF condition and effective time step.

% Calcul de la CFL
vp = eigs(K,M,1) ;
cfl = 2 / sqrt(vp) ;
 
% Calcul du pas de temps.
cfl_factor = 0.95;
dt = cfl_factor * cfl;


%% Computing initial conditions.

interpU0 = zeros(Nbpt, 1);
interpintU1 = zeros(Nbpt, 1);
for i = 1:Nbpt
    interpU0(i) = exp(-50*((Coorneu(i,1)-3)^2+(Coorneu(i,2)-1)^2)) ;
    interpintU1(i) = 0 ; % Petit détour pédagogique
end
interpU1 = dt*interpintU1 + (eye(Nbpt,'like',K)-M\K*dt^2 /2)*interpU0 ;

%% Propagating.

Tmax = 3.0;
niter = floor(Tmax / dt) + 1;
solverType = 'Cond' ; % 'Cholesky' ; 'Inv' ; 'Cond'
tic
if strcmp(solverType, 'Inv')
    [Us, Kinetic, Potential, Times] = propage(M, K, interpU0, interpU1, dt, niter);
elseif strcmp(solverType, 'Cholesky')
    [Us, Kinetic, Potential, Times] = propage_cholesky(M, K, interpU0, interpU1, dt, niter);
elseif strcmp(solverType, 'Cond')
    [Us, Kinetic, Potential, Times] = propage_cond(M, K, interpU0, interpU1, dt, niter);
end
toc
%% Plots of energy and solution.

figure
hold on;
plot(Times, Potential)
plot(Times, Kinetic)
plot(Times, Kinetic + Potential)
xlim([min(Times) max(Times)])
xlabel('Time')
ylabel('Energy')
legend('P', 'K', 'E', 'Location', 'SouthEast')

affiche(Us(:,end), Numtri, Coorneu);
afficheSigma(Numtri, Reftri, Coorneu);


%% Interpolation at point.

CoordsInterpPnts = [3, 1; 4, 1];
NbInterpPnts = size(CoordsInterpPnts, 1);

% Calcul de la matrice des coefficients d'interpolation.
interpolationOp = interpTriP1(Coorneu, Numtri, CoordsInterpPnts);
s1 = zeros(niter+2,1) ;
s2 = zeros(niter+2,1) ;
for i = 1:niter+2
    s1(i) =  interpolationOp(1,:)*Us(:, i);
    s2(i) =  interpolationOp(2,:)*Us(:, i);
end
figure;
plot(Times,s1(3:end))
hold on
plot(Times,s2(3:end))
xlabel('Time')
ylabel('valeurs en ces points')
legend('(3,1)', '(4,1)', 'Location', 'SouthEast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021
