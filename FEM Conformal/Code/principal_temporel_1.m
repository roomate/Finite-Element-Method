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

meshFilePath = "C:\cygwin64\home\hugon\ANN201\Routines\GeomRect\0.04\geomRect_0.04.msh";
massType ='Cond'; %'Exacte'%;    

[Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri]=lecture_msh(meshFilePath);

if strcmp(massType, 'Exacte')
    [M, K] = assembleMK(Coorneu, Refneu, Numtri, Reftri);
elseif strcmp(massType, 'Cond')
    [M, K] = assembleMCondK(Coorneu, Refneu, Numtri, Reftri);
end

%% Computing CLF condition and effective time step.

% Calcul de la CFL
lambda = eigs(K,M,1);
cfl = 2/sqrt(abs(lambda));
 
% Calcul du pas de temps.
cfl_factor = 0.95;
dt = cfl_factor * cfl;


%% Computing initial conditions.

interpU0 = zeros(Nbpt, 1);
interpU1int = zeros(Nbpt, 1);
for i = 1:Nbpt
    S = Coorneu(i,:);
    x = S(1);
    y = S(2);
    interpU0(i) = exp(-50*((x-3)^2 + (y-1)^2));
end


%% Propagating.

Tmax = 3.0;
niter = floor(Tmax / dt) + 1;
solverType =  'Cond';%'Inv'% ; % 'Cholesky'%;

if strcmp(solverType, 'Inv')
    [Us, Kinetic, Potential, Times] = propage(M, K, interpU0, interpU1, dt, niter);
elseif strcmp(solverType, 'Cholesky')
    [Us, Kinetic, Potential, Times] = propage_cholesky(M, K, interpU0, interpU1, dt, niter);
elseif strcmp(solverType, 'Cond')
    [Us, Kinetic, Potential, Times] = propage_cond(M, K, interpU0, interpU1, dt, niter);
end



%% Plots of energy and solution.

figure
hold on
hold on;
plot(Times, Potential)
plot(Times, Kinetic)
plot(Times, Kinetic + Potential)
xlim([min(Times) max(Times)])
xlabel('Time')
ylabel('Energy')
legend('P', 'K', 'E', 'Location', 'SouthEast')

% view(2);
% v = VideoWriter('Film_Exo_3_Q5.avi','Uncompressed AVI');
% v.FrameRate = 2;
% open(v);
% for i=1:5:niter+2
%     affiche(Us(:, i), Numtri, Coorneu,"Propagation d'une onde avec sigma_1 = 4 et sigma_2 = 9");
%     saveas(gcf,['Image_', num2str(i)],'png');
%     A = imread(['Image_', num2str(i)],'png');
%     writeVideo(v,A);
%     hold on;
%     pause(0.1);
% end
% close(v);
% %afficheSigma(Numtri, Reftri, Coorneu);
% winopen('Film_Exo_3_Q5.avi');

%% Interpolation at point.

CoordsInterpPnts = [3, 1; 4, 1];
NbInterpPnts = size(CoordsInterpPnts, 1);

% Calcul de la matrice des coéfficients d'interpolation.
interpolationOp = interpTriP1(Coorneu, Numtri, CoordsInterpPnts);

InterpU = zeros(NbInterpPnts,niter + 2);

for  i = 1 : niter+2
    InterpU(:,i) = interpolationOp * Us(:,i);
end
% figure
% plot(Times,InterpU(1,:));
% hold on;
% plot(Times, InterpU(2,:));
% xlabel("Time");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021
