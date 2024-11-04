[Nbpt1, Nbtri1, Coorneu1, Refneu1, Numtri1, Reftri1] = lecture_msh(meshFilePath(7));
interpU0 = zeros(Nbpt1,1);
interpU1 = zeros(Nbpt1,1);
for i=1:Nbpt
    S = Coorneu1(i,:);
    x = S(1);
    y = S(2);
    interpU0(i) = exp(-50*((x-3)^2 + (y-1)^2));
end

[MCond, K] = assembleMCondK(Coorneu1, Refneu1, Numtri1, Reftri1);
dt = 0.01; %On fait attention que la condition CFL soit vérifié
niter = 3/dt;

[Us, Kinetic, Potential, Times] = propage_cond(MCond, K, interpU0, interpU1, dt, niter);

plot(Times,Kinetic);
plot(Times,Potential);
plot(Times,Kinetic + Potential);