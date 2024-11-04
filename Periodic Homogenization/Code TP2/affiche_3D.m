function affiche_3D(UU,Numtri,Coorneu, titre)
Nbpt = size(Coorneu,1);
Nbtri = size(Numtri,1);
for i = 1:Nbtri
    sommet1 = Numtri(i,1);
    sommet2 = Numtri(i,2);
    sommet3 = Numtri(i,3);
    

   U(1) =  UU(sommet1);
   U(2) =  UU(sommet2);
   U(3) =  UU(sommet3);
   
   X = [Coorneu(sommet1,1) Coorneu(sommet2,1) Coorneu(sommet3,1)];
   Y = [Coorneu(sommet1,2) Coorneu(sommet2,2) Coorneu(sommet3,2)];
   Z = [U(1), U(2), U(3)];

   patch(X,Y,Z, [0, 3/5, 1]);
   view(3);
end
title(titre);
    
    
    