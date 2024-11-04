function [UU0_partial_d] = deriv_U0(UU0,Coorneu,Numtri,Nbtri)


UU0_partial_d = zeros(size(UU0,1),2);
for l=1:Nbtri
  % Coordonnees des sommets du triangle
  S1=[Coorneu(Numtri(l,1),:)];
  S2=[Coorneu(Numtri(l,2),:)];
  S3=[Coorneu(Numtri(l,3),:)];
  % calcul des matrices elementaires du triangle l 
  A = mat_A(S1,S2,S3); 
  n =  transpose(transpose(A)\[-UU0(Numtri(l,1)) + UU0(Numtri(l,2)); -UU0(Numtri(l,1)) + UU0(Numtri(l,3))]);
  UU0_partial_d(Numtri(l,1),:) = UU0_partial_d(Numtri(l,1),:) + n;
  UU0_partial_d(Numtri(l,2),:) = UU0_partial_d(Numtri(l,2),:) + n;
  UU0_partial_d(Numtri(l,3),:) = UU0_partial_d(Numtri(l,3),:) + n;
end % for l




