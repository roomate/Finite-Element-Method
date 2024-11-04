function val = mat_A(x,y,eps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_A :
% Evaluation de la matrice A .
%
% SYNOPSIS val = mat_A(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la matrice sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A COMPLETER
%TP2 : 1er exemple
%val = [1 0;0 1];

%TP2 : 2nd exemple (i)
%val = [1 0;0 2];

% %TP2 : 3rd exemple (ii)
%val = [2 + sin(2*pi.*x/eps) 0;0 4];
% 
% %TP2 : 4th ex/epsemple (iii)
%val = [2 + sin(2*pi.*x/eps) 0;0 4 + sin(2*pi.*x/eps)];
% 
% %TP2 : 5ème exemple (iv)
val = [1 0;0 1]*(2+sin(2*pi.*x/eps))*(4 + sin(2*pi.*y/eps));

%First part
%val = sin(2*pi.*x/eps).*sin(2*pi.*y/eps)+2;

%Second part
%val = sin(2*pi*x/eps) + 2;

%val = [2+sin(2*pi*x/(10*eps)) 0;0 2+sin(2*pi*x/(10*eps))];

%val = [2+sin(2*pi*x/(eps*100)) 0;0 2+sin(2*pi*x/(eps*100))];
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
