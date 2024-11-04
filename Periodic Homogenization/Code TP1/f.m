function val = f(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A COMPLETER
%Conditon de Neumann : A scalaire
%val = (1 + 5*pi*pi)*cos(pi*x).*cos(2*pi*y);

%val = cos(cos(x)).*sin(sin(y)); 

%Condition de Neumann et périodique: A variable 
val = cos(pi.*x).*cos(2*pi.*y) + pi*pi*(cos(2*pi.*y).*sin(2*pi.*y).*(2*cos(2*pi.*x).*sin(pi.*x) + 4*sin(2*pi.*x).*cos(pi.*x))+5*(2+sin(2*pi.*x).*sin(2*pi.*y)).*cos(pi.*x).*cos(2*pi.*y));

%Condition de Dirichlet : A variable
%val = sin(pi.*x).*sin(2*pi.*y) - pi*pi*(-5*sin(pi.*x).*sin(2*pi.*y).*A(x,y) + 2*cos(2*pi.*x).*cos(pi.*x).*sin(2*pi.*y));

%Condition de Dirichlet : A scalaire
%val = (1 + 5*pi*pi)*sin(pi*x).*sin(2*pi*y);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
