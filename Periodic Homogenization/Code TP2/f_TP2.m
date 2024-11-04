function val = f_TP2(x,y)
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

%TP2 pour A=Id et u = sin(pi*x)*sin(pi*y)
%val = (1 + 2*pi*pi)*sin(pi.*x).*sin(pi.*y);

%TP2 pour (i) et u =  sin(pi*x)*sin(pi*y)
%AA_eff = [1 0; 0 2]
%val = 3*pi*pi*sin(pi.*x).*sin(pi.*y);

%TP2 pour (ii) et u =  sin(pi*x)*sin(pi*y)
%AA_eff = [sqrt(3) 0; 0 4]
%val = -pi*pi*(2*cos(2*pi.*x).*cos(pi.*x).*sin(pi.*y) - (2 + sin(2*pi.*x)).*sin(pi.*x).*sin(pi.*y) - 4*sin(pi.*x).*sin(pi.*y));

%TP2 pour (iii) et u =  sin(pi*x)*sin(pi*y)
%AA_eff = [sqrt(3) 0; 0  4]
%val = -pi*pi*(2*cos(2*pi.*x).*cos(pi.*x).*sin(pi.*y) - (2 + sin(2*pi.*x)).*sin(pi.*x).*sin(pi.*y) - (4 + sin(2*pi.*x)).*sin(pi.*x).*sin(pi.*y));

%TP2 pour (iv) et u =  sin(pi*x)*sin(pi*y)
%AA_eff = [4*sqrt(3) 0; 0 2*sqrt(15)
val = -pi*pi*(2*cos(2*pi.*x).*(4+sin(2*pi.*y)).*cos(pi.*x).*sin(pi.*y) + 2*(2+sin(2*pi.*x)).*cos(2*pi.*y).*sin(pi.*y).*cos(pi.*y) - 2*(2 + sin(2*pi.*x)).*(4 + sin(2*pi.*y)).*sin(pi.*x).*sin(pi.*y));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
