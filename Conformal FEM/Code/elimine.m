function [tilde_AA,tilde_LL] = elimine(AA, LL, Refneu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elimine :
% Routine qui fait la pseudo-elimination pour les noeuds avec Refneu(i) = 1
% Conditions aux limites Dirichlet homogenes.
%
% SYNOPSIS elimine(AA,LL,Refneu)
%          
% INPUT * AA, LL, Refneu: La matrice et le second membre associes au
% probleme sans elimination et references des noeuds sur la frontiere
% Dirichlet.
%
% OUTPUT - tilde_AA, tilde_LL: On rend la matrice et le second membre une fois la 
% pseudo elimination realisee.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% On copie la matrice AA et le second membre dans AAE et FFE.
tilde_AA = AA;
tilde_LL = LL;

% Nombre de sommets dans la discretisation.
Nbpt = size(tilde_AA,1);

% Boucle sur les sommets.
for l=1:Nbpt   
    % Si le noeud est Dirichlet on fait la pseudo elimination.
    if Refneu(l)==1 || Refneu(l)==2
        tilde_LL(l)=0;
        for i=1:Npt
            if i==l
                tilde_AA(i,i)=1;
            else
                tilde_AA(i,l)=0;
                tilde_AA(l,i)=0;
            end
        end
    end     
                       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021