function [A_eff] = A_effi(W_x,W_y,x,y,KK)

A_eff = zeros(2);
A_eff(1,1) = transpose(W_x)*KK*W_x + transpose(W_x)*KK*x + transpose(x)*KK*W_x + transpose(x)*KK*x;
A_eff(1,2) = transpose(W_y)*KK*W_x + transpose(W_y)*KK*x + transpose(y)*KK*W_x + transpose(y)*KK*x;
A_eff(2,1) = transpose(W_x)*KK*W_y + transpose(W_x)*KK*y + transpose(x)*KK*W_y + transpose(x)*KK*y;
A_eff(2,2) = transpose(W_y)*KK*W_y + transpose(W_y)*KK*y + transpose(y)*KK*W_y + transpose(y)*KK*y;