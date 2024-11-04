function a=produit_scalaire(x,y)
a = 0;
for i=1:size(x,1)
    a = a + x(i)*y(i);
end
    