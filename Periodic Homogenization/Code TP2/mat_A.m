function [A] = mat_A(S1,S2,S3)

x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

A = [x2-x1, x3-x1; y2-y1, y3-y1]; 