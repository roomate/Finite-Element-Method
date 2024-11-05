# 🧱 Finite Elements Method

The Finite Elements Method (FEM for short) is a modern numerical tool to solve Partial Differential Equations. Usually
presented with their finite-difference and finite volume counterpart, this method remains indispensable for whoever 
wishes to compute numerical solutions in a competitive time. The FEM is part of the class of Direct Methods, in contrast with
the multi-grid method for example, who belongs to the class of relaxation methods.

## Approximations

The FEM usually relies on three fundamental approximations. 

The first one is the meshing of the domain of interest. In the case of a non-polyhedric domain, the meshing can be inexact and it adds a supplementary
approximation. The control parameter is the meshing size h.

The second one is the use of a finite dimensional space. In the case of Lagrange polynomials, the control parameters is the degree of the 
said polynomials.

The third and last one is the quadrature formula to compute integrals in Mass and Stiffness matrix. This method is usually necessary for 
polynomials of degree greater than one. If the quadrature is chosen wisely, it can lead to interesting features of the Mass matrix such as mass lumping.