function [ T ] = rotation_matrix3d(delta,gamma,beta )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
psi=delta;
phi=gamma;
theta=beta;


A1=cos(theta)*cos(psi);
A2=cos(theta)*sin(psi);
A3=-sin(theta);

B1=-cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
B2=cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi);
B3=sin(phi)*cos(theta);

C1=sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
C2=-sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);
C3=cos(phi)*cos(theta);

T=[A1 B1 C1 0;A2 B2 C2 0;A3 B3 C3 0;0 0 0 1];


end

