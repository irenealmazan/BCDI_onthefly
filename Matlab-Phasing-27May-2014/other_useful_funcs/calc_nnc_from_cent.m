function [ nnx ] = calc_nnc_from_cent(center,arr_x,fin_x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nnx(1)=center-fin_x/2;

nnx(2)=center+fin_x/2;


nnx(2)=arr_x-nnx(2);

disp(nnx)

end

