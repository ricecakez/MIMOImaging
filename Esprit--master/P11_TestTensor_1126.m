clear;clc;close all;

A(: , : ,1)=[1,0,0;1,0,0;0,0,0];
A(: , : ,2)=[0,0,0;0,1,0;0,0,0];
A(: , : ,3)=[0,0,0;0,0,0;0,0,1];
A = tensor(A);
A1=tenmat(A,1);
A2=tenmat(A,2);
A3=tenmat(A,3);
[U1,V1,W1]=svd(A1.data);
 [U2,V2,W2]=svd(A2.data);
 [U3,V3,W3]=svd(A3.data);
 U1(:,3)=[];
 S=ttm(A,{U1',U2',U3'});
 A1=ttm(S,{U1,U2,U3})