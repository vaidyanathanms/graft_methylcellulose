clc;
clear;
close all;
format long;
rHex=1;
% Use hexagons as markers
% vertHexX=[-sqrt(3)/2*rHex   0       +sqrt(3)/2*rHex +sqrt(3)/2*rHex     0       -sqrt(3)/2*rHex     -sqrt(3)/2*rHex ] ;
% vertHexY=[-rHex/2           -rHex   -rHex/2         +rHex/2             +rHex   +rHex/2             -rHex/2         ] ;
% scale = 0.5
% plotCustMark(sort(rand(20,1)*50),rand(20,1)*50,vertHexX,vertHexY,0.5)
% axis equal
% box on
% grid on

% Use crosses as markers
% vertHexX=[-0.5 +0.5 +0.5 +1.5 +1.5 +0.5 +0.5 -0.5 -0.5 -1.5 -1.5 -0.5 -0.5] ;
% vertHexY=[-3.5 -3.5 -0.5 -0.5 +0.5 +0.5	+1.5 +1.5 +0.5 +0.5 -0.5 -0.5 -3.5] ;
% scale = 0.2

%use tdt
h1 = figure;
scaleval = 1;
verttdtX1=[0 -0.5 0];
verttdtY1=[0.5 0 -0.5];

verttdtX2=[0 1 0];
verttdtY2=[0.5 -0.5 -0.5];


mf1 = 'r';
me1 = 'r';
mf2 = 'g';
me2 = 'g';

myCustMarker(0,0,verttdtX1,verttdtY1,verttdtX2,verttdtY2,scaleval,mf1,me1,mf2,me2,h1)
axis equal
box on
grid off

%use dtsh
h2 = figure;
scaleval = 1;

vert_dtsh_X1=[0 -0.5 0];
vert_dtsh_Y1=[0.5 -0.5 -0.5];
vert_dtsh_X2=[0 0.5 0.5 0];
vert_dtsh_Y2=[0.5 0.5 -0.5 -0.5];

myCustMarker(0,0,vert_dtsh_X1,vert_dtsh_Y1,vert_dtsh_X2,vert_dtsh_Y2,scaleval,mf1,me1,mf2,me2,h2)
axis equal
box on
grid off



%use dthex
h3 = figure;
scaleval = 0.5;

vert_thex_X1=[0 -0.5 0];
vert_thex_Y1=[0.5 0 -0.5];
vert_thex_X2=[0 0.25 0.5 0.25 0.5 0.25 0];
vert_thex_Y2=[0.5 0.25 0.25 0 -0.25 -0.25 -0.5];

x = [1 2 3];
y = x.^2;
myCustMarker(x,y,vert_thex_X1,vert_thex_Y1,vert_thex_X2,vert_thex_Y2,scaleval,mf1,me1,mf2,me2,h3)
axis equal
box on
grid off

%plot two at the same time
y1 = 2+x.^3;
h4 = figure;
myCustMarker(x,y1,vert_dtsh_X1,vert_dtsh_Y1,vert_dtsh_X2,vert_dtsh_Y2,scaleval,mf1,me1,mf2,me2,h4)
myCustMarker(x,y,vert_thex_X1,vert_thex_Y1,vert_thex_X2,vert_thex_Y2,scaleval,mf1,me1,mf2,me2,h4)
