clear all
clc
%% Linear test case
xy = [0 1; 0 0; 1 0; 0 1];
spline_order_vett = [2 2; 2 3; 2 4];
X=xy(:,1); Y=xy(:,2);
[Sx,Sy]=Sx_Sy_preparation(X,Y,1,spline_order_vett);
[xyw,momsres,Z,Zin] = splcub(4,Sx,Sy,1);
nodes_x = xyw(:,1);
nodes_y = xyw(:,2);
weights = xyw(:, 3);

for n = 0:4
    true_sol = 1 / ((n+1)*(n+2));
    quad = weights' * (nodes_x .^ n);
    true_sol-quad;
end

%% Quadratic test case
clc; clear all;
xy = [0 1; 0 0; 1 0; sqrt(2)/2 sqrt(2)/2; 0 1];
spline_order_vett = [2 2; 2 3; 3 5];
X=xy(:,1); Y=xy(:,2);
[Sx,Sy]=Sx_Sy_preparation(X,Y,1,spline_order_vett);
[xyw,momsres,Z,Zin] = splcub(5,Sx,Sy,1);
nodes_x = xyw(:,1);
nodes_y = xyw(:,2);
weights = xyw(:, 3);

%true_sol = 1 / ((n+1)*(n+2));
quad = weights' * (nodes_x.^2)
%true_sol-quad

%%
polygon_sides(:,1) = [1 3/4  1/2 1/4  0   -1/4    -1/2     -3/4     -1 -0.8    -0.6    -0.4   -0.2     0     1/4     1/2        3/4    1];
polygon_sides(:,2) = [0 9/16 1/4 1/16 1   (1/4)^3  (1/2)^3  (3/4)^3  0  -0.2^5 -0.4^5  -0.6^5 -0.8^5  -1   -(3/4)^4 -(1/2)^4  -(1/4)^4 0];

spline_order_vett=[3, 3;
                   3, 5;
                   5, 9;
                   6, 14;
                   5, 18];
X=polygon_sides(:,1); Y=polygon_sides(:,2);
[Sx,Sy]=Sx_Sy_preparation(X,Y,1,spline_order_vett);
[xyw,momsres,Z,Zin] = splcub(5,Sx,Sy,1);
nodes_x = xyw(:,1);
nodes_y = xyw(:,2);
weights = xyw(:, 3);

weights' * (nodes_x.^0) - 0.95
weights' * (nodes_x.^1) - 17/315
weights' * (nodes_x.^2) - 557/2520
weights' * (nodes_x.^3) - 101/5040
weights' * (nodes_x.^4) - 79/630

weights' * (nodes_y) - (0.048958333333336+0.034375+41/756-1/22-1/18)
weights' * (nodes_x.*nodes_y) - (9/256+3/1280-317/15120+1/264-1/180)
weights' * (nodes_x.^2.*nodes_y) - (0.02578125+0.00040922619047619+43/2835-1/1716-1/990)
weights' * (nodes_x.^3.*nodes_y) - (0.019287109375+0.00013253348214286-(89/7560)+(1/8008)-1/3960)

weights' * (nodes_y.^2) - (0.016145833333333+0.015997023809524+37/1890+1/48+1/39)
weights' * (nodes_x.*nodes_y.^2) - (0.011555989583333 + 0.00072079613095238 - 3919/748440 - 1/816+1/546)
weights' * (nodes_x.^2.*nodes_y.^2) - (0.0084129050925926+0.00007750496031746+277/74844+1/7344+1/4095)

weights' * (nodes_y.^3) - (0.0061740451388889+0.00927734375+11549/1216215-1/84-1/68)
weights' * (nodes_x.*nodes_y.^3) - (0.0044135199652778+0.0003173828125-(6113/3891888)+(1/1848)-1/1224)

weights' * (nodes_y.^4) - (0.002569641729798+0.0060635653409091+8363/1496880+1/130+1/105)

theta = rand * 2 * pi;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
rpolygon_sides = polygon_sides * R';
rX=rpolygon_sides(:,1); rY=rpolygon_sides(:,2);
[Sx,Sy]=Sx_Sy_preparation(rX,rY,1,spline_order_vett);
[xyw,momsres,Z,Zin] = splcub(5,Sx,Sy,1);
rnodes_x = xyw(:,1);
rnodes_y = xyw(:,2);
rweights = xyw(:, 3);

%true_sol = 1 / ((n+1)*(n+2));
rweights' * (rnodes_x.^0) - 0.95;

%%
x1 = X(1:3);
vandermonde1 = [x1.^0 x1.^1 x1.^2 ];
coeff1 = vandermonde1\Y(1:3)

x2 = X(3:5);
vandermonde2 = [x2.^0 x2.^1 x2.^2 ];
coeff2 = vandermonde2\Y(3:5)

x3 = X(5:9);
vandermonde3 = [x3.^0 x3.^1 x3.^2 x3.^3 x3.^4];
coeff3 = vandermonde3\Y(5:9)

x4 = X(9:14);
vandermonde4 = [x4.^0 x4.^1 x4.^2 x4.^3 x4.^4 x4.^5];
coeff4 = vandermonde4\Y(9:14)

x5 = X(14:18);
vandermonde5 = [x5.^0 x5.^1 x5.^2 x5.^3 x5.^4];
coeff5 = vandermonde5\Y(14:18)
%%
function [Sx,Sy,Nsub]=Sx_Sy_preparation(X,Y,Nsub,spline_parms)

% OBJECT:
%
%
% INPUT:
% X,Y:
% Nsub:
% spline_parms:
% SPLtypestring:
%
% OUTPUT:
% Sx,Sy:
% Nsub:
%

% check input variables.
L=size(spline_parms,1);

for i=1:L % define spline
    
    if i == 1
        imin=1;
    else
        imin=spline_parms(i-1,2);
    end
    imax=spline_parms(i,2);
    
    tL=imin:imax;
    xL=X(imin:imax);
    yL=Y(imin:imax);
    
    [SxL,SyL]=compute_parametric_spline(tL,xL,yL,...
        spline_parms(i,1));
    
    Sx(i)=SxL;
    Sy(i)=SyL;
    
end
end

function [ppx,ppy]=compute_parametric_spline(s,x,y,spline_order)

% OBJECT:
% compute parametric spline relavant parameters "ppx", "ppy" so that a
% point at the boundary of the  domain has coordinates (x(s),y(s))
%
% INPUT:
% s: parameter data.
% x: determine spline x(s) interpolating (s,x)
% y: determine spline y(s) interpolating (s,y)
% spline_order: spline order (i.e. degree + 1)
% SPLtypestring: string with the spline type i.e.
%             'complete'   : match endslopes (as given in VALCONDS, with
%                     default as under *default*).
%             'not-a-knot' : make spline C^3 across first and last interior
%                     break (ignoring VALCONDS if given).
%             'periodic'   : match first and second derivatives at first
%                     data point with those at last data point (ignoring
%                     VALCONDS if given).
%             'second'     : match end second derivatives (as given in
%                    VALCONDS, with default [0 0], i.e., as in variational).
%             'variational': set end second derivatives equal to zero
%                     (ignoring VALCONDS if given).
%
% OUTPUT:
% ppx: spline x(t) data
% ppy: spline y(t) data


ppx=spapi(spline_order,s,x);
ppy=spapi(spline_order,s,y);

ppx=fn2fm(ppx,'pp');
ppy=fn2fm(ppy,'pp');
end
