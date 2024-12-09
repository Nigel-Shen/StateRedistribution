
function [xyw,res,Z,Zin,cmom,bbox] = splcub(n,Sx,Sy,extraction_type)

% OBJECT:
% it computes an algebraic formula with positive weights and interior nodes
% for degree of exactness n with respect to the Lebesgue measure in a
% Jordan spline polygon, whose boundary is given by the counterclockwise
% concatened spline arcs (Sx,Sy)
%
% INPUT:
% n: polynomial degree
% Sx,Sy: arrays of spline structures; (SX(i),Sy(i)) is the i-th arc
% the arcs must be counterclockwise concatenated forming a Jordan curve
% Sx(i).breaks(end)=Sx(i+1).breaks(1),Sy(i).breaks(end)=Sy(i+1).breaks(1)
% i=1,...,end, with end+1:=1
% extraction_type: 1: lsqnonneg, 2: fast lawson-hanson 3: backslash
%
% OUTPUT:
% xyw: 3-column array of nodes coords (cols 1-2) and weights (col 3)
% res: norm 2 of compressed moments vs moments
% Zin: internal points
% Z  : all analysed points
% cmom: Chebyshev-Vandermonde moments
% bbox: 
%
% DATA:
%
% built: April 2019
% modified: June 04, 2020
%
% CODE VERSION:
%
% v: 11.0.1 (conditioning issues)
% v: 11.0.2 (some minor bugs for low degrees)

% ....................... troubleshooting .................................
if nargin < 4
    extraction_type=1;
end

% ......................... initialize ....................................
alpha=1.5; beta=1.5; % discretization parameters

k=floor(max(n,5)^(alpha)); % we use max to avoid too many iterations for 
                           % low "n".

% Note: on the setting of restol: 10^(-14) seems fine, 10^(-15) is often  
% not achieved.                           
restol=1*10^(-14);

Z=[];   % all points
Zin=[]; % points in the domain
V=[];   % Vandermonde matrix

Lmax=6; % maximum number of sublevels to analyse.

% ......................... exact moments .................................

[cmom,bbox]=chebmom(n,Sx,Sy); % compute Chebyshev-Vandermonde moments

for i=1:Lmax
    
    % ....................... define local set ............................
    x=linspace(bbox(1),bbox(2),k);
    y=linspace(bbox(3),bbox(4),k);
    [u,v]=meshgrid(x,y);
    ZL=[u(:) v(:)];
    
    % .......................... indomain .................................
    
    [inside,onL]=incurvpolygon(ZL,Sx,Sy);
    ind=(inside==1 & isnan(onL)==1);
    
    % .......................... storage ..................................
    
    X=ZL(ind,:);    % points inside the domain
    
    if length(X) == 0
        % mesh refinement in case of failures (no points in domain)
        k=floor(beta*(k+1));
    else
        Zin=[Zin; X];   % all points inside the domain
        VL=cvand(n,X,bbox);
        V=[V; VL];      % Vandermonde of all points inside the domain
        Z=[Z; ZL];      % all points: inside/outside domain
        [Q,R]=qr(V,0);
        cmom1=R'\cmom;
        
        % .......................... compression ..........................
        
        switch extraction_type
            case 1
                % ... w=lsqnonneg(V',cmom);
                w=lsqnonneg(Q',cmom1);
            case 2
                % ... w = lawsonhanson(V', cmom);
                w = lawsonhanson(Q', cmom1);
            otherwise
                % ... w=V'\cmom;
                w=Q'\cmom1;
        end
        
        
        ind1=find(abs(w)>0); 
        w=w(ind1);
        
        % ......................... error analysis ........................
        % res=norm(Q(ind1,:)'*w-cmom1); 
        res=norm(V(ind1,:)'*w-cmom); 
        
        if res>restol | length(w) == 0
            k=floor(beta*k);
        else
            xyw=[Zin(ind1,:) w]; 
            return;
        end
        
    end
    
end

warning('\n \t Compression not completed. Moment error: %1.3e',...
    res);
xyw=[Zin(ind1,:) w];





















%==========================================================================
% ADDITIONAL ROUTINES
%==========================================================================


%==========================================================================
% cvand
%==========================================================================

function V = cvand(n,pts,bbox)

% OBJECT:
% computes by recurrence the Chebyshev-Vandermonde matrix on a 2d
% arbitrarily located mesh, in a total-degree product Chebyshev basis
% of a given rectangle
%
% INPUT:
% n: polynomial degree
% pts: 2-column array of mesh point coordinates
% bbox: [bbox(1),bbox(2)] x [bbox(3),bbox(4)] bounding box for the mesh
%
% OUTPUT:
% V: Chebyshev-Vandermonde matrix
%
% DATA:
%
% built: March 2019
% check: May 2, 2019


% default rectangle containing the mesh if not passed
if isempty(bbox)
    bbox=[min(pts(:,1)) max(pts(:,1)) min(pts(:,2)) max(pts(:,2))];
end

a=bbox(1);b=bbox(2);c=bbox(3);d=bbox(4);

Vx=chebpolys(n,(2*pts(:,1)-b-a)/(b-a));
Vy=chebpolys(n,(2*pts(:,2)-d-c)/(d-c));

k=0;
for i=0:n
    for j=0:n-i
        k=k+1;
        V(:,k)=Vx(:,i+1).*Vy(:,j+1);
    end
end




%==========================================================================
% chebpolys
%==========================================================================

function T=chebpolys(deg,x)

% OBJECT:
% computes the Chebyshev-Vandermonde matrix on the real line by recurrence
%
% INPUT:
% deg = maximum polynomial degree
% x = 1-column array of abscissas
%
% OUTPUT
% T: Chebyshev-Vandermonde matrix at x, T(i,j+1)=T_j(x_i), j=0,...,deg
%
% DATA:
%
% built: March 2019
% check: May 2, 2019

% inizialization
T=zeros(length(x),deg+1);
t0=ones(length(x),1);
T(:,1)=t0;
t1=x;
T(:,2)=t1;

% 3-term recurrence
for j=2:deg
    t2=2*x.*t1-t0;
    T(:,j+1)=t2;
    t0=t1;
    t1=t2;
end





%==========================================================================
% chebmom
%==========================================================================

function [cmom,bbox] = chebmom(n,Sx,Sy)

% OBJECT:
% computes the moments up to degree m of a total-degree product Chebyshev
% basis with respect to the Lebesgue measure in a Jordan spline polygon,
% whose boundary is given by the counterclockwise concatened spline
% arcs (Sx,Sy)
%
% INPUT:
% n: polynomial degree
% Sx,Sy: arrays of spline structures; (SX(i),Sy(i)) is the i-th arc
% the arcs must be counterclockwise concatenated forming a Jordan curve
% Sx(i).breaks(end)=Sx(i+1).breaks(1),Sy(i).breaks(end)=Sy(i+1).breaks(1)
% i=1,...,end, with end+1:=1

% OUTPUT:
% cmom: array of bivariate Chebyshev moments
% bbox: bounding box for the spline polygon
%
% DATA:
% built: March 2019
% check: May 2, 2019
% modified: June 4, 2020

% Note. From this routine, we compute the approximation of a rectangle bbox
% that is a bounding box for the spline polygon. For low "n" such rectangle
% may not be good enough. Thus, in this case, we increase a little the 
% degree so to be acceptable. 

xyw=lineint(n,Sx,Sy); 

bbox=[min(xyw(:,1)) max(xyw(:,1)) min(xyw(:,2)) max(xyw(:,2))];
intV=intcvand(n,xyw(:,1:2),bbox); 
cmom=intV'*xyw(:,3);




%==========================================================================
% lineint
%==========================================================================

function xyw = lineint(m,Sx,Sy)

% OBJECT:
% computes weights and nodes for the line integral on concatenated 
% spline arcs; the formula is exact on bivariate polynomials up to deg m
%
% INPUT:
% m = polynomial degree to be integrated, via its x-primitive, on the
%     spline-curvilinear domain.
% Sx,Sy: arrays of spline structures; (SX(i),Sy(i)) is the i-th arc
% the arcs must be concatenated
% Sx(i).breaks(end)=Sx(i+1).breaks(1),Sy(i).breaks(end)=Sy(i+1).breaks(1)
%
% OUTPUT:
% xyw: 3-column array of nodes coords (cols 1-2) and weights (col 3)
%
% DATA:
% built: March 2019
% check: May 2, 2019

xyw=[];
for i=1:length(Sx)
    
    ord_spline=Sx(i).order; % degree is the spline order minus 1.
    
    % Gauss-Legendre nodes in [-1,1] and corresponding weights on the side.
    k=ceil(((ord_spline-1)*(m+2))/2)+2; ab=r_jacobi(k,0,0); xw=gauss(k,ab);
    t=xw(:,1); w=xw(:,2);
    
    a=Sx(i).breaks(1:end-1); b=Sx(i).breaks(2:end);
    alpha=(b-a)/2; beta=(b+a)/2;
    dSy(i)=fnder(Sy(i));
    for j=1:length(a)
        nodes=alpha(j)*t+beta(j);
        wloc=w*alpha(j);
        xyw=[xyw;[ppval(Sx(i),nodes) ppval(Sy(i),nodes) ...
            wloc.*ppval(dSy(i),nodes)]];
    end
end





%==========================================================================
% intcvand
%==========================================================================

function intV = intcvand(n,pts,bbox)

% computes by recurrence an x-primitive Chebyshev-Vandermonde matrix on a
% 2d arbitrarily located mesh, in a total-degree product Chebyshev basis
% of a given rectangle

% March 2019

% INPUT:
% n: polynomial degree
% pts: 2-column array of mesh point coordinates
% bbox: [bbox(1),bbox(2)] x [bbox(3),bbox(4)] bounding box for the mesh

% OUTPUT:
% intV: x-primitive of the Chebyshev-Vandermonde matrix


% default rectangle containing the mesh if not passed

if isempty(bbox)
    bbox=[min(pts(:,1)) max(pts(:,1)) min(pts(:,2)) max(pts(:,2))];
end

a=bbox(1);b=bbox(2);c=bbox(3);d=bbox(4);

Vx=chebpolys(n+1,(2*pts(:,1)-b-a)/(b-a));
Vy=chebpolys(n,(2*pts(:,2)-d-c)/(d-c));

k=0;
for i=0:n
    for j=0:n-i
        k=k+1;
        
        if i==0
            intV(:,k)=pts(:,1).*Vy(:,j+1);
        end
        
        if i==1
            intV(:,k)=0.25*(b-a)*(((2*pts(:,1)-b-a)/(b-a)).^2).*Vy(:,j+1);
        end
        
        if i>1
            intV(:,k)=0.25*(b-a)*(Vx(:,i+2).*Vy(:,j+1)/(i+1)-...
                Vx(:,i).*Vy(:,j+1)/(i-1));
        end
        
    end
end


%==========================================================================
% r_jacobi
%==========================================================================

function ab=r_jacobi(N,a,b)
%R_JACOBI Recurrence coefficients for monic Jacobi polynomials.
%   AB=R_JACOBI(N,A,B) generates the Nx2 array AB of the first
%   N recurrence coefficients for the monic Jacobi polynomials
%   orthogonal on [-1,1] relative to the weight function
%   w(x)=(1-x)^A*(1+x)^B. The call AB=R_JACOBI(N,A) is the same
%   as AB=R_JACOBI(N,A,A) and AB=R_JACOBI(N) the same as
%   AB=R_JACOBI(N,0,0).
%
%   Supplied by Dirk Laurie, 6-22-1998; edited by Walter
%   Gautschi, 4-4-2002.
%   Edited by Walter Leopardi 10-22-2006.

if nargin<2, a=0; end
if nargin<3, b=a; end
if((N<=0)|(a<=-1)|(b<=-1)) error('parameter(s) out of range'), end
nu=(b-a)/(a+b+2);
if a+b+2 > 128
    mu=exp((a+b+1)*log(2)+((gammaln(a+1)+gammaln(b+1))-gammaln(a+b+2)));
else
    mu=2^(a+b+1)*((gamma(a+1)*gamma(b+1))/gamma(a+b+2));
end
if N==1, ab=[nu mu]; return, end
N=N-1; n=1:N; nab=2*n+a+b;
A=[nu (b^2-a^2)*ones(1,N)./(nab.*(nab+2))];
n=2:N; nab=nab(n);
B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
ab=[A' [mu; B1; B']];





%==========================================================================
% gauss
%==========================================================================

function xw=gauss(N,ab)
%GAUSS Gauss quadrature rule.
%   GAUSS(N,AB) generates the Nx2 array XW of Gauss quadrature
%   nodes and weights for a given weight function W. The nodes,
%   in increasing order, are placed into the first column of XW,
%   and the corresponding weights into the second column. The
%   weight function W is specified by the Nx2 input array AB
%   of recurrence coefficients for the polynomials orthogonal
%   with respect to the weight function W.

N0=size(ab,1); if N0<N, error('input array ab too short'), end
J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end
for n=2:N
    J(n,n-1)=sqrt(ab(n,2));
    J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);
[D,I]=sort(diag(D));
V=V(:,I);
xw=[D ab(1,2)*V(1,:)'.^2];











