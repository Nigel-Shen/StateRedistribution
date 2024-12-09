
function [in,on,Sx,Sy,boxVx,turning_points_Sx,turning_points_Sy,...
    singular_junctions_Sx]=incurvpolygon(P,X,Y,Nsub,spline_parms,...
    SPLtypestring,boxVx,singular_points)

%
% OBJECT:
% This routine determines if a sequence of points is inside a curvilinear
% polygon whose boundary is defined by a sequence of parametric splines,
% possibly of different order. Domain of such kind include the classical
% polygons, for which there already exists the "inpolygon" function.
%
% INPUT:
% P: a N x 2 matrix of N bivariate points for which we want to test if they
%    are inside the polygon.
% X,Y: can be
%   1. a sequence of vertices whose dimension, in each case, is M x 1.
%   2. a vector of splines, say "Sx", "Sy" whose boundary is defined by
%      as the union of the arcs ((Sx(i))(t),((Sy(i))(t)) with "t" in
%      "[t(i),t(i+1)]", so that the domain has no self intersections and
%    ((Sx(1))(t(1)),((Sy(1))(t(1)))=((Sx(end))(t(end)),((Sy(end))(t(end)))
%      i.e. the boundary is a closed curve.
%      It is supposed that "t(i) < t(i+1)".
% spline_parms: this variable is useful only when X,Y are vectors and can
%      be dropped otherwise, for instance setting "spline_parms=[]";
%      spline_parms is a vector whose i-th row is of the form
%                        [order(i) index(i)]
%      In this case the i-th arc will be defined by the splines
%      "Sx", "Sy" of  order "order(i)" so that the curve
%                    "((Sx(i))(t),((Sy(i))(t))"
%      with "t" in "[t(i),t(i+1)]" will interpolate all the points
%                         (X(t(k)),Y(t((k)))
%      where "k" ranges in the interval [index(i-1), index(i)] (where as
%      default index(0)=1).
% SPLtypestring: extra parameter defining the particular cubic spline (e.g
%      SPLtypestring='not-a-knot').
% boxVx: indomain boxes (if available from previous usage of the function);
%      if not available do not declare the variable or set "boxVx=[]".
% singular_points: singular points (if available from previous usage of
%      the function); if not available do not declare the variable or set
%      "singular_points=[]".
%
% OUTPUT:
% in: vector with the cardinality of P, where 
%    * in(k)=0 means that P(k,:) is not in the domain, 
%    * in(k)=1 means that P(k,:) is in the domain, 
% on: vector with the cardinality of P, where 
%    * on(k)=NaN means that P(k,:) is not on the boundary, 
%    * on(k)=is a number,  means that P(k,:) is "close" to the boundary,
% Sx,Sy: parametric splines describing the boundary in the x and y
%    variable.
% boxVx: monotone boxes used in the algorithm
% turning_points_Sx: points where the derivative of Sx is null.
% turning_points_Sy: points where the derivative of Sy is null.
% singular_junctions_Sx: abscissae where is difficult to determine the
%    crossing.
%
% NOTE:
% If "boxVx" and "singular_points" are not empty we do not need to
% recompute them and we proceed only with determining points in the domain.
% The typical situation is that this routine has already been applied
% before and for some reasons needs to be called again. In this case the
% splines that describe the border are available and passed via X and Y.
% Notice that in this situation the variables "turning_points_Sx",
% "turning_points_Sy", "singular_junctions_Sx" were already at hand and
% thus one should use a call of the form
%
% [in,on]=incurvpolygon_v10c(P,X,Y,[],spline_parms,SPLtypestring,boxVx,...
%    singular_points)
%
% ADDITIONAL DATA:
%
% built: april 2019
% modified: May 29, 2019
%
% CODE VERSION:
%
% v: 11.0.0 (introduction of not simply connected domains)

% ...................... troubleshooting/setting ..........................

if nargin < 4
    Nsub=[];
end

if nargin < 7
    boxVx=[];
end

if nargin < 8
    singular_points=[];
end

% *************************** special case ********************************

% Note:
% In this case it is supposed that boxes (as well as splines) are already
% available.

tol=10^(-11); % boundary tolerance.

% case of further application of the routine
if isempty(boxVx) == 0 & isempty(singular_points) == 0
    Sx=X; Sy=Y;
    if size(P,1) > 0
        [in,on]=indomain(P,Sx,Sy,boxVx,singular_points(:,1),tol);
    else
        in=[];
        on=[];
    end
    
    % these variables were already at hand from previous calls
    turning_points_Sx=[];
    turning_points_Sy=[];
    singular_junctions_Sx=[];
    
    return
end


% *************************** general case ********************************

% .................. check if X,Y are splines or scalars ..................

if isstruct(X(1)) == 0 % X is a vector of numbers
    if nargin < 5, spline_parms=[2 length(X)]; end
    if nargin < 6, SPLtypestring='not-a-knot'; end
    [Sx,Sy,Nsub]=Sx_Sy_preparation(X,Y,Nsub,spline_parms,SPLtypestring);
else
    % X is a spline
    Sx=X; Sy=Y;
    for ii=1:length(Sx)
        SxL=Sx(ii);
        LV(ii)=length(SxL.breaks);
    end
    L=sum(LV); Nsub=ceil(210/L);
end


% ........................ compute boxes ..................................
[boxVx,turning_points_Sx,turning_points_Sy]=spline_boxer(Sx,Sy,Nsub);

% ........................ singular points ................................

singular_junctions_Sx=singular_vertices_v8a(Sx,Sy);
singular_points=define_singpts(turning_points_Sx,singular_junctions_Sx);

% ........................... indomain ....................................

[in,on]=indomain(P,Sx,Sy,boxVx,singular_points(:,1),tol);

% ................. indomain for points in singular zones .................

index_sz=find(isnan(in) == 1);
for ii=1:length(index_sz)
    windn = winding_algorithm(P,Sx,Sy);
    in(index_sz)=rem(windn,2);
end

if min(size(index_sz)) > 1
    fprintf('\n \t *** used winding algorithm %3.0f times',...
        length(index_sz));
end





%--------------------------------------------------------------------------
% Attached functions.
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Sx_Sy_preparation
%--------------------------------------------------------------------------

function [Sx,Sy,Nsub]=Sx_Sy_preparation(X,Y,Nsub,spline_parms,SPLtypestring)

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
        spline_parms(i,1),SPLtypestring);
    
    Sx(i)=SxL;
    Sy(i)=SyL;
    
end




%--------------------------------------------------------------------------
% define_singpts
%--------------------------------------------------------------------------

function singular_points=define_singpts(turning_points_Sx,...
    singular_junctions_Sx)

% OBJECT:
%
%
% INPUT:
% turning_points_Sx:
% singular_junctions_Sx:
%
% OUTPUT:
% singular_points:
%

if isempty(turning_points_Sx) == 0
    singular_points=[turning_points_Sx(:,1:2)];
else
    singular_points=[];
end

if isempty(singular_junctions_Sx) == 0
    singular_points=[singular_points; singular_junctions_Sx(:,1:2)];
end





%--------------------------------------------------------------------------
% indomain (and subroutines)
%--------------------------------------------------------------------------

function [in,on]=indomain(P,Sx,Sy,boxV,singptsX,tol)

% OBJECT:
% indomain_v6 True for points inside or on a spline curvilinear region.
%
% INPUT:
% P: points to test, defined as an N x 2 matrix.
% Sx: M x 1 vector, whose i-th component is a spline defining the "x=x(t)"
%     parametric equation that defines the boundary.
% Sy: M x 1 vector, whose i-th component is a spline defining the "y=y(t)"
%     parametric equation that defines the boundary.
% boxV: S x k vector, whose i-th component is a spline defines the
%      properties of the i-th box of the boundary, i.e.
%      *           [xm,xM,ym,yM,hardbox,ispline,iblock]
%      where [xm,xM,ym,yM] defines the boundary of the box, i.e. vertices
%                   X=[xm xM xM xm]   Y=[ym ym yM yM]
%      * "hardbox" says if it is a box with changes in the derivatives of
%      x(y), y(t),
%      * ispline says what spline index is defined in the box, i.e.
%        x(t)=[Sx(ispline)](t), y(t)=[Sy(ispline)](t)
%      * iblock says what spline block index is defined in the box, i.e.
%        x(t)=[Sx(ispline)](t), y(t)=[Sy(ispline)](t) with
%        "t in [t1,t2]", where "t1=(Sx(ispline)).breaks(iblock)" and
%        "t2=(Sx(ispline)).breaks(iblock+1)".
% singptsX: abscissae of singular points.
% tol: tolerance in detecting the point in the box by polynomial solver.
%
% OUTPUT:
% in: the i-th component is 1 if P(i,:) is inside the domain,
%                           0 if not inside,
%                           NaN in case of doubts and other winding
%                           algorithms must be used,
% on: the i-th component is a number if P(i,:) is "close" to the boundary
%     of the domain, NaN if not on the domain.
%
% ROUTINES USED IN THIS FUNCTION:
% 1. spline_boxer (attached)
% 2. box_analysis (attached)
% 3. singularity_checks (attached)
% plus all the functions called by these functions.
%
% Version: 8a1.
% * November 26, 2019 (bug dealing with very peculiar boxes)

% ................. troubleshoouting and defaults .........................

if nargin < 4
    [boxV,turning_points_Sx,boxVy,turning_points_Sy]=spline_boxer(Sx,Sy);
end

if nargin < 5, singptsX=[]; end
if nargin < 6, tol=10^(-11); end

% ......................... general study .................................

X=(P(:,1))'; Y=(P(:,2))';
a=boxV(:,1); b=boxV(:,2); c=boxV(:,3); d=boxV(:,4);
iboxes=(X > a & X <=b & Y >= c); 
itestL=(iboxes & Y <= d);
iregL=(iboxes & Y > d);
crossings=sum(iregL,1);


% ..................... study points inside boxes .........................

further_analysis_index=find(sum(itestL,1) > 0);

on=NaN*ones(size(X));
for jj=1:length(further_analysis_index)
    Pindex=further_analysis_index(jj);
    PL=P(Pindex,:);
    itestPL0=itestL(:,Pindex);
    index_testPL=find(itestPL0 > 0);
    [final_value0,on(Pindex)]=box_analysis(PL,index_testPL,boxV,Sx,Sy,tol);
    crossings(Pindex)=crossings(Pindex)+final_value0;
end
in=rem(crossings,2);

% ............... study points close to singularities .....................

if isempty(singptsX) == 0
    iflag=singularity_checks(P(:,1),singptsX);
    iflag0=find(isfinite(on(iflag)) == 0);
    iflagF=iflag(iflag0);
    on(iflagF)=-Inf*ones(length(iflagF),1);
    in(iflagF)=on(iflagF);
end

% ................... safe exit for boundary points .......................

% Though it is not necessary, we double check that points on the boundary
% cannot also be addressed as inside the domain.

ion=find(isfinite(on) == 1);
in(ion)=0;





%--------------------------------------------------------------------------
% box_analysis
%--------------------------------------------------------------------------

function [final_value,flag]=box_analysis(P,itest,boxV,Sx,Sy,tol)

% INPUT:
% P: point to test, defined as an 1 x 2 matrix.
% itest: index of the box in the set described by "boxV", to be studied.
% Sx: M x 1 vector, whose i-th component is a spline defining the "x=x(t)"
%     parametric equation that defines the boundary.
% Sy: M x 1 vector, whose i-th component is a spline defining the "y=y(t)"
%     parametric equation that defines the boundary.
% boxV: S x k vector, whose i-th component is a spline defines the
%      properties of the i-th box of the boundary, i.e.
%      *           [xm xM ym yM tm tM isp ibl der_signx der_signy];
%      where [xm,xM,ym,yM] defines the boundary of the box, i.e. vertices
%                   X=[xm xM xM xm]   Y=[ym ym yM yM]
%      * [tm,tM] t variable extremes
%      * ispline says what spline index is defined in the box, i.e.
%        x(t)=[Sx(ispline)](t), y(t)=[Sy(ispline)](t)
%      * iblock says what spline block index is defined in the box, i.e.
%        x(t)=[Sx(ispline)](t), y(t)=[Sy(ispline)](t) with
%        "t in [t1,t2]", where "t1=(Sx(ispline)).breaks(iblock)" and
%        "t2=(Sx(ispline)).breaks(iblock+1)".
%      * der_signx der_signy: signs of the derivative of the spline in the
%      box
% tol: tolerance in detecting the point in the box by polynomial solver.
%
% OUTPUT:
% final_value: crossing value, odd if in even if outside the curv. domain.
% flag:
%      * Nan    : point is "safely" outside or inside the domain.
%      * finite : point is close to the boundary
%      * -Inf   : problems with the determination inside/on/outside.

if nargin < 5
    singptsX=[];
end

if nargin < 6
    tol=10^(-11);
end

x=P(1); y=P(2);

flag=NaN;

final_value=0;

%--------------------------------------------------------------------------
% 3. Further analysis on special boxes (hard boxes or point in a box).
%--------------------------------------------------------------------------

if length(itest) > 0
    
    for ii=1:length(itest)
        
        % "itestL" indices to test
        itestL=itest(ii);
        isplL=boxV(itestL,7);
        icoefsL=boxV(itestL,8);
        SxL=Sx(isplL);
        SyL=Sy(isplL);
        
        % ................ SPLINES OF ORDER > 2 ................
        
        if SxL.order > 2
            
            SxL_coefs=SxL.coefs(icoefsL,:);
            SyL_coefs=SyL.coefs(icoefsL,:);
            a=boxV(itestL,5);
            b=boxV(itestL,6);
            
            t0=SxL.breaks(icoefsL);
            t1=SxL.breaks(icoefsL+1);
            
            [local_value,flag]=spline_over(P,SxL_coefs,SyL_coefs,a,b,...
                t0,t1,tol);
            
            % Technical note:
            % point is on the boundary: here local value is 0 but some
            % other temporary "final_value" could have stored before, so we
            % set definitively "final_value=0" since point is approx. on the
            % boundary more than inside the domain.
            %
            % At the end of the routine "flag" is a number, i.e. on the
            % boundary, or NaN (point has been well analysed in the
            % crossing).
            
            if isfinite(flag) == 1
                final_value=0; % point is on the boundary and
                % not inside the domain
                return;
            end
            
        else
            
            % ................ SPLINES OF ORDER = 2 ................
            
            nodes=SxL.breaks;
            a=boxV(itestL,5);
            b=boxV(itestL,6);
            
            P1=[ppval(SxL,a) ppval(SyL,a)];
            P2=[ppval(SxL,b) ppval(SyL,b)];
            
            % ... ... ... sides are not vertical ... ... ...
            
            if abs(P1(1)-P2(1)) > 0
                
                x0=P1(1); y0=P1(2);
                x1=P2(1); y1=P2(2);
                m=(y1-y0)/(x1-x0);
                yy=y0+m*(x-x0);
                
                % point on the boundary
                if abs(yy-y) <= tol
                    final_value=0;
                    flag=abs(yy-y);
                    return;
                else % point is not on the boundary
                    if yy < y
                        local_value=1; % crossing
                    else
                        local_value=0; % not crossing
                    end
                end
                
                
            else
                % In this case the point has coordinates (x,y) with
                % x=P1(1)=P2(1). It is a difficult istance and must be
                % analyzed by winding if point is not on the side.
                % Winding requirement is communicated by "flag=-Inf".
                
                min_y=min(P1(2),P2(2));
                max_y=max(P1(2),P2(2));
                
                % point on a vertical side
                if (y >= min_y) & (y <= max_y)
                    flag=0;
                    final_value=0;
                    return;
                else
                    flag=-Inf;
                    final_value=NaN;
                    return;
                end
                
            end
            
            
        end
        
        if isnan(flag) == 1
            final_value=final_value+local_value;
        else % points has some issues
            final_value=0;
            return;
        end
        
    end
end




%--------------------------------------------------------------------------
% spline_over
%--------------------------------------------------------------------------

function [val,flag]=spline_over(P,xtv,ytv,a,b,t0,t1,tol)

% INPUT:
% P   : point to test, 1 x 2 vector. It is supposed that the point is in
%       the rectangle with "SxL([t0,t1]) x SyL([t0,t1])".
% xtv : x=x(t) "local" coefficients of a spline "SxL".
% ytv : y=y(t) "local" coefficients of a spline "SyL".
% a,b : spline "SxL","SyL" extrema (subsequent break points of the SxL,SyL).
% t0,t1: the interval [t0,t1] is a subset of [a,b] (the portion of the
% boundary to analyse has the form (SxL(t),SyL(t)) with "t in [t0,t1]".
% tol : tolerance of the boundary (to say how much P is close to the
%      boundary).
%
% OUTPUT:
% val: how many times a straight line with x=P(1) intersects the spline
%      boundary.
% flag:
%     * NaN           : normal exit, point not on the boundary.
%     * finite number : parameter that says how close is to the border; in
%     this case the point is accepted as on the boundary since its distance
%     from the boundary is below a certain threshold (given by "tol").
%     * -Inf: the point P has abscissa not in the range of "SxL([t0,t1])".
%
% Checked: May 1, 2019.


% ..................... setting missing defaults ..........................

if nargin < 6
    tol=10^(-12);
end


% ......................... initial settings ..............................

val=0;
flag=NaN;
x=P(1); y=P(2);
s0=a-t0; s1=b-t0;

% ...................... compute intersections ............................

xtv0=xtv;
xtv0(end)=xtv0(end)-x;

% ..... special cases .....

if norm(xtv0(1:end-1)) == 0
    % All the points are solution.
    % Geometrical interpretation: vertical line x=P(1). The box is not a
    % rectangle but a vertical segment.
    if xtv0(end) == 0
        val=0;
        flag=0;
        % No solution to the problem.
        % Geometrical interpretation: the point P has abscissa not in the
        %                             range of "SxL([t0,t1])"
    else
        val=NaN;
        flag=-Inf;
    end
    return;
end


% ..... normal case .....
tR=myroots(xtv0,s0,s1);


% check relevant boxes
for j=1:length(tR)
    
    tRL=tR(j);
    yL=polyval(ytv,tRL);
    
    % point "safely" over the boundary
    if yL < y+tol
        val=val+1;
    end
    
    % point "safely" close to the boundary
    if abs(yL-y) <= tol
        val=0;
        flag=abs(yL-y);
        return;
    end
    
end





%--------------------------------------------------------------------------
% singularity_checks.
%--------------------------------------------------------------------------

function iflag=singularity_checks(P,singular_pointsX)

% OBJECT:
% This routine analyses points of the set P that are "close" to
% singularities. These points devote a special analysis in the crossing
% algorithm.
%
% INPUT:
% P: matrix of size M x 2 (coordinates of the points)
% singular points: abscissas of the domain in which problems may arise.
%
% OUTPUT:
% iflag: vector of indices determining the indices of P where special
%    algorithms are required.
%
% Checked: May 1, 2019.

if isempty(singular_pointsX) == 0
    toll=10^(-12);
    iflag=find(min(abs(singular_pointsX-(P(:,1))'),[],1) < toll);
else
    iflag=[];
end






%--------------------------------------------------------------------------
% myroots.
%--------------------------------------------------------------------------

function [sols,flag]=myroots(p,t1,t2)

% OBJECT:
% Solves polynomial equation
%                     p(1)*x^n+...+p(n+1)=0
% determining the solutions in the interval (t0,t1).
%
% INPUT:
% p: vector determining the polynomial "p(1)*x^n+...+p(n+1)=0".
% t0,t1: the solutions are founf in the interval "[t0,t1]".
%
% OUTPUT:
% sols: vector of solutions.
% flag: 0: finite complex solutions,
%       Inf: infinite complex solutions (problem p=[0 ... 0]).
%       NaN: no complex solutions (problem p=[0 ... 0 c] with "c" not 0).
%
% Checked: May 2, 2019.

% .................... troubleshooting and defaults .......................

if nargin < 3
    t2=+inf;
end

if nargin < 2
    t1=-inf;
end

% ........................ computing roots ................................

sols=roots(p);

if isempty(sols) == 0
    % ..... In this case the equation has complex roots .....
    isols=imag(sols); rsols=real(sols);
    isols_in=find(abs(isols)< 10^(-12) & rsols >= t1 & sols <= t2);
    sols=rsols(isols_in);
    flag=0;
else
    % ... In this case the equation has infinite or no complex roots ......
    if p(end) == 0 % infinite complex roots
        flag=Inf;
    else
        flag=NaN; % no complex roots
    end
end







%--------------------------------------------------------------------------
% singular_vertices_v8a.
%--------------------------------------------------------------------------

function singular_junctions=singular_vertices_v8a(Sx,Sy)

% OBJECT:
%
%
% INPUT:
%
%
% OUTPUT:
%
%

singular_junctions=[];
derivatives_at_endpoints=[];

for ii=1:length(Sx)
    SxL=Sx(ii);
    SyL=Sy(ii);
    nodes=SxL.breaks;
    derivatives_at_endpointsL=[];
    switch SxL.order
        case 2
            
            Nnodes=length(nodes);
            t0=nodes(1:end-1);
            t1=nodes(2:end);
            SxL_nodes=ppval(SxL,nodes);
            x0=SxL_nodes(1:end-1);
            x1=SxL_nodes(2:end);
            m=(x1-x0)./(t1-t0);
            
            % find turning points (can be improved)
            if length(m) > 1
                
                mend=m(2:end); minit=m(1:end-1);
                ichangeder=find(mend.*minit <= 0);
                inodes_ch=ichangeder+1;
                singular_pointsT=nodes(inodes_ch);
                % singular_pointsX=ppval(SxL,singular_pointsT);
                singular_pointsX=SxL_nodes(inodes_ch);
                singular_pointsY=ppval(SyL,singular_pointsT);
                singular_junctions=[singular_junctions; singular_pointsX' ...
                    singular_pointsY' singular_pointsT'];
                
            end
            
            derivatives_at_endpointsL=[m(1) t0(1); ...
                m(end) t1(end)];
            
        otherwise
            
            
            Sx1L=fnder(SxL);
            Sx1nodes=ppval(Sx1L,nodes);
            nodes_in=nodes(2:end-1);
            Sx1nodes_in=Sx1nodes(2:end-1);
            
            inullders=find(abs(Sx1nodes_in) <= 10^(-14));
            
            if min(size(inullders)) > 0
                singular_pointsT=nodes_in(inullders);
                singular_pointsX=ppval(SxL,singular_pointsT);
                singular_pointsY=ppval(SyL,singular_pointsT);
                singular_junctions=[singular_junctions; singular_pointsX' ...
                    singular_pointsY' singular_pointsT'];
            end
            
            Dinit=Sx1nodes(1);
            Dend=Sx1nodes(end);
            derivatives_at_endpointsL=[Dinit nodes(1); Dend nodes(end)];
            
    end
    derivatives_at_endpoints=[derivatives_at_endpoints; ...
        derivatives_at_endpointsL];
end


%derivatives_at_endpoints

% changes at spline junctions.
DeL=derivatives_at_endpoints(2:2:end-2,1); % from left
DeR=derivatives_at_endpoints(3:2:end-1,1); % from right
ichangeder=find(DeL.*DeR <= 0);
LL=length(ichangeder);

if LL >= 1
    for kk=1:LL
        SxL=Sx(ichangeder(kk));
        SyL=Sy(ichangeder(kk));
        nodesL=SxL.breaks;
        singular_pointsTJ(kk,1)=nodesL(end);
        singular_pointsXJ(kk,1)=ppval(SxL,singular_pointsTJ(kk,1));
        singular_pointsYJ(kk,1)=ppval(SyL,singular_pointsTJ(kk,1));
    end
else
    singular_pointsXJ=[]; singular_pointsYJ=[]; singular_pointsTJ=[];
end


% .... final and initial point ....
temp=derivatives_at_endpoints(1,1)* derivatives_at_endpoints(end,1);

ichangeder1=find(temp <= 0);
if isempty(ichangeder1) == 0
    SxL=Sx(1);
    SyL=Sy(1);
    nodesL=SxL.breaks;
    singular_pointsTJ(end+1,1)=nodesL(1);
    singular_pointsXJ(end+1,1)=ppval(SxL,nodesL(1));
    singular_pointsYJ(end+1,1)=ppval(SyL,nodesL(1));
end

singular_junctions=[singular_junctions;
    singular_pointsXJ singular_pointsYJ singular_pointsTJ];




%--------------------------------------------------------------------------
% spline_boxer (and subroutines)
%--------------------------------------------------------------------------

function [boxVx,turning_points_Sx,turning_points_Sy]=spline_boxer(...
    Sx,Sy,Nsub)

% OBJECT:
% Given a sequence of splines "Sx", "Sy" representing the boundary of
% the domain, this routine define its boxes.
% If "t(k)" and "t(k+1)" are two subsequent breaks of some "Sx", it builds
% a rectangle (i.e. a box) [x0,x1] x [y0,y1] , where
%    x0=min(Sx(t))   x1=max(Sx(t))  y0=min(Sy(t))   y1=max(Sy(t))
% It collects all this boxes in "boxVx".
%
% Turning points "turning_points_Sx" are points where the derivative of
% "Sx" in "t(k)" and "t(k+1)" is null.
%
% The same applies to "boxVy", "turning_points_Sy" for the spline "Sy".
%
% INPUT:
%
%
% OUTPUT:
%
%
% Version: 8a.

if nargin < 3
    Nsub=10;
end

boxVx=[];
turning_points_Sx=[]; % N x 3 (turning points coor. (x,y,t), in this order)
turning_points_Sy=[];

for isp=1:length(Sx)
    
    SxL=Sx(isp);
    SyL=Sy(isp);
    
    [box_VxL,turning_points_SxL,turning_points_SyL]=spline_boxer0(...
        SxL,SyL,isp,Nsub);
    
    boxVx=[boxVx; box_VxL];
    turning_points_Sx=[turning_points_Sx; turning_points_SxL];
    turning_points_Sy=[turning_points_Sy; turning_points_SyL];
    
end




function [box_VxL,turning_points_SxL,turning_points_SyL]=spline_boxer0(...
    SxL,SyL,isp,Nsub)


% OBJECT:
%
%
% INPUT:
%
%
% OUTPUT:
%
%

% Technical note.
% spline derivative is stored as a polynomial "q'(t)=p'(t+a)" where "a"
% is the initial point of a block.

if SxL.order == 2
    box_VxL=spline_boxer_order2(SxL,SyL,isp,Nsub);
    turning_points_SxL=[]; turning_points_SyL=[];
else
    [box_VxL,turning_points_SxL,turning_points_SyL]=...
        spline_boxer_ordernot2(SxL,SyL,isp,Nsub);
end




%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

function box_VxL=spline_boxer_order2(SxL,SyL,isp,Nsub)

% OBJECT:
%
%
% INPUT:
%
%
% OUTPUT:
%
%

derivatives_infos=0;

nodes=SxL.breaks;
box_VxL=[];

% analysis of the ii-th spline component.
for ibl=1:length(nodes)-1
    
    tm=nodes(ibl); tM=nodes(ibl+1);
    tt=linspace(tm,tM,Nsub+1); tt=tt';
    
    % spline is stored as a polynomial "q(t)=p(t+a)".
    xvalues=fnval(SxL,tt);
    xvaluesM=[xvalues(1:end-1) xvalues(2:end)];
    xm=min(xvaluesM,[],2); xM=max(xvaluesM,[],2);
    
    yvalues=fnval(SyL,tt);
    yvaluesM=[yvalues(1:end-1) yvalues(2:end)];
    ym=min(yvaluesM,[],2); yM=max(yvaluesM,[],2);
    
    ttm=tt(1:end-1); ttM=tt(2:end);
    
    ss=ones(size(xm));
    
    
    % technical:
    % if SxL.coeffs(ibl,1)=[a b] then x(t)=a*t+b with t in [0,tM-tm]
    if derivatives_infos == 1
        der_signx=sign(SxL.coefs(ibl,1)); der_signy=sign(SyL.coefs(ibl,1));
        boxL=[xm xM ym yM ttm ttM isp*ss ibl*ss der_signx*ss der_signy*ss];
    else
        boxL=[xm xM ym yM ttm ttM isp*ss ibl*ss];
    end
    
    box_VxL=[box_VxL; boxL];
    
end




%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

function [box_VxL,turning_points_SxL,turning_points_SyL]=...
    spline_boxer_ordernot2(SxL,SyL,isp,Nsub)

% OBJECT.
% This function builds the "boxes" for a piece of spline of order not 2.
% Each block, determined by two successive breaks and possible turning
% points in the "x" direction (i.e. in which it is null the x derivative),
% is subdivided into "Nsub" intervals.
%
% INPUTS:
% SxL, SyL: splines determining locally the boundary of the spline.
% isp: the splines SxL, SyL are usually a component of a vector of spline
%      say Sx, Sy; in this case "SxL=Sx(isp)", "SyL=Sy(isp)".
% Nsub: number of subdvisions of each subinterval determined by two
%       successive breaks and x-turning points.
%
% OUTPUTS:
% box_VxL: each box contains the variables:
%           boxL=[preboxVxL isp ibl der_signx der_signy];
%          where "prebox" is a row vector of the form "[a b c d]", that
%          determines the boundary of the box, whose vertices are
%          (a,c), (b,c), (b,d), (a,d), "isp" is the input variable
%          described above, "ibl" determines that the "t" variable is
%          included in an interval "[t0,t1]" that is a subset of
%          "[T(ibl),T(ibl+1)]" where "T(k)" are the ordered spline breaks;
%          in each box the sign of the "x", "y" derivatives is constant and
%          their sign is stored in "der_signx", "der_signy".
% turning_points_SxL: (x,y) points of the interval "[T(1),T(end)]" where
%          "T(k)" are the ordered "SxL" spline breaks in which the
%          derivative of "SxL" is null.
% turning_points_SyL: (x,y) points of the interval "[T(1),T(end)]" where
%          "T(k)" are the ordered "SyL" spline breaks in which the
%          derivative of "SyL" is null.

derivatives_infos=0;

SxL1=fnder(SxL); SxL1_coeffs=SxL1.coefs; % derivative "x" direction
SyL1=fnder(SyL); SyL1_coeffs=SyL1.coefs; % derivative "y" direction

nodes=SxL.breaks;

% initialisation
box_VxL=[]; turning_points_SxL=[]; turning_points_SyL=[];

% analysis of the ibl-th spline component.

for ibl=1:length(nodes)-1
    
    % .................. retrieving data ..................
    
    a=nodes(ibl); b=nodes(ibl+1);
    
    % spline is stored as a polynomial "q(t)=p(t+a)".
    t1=0; t2=b-a;
    
    % spline coefficients
    
    SxL_coeffsL=SxL.coefs(ibl,:); SyL_coeffsL=SyL.coefs(ibl,:);
    SxL1_coeffsL=SxL1_coeffs(ibl,:); SyL1_coeffsL=SyL1_coeffs(ibl,:);
    
    % .................. turning points analysis ..................
    
    turning_points_xL=turning_points_analysis(SxL_coeffsL,...
        SyL_coeffsL,SxL1_coeffsL,a,t1,t2);
    
    turning_points_yL=turning_points_analysis(SxL_coeffsL,...
        SyL_coeffsL,SyL1_coeffsL,a,t1,t2);
    
    turning_points_SxL=[turning_points_SxL; turning_points_xL];
    turning_points_SyL=[turning_points_SyL; turning_points_yL];
    
    % ..................  making monotone boxes ..................
    
    % determining box (or boxes) in intervals relative to "a", "b" and
    % "turning points", so that in each box the spline is a monotone
    % function.
    
    t12=linspace(t1,t2,Nsub+1); t12=t12';
    SxLvalues_ab=polyval(SxL_coeffsL,t12);
    SyLvalues_ab=polyval(SyL_coeffsL,t12);
    SxyLvalues_ab_ext=[SxLvalues_ab SyLvalues_ab a+t12];
    
    SxLvalues=[SxyLvalues_ab_ext; turning_points_xL; turning_points_yL];
    
    % prebox [xm xM ym yM tm tM] (can be matrix, depending on the sizes of
    % turning_points_xL, turning_points_yL.
    prebox_VxL0=sortrows(SxLvalues,3);
    
    % derivatives in each box.
    tt=prebox_VxL0(:,3);
    
    % boxes coordinates
    Xv=prebox_VxL0(:,1); Yv=prebox_VxL0(:,2); Tv=prebox_VxL0(:,3);
    Xm=min(Xv(1:end-1),Xv(2:end)); XM=max(Xv(1:end-1),Xv(2:end));
    Ym=min(Yv(1:end-1),Yv(2:end)); YM=max(Yv(1:end-1),Yv(2:end));
    preboxVxL=[Xm XM Ym YM Tv(1:end-1) Tv(2:end)];
    
    
    if derivatives_infos == 1
        [der_signx,der_signy]=derivatives_sign_monotone_box(...
            SxL1_coeffsL,SyL1_coeffsL,tt);
        ss=ones(size(der_signx));
        boxL=[preboxVxL isp*ss ibl*ss der_signx der_signy];
    else
        ss=ones(size(preboxVxL,1),1);
        boxL=[preboxVxL isp*ss ibl*ss];
    end
    
    box_VxL=[box_VxL; boxL];
    
end




%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

function turning_pointsL=turning_points_analysis(SxL_coeffsL,...
    SyL_coeffsL,SxL1_coeffsL,a,t1,t2)

% OBJECT:
%
%
% INPUT:
%
%
% OUTPUT:
%
%

[tR,flag]=myroots(SxL1_coeffsL,t1,t2);

if flag == 0
    if length(tR) > 0
        SxLvalues_tR=polyval(SxL_coeffsL,tR);
        SyLvalues_tR=polyval(SyL_coeffsL,tR);
        StLvalues_tR=a+tR;
        turning_pointsL=[SxLvalues_tR SyLvalues_tR StLvalues_tR];
    else
        tR=[t1; t2];
        SxLvalues_tR=polyval(SxL_coeffsL,tR);
        SyLvalues_tR=polyval(SyL_coeffsL,tR);
        StLvalues_tR=a+tR;
        turning_pointsL=[SxLvalues_tR SyLvalues_tR StLvalues_tR];
    end
else
    if flag == Inf % Analysis: infinite turning points.
        % Geometrical view: vertical segment.
        tR=[t1; t2];
        SxLvalues_tR=polyval(SxL_coeffsL,tR);
        SyLvalues_tR=polyval(SyL_coeffsL,tR);
        StLvalues_tR=a+tR;
        turning_pointsL=[SxLvalues_tR SyLvalues_tR StLvalues_tR];
    else % Analysis: no turning points.
        turning_pointsL=[];
    end
end



%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

function [der_signx,der_signy]=derivatives_sign_monotone_box(...
    SxL1_coeffsL,SyL1_coeffsL,tt)

% OBJECT:
%
%
% INPUT:
%
%
% OUTPUT:
%
%

ttmed=(tt(2:end)-tt(1:end-1))/2;
der_signx=sign(polyval(SxL1_coeffsL,ttmed));
der_signy=sign(polyval(SyL1_coeffsL,ttmed));




%--------------------------------------------------------------------------
% compute_parametric_spline
%--------------------------------------------------------------------------

function [ppx,ppy]=compute_parametric_spline(s,x,y,spline_order,...
    SPLtypestring)

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


switch spline_order
    case 4
        
        % CUBIC SPLINES BY CSAPE. AS DEFAULT WE USE PERIODIC CUBIC SPLINES.
        % Derivatives parameters are computed as well.
        ppx=csape(s,x,SPLtypestring);
        ppy=csape(s,y,SPLtypestring);
        
    otherwise
        
        ppx=spapi(spline_order,s,x);
        ppy=spapi(spline_order,s,y);
        
end

if (ppx.form =='B-')
    ppx=fn2fm(ppx,'pp');
    ppy=fn2fm(ppy,'pp');
end








%--------------------------------------------------------------------------
% winding_algorithm
%--------------------------------------------------------------------------

function windn = winding_algorithm(P,Sx,Sy)

% OBJECT:
% computes the winding number of P with respect to the counterclockwise
% concatened cubic spline arcs (Sx,Sy) forming a Jordan spline polygon,
% by contour integration
%
% INPUT:
% P: 2-column array of planar points
% Sx,Sy: arrays of cubic spline structures; (SX(i),Sy(i)) is the i-th arc
% the arcs must be counterclockwise concatenated forming a Jordan curve
% Sx(i).breaks(end)=Sx(i+1).breaks(1),Sy(i).breaks(end)=Sy(i+1).breaks(1)
% i=1,...,end, with end+1:=1
%
% OUTPUT:
% windn: winding numbers of P with respect to the Jordan curve (Sx,Sy)
% windn(j)==1 point is inside, windn(j)==0 point is outside
%
% DATA:
% built: April 2019
% checked: May 2, 2019
%

windn=zeros(1,length(P(:,1)));

xw=GL_1000;

for i=1:length(Sx)
    
    a=Sx(i).breaks(1);b=Sx(i).breaks(end);
    
    dSx(i)=fnder(Sx(i));
    dSy(i)=fnder(Sy(i));
    
    x=@(s) ppval(Sx(i),s);
    dx=@(s) ppval(dSx(i),s);
    y=@(s) ppval(Sy(i),s);
    dy=@(s) ppval(dSy(i),s);
    
    f=@(s,j) (dy(s).*(x(s)-P(j,1))-dx(s).*(y(s)-P(j,2)))./...
        ((x(s)-P(j,1)).^2+(y(s)-P(j,2)).^2);
    
    nodes=xw(:,1)*(b-a)/2+(b+a)/2;
    weights=xw(:,2)*(b-a)/2;
    
    [u,v]=meshgrid(nodes,1:length(P(:,1)));
    fval=f(u(:),v(:));
    FV=reshape(fval,length(P(:,1)),length(nodes));
    int=FV*weights;
    
    windn=windn+int';
    
end

den=Sx(end).breaks(end)-Sx(1).breaks(1);
% windn=round(windn/(2*pi));
windn=round(windn/den);




%==========================================================================
% attached function: GL_1000
%==========================================================================

function xw=GL_1000

x=[-9.999971112980756e-01
    -9.999847796329174e-01
    -9.999625941483602e-01
    -9.999305501355009e-01
    -9.998886473067012e-01
    -9.998368859309700e-01
    -9.997752664706340e-01
    -9.997037895136229e-01
    -9.996224557554706e-01
    -9.995312659933240e-01
    -9.994302211236081e-01
    -9.993193221410008e-01
    -9.991985701379369e-01
    -9.990679663043477e-01
    -9.989275119275132e-01
    -9.987772083919715e-01
    -9.986170571794587e-01
    -9.984470598688657e-01
    -9.982672181362039e-01
    -9.980775337545768e-01
    -9.978780085941545e-01
    -9.976686446221500e-01
    -9.974494439027950e-01
    -9.972204085973179e-01
    -9.969815409639196e-01
    -9.967328433577501e-01
    -9.964743182308847e-01
    -9.962059681322979e-01
    -9.959277957078384e-01
    -9.956398037002028e-01
    -9.953419949489071e-01
    -9.950343723902594e-01
    -9.947169390573305e-01
    -9.943896980799236e-01
    -9.940526526845432e-01
    -9.937058061943637e-01
    -9.933491620291958e-01
    -9.929827237054533e-01
    -9.926064948361183e-01
    -9.922204791307051e-01
    -9.918246803952242e-01
    -9.914191025321440e-01
    -9.910037495403531e-01
    -9.905786255151200e-01
    -9.901437346480535e-01
    -9.896990812270609e-01
    -9.892446696363055e-01
    -9.887805043561644e-01
    -9.883065899631827e-01
    -9.878229311300296e-01
    -9.873295326254522e-01
    -9.868263993142279e-01
    -9.863135361571168e-01
    -9.857909482108129e-01
    -9.852586406278939e-01
    -9.847166186567707e-01
    -9.841648876416356e-01
    -9.836034530224089e-01
    -9.830323203346869e-01
    -9.824514952096856e-01
    -9.818609833741858e-01
    -9.812607906504773e-01
    -9.806509229563004e-01
    -9.800313863047881e-01
    -9.794021868044072e-01
    -9.787633306588969e-01
    -9.781148241672090e-01
    -9.774566737234449e-01
    -9.767888858167928e-01
    -9.761114670314638e-01
    -9.754244240466269e-01
    -9.747277636363431e-01
    -9.740214926694986e-01
    -9.733056181097374e-01
    -9.725801470153922e-01
    -9.718450865394149e-01
    -9.711004439293064e-01
    -9.703462265270447e-01
    -9.695824417690128e-01
    -9.688090971859251e-01
    -9.680262004027537e-01
    -9.672337591386525e-01
    -9.664317812068817e-01
    -9.656202745147303e-01
    -9.647992470634386e-01
    -9.639687069481189e-01
    -9.631286623576756e-01
    -9.622791215747253e-01
    -9.614200929755139e-01
    -9.605515850298352e-01
    -9.596736063009463e-01
    -9.587861654454841e-01
    -9.578892712133795e-01
    -9.569829324477710e-01
    -9.560671580849179e-01
    -9.551419571541118e-01
    -9.542073387775878e-01
    -9.532633121704347e-01
    -9.523098866405036e-01
    -9.513470715883170e-01
    -9.503748765069749e-01
    -9.493933109820625e-01
    -9.484023846915548e-01
    -9.474021074057215e-01
    -9.463924889870307e-01
    -9.453735393900514e-01
    -9.443452686613558e-01
    -9.433076869394198e-01
    -9.422608044545232e-01
    -9.412046315286491e-01
    -9.401391785753814e-01
    -9.390644560998033e-01
    -9.379804746983922e-01
    -9.368872450589167e-01
    -9.357847779603303e-01
    -9.346730842726655e-01
    -9.335521749569263e-01
    -9.324220610649806e-01
    -9.312827537394508e-01
    -9.301342642136041e-01
    -9.289766038112420e-01
    -9.278097839465882e-01
    -9.266338161241762e-01
    -9.254487119387360e-01
    -9.242544830750798e-01
    -9.230511413079867e-01
    -9.218386985020864e-01
    -9.206171666117428e-01
    -9.193865576809351e-01
    -9.181468838431406e-01
    -9.168981573212134e-01
    -9.156403904272650e-01
    -9.143735955625426e-01
    -9.130977852173067e-01
    -9.118129719707079e-01
    -9.105191684906634e-01
    -9.092163875337312e-01
    -9.079046419449854e-01
    -9.065839446578886e-01
    -9.052543086941648e-01
    -9.039157471636710e-01
    -9.025682732642681e-01
    -9.012119002816904e-01
    -8.998466415894147e-01
    -8.984725106485290e-01
    -8.970895210075985e-01
    -8.956976863025338e-01
    -8.942970202564546e-01
    -8.928875366795559e-01
    -8.914692494689709e-01
    -8.900421726086343e-01
    -8.886063201691446e-01
    -8.871617063076249e-01
    -8.857083452675839e-01
    -8.842462513787749e-01
    -8.827754390570548e-01
    -8.812959228042421e-01
    -8.798077172079737e-01
    -8.783108369415609e-01
    -8.768052967638450e-01
    -8.752911115190520e-01
    -8.737682961366459e-01
    -8.722368656311812e-01
    -8.706968351021556e-01
    -8.691482197338604e-01
    -8.675910347952316e-01
    -8.660252956396985e-01
    -8.644510177050330e-01
    -8.628682165131966e-01
    -8.612769076701884e-01
    -8.596771068658904e-01
    -8.580688298739131e-01
    -8.564520925514401e-01
    -8.548269108390713e-01
    -8.531933007606665e-01
    -8.515512784231865e-01
    -8.499008600165348e-01
    -8.482420618133982e-01
    -8.465749001690859e-01
    -8.448993915213685e-01
    -8.432155523903155e-01
    -8.415233993781333e-01
    -8.398229491690001e-01
    -8.381142185289034e-01
    -8.363972243054727e-01
    -8.346719834278143e-01
    -8.329385129063448e-01
    -8.311968298326226e-01
    -8.294469513791792e-01
    -8.276888947993514e-01
    -8.259226774271093e-01
    -8.241483166768867e-01
    -8.223658300434088e-01
    -8.205752351015196e-01
    -8.187765495060093e-01
    -8.169697909914395e-01
    -8.151549773719688e-01
    -8.133321265411768e-01
    -8.115012564718881e-01
    -8.096623852159945e-01
    -8.078155309042779e-01
    -8.059607117462307e-01
    -8.040979460298764e-01
    -8.022272521215899e-01
    -8.003486484659154e-01
    -7.984621535853855e-01
    -7.965677860803383e-01
    -7.946655646287335e-01
    -7.927555079859688e-01
    -7.908376349846947e-01
    -7.889119645346292e-01
    -7.869785156223710e-01
    -7.850373073112119e-01
    -7.830883587409498e-01
    -7.811316891276996e-01
    -7.791673177637031e-01
    -7.771952640171397e-01
    -7.752155473319350e-01
    -7.732281872275689e-01
    -7.712332032988839e-01
    -7.692306152158909e-01
    -7.672204427235757e-01
    -7.652027056417049e-01
    -7.631774238646295e-01
    -7.611446173610890e-01
    -7.591043061740155e-01
    -7.570565104203343e-01
    -7.550012502907671e-01
    -7.529385460496325e-01
    -7.508684180346457e-01
    -7.487908866567183e-01
    -7.467059723997576e-01
    -7.446136958204634e-01
    -7.425140775481267e-01
    -7.404071382844252e-01
    -7.382928988032200e-01
    -7.361713799503499e-01
    -7.340426026434269e-01
    -7.319065878716291e-01
    -7.297633566954940e-01
    -7.276129302467115e-01
    -7.254553297279143e-01
    -7.232905764124696e-01
    -7.211186916442698e-01
    -7.189396968375211e-01
    -7.167536134765330e-01
    -7.145604631155059e-01
    -7.123602673783193e-01
    -7.101530479583182e-01
    -7.079388266180990e-01
    -7.057176251892955e-01
    -7.034894655723630e-01
    -7.012543697363630e-01
    -6.990123597187459e-01
    -6.967634576251345e-01
    -6.945076856291054e-01
    -6.922450659719704e-01
    -6.899756209625583e-01
    -6.876993729769928e-01
    -6.854163444584741e-01
    -6.831265579170562e-01
    -6.808300359294256e-01
    -6.785268011386784e-01
    -6.762168762540970e-01
    -6.739002840509267e-01
    -6.715770473701509e-01
    -6.692471891182650e-01
    -6.669107322670522e-01
    -6.645676998533557e-01
    -6.622181149788520e-01
    -6.598620008098232e-01
    -6.574993805769284e-01
    -6.551302775749750e-01
    -6.527547151626888e-01
    -6.503727167624830e-01
    -6.479843058602290e-01
    -6.455895060050229e-01
    -6.431883408089545e-01
    -6.407808339468741e-01
    -6.383670091561594e-01
    -6.359468902364808e-01
    -6.335205010495674e-01
    -6.310878655189712e-01
    -6.286490076298320e-01
    -6.262039514286397e-01
    -6.237527210229985e-01
    -6.212953405813885e-01
    -6.188318343329273e-01
    -6.163622265671316e-01
    -6.138865416336770e-01
    -6.114048039421588e-01
    -6.089170379618506e-01
    -6.064232682214638e-01
    -6.039235193089045e-01
    -6.014178158710324e-01
    -5.989061826134173e-01
    -5.963886443000951e-01
    -5.938652257533241e-01
    -5.913359518533402e-01
    -5.888008475381117e-01
    -5.862599378030930e-01
    -5.837132477009781e-01
    -5.811608023414545e-01
    -5.786026268909548e-01
    -5.760387465724085e-01
    -5.734691866649939e-01
    -5.708939725038883e-01
    -5.683131294800189e-01
    -5.657266830398111e-01
    -5.631346586849395e-01
    -5.605370819720750e-01
    -5.579339785126332e-01
    -5.553253739725217e-01
    -5.527112940718880e-01
    -5.500917645848646e-01
    -5.474668113393159e-01
    -5.448364602165826e-01
    -5.422007371512277e-01
    -5.395596681307794e-01
    -5.369132791954765e-01
    -5.342615964380096e-01
    -5.316046460032660e-01
    -5.289424540880701e-01
    -5.262750469409269e-01
    -5.236024508617610e-01
    -5.209246922016596e-01
    -5.182417973626101e-01
    -5.155537927972434e-01
    -5.128607050085692e-01
    -5.101625605497169e-01
    -5.074593860236734e-01
    -5.047512080830209e-01
    -5.020380534296736e-01
    -4.993199488146152e-01
    -4.965969210376339e-01
    -4.938689969470604e-01
    -4.911362034395000e-01
    -4.883985674595709e-01
    -4.856561159996345e-01
    -4.829088760995346e-01
    -4.801568748463252e-01
    -4.774001393740068e-01
    -4.746386968632583e-01
    -4.718725745411689e-01
    -4.691017996809685e-01
    -4.663263996017608e-01
    -4.635464016682530e-01
    -4.607618332904848e-01
    -4.579727219235606e-01
    -4.551790950673760e-01
    -4.523809802663499e-01
    -4.495784051091496e-01
    -4.467713972284212e-01
    -4.439599843005154e-01
    -4.411441940452169e-01
    -4.383240542254689e-01
    -4.354995926470994e-01
    -4.326708371585488e-01
    -4.298378156505943e-01
    -4.270005560560748e-01
    -4.241590863496144e-01
    -4.213134345473503e-01
    -4.184636287066516e-01
    -4.156096969258460e-01
    -4.127516673439428e-01
    -4.098895681403530e-01
    -4.070234275346142e-01
    -4.041532737861100e-01
    -4.012791351937939e-01
    -3.984010400959077e-01
    -3.955190168697040e-01
    -3.926330939311652e-01
    -3.897432997347245e-01
    -3.868496627729833e-01
    -3.839522115764340e-01
    -3.810509747131738e-01
    -3.781459807886273e-01
    -3.752372584452621e-01
    -3.723248363623067e-01
    -3.694087432554689e-01
    -3.664890078766510e-01
    -3.635656590136674e-01
    -3.606387254899606e-01
    -3.577082361643171e-01
    -3.547742199305824e-01
    -3.518367057173760e-01
    -3.488957224878070e-01
    -3.459512992391884e-01
    -3.430034650027503e-01
    -3.400522488433541e-01
    -3.370976798592069e-01
    -3.341397871815732e-01
    -3.311785999744878e-01
    -3.282141474344693e-01
    -3.252464587902320e-01
    -3.222755633023960e-01
    -3.193014902632018e-01
    -3.163242689962176e-01
    -3.133439288560541e-01
    -3.103604992280731e-01
    -3.073740095280963e-01
    -3.043844892021185e-01
    -3.013919677260152e-01
    -2.983964746052525e-01
    -2.953980393745963e-01
    -2.923966915978200e-01
    -2.893924608674150e-01
    -2.863853768042977e-01
    -2.833754690575169e-01
    -2.803627673039623e-01
    -2.773473012480731e-01
    -2.743291006215421e-01
    -2.713081951830252e-01
    -2.682846147178464e-01
    -2.652583890377055e-01
    -2.622295479803826e-01
    -2.591981214094458e-01
    -2.561641392139540e-01
    -2.531276313081655e-01
    -2.500886276312411e-01
    -2.470471581469482e-01
    -2.440032528433684e-01
    -2.409569417325970e-01
    -2.379082548504527e-01
    -2.348572222561768e-01
    -2.318038740321401e-01
    -2.287482402835438e-01
    -2.256903511381238e-01
    -2.226302367458543e-01
    -2.195679272786493e-01
    -2.165034529300667e-01
    -2.134368439150074e-01
    -2.103681304694220e-01
    -2.072973428500084e-01
    -2.042245113339161e-01
    -2.011496662184466e-01
    -1.980728378207556e-01
    -1.949940564775524e-01
    -1.919133525448032e-01
    -1.888307563974286e-01
    -1.857462984290076e-01
    -1.826600090514756e-01
    -1.795719186948246e-01
    -1.764820578068050e-01
    -1.733904568526228e-01
    -1.702971463146424e-01
    -1.672021566920822e-01
    -1.641055185007176e-01
    -1.610072622725773e-01
    -1.579074185556439e-01
    -1.548060179135523e-01
    -1.517030909252879e-01
    -1.485986681848864e-01
    -1.454927803011297e-01
    -1.423854578972470e-01
    -1.392767316106107e-01
    -1.361666320924357e-01
    -1.330551900074751e-01
    -1.299424360337224e-01
    -1.268284008621025e-01
    -1.237131151961749e-01
    -1.205966097518273e-01
    -1.174789152569755e-01
    -1.143600624512574e-01
    -1.112400820857331e-01
    -1.081190049225787e-01
    -1.049968617347854e-01
    -1.018736833058549e-01
    -9.874950042949604e-02
    -9.562434390932102e-02
    -9.249824455854230e-02
    -8.937123319966822e-02
    -8.624334066419928e-02
    -8.311459779232394e-02
    -7.998503543261501e-02
    -7.685468444172568e-02
    -7.372357568408372e-02
    -7.059174003158945e-02
    -6.745920836330904e-02
    -6.432601156517212e-02
    -6.119218052966611e-02
    -5.805774615553117e-02
    -5.492273934745737e-02
    -5.178719101577792e-02
    -4.865113207616557e-02
    -4.551459344932774e-02
    -4.237760606070101e-02
    -3.924020084014754e-02
    -3.610240872164740e-02
    -3.296426064299765e-02
    -2.982578754550276e-02
    -2.668702037367428e-02
    -2.354799007492080e-02
    -2.040872759924663e-02
    -1.726926389894526e-02
    -1.412962992829401e-02
    -1.098985664324877e-02
    -7.849975001139075e-03
    -4.710015960364133e-03
    -1.570010480083192e-03
    1.570010480083192e-03
    4.710015960364133e-03
    7.849975001139075e-03
    1.098985664324877e-02
    1.412962992829401e-02
    1.726926389894526e-02
    2.040872759924663e-02
    2.354799007492080e-02
    2.668702037367428e-02
    2.982578754550276e-02
    3.296426064299765e-02
    3.610240872164740e-02
    3.924020084014754e-02
    4.237760606070101e-02
    4.551459344932774e-02
    4.865113207616557e-02
    5.178719101577792e-02
    5.492273934745737e-02
    5.805774615553117e-02
    6.119218052966611e-02
    6.432601156517212e-02
    6.745920836330904e-02
    7.059174003158945e-02
    7.372357568408372e-02
    7.685468444172568e-02
    7.998503543261501e-02
    8.311459779232394e-02
    8.624334066419928e-02
    8.937123319966822e-02
    9.249824455854230e-02
    9.562434390932102e-02
    9.874950042949604e-02
    1.018736833058549e-01
    1.049968617347854e-01
    1.081190049225787e-01
    1.112400820857331e-01
    1.143600624512574e-01
    1.174789152569755e-01
    1.205966097518273e-01
    1.237131151961749e-01
    1.268284008621025e-01
    1.299424360337224e-01
    1.330551900074751e-01
    1.361666320924357e-01
    1.392767316106107e-01
    1.423854578972470e-01
    1.454927803011297e-01
    1.485986681848864e-01
    1.517030909252879e-01
    1.548060179135523e-01
    1.579074185556439e-01
    1.610072622725773e-01
    1.641055185007176e-01
    1.672021566920822e-01
    1.702971463146424e-01
    1.733904568526228e-01
    1.764820578068050e-01
    1.795719186948246e-01
    1.826600090514756e-01
    1.857462984290076e-01
    1.888307563974286e-01
    1.919133525448032e-01
    1.949940564775524e-01
    1.980728378207556e-01
    2.011496662184466e-01
    2.042245113339161e-01
    2.072973428500084e-01
    2.103681304694220e-01
    2.134368439150074e-01
    2.165034529300667e-01
    2.195679272786493e-01
    2.226302367458543e-01
    2.256903511381238e-01
    2.287482402835438e-01
    2.318038740321401e-01
    2.348572222561768e-01
    2.379082548504527e-01
    2.409569417325970e-01
    2.440032528433684e-01
    2.470471581469482e-01
    2.500886276312411e-01
    2.531276313081655e-01
    2.561641392139540e-01
    2.591981214094458e-01
    2.622295479803826e-01
    2.652583890377055e-01
    2.682846147178464e-01
    2.713081951830252e-01
    2.743291006215421e-01
    2.773473012480731e-01
    2.803627673039623e-01
    2.833754690575169e-01
    2.863853768042977e-01
    2.893924608674150e-01
    2.923966915978200e-01
    2.953980393745963e-01
    2.983964746052525e-01
    3.013919677260152e-01
    3.043844892021185e-01
    3.073740095280963e-01
    3.103604992280731e-01
    3.133439288560541e-01
    3.163242689962176e-01
    3.193014902632018e-01
    3.222755633023960e-01
    3.252464587902320e-01
    3.282141474344693e-01
    3.311785999744878e-01
    3.341397871815732e-01
    3.370976798592069e-01
    3.400522488433541e-01
    3.430034650027503e-01
    3.459512992391884e-01
    3.488957224878070e-01
    3.518367057173760e-01
    3.547742199305824e-01
    3.577082361643171e-01
    3.606387254899606e-01
    3.635656590136674e-01
    3.664890078766510e-01
    3.694087432554689e-01
    3.723248363623067e-01
    3.752372584452621e-01
    3.781459807886273e-01
    3.810509747131738e-01
    3.839522115764340e-01
    3.868496627729833e-01
    3.897432997347245e-01
    3.926330939311652e-01
    3.955190168697040e-01
    3.984010400959077e-01
    4.012791351937939e-01
    4.041532737861100e-01
    4.070234275346142e-01
    4.098895681403530e-01
    4.127516673439428e-01
    4.156096969258460e-01
    4.184636287066516e-01
    4.213134345473503e-01
    4.241590863496144e-01
    4.270005560560748e-01
    4.298378156505943e-01
    4.326708371585488e-01
    4.354995926470994e-01
    4.383240542254689e-01
    4.411441940452169e-01
    4.439599843005154e-01
    4.467713972284212e-01
    4.495784051091496e-01
    4.523809802663499e-01
    4.551790950673760e-01
    4.579727219235606e-01
    4.607618332904848e-01
    4.635464016682530e-01
    4.663263996017608e-01
    4.691017996809685e-01
    4.718725745411689e-01
    4.746386968632583e-01
    4.774001393740068e-01
    4.801568748463252e-01
    4.829088760995346e-01
    4.856561159996345e-01
    4.883985674595709e-01
    4.911362034395000e-01
    4.938689969470604e-01
    4.965969210376339e-01
    4.993199488146152e-01
    5.020380534296736e-01
    5.047512080830209e-01
    5.074593860236734e-01
    5.101625605497169e-01
    5.128607050085692e-01
    5.155537927972434e-01
    5.182417973626101e-01
    5.209246922016596e-01
    5.236024508617610e-01
    5.262750469409269e-01
    5.289424540880701e-01
    5.316046460032660e-01
    5.342615964380096e-01
    5.369132791954765e-01
    5.395596681307794e-01
    5.422007371512277e-01
    5.448364602165826e-01
    5.474668113393159e-01
    5.500917645848646e-01
    5.527112940718880e-01
    5.553253739725217e-01
    5.579339785126332e-01
    5.605370819720750e-01
    5.631346586849395e-01
    5.657266830398111e-01
    5.683131294800189e-01
    5.708939725038883e-01
    5.734691866649939e-01
    5.760387465724085e-01
    5.786026268909548e-01
    5.811608023414545e-01
    5.837132477009781e-01
    5.862599378030930e-01
    5.888008475381117e-01
    5.913359518533402e-01
    5.938652257533241e-01
    5.963886443000951e-01
    5.989061826134173e-01
    6.014178158710324e-01
    6.039235193089045e-01
    6.064232682214638e-01
    6.089170379618506e-01
    6.114048039421588e-01
    6.138865416336770e-01
    6.163622265671316e-01
    6.188318343329273e-01
    6.212953405813885e-01
    6.237527210229985e-01
    6.262039514286397e-01
    6.286490076298320e-01
    6.310878655189712e-01
    6.335205010495674e-01
    6.359468902364808e-01
    6.383670091561594e-01
    6.407808339468741e-01
    6.431883408089545e-01
    6.455895060050229e-01
    6.479843058602290e-01
    6.503727167624830e-01
    6.527547151626888e-01
    6.551302775749750e-01
    6.574993805769284e-01
    6.598620008098232e-01
    6.622181149788520e-01
    6.645676998533557e-01
    6.669107322670522e-01
    6.692471891182650e-01
    6.715770473701509e-01
    6.739002840509267e-01
    6.762168762540970e-01
    6.785268011386784e-01
    6.808300359294256e-01
    6.831265579170562e-01
    6.854163444584741e-01
    6.876993729769928e-01
    6.899756209625583e-01
    6.922450659719704e-01
    6.945076856291054e-01
    6.967634576251345e-01
    6.990123597187459e-01
    7.012543697363630e-01
    7.034894655723630e-01
    7.057176251892955e-01
    7.079388266180990e-01
    7.101530479583182e-01
    7.123602673783193e-01
    7.145604631155059e-01
    7.167536134765330e-01
    7.189396968375211e-01
    7.211186916442698e-01
    7.232905764124696e-01
    7.254553297279143e-01
    7.276129302467115e-01
    7.297633566954940e-01
    7.319065878716291e-01
    7.340426026434269e-01
    7.361713799503499e-01
    7.382928988032200e-01
    7.404071382844252e-01
    7.425140775481267e-01
    7.446136958204634e-01
    7.467059723997576e-01
    7.487908866567183e-01
    7.508684180346457e-01
    7.529385460496325e-01
    7.550012502907671e-01
    7.570565104203343e-01
    7.591043061740155e-01
    7.611446173610890e-01
    7.631774238646295e-01
    7.652027056417049e-01
    7.672204427235757e-01
    7.692306152158909e-01
    7.712332032988839e-01
    7.732281872275689e-01
    7.752155473319350e-01
    7.771952640171397e-01
    7.791673177637031e-01
    7.811316891276996e-01
    7.830883587409498e-01
    7.850373073112119e-01
    7.869785156223710e-01
    7.889119645346292e-01
    7.908376349846947e-01
    7.927555079859688e-01
    7.946655646287335e-01
    7.965677860803383e-01
    7.984621535853855e-01
    8.003486484659154e-01
    8.022272521215899e-01
    8.040979460298764e-01
    8.059607117462307e-01
    8.078155309042779e-01
    8.096623852159945e-01
    8.115012564718881e-01
    8.133321265411768e-01
    8.151549773719688e-01
    8.169697909914395e-01
    8.187765495060093e-01
    8.205752351015196e-01
    8.223658300434088e-01
    8.241483166768867e-01
    8.259226774271093e-01
    8.276888947993514e-01
    8.294469513791792e-01
    8.311968298326226e-01
    8.329385129063448e-01
    8.346719834278143e-01
    8.363972243054727e-01
    8.381142185289034e-01
    8.398229491690001e-01
    8.415233993781333e-01
    8.432155523903155e-01
    8.448993915213685e-01
    8.465749001690859e-01
    8.482420618133982e-01
    8.499008600165348e-01
    8.515512784231865e-01
    8.531933007606665e-01
    8.548269108390713e-01
    8.564520925514401e-01
    8.580688298739131e-01
    8.596771068658904e-01
    8.612769076701884e-01
    8.628682165131966e-01
    8.644510177050330e-01
    8.660252956396985e-01
    8.675910347952316e-01
    8.691482197338604e-01
    8.706968351021556e-01
    8.722368656311812e-01
    8.737682961366459e-01
    8.752911115190520e-01
    8.768052967638450e-01
    8.783108369415609e-01
    8.798077172079737e-01
    8.812959228042421e-01
    8.827754390570548e-01
    8.842462513787749e-01
    8.857083452675839e-01
    8.871617063076249e-01
    8.886063201691446e-01
    8.900421726086343e-01
    8.914692494689709e-01
    8.928875366795559e-01
    8.942970202564546e-01
    8.956976863025338e-01
    8.970895210075985e-01
    8.984725106485290e-01
    8.998466415894147e-01
    9.012119002816904e-01
    9.025682732642681e-01
    9.039157471636710e-01
    9.052543086941648e-01
    9.065839446578886e-01
    9.079046419449854e-01
    9.092163875337312e-01
    9.105191684906634e-01
    9.118129719707079e-01
    9.130977852173067e-01
    9.143735955625426e-01
    9.156403904272650e-01
    9.168981573212134e-01
    9.181468838431406e-01
    9.193865576809351e-01
    9.206171666117428e-01
    9.218386985020864e-01
    9.230511413079867e-01
    9.242544830750798e-01
    9.254487119387360e-01
    9.266338161241762e-01
    9.278097839465882e-01
    9.289766038112420e-01
    9.301342642136041e-01
    9.312827537394508e-01
    9.324220610649806e-01
    9.335521749569263e-01
    9.346730842726655e-01
    9.357847779603303e-01
    9.368872450589167e-01
    9.379804746983922e-01
    9.390644560998033e-01
    9.401391785753814e-01
    9.412046315286491e-01
    9.422608044545232e-01
    9.433076869394198e-01
    9.443452686613558e-01
    9.453735393900514e-01
    9.463924889870307e-01
    9.474021074057215e-01
    9.484023846915548e-01
    9.493933109820625e-01
    9.503748765069749e-01
    9.513470715883170e-01
    9.523098866405036e-01
    9.532633121704347e-01
    9.542073387775878e-01
    9.551419571541118e-01
    9.560671580849179e-01
    9.569829324477710e-01
    9.578892712133795e-01
    9.587861654454841e-01
    9.596736063009463e-01
    9.605515850298352e-01
    9.614200929755139e-01
    9.622791215747253e-01
    9.631286623576756e-01
    9.639687069481189e-01
    9.647992470634386e-01
    9.656202745147303e-01
    9.664317812068817e-01
    9.672337591386525e-01
    9.680262004027537e-01
    9.688090971859251e-01
    9.695824417690128e-01
    9.703462265270447e-01
    9.711004439293064e-01
    9.718450865394149e-01
    9.725801470153922e-01
    9.733056181097374e-01
    9.740214926694986e-01
    9.747277636363431e-01
    9.754244240466269e-01
    9.761114670314638e-01
    9.767888858167928e-01
    9.774566737234449e-01
    9.781148241672090e-01
    9.787633306588969e-01
    9.794021868044072e-01
    9.800313863047881e-01
    9.806509229563004e-01
    9.812607906504773e-01
    9.818609833741858e-01
    9.824514952096856e-01
    9.830323203346869e-01
    9.836034530224089e-01
    9.841648876416356e-01
    9.847166186567707e-01
    9.852586406278939e-01
    9.857909482108129e-01
    9.863135361571168e-01
    9.868263993142279e-01
    9.873295326254522e-01
    9.878229311300296e-01
    9.883065899631827e-01
    9.887805043561644e-01
    9.892446696363055e-01
    9.896990812270609e-01
    9.901437346480535e-01
    9.905786255151200e-01
    9.910037495403531e-01
    9.914191025321440e-01
    9.918246803952242e-01
    9.922204791307051e-01
    9.926064948361183e-01
    9.929827237054533e-01
    9.933491620291958e-01
    9.937058061943637e-01
    9.940526526845432e-01
    9.943896980799236e-01
    9.947169390573305e-01
    9.950343723902594e-01
    9.953419949489071e-01
    9.956398037002028e-01
    9.959277957078384e-01
    9.962059681322979e-01
    9.964743182308847e-01
    9.967328433577501e-01
    9.969815409639196e-01
    9.972204085973179e-01
    9.974494439027950e-01
    9.976686446221500e-01
    9.978780085941545e-01
    9.980775337545768e-01
    9.982672181362039e-01
    9.984470598688657e-01
    9.986170571794587e-01
    9.987772083919715e-01
    9.989275119275132e-01
    9.990679663043477e-01
    9.991985701379369e-01
    9.993193221410008e-01
    9.994302211236081e-01
    9.995312659933240e-01
    9.996224557554706e-01
    9.997037895136229e-01
    9.997752664706340e-01
    9.998368859309700e-01
    9.998886473067012e-01
    9.999305501355009e-01
    9.999625941483602e-01
    9.999847796329174e-01
    9.999971112980756e-01];



w=[7.413338416432071e-06
    1.725676977373923e-05
    2.711460656520586e-05
    3.697344200643550e-05
    4.683216706971275e-05
    5.669050651151729e-05
    6.654831593030788e-05
    7.640548208416073e-05
    8.626190132806909e-05
    9.611747354547055e-05
    1.059721000990171e-04
    1.158256830390418e-04
    1.256781247645614e-04
    1.355293278657728e-04
    1.453791950459568e-04
    1.552276290807009e-04
    1.650745327957399e-04
    1.749198090545718e-04
    1.847633607514339e-04
    1.946050908073345e-04
    2.044449021678837e-04
    2.142826978022177e-04
    2.241183807026033e-04
    2.339518538844774e-04
    2.437830203867700e-04
    2.536117832724149e-04
    2.634380456289895e-04
    2.732617105694400e-04
    2.830826812328697e-04
    2.929008607853699e-04
    3.027161524208810e-04
    3.125284593620774e-04
    3.223376848612676e-04
    3.321437322013066e-04
    3.419465046965173e-04
    3.517459056936195e-04
    3.615418385726617e-04
    3.713342067479596e-04
    3.811229136690341e-04
    3.909078628215541e-04
    4.006889577282802e-04
    4.104661019500086e-04
    4.202391990865167e-04
    4.300081527775093e-04
    4.397728667035656e-04
    4.495332445870837e-04
    4.592891901932303e-04
    4.690406073308848e-04
    4.787873998535875e-04
    4.885294716604855e-04
    4.982667266972789e-04
    5.079990689571673e-04
    5.177264024817950e-04
    5.274486313621963e-04
    5.371656597397405e-04
    5.468773918070770e-04
    5.565837318090788e-04
    5.662845840437867e-04
    5.759798528633512e-04
    5.856694426749778e-04
    5.953532579418668e-04
    6.050312031841557e-04
    6.147031829798617e-04
    6.243691019658201e-04
    6.340288648386262e-04
    6.436823763555741e-04
    6.533295413355955e-04
    6.629702646601990e-04
    6.726044512744057e-04
    6.822320061876894e-04
    6.918528344749099e-04
    7.014668412772516e-04
    7.110739318031568e-04
    7.206740113292620e-04
    7.302669852013290e-04
    7.398527588351816e-04
    7.494312377176357e-04
    7.590023274074321e-04
    7.685659335361670e-04
    7.781219618092233e-04
    7.876703180066994e-04
    7.972109079843394e-04
    8.067436376744594e-04
    8.162684130868772e-04
    8.257851403098364e-04
    8.352937255109361e-04
    8.447940749380504e-04
    8.542860949202592e-04
    8.637696918687673e-04
    8.732447722778282e-04
    8.827112427256674e-04
    8.921690098754021e-04
    9.016179804759605e-04
    9.110580613630041e-04
    9.204891594598446e-04
    9.299111817783607e-04
    9.393240354199172e-04
    9.487276275762778e-04
    9.581218655305246e-04
    9.675066566579679e-04
    9.768819084270612e-04
    9.862475284003150e-04
    9.956034242352054e-04
    1.004949503685087e-03
    1.014285674600100e-03
    1.023611844928081e-03
    1.032927922715472e-03
    1.042233816108219e-03
    1.051529433352690e-03
    1.060814682796570e-03
    1.070089472889766e-03
    1.079353712185315e-03
    1.088607309340280e-03
    1.097850173116653e-03
    1.107082212382255e-03
    1.116303336111633e-03
    1.125513453386959e-03
    1.134712473398925e-03
    1.143900305447640e-03
    1.153076858943521e-03
    1.162242043408192e-03
    1.171395768475371e-03
    1.180537943891762e-03
    1.189668479517946e-03
    1.198787285329270e-03
    1.207894271416733e-03
    1.216989347987873e-03
    1.226072425367654e-03
    1.235143413999347e-03
    1.244202224445417e-03
    1.253248767388401e-03
    1.262282953631792e-03
    1.271304694100914e-03
    1.280313899843806e-03
    1.289310482032094e-03
    1.298294351961871e-03
    1.307265421054567e-03
    1.316223600857826e-03
    1.325168803046379e-03
    1.334100939422908e-03
    1.343019921918926e-03
    1.351925662595635e-03
    1.360818073644801e-03
    1.369697067389615e-03
    1.378562556285561e-03
    1.387414452921273e-03
    1.396252670019407e-03
    1.405077120437490e-03
    1.413887717168789e-03
    1.422684373343161e-03
    1.431467002227917e-03
    1.440235517228671e-03
    1.448989831890196e-03
    1.457729859897278e-03
    1.466455515075564e-03
    1.475166711392415e-03
    1.483863362957751e-03
    1.492545384024901e-03
    1.501212688991444e-03
    1.509865192400058e-03
    1.518502808939362e-03
    1.527125453444752e-03
    1.535733040899247e-03
    1.544325486434322e-03
    1.552902705330750e-03
    1.561464613019433e-03
    1.570011125082237e-03
    1.578542157252828e-03
    1.587057625417495e-03
    1.595557445615986e-03
    1.604041534042337e-03
    1.612509807045688e-03
    1.620962181131121e-03
    1.629398572960475e-03
    1.637818899353169e-03
    1.646223077287026e-03
    1.654611023899084e-03
    1.662982656486419e-03
    1.671337892506963e-03
    1.679676649580308e-03
    1.687998845488527e-03
    1.696304398176983e-03
    1.704593225755132e-03
    1.712865246497340e-03
    1.721120378843683e-03
    1.729358541400747e-03
    1.737579652942445e-03
    1.745783632410800e-03
    1.753970398916756e-03
    1.762139871740975e-03
    1.770291970334624e-03
    1.778426614320179e-03
    1.786543723492215e-03
    1.794643217818192e-03
    1.802725017439252e-03
    1.810789042670997e-03
    1.818835214004284e-03
    1.826863452106004e-03
    1.834873677819862e-03
    1.842865812167164e-03
    1.850839776347589e-03
    1.858795491739974e-03
    1.866732879903076e-03
    1.874651862576361e-03
    1.882552361680765e-03
    1.890434299319467e-03
    1.898297597778659e-03
    1.906142179528309e-03
    1.913967967222926e-03
    1.921774883702324e-03
    1.929562851992384e-03
    1.937331795305809e-03
    1.945081637042882e-03
    1.952812300792226e-03
    1.960523710331551e-03
    1.968215789628411e-03
    1.975888462840948e-03
    1.983541654318645e-03
    1.991175288603067e-03
    1.998789290428615e-03
    2.006383584723250e-03
    2.013958096609253e-03
    2.021512751403949e-03
    2.029047474620451e-03
    2.036562191968393e-03
    2.044056829354659e-03
    2.051531312884117e-03
    2.058985568860349e-03
    2.066419523786370e-03
    2.073833104365362e-03
    2.081226237501393e-03
    2.088598850300136e-03
    2.095950870069590e-03
    2.103282224320793e-03
    2.110592840768541e-03
    2.117882647332103e-03
    2.125151572135922e-03
    2.132399543510334e-03
    2.139626489992268e-03
    2.146832340325951e-03
    2.154017023463618e-03
    2.161180468566204e-03
    2.168322605004043e-03
    2.175443362357568e-03
    2.182542670418008e-03
    2.189620459188072e-03
    2.196676658882644e-03
    2.203711199929472e-03
    2.210724012969852e-03
    2.217715028859311e-03
    2.224684178668293e-03
    2.231631393682832e-03
    2.238556605405238e-03
    2.245459745554762e-03
    2.252340746068279e-03
    2.259199539100956e-03
    2.266036057026913e-03
    2.272850232439903e-03
    2.279641998153969e-03
    2.286411287204108e-03
    2.293158032846927e-03
    2.299882168561310e-03
    2.306583628049067e-03
    2.313262345235589e-03
    2.319918254270500e-03
    2.326551289528307e-03
    2.333161385609047e-03
    2.339748477338928e-03
    2.346312499770981e-03
    2.352853388185686e-03
    2.359371078091624e-03
    2.365865505226107e-03
    2.372336605555810e-03
    2.378784315277403e-03
    2.385208570818183e-03
    2.391609308836700e-03
    2.397986466223378e-03
    2.404339980101141e-03
    2.410669787826031e-03
    2.416975826987830e-03
    2.423258035410665e-03
    2.429516351153631e-03
    2.435750712511404e-03
    2.441961058014833e-03
    2.448147326431567e-03
    2.454309456766642e-03
    2.460447388263092e-03
    2.466561060402543e-03
    2.472650412905817e-03
    2.478715385733515e-03
    2.484755919086618e-03
    2.490771953407074e-03
    2.496763429378383e-03
    2.502730287926187e-03
    2.508672470218845e-03
    2.514589917668023e-03
    2.520482571929258e-03
    2.526350374902550e-03
    2.532193268732920e-03
    2.538011195810989e-03
    2.543804098773543e-03
    2.549571920504098e-03
    2.555314604133467e-03
    2.561032093040317e-03
    2.566724330851727e-03
    2.572391261443745e-03
    2.578032828941942e-03
    2.583648977721962e-03
    2.589239652410075e-03
    2.594804797883713e-03
    2.600344359272021e-03
    2.605858281956397e-03
    2.611346511571029e-03
    2.616808994003433e-03
    2.622245675394985e-03
    2.627656502141452e-03
    2.633041420893520e-03
    2.638400378557324e-03
    2.643733322294964e-03
    2.649040199525035e-03
    2.654320957923137e-03
    2.659575545422398e-03
    2.664803910213981e-03
    2.670006000747599e-03
    2.675181765732026e-03
    2.680331154135593e-03
    2.685454115186698e-03
    2.690550598374310e-03
    2.695620553448458e-03
    2.700663930420736e-03
    2.705680679564785e-03
    2.710670751416793e-03
    2.715634096775979e-03
    2.720570666705080e-03
    2.725480412530825e-03
    2.730363285844427e-03
    2.735219238502055e-03
    2.740048222625306e-03
    2.744850190601680e-03
    2.749625095085050e-03
    2.754372888996129e-03
    2.759093525522930e-03
    2.763786958121232e-03
    2.768453140515037e-03
    2.773092026697030e-03
    2.777703570929024e-03
    2.782287727742421e-03
    2.786844451938652e-03
    2.791373698589630e-03
    2.795875423038185e-03
    2.800349580898514e-03
    2.804796128056608e-03
    2.809215020670695e-03
    2.813606215171669e-03
    2.817969668263520e-03
    2.822305336923760e-03
    2.826613178403850e-03
    2.830893150229617e-03
    2.835145210201678e-03
    2.839369316395851e-03
    2.843565427163576e-03
    2.847733501132313e-03
    2.851873497205961e-03
    2.855985374565260e-03
    2.860069092668195e-03
    2.864124611250388e-03
    2.868151890325506e-03
    2.872150890185647e-03
    2.876121571401736e-03
    2.880063894823915e-03
    2.883977821581920e-03
    2.887863313085473e-03
    2.891720331024665e-03
    2.895548837370322e-03
    2.899348794374390e-03
    2.903120164570302e-03
    2.906862910773354e-03
    2.910576996081061e-03
    2.914262383873530e-03
    2.917919037813817e-03
    2.921546921848289e-03
    2.925146000206969e-03
    2.928716237403907e-03
    2.932257598237510e-03
    2.935770047790903e-03
    2.939253551432268e-03
    2.942708074815184e-03
    2.946133583878970e-03
    2.949530044849017e-03
    2.952897424237124e-03
    2.956235688841825e-03
    2.959544805748718e-03
    2.962824742330790e-03
    2.966075466248741e-03
    2.969296945451293e-03
    2.972489148175522e-03
    2.975652042947158e-03
    2.978785598580897e-03
    2.981889784180712e-03
    2.984964569140165e-03
    2.988009923142688e-03
    2.991025816161902e-03
    2.994012218461904e-03
    2.996969100597564e-03
    2.999896433414806e-03
    3.002794188050909e-03
    3.005662335934783e-03
    3.008500848787250e-03
    3.011309698621333e-03
    3.014088857742515e-03
    3.016838298749028e-03
    3.019557994532115e-03
    3.022247918276299e-03
    3.024908043459645e-03
    3.027538343854029e-03
    3.030138793525386e-03
    3.032709366833974e-03
    3.035250038434625e-03
    3.037760783276995e-03
    3.040241576605799e-03
    3.042692393961083e-03
    3.045113211178437e-03
    3.047504004389244e-03
    3.049864750020924e-03
    3.052195424797148e-03
    3.054496005738085e-03
    3.056766470160619e-03
    3.059006795678572e-03
    3.061216960202934e-03
    3.063396941942065e-03
    3.065546719401931e-03
    3.067666271386294e-03
    3.069755576996934e-03
    3.071814615633858e-03
    3.073843366995487e-03
    3.075841811078876e-03
    3.077809928179896e-03
    3.079747698893435e-03
    3.081655104113591e-03
    3.083532125033856e-03
    3.085378743147301e-03
    3.087194940246763e-03
    3.088980698425022e-03
    3.090736000074978e-03
    3.092460827889820e-03
    3.094155164863208e-03
    3.095818994289429e-03
    3.097452299763566e-03
    3.099055065181663e-03
    3.100627274740879e-03
    3.102168912939644e-03
    3.103679964577818e-03
    3.105160414756831e-03
    3.106610248879840e-03
    3.108029452651865e-03
    3.109418012079937e-03
    3.110775913473224e-03
    3.112103143443181e-03
    3.113399688903676e-03
    3.114665537071113e-03
    3.115900675464568e-03
    3.117105091905902e-03
    3.118278774519892e-03
    3.119421711734338e-03
    3.120533892280183e-03
    3.121615305191622e-03
    3.122665939806213e-03
    3.123685785764977e-03
    3.124674833012502e-03
    3.125633071797050e-03
    3.126560492670637e-03
    3.127457086489142e-03
    3.128322844412387e-03
    3.129157757904231e-03
    3.129961818732646e-03
    3.130735018969810e-03
    3.131477350992170e-03
    3.132188807480534e-03
    3.132869381420126e-03
    3.133519066100670e-03
    3.134137855116449e-03
    3.134725742366366e-03
    3.135282722054010e-03
    3.135808788687708e-03
    3.136303937080585e-03
    3.136768162350607e-03
    3.137201459920638e-03
    3.137603825518479e-03
    3.137975255176912e-03
    3.138315745233739e-03
    3.138625292331819e-03
    3.138903893419101e-03
    3.139151545748650e-03
    3.139368246878682e-03
    3.139553994672580e-03
    3.139708787298920e-03
    3.139832623231489e-03
    3.139925501249297e-03
    3.139987420436592e-03
    3.140018380182867e-03
    3.140018380182867e-03
    3.139987420436592e-03
    3.139925501249297e-03
    3.139832623231489e-03
    3.139708787298920e-03
    3.139553994672580e-03
    3.139368246878682e-03
    3.139151545748650e-03
    3.138903893419101e-03
    3.138625292331819e-03
    3.138315745233739e-03
    3.137975255176912e-03
    3.137603825518479e-03
    3.137201459920638e-03
    3.136768162350607e-03
    3.136303937080585e-03
    3.135808788687708e-03
    3.135282722054010e-03
    3.134725742366366e-03
    3.134137855116449e-03
    3.133519066100670e-03
    3.132869381420126e-03
    3.132188807480534e-03
    3.131477350992170e-03
    3.130735018969810e-03
    3.129961818732646e-03
    3.129157757904231e-03
    3.128322844412387e-03
    3.127457086489142e-03
    3.126560492670637e-03
    3.125633071797050e-03
    3.124674833012502e-03
    3.123685785764977e-03
    3.122665939806213e-03
    3.121615305191622e-03
    3.120533892280183e-03
    3.119421711734338e-03
    3.118278774519892e-03
    3.117105091905902e-03
    3.115900675464568e-03
    3.114665537071113e-03
    3.113399688903676e-03
    3.112103143443181e-03
    3.110775913473224e-03
    3.109418012079937e-03
    3.108029452651865e-03
    3.106610248879840e-03
    3.105160414756831e-03
    3.103679964577818e-03
    3.102168912939644e-03
    3.100627274740879e-03
    3.099055065181663e-03
    3.097452299763566e-03
    3.095818994289429e-03
    3.094155164863208e-03
    3.092460827889820e-03
    3.090736000074978e-03
    3.088980698425022e-03
    3.087194940246763e-03
    3.085378743147301e-03
    3.083532125033856e-03
    3.081655104113591e-03
    3.079747698893435e-03
    3.077809928179896e-03
    3.075841811078876e-03
    3.073843366995487e-03
    3.071814615633858e-03
    3.069755576996934e-03
    3.067666271386294e-03
    3.065546719401931e-03
    3.063396941942065e-03
    3.061216960202934e-03
    3.059006795678572e-03
    3.056766470160619e-03
    3.054496005738085e-03
    3.052195424797148e-03
    3.049864750020924e-03
    3.047504004389244e-03
    3.045113211178437e-03
    3.042692393961083e-03
    3.040241576605799e-03
    3.037760783276995e-03
    3.035250038434625e-03
    3.032709366833974e-03
    3.030138793525386e-03
    3.027538343854029e-03
    3.024908043459645e-03
    3.022247918276299e-03
    3.019557994532115e-03
    3.016838298749028e-03
    3.014088857742515e-03
    3.011309698621333e-03
    3.008500848787250e-03
    3.005662335934783e-03
    3.002794188050909e-03
    2.999896433414806e-03
    2.996969100597564e-03
    2.994012218461904e-03
    2.991025816161902e-03
    2.988009923142688e-03
    2.984964569140165e-03
    2.981889784180712e-03
    2.978785598580897e-03
    2.975652042947158e-03
    2.972489148175522e-03
    2.969296945451293e-03
    2.966075466248741e-03
    2.962824742330790e-03
    2.959544805748718e-03
    2.956235688841825e-03
    2.952897424237124e-03
    2.949530044849017e-03
    2.946133583878970e-03
    2.942708074815184e-03
    2.939253551432268e-03
    2.935770047790903e-03
    2.932257598237510e-03
    2.928716237403907e-03
    2.925146000206969e-03
    2.921546921848289e-03
    2.917919037813817e-03
    2.914262383873530e-03
    2.910576996081061e-03
    2.906862910773354e-03
    2.903120164570302e-03
    2.899348794374390e-03
    2.895548837370322e-03
    2.891720331024665e-03
    2.887863313085473e-03
    2.883977821581920e-03
    2.880063894823915e-03
    2.876121571401736e-03
    2.872150890185647e-03
    2.868151890325506e-03
    2.864124611250388e-03
    2.860069092668195e-03
    2.855985374565260e-03
    2.851873497205961e-03
    2.847733501132313e-03
    2.843565427163576e-03
    2.839369316395851e-03
    2.835145210201678e-03
    2.830893150229617e-03
    2.826613178403850e-03
    2.822305336923760e-03
    2.817969668263520e-03
    2.813606215171669e-03
    2.809215020670695e-03
    2.804796128056608e-03
    2.800349580898514e-03
    2.795875423038185e-03
    2.791373698589630e-03
    2.786844451938652e-03
    2.782287727742421e-03
    2.777703570929024e-03
    2.773092026697030e-03
    2.768453140515037e-03
    2.763786958121232e-03
    2.759093525522930e-03
    2.754372888996129e-03
    2.749625095085050e-03
    2.744850190601680e-03
    2.740048222625306e-03
    2.735219238502055e-03
    2.730363285844427e-03
    2.725480412530825e-03
    2.720570666705080e-03
    2.715634096775979e-03
    2.710670751416793e-03
    2.705680679564785e-03
    2.700663930420736e-03
    2.695620553448458e-03
    2.690550598374310e-03
    2.685454115186698e-03
    2.680331154135593e-03
    2.675181765732026e-03
    2.670006000747599e-03
    2.664803910213981e-03
    2.659575545422398e-03
    2.654320957923137e-03
    2.649040199525035e-03
    2.643733322294964e-03
    2.638400378557324e-03
    2.633041420893520e-03
    2.627656502141452e-03
    2.622245675394985e-03
    2.616808994003433e-03
    2.611346511571029e-03
    2.605858281956397e-03
    2.600344359272021e-03
    2.594804797883713e-03
    2.589239652410075e-03
    2.583648977721962e-03
    2.578032828941942e-03
    2.572391261443745e-03
    2.566724330851727e-03
    2.561032093040317e-03
    2.555314604133467e-03
    2.549571920504098e-03
    2.543804098773543e-03
    2.538011195810989e-03
    2.532193268732920e-03
    2.526350374902550e-03
    2.520482571929258e-03
    2.514589917668023e-03
    2.508672470218845e-03
    2.502730287926187e-03
    2.496763429378383e-03
    2.490771953407074e-03
    2.484755919086618e-03
    2.478715385733515e-03
    2.472650412905817e-03
    2.466561060402543e-03
    2.460447388263092e-03
    2.454309456766642e-03
    2.448147326431567e-03
    2.441961058014833e-03
    2.435750712511404e-03
    2.429516351153631e-03
    2.423258035410665e-03
    2.416975826987830e-03
    2.410669787826031e-03
    2.404339980101141e-03
    2.397986466223378e-03
    2.391609308836700e-03
    2.385208570818183e-03
    2.378784315277403e-03
    2.372336605555810e-03
    2.365865505226107e-03
    2.359371078091624e-03
    2.352853388185686e-03
    2.346312499770981e-03
    2.339748477338928e-03
    2.333161385609047e-03
    2.326551289528307e-03
    2.319918254270500e-03
    2.313262345235589e-03
    2.306583628049067e-03
    2.299882168561310e-03
    2.293158032846927e-03
    2.286411287204108e-03
    2.279641998153969e-03
    2.272850232439903e-03
    2.266036057026913e-03
    2.259199539100956e-03
    2.252340746068279e-03
    2.245459745554762e-03
    2.238556605405238e-03
    2.231631393682832e-03
    2.224684178668293e-03
    2.217715028859311e-03
    2.210724012969852e-03
    2.203711199929472e-03
    2.196676658882644e-03
    2.189620459188072e-03
    2.182542670418008e-03
    2.175443362357568e-03
    2.168322605004043e-03
    2.161180468566204e-03
    2.154017023463618e-03
    2.146832340325951e-03
    2.139626489992268e-03
    2.132399543510334e-03
    2.125151572135922e-03
    2.117882647332103e-03
    2.110592840768541e-03
    2.103282224320793e-03
    2.095950870069590e-03
    2.088598850300136e-03
    2.081226237501393e-03
    2.073833104365362e-03
    2.066419523786370e-03
    2.058985568860349e-03
    2.051531312884117e-03
    2.044056829354659e-03
    2.036562191968393e-03
    2.029047474620451e-03
    2.021512751403949e-03
    2.013958096609253e-03
    2.006383584723250e-03
    1.998789290428615e-03
    1.991175288603067e-03
    1.983541654318645e-03
    1.975888462840948e-03
    1.968215789628411e-03
    1.960523710331551e-03
    1.952812300792226e-03
    1.945081637042882e-03
    1.937331795305809e-03
    1.929562851992384e-03
    1.921774883702324e-03
    1.913967967222926e-03
    1.906142179528309e-03
    1.898297597778659e-03
    1.890434299319467e-03
    1.882552361680765e-03
    1.874651862576361e-03
    1.866732879903076e-03
    1.858795491739974e-03
    1.850839776347589e-03
    1.842865812167164e-03
    1.834873677819862e-03
    1.826863452106004e-03
    1.818835214004284e-03
    1.810789042670997e-03
    1.802725017439252e-03
    1.794643217818192e-03
    1.786543723492215e-03
    1.778426614320179e-03
    1.770291970334624e-03
    1.762139871740975e-03
    1.753970398916756e-03
    1.745783632410800e-03
    1.737579652942445e-03
    1.729358541400747e-03
    1.721120378843683e-03
    1.712865246497340e-03
    1.704593225755132e-03
    1.696304398176983e-03
    1.687998845488527e-03
    1.679676649580308e-03
    1.671337892506963e-03
    1.662982656486419e-03
    1.654611023899084e-03
    1.646223077287026e-03
    1.637818899353169e-03
    1.629398572960475e-03
    1.620962181131121e-03
    1.612509807045688e-03
    1.604041534042337e-03
    1.595557445615986e-03
    1.587057625417495e-03
    1.578542157252828e-03
    1.570011125082237e-03
    1.561464613019433e-03
    1.552902705330750e-03
    1.544325486434322e-03
    1.535733040899247e-03
    1.527125453444752e-03
    1.518502808939362e-03
    1.509865192400058e-03
    1.501212688991444e-03
    1.492545384024901e-03
    1.483863362957751e-03
    1.475166711392415e-03
    1.466455515075564e-03
    1.457729859897278e-03
    1.448989831890196e-03
    1.440235517228671e-03
    1.431467002227917e-03
    1.422684373343161e-03
    1.413887717168789e-03
    1.405077120437490e-03
    1.396252670019407e-03
    1.387414452921273e-03
    1.378562556285561e-03
    1.369697067389615e-03
    1.360818073644801e-03
    1.351925662595635e-03
    1.343019921918926e-03
    1.334100939422908e-03
    1.325168803046379e-03
    1.316223600857826e-03
    1.307265421054567e-03
    1.298294351961871e-03
    1.289310482032094e-03
    1.280313899843806e-03
    1.271304694100914e-03
    1.262282953631792e-03
    1.253248767388401e-03
    1.244202224445417e-03
    1.235143413999347e-03
    1.226072425367654e-03
    1.216989347987873e-03
    1.207894271416733e-03
    1.198787285329270e-03
    1.189668479517946e-03
    1.180537943891762e-03
    1.171395768475371e-03
    1.162242043408192e-03
    1.153076858943521e-03
    1.143900305447640e-03
    1.134712473398925e-03
    1.125513453386959e-03
    1.116303336111633e-03
    1.107082212382255e-03
    1.097850173116653e-03
    1.088607309340280e-03
    1.079353712185315e-03
    1.070089472889766e-03
    1.060814682796570e-03
    1.051529433352690e-03
    1.042233816108219e-03
    1.032927922715472e-03
    1.023611844928081e-03
    1.014285674600100e-03
    1.004949503685087e-03
    9.956034242352054e-04
    9.862475284003150e-04
    9.768819084270612e-04
    9.675066566579679e-04
    9.581218655305246e-04
    9.487276275762778e-04
    9.393240354199172e-04
    9.299111817783607e-04
    9.204891594598446e-04
    9.110580613630041e-04
    9.016179804759605e-04
    8.921690098754021e-04
    8.827112427256674e-04
    8.732447722778282e-04
    8.637696918687673e-04
    8.542860949202592e-04
    8.447940749380504e-04
    8.352937255109361e-04
    8.257851403098364e-04
    8.162684130868772e-04
    8.067436376744594e-04
    7.972109079843394e-04
    7.876703180066994e-04
    7.781219618092233e-04
    7.685659335361670e-04
    7.590023274074321e-04
    7.494312377176357e-04
    7.398527588351816e-04
    7.302669852013290e-04
    7.206740113292620e-04
    7.110739318031568e-04
    7.014668412772516e-04
    6.918528344749099e-04
    6.822320061876894e-04
    6.726044512744057e-04
    6.629702646601990e-04
    6.533295413355955e-04
    6.436823763555741e-04
    6.340288648386262e-04
    6.243691019658201e-04
    6.147031829798617e-04
    6.050312031841557e-04
    5.953532579418668e-04
    5.856694426749778e-04
    5.759798528633512e-04
    5.662845840437867e-04
    5.565837318090788e-04
    5.468773918070770e-04
    5.371656597397405e-04
    5.274486313621963e-04
    5.177264024817950e-04
    5.079990689571673e-04
    4.982667266972789e-04
    4.885294716604855e-04
    4.787873998535875e-04
    4.690406073308848e-04
    4.592891901932303e-04
    4.495332445870837e-04
    4.397728667035656e-04
    4.300081527775093e-04
    4.202391990865167e-04
    4.104661019500086e-04
    4.006889577282802e-04
    3.909078628215541e-04
    3.811229136690341e-04
    3.713342067479596e-04
    3.615418385726617e-04
    3.517459056936195e-04
    3.419465046965173e-04
    3.321437322013066e-04
    3.223376848612676e-04
    3.125284593620774e-04
    3.027161524208810e-04
    2.929008607853699e-04
    2.830826812328697e-04
    2.732617105694400e-04
    2.634380456289895e-04
    2.536117832724149e-04
    2.437830203867700e-04
    2.339518538844774e-04
    2.241183807026033e-04
    2.142826978022177e-04
    2.044449021678837e-04
    1.946050908073345e-04
    1.847633607514339e-04
    1.749198090545718e-04
    1.650745327957399e-04
    1.552276290807009e-04
    1.453791950459568e-04
    1.355293278657728e-04
    1.256781247645614e-04
    1.158256830390418e-04
    1.059721000990171e-04
    9.611747354547055e-05
    8.626190132806909e-05
    7.640548208416073e-05
    6.654831593030788e-05
    5.669050651151729e-05
    4.683216706971275e-05
    3.697344200643550e-05
    2.711460656520586e-05
    1.725676977373923e-05
    7.413338416432071e-06];

xw=[x w];



















