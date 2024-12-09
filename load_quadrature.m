%%
clear all
clc

%%
P = 5;
filename = "channel_11_11_q5";

%%
V = fopen(filename+'_v.txt', 'r');
formatSpec = '%f';
tline = fgetl(V);
N_V = sscanf(tline, '%i');
vertices = reshape(fscanf(V,formatSpec), 2, [])';

%%
C = fopen(filename+'_c.txt', 'r');
tline = fgetl(C); 
N_C = sscanf(tline, '%i');
cells = {};

for i = 1:N_C
    tline = fgetl(C);
    cells{end+1} = sscanf(tline, '%i');
end

%%
F = fopen(filename+'_f.txt', 'r');
tline = fgetl(F);
N_F = sscanf(tline, '%i');
faces = {};

for i = 1:N_F
    tline = fgetl(F);
    faces{end+1} = sscanf(tline, '%i');
end

%%

for i = 1:N_C
    polygon_faces = {}; % the indices of vertices of each face of the polygon
    f_idx = cells{i}(2:end); % face index
    
    for idx = 1:length(f_idx)
        face = faces{f_idx(idx)+1}(2:end)+1; % the indeces of the vertices
        if ~isempty(polygon_faces)
            if ismember(polygon_faces{end}(1), face)
                polygon_faces{end} = flip(polygon_faces{end});
            end
        end
        polygon_faces{end+1} = face;
    end
    if ismember(polygon_faces{end}(1), polygon_faces{1})
        polygon_faces{end} = flip(polygon_faces{end});
    end
    
    polygon_vertices = []; % the vertices of the polygon, including control points
    spline_order_vett = [];
    
    v_idx = 1;
    for idx = 1:length(f_idx)
        face = polygon_faces{idx};
        v_idx = v_idx + length(face) - 1;
        for j = 1:length(face)-1
            polygon_vertices = [polygon_vertices; vertices(face(j), :)];
        end
        spline_order_vett = [spline_order_vett; length(face), v_idx];
    end
    polygon_vertices = [polygon_vertices; polygon_vertices(1,:)];
    
    X=polygon_vertices(:,1); Y=polygon_vertices(:,2);
    xmin = min(X);
    xmax = max(X);
    ymin = min(Y);
    ymax = max(Y);
    X = (X-xmin) / (xmax-xmin);
    Y = (Y-ymin) / (ymax-ymin);
    [Sx,Sy]=Sx_Sy_preparation(X,Y,1,spline_order_vett);
    [xyw,momsres,Z,Zin] = splcub(P*2,Sx,Sy,1); % need to evaluate polynomial up to degree 2P
    xyw(:, 1) = xyw(:, 1) * (xmax-xmin) + xmin;
    xyw(:, 2) = xyw(:, 2) * (ymax-ymin) + ymin;
    xyw(:, 3) = xyw(:, 3) * (xmax-xmin) * (ymax-ymin);
    dlmwrite(filename+'_q.txt', xyw','Delimiter',' ', '-append', 'precision', 16)
end
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
