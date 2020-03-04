%Creates a regularly spaced grid and interpolates the data points close to
%it, so it can be display as a rasterized image
%q_file: File created by rsm_map.m
%xstep: grid step size in x
%ystep: grid step size in y
%out(optional): output file (otherwise q_file_gridded.h5)
%xbottom,ybottom,xtop,ytop: definition of a rectangle around region of
%interest. Allows to reduce the number of data points, if only a specific
%region is interesting

function create_gridded_rsm_view(q_file,meshsize,xaxis,yaxis,xstep,ystep,xbottom,ybottom,xtop,ytop)
q_hkl=h5read(q_file,'/QCOORDS');
hkl=q_hkl(:,1:3);

[r,~] = size(hkl);

X=zeros(r,1);
Y=zeros(r,1);

xaxis=xaxis/norm(xaxis);
yaxis=yaxis/norm(yaxis);

for i=1:r
    X(i) = dot(hkl(i,1:3),xaxis);
    Y(i) = dot(hkl(i,1:3),yaxis);
end

Z=q_hkl(:,4);

%Good for fource xstep: 3e-4, ystep: 3e-2
%Guess a good xstep, ystep
if(nargin < 3)
    Xu=unique(X);
    
    Xdiff=zeros(numel(Xu),1);
    for i=2:length(Xu)
        Xdiff(i)=(Xu(i)-Xu(i-1));
    end
    
    %Resolution too high, if mean difference is taken
    xstep=mean(Xdiff)*1000;
    
    Yu=unique(Y);
    
    Ydiff=zeros(numel(Xu),1);
    for i=2:length(Yu)
        Ydiff(i)=(Yu(i)-Yu(i-1));
    end
    
    ystep=mean(Ydiff)*1000;
end

%If rectangle is defined, skip points outside of it
if(nargin > 6)
    if(xbottom~=-1)
        XWithinRectIdx=find((xbottom<=X) & (X<=xtop));
        YWithinRectIdx=find((ybottom<=Y) & (Y<=ytop));
        
        XYWithinRectIdx=intersect(XWithinRectIdx,YWithinRectIdx);
        
        X=X(XYWithinRectIdx);
        Y=Y(XYWithinRectIdx);
        Z=Z(XYWithinRectIdx);
    end
end

%First treat points as multiplicates of xstep / ystep to make calculation
%faster and easier
% Look for closest grid points around one data point
leftX=floor(X/xstep);
bottomY=floor(Y/ystep);

if(meshsize==0)
    ar_size=1;
elseif(meshsize==1)
    ar_size=9;
elseif(meshsize==2)
    ar_size=16;
end
GridPoints=zeros(ar_size,numel(leftX),2);

i=1;

if(meshsize==0)
    GriddedRSMPoints=[leftX bottomY];
else
    for x=-meshsize:meshsize
        for y=-meshsize:meshsize
            add=[leftX+x bottomY+y];
            GridPoints(i,:,:)=add;
            i=i+1;
        end
    end
    
    GriddedRSMPoints = [];
    for a=1:ar_size
        add=squeeze(GridPoints(a,:,:));
        GriddedRSMPoints = [GriddedRSMPoints;add];
    end
end

%Skip duplicates
GriddedRSMPoints=unique(GriddedRSMPoints,'rows');

xgridmin=min(GriddedRSMPoints(:,1));
ygridmin=min(GriddedRSMPoints(:,2));
xgridmax=max(GriddedRSMPoints(:,1));
ygridmax=max(GriddedRSMPoints(:,2));

fprintf('Gridding %d points... \n',numel(GriddedRSMPoints(:,1)));

if(numel(GriddedRSMPoints(:,1)) > 800000)
    fprintf('Too many points, aborting ...');
    return;
end
tic
F = scatteredInterpolant(X,Y,Z);
fprintf('Interpolation function generated \n');
GridValues = F(GriddedRSMPoints(:,1)*xstep,GriddedRSMPoints(:,2)*ystep);
et=toc;
rsm_grid=zeros((xgridmax-xgridmin+1),(ygridmax-ygridmin+1));

grid_idx=[(GriddedRSMPoints(:,1)-xgridmin) (GriddedRSMPoints(:,2)-ygridmin)];
grid_idx(:,1)=grid_idx(:,1)+1;
grid_idx(:,2)=grid_idx(:,2)+1;

[n_grid_points,~]=size(grid_idx);

%Interpolated points at sites where there was originally no data point
%closeby produce weird results and often negative (unphysical) results. Clean up by removing negative values at the
%end
for z=1:n_grid_points
    val = GridValues(z);
    if(val < 0)
        rsm_grid(grid_idx(z,1),grid_idx(z,2))=0;
    else
        rsm_grid(grid_idx(z,1),grid_idx(z,2))=val;
    end
end

fprintf('Gridding finished, took %f s\n',et);

XAxis=(xgridmin:xgridmax)*xstep;
YAxis=(ygridmin:ygridmax)*ystep;


file=strcat(q_file,'_gridded.h5');

%Delete file, if existing
if(exist(file,'file'))
    delete(file);
end


d=log10(rsm_grid);
hf=figure();
imagesc(XAxis,YAxis,d');
set(gca,'YDir','normal')
grid on

% Use maximum compression
[x,y] = size(rsm_grid);
h5create(file,'/RSM',[x y],'ChunkSize',[x y],'Deflate',9);
h5write(file, '/RSM', rsm_grid);

h5create(file,'/XAxis',numel(XAxis));
h5write(file, '/XAxis', XAxis);
h5create(file,'/YAxis',numel(YAxis));
h5write(file, '/YAxis', YAxis);

disp('RSM plotted\n');
end