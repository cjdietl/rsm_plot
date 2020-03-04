% This Function takes a spc file and creates
% an HDF output file with the corresponding reciprocal space coordinates.
% spc_file: spc_file to be read out
% spc_path: Path to the folder of the .myt files
% d_mythen: sample detector distance in mm
% myt_max: The Program will not read .myt files with a Mythen_max less than
% myt_max Only channels with intensities greater than myt_max will be added

function [out_file,q_hkl]=rsm_map_fourc(spc_file,spc_path,scan,d_mythen,myt_max)
tic

fprintf('Start reading SPC file ...');
[ub,scandata,ch_ar,~,~,lambda,mean_ch] = specread(spc_file,spc_path,scan,myt_max);
et=toc;
fprintf('Finished reading, took %f s \n',et);

inv_ub=inv(ub);

inv_eta=eye(3,3);
inv_phi=eye(3,3);
inv_chi=eye(3,3);

init=false;

[n_scans,~] = size(scandata);

%ch_mean = 640;
dg2rad=pi/180.0;

tic
fprintf('(Takagi) Converting angles to reciprocal space ... ');
clr_l=0;
q_entry=1;
old_progress=0;
% The delta, eta and phi circles all have a left-handed rotation sense. The
% delta circle is at zero when the detector arm is along the positive y
% direction. Eta is zero when the chi circle is in the x-z plane, with the
% chi axis along the y axis. The chi rotation is right-handed. The zero of
% chi puts the phi rotation axis along the z direction with the "surface
% normal" of the phi table pointed in the positive z direction. The zero of
% the phi circle is arbitrary.
% http://www.certif.com/spec_help/psic.html
 
for i=1:n_scans
    if(isempty(ch_ar{i,1}))
        continue;
    end
    
    %First available channel array sets the number of channels
    if(~init)
        n_ch = numel(ch_ar{i,1}{1,1});
        q_hkl=zeros(n_scans*n_ch,4,'double');
        init=true;
    end
    
    eta = scandata(i,2)*dg2rad;
    c=cos(-eta);
    s=sin(-eta);
    inv_eta(1,1)=c;
    inv_eta(1,2)=s;
    inv_eta(2,1)=-s;
    inv_eta(2,2)=c;
    
    phi = scandata(i,4)*dg2rad;
    
    if(phi~=0)
        c=cos(-phi);
        s=sin(-phi);
        inv_phi(1,1)=c;
        inv_phi(1,2)=s;
        inv_phi(2,1)=-s;
        inv_phi(2,2)=c;
    end
    
    chi = scandata(i,3)*dg2rad;
    c=cos(-chi);
    s=sin(-chi);
    inv_chi(1,1)=c;
    inv_chi(1,3)=s;
    inv_chi(3,1)=-s;
    inv_chi(3,3)=c;
    
    sec=scandata(i,5);
    
    % Treat the detector as a circular shaped ('channel per degree
    % approx').
    %Normalize to time
    
    for ii=1:n_ch
        %Look only at pixels > myt_max
        if(double(ch_ar{i,1}{1,2}(ii))>myt_max)
            delta = scandata(i,1)*dg2rad-(double((ch_ar{i,1}{1,1}(ii)-mean_ch))*0.05/d_mythen);
            [h, k, l]=getQ(inv_ub,inv_phi,inv_chi,inv_eta,delta,lambda);
            q_hkl(q_entry,:)=[h k l double(ch_ar{i,1}{1,2}(ii))/sec];
            q_entry = q_entry+1;
        else
            continue;
        end
        
        progress=round((n_ch*(i-1)+ii)/(n_scans*n_ch)*100);
        if(old_progress~=progress)
            for a=1:clr_l
                fprintf('\b');
            end
            str=strcat(num2str(progress),' percent');
            clr_l=length(str);
            fprintf(str);
        end
        old_progress=progress;
    end
end

%Overwrite Line
for a=1:clr_l
    fprintf('\b');
end
et=toc;
fprintf('done, took %f s\n',et);

rs_file=sprintf('%s_%d_q_coords.h5',spc_file,scan);
%rs_file='q_coords.h5';

%Delete file, if existing
if(exist(rs_file,'file'))
    delete(rs_file);
end

h5create(rs_file,'/QCOORDS',size(q_hkl(1:q_entry-1,:)),'ChunkSize',size(q_hkl(1:q_entry-1,:)),'Deflate',9);
h5write(rs_file,'/QCOORDS',q_hkl(1:q_entry-1,:));

fprintf('Q coordinates saved to %s\n',rs_file);
out_file=rs_file;
end

% Convert Angles to reciprocal space. Formalisms follows the geometry from
% the paper (same as SPEC)
% You, H. (1999). Journal of Applied Crystallography, 32(4), 614–623. doi:10.1107/S0021889899001223
function [h,k,l]=getQ(inv_ub,inv_phi,inv_chi,inv_eta,delta,lambda)
%Incoming Beam
kappa=zeros(3);
kappa(2)=2*pi/lambda;

delta_m=eye(3,3);
c=cos(delta);
s=sin(delta);
delta_m(1,1)=c;
delta_m(1,2)=s;
delta_m(2,1)=-s;
delta_m(2,2)=c;
q=inv_ub*inv_phi*inv_chi*inv_eta*(delta_m*kappa-kappa);
h=q(1,1);
k=q(2,1);
l=q(3,1);
end

function [ub,data,ch_ar,alpha,lat,lambda,mean_ch]=specread(spc_file,myt_path,scanno,myt_max)
filename=[myt_path spc_file];
disp(filename);
fid=fopen(filename);
F=fread(fid);FS=char(F');
fposS=strfind(FS,['#S ' num2str(scanno)]);
rest = FS(fposS(1):end);
fposS_rest=strfind(rest,'#S');

% Last Scan
if(numel(fposS_rest)==1)
    scandata=rest(1:end-1);
else
    scandata=rest(1:fposS_rest(2));
end

fposL=strfind(scandata,'#L');

hdr = strsplit(scandata(1:fposL(1)),'#');
hdr = hdr(2:end);

alpha=zeros(3);
lat=zeros(3);

lat_str=hdr{5};
a=sscanf(lat_str,'G1 %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
a=a';
lat(1:3)=a(1:3);
alpha(1:3)=a(3:5);
lambda=a(end);

om_str=hdr{6};
a=sscanf(om_str,'G3 %f %f %f %f %f %f %f %f %f');
ub=zeros(3,3);

a=a';
ub(1,1:3)=a(1:3);
ub(2,1:3)=a(4:6);
ub(3,1:3)=a(7:9);

data=scandata(fposL(1):end);
hdrlen=length(hdr);
for i=1:hdrlen
    if(strcmp(hdr{i}(1),'N') == 1)
        cols_str=hdr{i};
        break;
    end
end

s = sscanf(cols_str,'N %i');

format='%s';
for i=1:s-1
    format = strcat(format,' %s');
end

data=strrep(data,'Two Theta','two_theta');
data=strrep(data,'Mythen Max','mythen_max');
data=strrep(data,'Mythen ROI','mythen_roi');
data=strrep(data,'Mythen total','mythen_total');
data=strrep(data,'Mythen Total','mythen_total');

%Skip '#L ' at the beginning
dataM=textscan(data(4:end),format);

for i=1:s
    if(strcmp(dataM{1,i}(1),'two_theta'))
        delta_idx = i;
    elseif(strcmp(dataM{1,i}(1),'Del'))
        delta_idx = i;
    elseif(strcmp(dataM{1,i}(1),'Theta'))
        eta_idx=i;
    elseif(strcmp(dataM{1,i}(1),'Th'))
        eta_idx=i;
    elseif(strcmp(dataM{1,i}(1),'Chi'))
        chi_idx = i;
    elseif(strcmp(dataM{1,i}(1),'Phi'))
        phi_idx = i;
    elseif(strcmp(dataM{1,i}(1),'Seconds'))
        time_idx = i;
    end
end

h_ar = str2double(dataM{1,1}(2:end));
n = numel(h_ar);

for a=n:-1:1
    if(~isnan(h_ar(a)))
        break;
    end
end

delta_ar = str2double(dataM{1,delta_idx}(2:a+1));
eta_ar = str2double(dataM{1,eta_idx}(2:a+1));
chi_ar = str2double(dataM{1,chi_idx}(2:a+1));
phi_ar = str2double(dataM{1,phi_idx}(2:a+1));
time_ar = str2double(dataM{1,time_idx}(2:a+1));

n=numel(delta_ar);

ch_ar{n}=[];
binned_ch_ar{n}=[];

fprintf('Reading Mythen Files ... ');
clr_l=0;

%For aborted scans, maximum in filename doesnt correspond to actual maximum
test=dir(strcat(myt_path,sprintf('%s_scan%d_%dof*',spc_file,scanno,1),'.myt'));
max_num=sscanf(test.name,strcat(sprintf('%s_scan%d_%dof',spc_file,scanno,1),'%d.myt'));

for i=1:n
    %Example: d1_34.spc_scan132_1843of2131.myt
    myt_file=sprintf('%s_scan%d_%dof%d',spc_file,scanno,i,max_num);
    myt_file=strcat(myt_path,myt_file,'.myt');
    myt_fid=fopen(myt_file);
  %  fprintf(i);
    ch_ar{i}=textscan(myt_fid,'%d %d');
    
    if(max(ch_ar{1,i}{1,1}) < myt_max)
        continue;
    end
    
    fclose(myt_fid);
    
    for a=1:clr_l
        fprintf('\b');
    end
    
   str=strcat(num2str(floor(i/n*100)),' percent');  
    %str=strcat(num2str(i),' percent');
    fprintf(str);
    clr_l=length(str);
end

mean_ch=641;

for a=1:clr_l
    fprintf('\b');
end
fprintf('done\n');
fclose(fid);

data=[delta_ar eta_ar chi_ar phi_ar time_ar];
ch_ar=ch_ar';

end

