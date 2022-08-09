%% Select the rocking curve

code_location = '/home/dzhigd/Software/nanomax_preprocessing';

% add class definitions and functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(fullfile(code_location,'functions'));
addpath(fullfile(code_location,'functions/openFunctions'));
addpath(fullfile(code_location,'classDef'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clc
clear;
close all;
scan_number = 442;

scan_path = sprintf('/home/dzhigd/work/projects/CsPbBr3_NC_BCDI_NanoMAX/data/sample0609_%d/scan_%06d_merlin.mat',scan_number,scan_number);
load(scan_path);

% Vertical coordinate on the detector is 1st
% Horizontal coordinate on the detector is 2nd

% Experiment parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correction angles
detector_delta_correction = 0.48;
detector_gamma_correction = 0;

% NanoMax convention:
% gamma - horizontal detector
% delta - vertical detector
% gonphi - rotation about vertical axis
% gontheta - rotation about horizontal axis
nanomax.photon_energy   = scan.motor_positions.energy;

if strcmp(scan.rocking_motor,'gonphi')
    nanomax.gonphi          = -scan.rocking_angles; % [deg] - can be a range of angles
    nanomax.gontheta        = -scan.motor_positions.gontheta; % [deg] - can be a range of angles
elseif strcmp(scan.rocking_motor,'gontheta')
    nanomax.gonphi          = scan.motor_positions.gonphi; % [deg] - can be a range of angles
    nanomax.gontheta        = scan.rocking_angles; % [deg] - can be a range of angles
end

nanomax.radius          = 0.9984; % [m]
nanomax.delta           = scan.motor_positions.delta+detector_delta_correction; % [deg] these angles are corrected with the sign respecting the rotation rules
nanomax.gamma           = scan.motor_positions.gamma+detector_gamma_correction; % [deg] 
nanomax.detector_pitch  = 55e-6; % [m]

nanomax.direct_beam     = round([251.7859,  250.4288]);       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants. They are needed for correct labeling of axes
h                       = 4.1357e-15;                                  % Plank's constant
c                       = 2.99792458e8;                                % Speed of light in vacuum

wavelength = h*c/nanomax.photon_energy;

k = 2*pi/wavelength; % wave vector

dq = k*2*atan(nanomax.detector_pitch/(2*nanomax.radius)); % q-space pitch at the detector plane

[hd,vd] = meshgrid(-size(scan.data,2)/2:(size(scan.data,2)/2-1),size(scan.data,1)/2:-1:-(size(scan.data,1)/2-1));

hd = (hd+size(scan.data,2)/2-nanomax.direct_beam(2)).*nanomax.detector_pitch;
vd = (vd+size(scan.data,1)/2-nanomax.direct_beam(1)).*nanomax.detector_pitch;
zd = ones(size(vd)).*nanomax.radius;

% figure;
% imagesc(hd(1,:),vd(:,1),log10(sum(scan.data,3)));

%% Data reduction

try 
    % TEST
%     scan.roi = [1,1,514,514];
    fprintf('The roi is used: %d %d %d %d\n',scan.roi);
catch
    f = figure;
    imagesc(log10(sum(scan.data,3)));axis image
    % set(gca,'FontSize',12);
    h = imrect;
    scan.roi = round(getPosition(h)); %[xmin ymin width height]
    save(scan_path,'scan');
    close(f);
end
% scan_path = sprintf('/home/dzhigd/work/projects/CsPbBr3_NC_BCDI_NanoMAX/data/sample0609_%d/sample0609_%d.mat',scan_number,scan_number);

data = scan.data(scan.roi(2):scan.roi(2)+scan.roi(4),scan.roi(1):scan.roi(1)+scan.roi(3),:,:,:);
hd = hd(scan.roi(2):scan.roi(2)+scan.roi(4),scan.roi(1):scan.roi(1)+scan.roi(3),:);
vd = vd(scan.roi(2):scan.roi(2)+scan.roi(4),scan.roi(1):scan.roi(1)+scan.roi(3),:);
zd = zd(scan.roi(2):scan.roi(2)+scan.roi(4),scan.roi(1):scan.roi(1)+scan.roi(3),:);

close;

d = [hd(:),vd(:),zd(:)]';

r = squeeze(sqrt(sum(d.^2,1)));

hq = k*(d(1,:)./r);
vq = k*(d(2,:)./r);
zq = k*(1-d(3,:)./r);

q = [hq;vq;zq];

% Check the coordinate mapping - real space m, direct scattering
% test = data(:,:,48,1,1);
% figure;
% scatter3(q(1,:),q(2,:),q(3,:),10,log10(test(:)),'filled','s');
% xlabel('qh');
% ylabel('qv');
% zlabel('qz');
% view(7.7549,53.9161);
% PROVED!

%% Estimation of the real space support
h_scale = 0.6e-6;
v_scale = 0.6e-6;
z_scale = 0.5e-6;

angular_step = abs(mean(diff(nanomax.gonphi)));
H = 2*3.14/5.8795e-10;

% Real space sampling
h_sampling = wavelength*nanomax.radius/(size(data,2)*nanomax.detector_pitch);
v_sampling = wavelength*nanomax.radius/(size(data,1)*nanomax.detector_pitch);
z_sampling = 360/(angular_step*size(data,3)*H);

h_N = round(h_scale/h_sampling);
v_N = round(h_scale/v_sampling);
z_N = round(h_scale/z_sampling);

support   = ones(h_N,v_N,z_N);
support   = padarray(support,round([(size(data,1)-h_N)/2,(size(data,2)-v_N)/2,(size(data,3)-z_N)/2]),'pre');
support   = padarray(support,round([(size(data,1)-h_N)/2-1,(size(data,2)-v_N)/2-1,(size(data,3)-z_N)/2-1]),'post');

%% Sample orientation matrix. Bounds the sample crystal with the laboratory frame
% Angles alpha beta gamma were manually adjusted so that known peaks 
% are exactly in their places

% X is horizontal, perp to the beam, Y is vertical

Rh = [1         0                    0; % detector rotation around horizintal axis 
      0         cosd(nanomax.delta) -sind(nanomax.delta);
      0         sind(nanomax.delta)  cosd(nanomax.delta)]; 

Rv = [cosd(nanomax.gamma)  0  sind(nanomax.gamma); % detector rotation around vertical axis 
      0                    1  0;
      -sind(nanomax.gamma) 0  cosd(nanomax.gamma)];

Rz = [cosd(0) -sind(0) 0; % detector rotation around beam axis 
      sind(0)  cosd(0)  0;
      0        0         1];

U = Rh*Rv*Rz;
    
qR = (U*q); % correct so far in real space

% Check the coordinate mapping - real space m, Bragg condition
% figure;
% scatter3(qR(1,:),qR(2,:),qR(3,:),10,test(:),'fill','s');
% xlabel('qh');
% ylabel('qv');
% zlabel('qz');
% view(7.7549,53.9161);
% PROVED!

%
%% Initial coordinate of ki
ki = [0,0,k]';

kf = U*ki;

Q = kf-ki;

% Lab coordinate system: accosiated with the ki
QLab(1,:) = (qR(1,:)+Q(1));
QLab(2,:) = (qR(2,:)+Q(2));
QLab(3,:) = (qR(3,:)+Q(3));

% Check the lab space mapping
% figure;
% scatter3(QLab(1,:),QLab(2,:),QLab(3,:),10,scan.data_average(:),'fill','s');
% xlabel('Qx');
% ylabel('Qy');
% zlabel('Qz');
% view(7.7549,53.9161);
% PROVED!

%%
% Small corrections to misalignment of the sample
% Here the rocking curve should be introduced
% alpha
% beta

% Gonphi correction
sample_alpha_correction = 0; % Qx+Qz 
sample_beta_correction = 1.19; % Qz should be positive
sample_gamma_correction = 0; % Qz
dphi = mean(diff(nanomax.gonphi));

for ii = 1:length(nanomax.gonphi)
    % Rotations to bring the q vector into sample coordinate system
    Rsh = [1         0                       0; % detector rotation around horizintal axis 
           0         cosd(nanomax.gontheta+sample_alpha_correction) -sind(nanomax.gontheta+sample_alpha_correction);
           0         sind(nanomax.gontheta+sample_alpha_correction)  cosd(nanomax.gontheta+sample_alpha_correction)]; 

    Rsv = [cosd(-nanomax.gamma/2 + dphi*(ii-1)+sample_beta_correction)  0  sind(-nanomax.gamma/2 + dphi*(ii-1)+sample_beta_correction); % detector rotation around vertical axis 
           0                                                     1  0;
          -sind(-nanomax.gamma/2 + dphi*(ii-1)+sample_beta_correction)  0  cosd(-nanomax.gamma/2 + dphi*(ii-1)+sample_beta_correction)];

    Rsz = [cosd(sample_gamma_correction) -sind(sample_gamma_correction) 0; 
           sind(sample_gamma_correction)  cosd(sample_gamma_correction)  0;
           0        0         1];

    Rs = Rsh*Rsv*Rsz; 

    % Sample coordinate system: accosiated with the ki
    QSample(:,:,ii) = Rs*QLab;
%     Qs = Rs*Q;
% 
%     modQSample = squeeze(sqrt(sum(Qs.^2,1)));
end

% scan.data_meta.QSample = QSample;

% Check the sample space mapping

% figure;
% for ii = 1:length(nanomax.gonphi)
%     Q1 = squeeze(QSample(1,:,ii));
%     Q2 = squeeze(QSample(2,:,ii));
%     Q3 = squeeze(QSample(3,:,ii));
%     t  = data(:,:,ii);
%     
%     scatter3(Q1(:).*1e-10,Q2(:).*1e-10,Q3(:).*1e-10,10,log10(t(:)),'fill','s');
%     xlabel('Qx');ylabel('Qy');zlabel('Qz');
% %     view(7.7549,53.9161);
%     drawnow
% %     waitforbuttonpress
%     hold on;
% end

% view(0,0);
% view(0,-90);
% PROVED!

scaleCoefficient = 1;

qx = squeeze(QSample(1,:,:));
qy = squeeze(QSample(2,:,:));
qz = squeeze(QSample(3,:,:));

dqX = scaleCoefficient*(2*pi*nanomax.detector_pitch/(nanomax.radius*wavelength));
dqY = scaleCoefficient*(2*pi*nanomax.detector_pitch/(nanomax.radius*wavelength));
dqZ = scaleCoefficient*(2*pi*nanomax.detector_pitch/(nanomax.radius*wavelength));

[Qx,Qy,Qz] = meshgrid(min(qx(:)):dqX:max(qx(:)), min(qy(:)):dqY:max(qy(:)),min(qz(:)):dqZ:max(qz(:)));

%% interpolate Q
scan.q_space.data = [];
tic;

F = TriScatteredInterp(qx(:),qy(:),qz(:),double(data(:))); 
lab_data = F(Qx,Qy,Qz);
lab_data(isnan(lab_data))=0;

scan.q_space.data(:,:,:) = abs(lab_data);

toc

%%
scan.q_space.qx = Qx;
scan.q_space.qy = Qy;
scan.q_space.qz = Qz;

save(scan_path,'scan');

%%
% 3D view
% imagery = squeeze(sum((scan.q_space.data),[4,5]));
% val = [0.0001]*max(imagery(:));
% figure('Units','normalized','Position', [0.1 0.1 0.300 0.600]);
% for ii = 1:numel(val)
%     alpha_val = 1-0.5*(ii-1);
%     isosurface(Qx(1,:,1).*1e-10,Qy(:,1,1).*1e-10,squeeze(Qz(1,1,:)).*1e-10,imagery,val(ii));
%     alpha(alpha_val); hold on;
% end
% xlabel('Qx, [A^{-1}]');
% ylabel('Qy, [A^{-1}]');
% zlabel('Qz, [A^{-1}]');
% set(gca,'FontSize',24);
% grid on;axis image;%colormap jet

%% Vis single points
pos = [3,3];
save_path_figures = sprintf('/home/dzhigd/work/projects/CsPbBr3_NC_BCDI_NanoMAX/data/sample0609_%d/figures',scan_number);
mkdir(save_path_figures);

% Reference value for 002 reflection, align the diffraction peak to it
H_CsPbBr3 = 2*3.14/5.8795; % A^-1

min_scale = 0.15;
handles.figHandle = figure('Units','normalized','Position', [0.1 0.1 0.400 0.300]);
imagery = log10(squeeze(sum(scan.q_space.data(:,:,:),3)));
subplot(1,3,1);imagesc(squeeze(Qx(1,:,1)).*1e-10,squeeze(Qy(:,1,1)).*1e-10,imagery,[min_scale*max(imagery(:)),max(imagery(:))]); axis image; title('Integrated intensity');hold on;
xline(H_CsPbBr3,'--','color','white','linewidth',3);
yline(0,'--','color','white','linewidth',3);
xline(H_CsPbBr3,'--','color','white');
xlabel('Qx, [A^{-1}]');
ylabel('Qy, [A^{-1}]');
set(gca,'FontSize',12);

imagery = log10(squeeze(sum(scan.q_space.data(:,:,:),1)));
subplot(1,3,2);imagesc(squeeze(Qz(1,1,:)).*1e-10,squeeze(Qx(1,:,1)).*1e-10,imagery,[min_scale*max(imagery(:)),max(imagery(:))]); axis image; title('Integrated intensity');
xlabel('Qz, [A^{-1}]');
ylabel('Qx, [A^{-1}]');
xline(0,'--','color','white','linewidth',3);
yline(H_CsPbBr3,'--','color','white','linewidth',3);
set(gca,'FontSize',12);

imagery = log10(squeeze(sum(scan.q_space.data(:,:,:),2)));
subplot(1,3,3);imagesc(squeeze(Qz(1,1,:)).*1e-10,squeeze(Qy(:,1,1)).*1e-10,imagery,[min_scale*max(imagery(:)),max(imagery(:))]); axis image; title('Integrated intensity');
colormap hot;
% xline(0,'--','color','white','linewidth',3);
% yline(0,'--','color','white','linewidth',3);
xlabel('Qz, [A^{-1}]');
ylabel('Qy, [A^{-1}]');
set(gca,'FontSize',12);

hFig = gcf;
            set(hFig,'Units','Inches');
            pos = get(hFig,'Position');
            set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
            print(hFig,fullfile(save_path_figures,'qspace_projection.pdf'),'-dpdf','-r0');
            print(hFig,fullfile(save_path_figures,'qspace_projection.png'),'-dpng','-r0');
%             close(hFig);

%% After BCDI
% H_CsPbBr3 = -2*3.14/5.8795; % A^-1
% 
% min_scale = 0.01;
% handles.figHandle = figure('Units','normalized','Position', [0.1 0.1 0.400 0.300]);
% imagery = log10(squeeze(sum(data,3)));
% subplot(1,3,1);imagesc(squeeze(qy),squeeze(qz),imagery,[min_scale*max(imagery(:)),max(imagery(:))]); axis image; title('Integrated intensity');hold on;
% % xline(H_CsPbBr3,'--','color','white','linewidth',3);
% yline(0,'--','color','white','linewidth',3);
% % xline(H_CsPbBr3,'--','color','white');
% xlabel('Qx, [A^{-1}]');
% ylabel('Qy, [A^{-1}]');
% set(gca,'FontSize',12);
% 
% imagery = log10(squeeze(sum(data,1)));
% subplot(1,3,2);imagesc(squeeze(qx),squeeze(qy),imagery,[min_scale*max(imagery(:)),max(imagery(:))]); axis image; title('Integrated intensity');
% xlabel('Qz, [A^{-1}]');
% ylabel('Qx, [A^{-1}]');
% xline(0,'--','color','white','linewidth',3);
% % yline(H_CsPbBr3,'--','color','white','linewidth',3);
% set(gca,'FontSize',12);
% 
% imagery = log10(squeeze(sum(data,2)));
% subplot(1,3,3);imagesc(squeeze(qx),squeeze(qz),imagery,[min_scale*max(imagery(:)),max(imagery(:))]); axis image; title('Integrated intensity');
% colormap hot;
% % xline(0,'--','color','white','linewidth',3);
% % yline(0,'--','color','white','linewidth',3);
% xlabel('Qz, [A^{-1}]');
% ylabel('Qy, [A^{-1}]');
% set(gca,'FontSize',12);
% 
% hFig = gcf;
%             set(hFig,'Units','Inches');
%             pos = get(hFig,'Position');
%             set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
