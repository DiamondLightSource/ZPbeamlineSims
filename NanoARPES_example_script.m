a.energy=250;           % Photon energy (eV)
a.el_h=35.42;           % Horizontal electron source size [sigma] (um)
a.el_v=5.706;           % Vertical electron source size [sigma] (um)
a.el_hdiv=3.980;        % Horizontal electron source divergence [sigma] (urad)
a.el_vdiv=1.402;        % Vertical electron source divergence [sigma] (urad)
a.L=5;                  % Length of the photon source (m)
a.order=1;              % Diffraction order (normally 1)
a.N=800;               % Line density of grating (lines/mm)
a.cff=2.25;             % cos(beta)/cos(alpha) [alpha/beta is Grating incidence/diffraction angle]
a.r1=26.162;            % Distance from the source to the collimating mirror (m)
a.r2=13.9;                % Distance from the collimating mirror to the exit slit (m)
a.r3=9.1;                % Distance from the focusing mirror to the exit slit (m)
a.CMa=87;               % Incidence angle for the Collimating Mirror (CM) relative to normal (deg)
a.FMa=87;               % Incfidence angle for the Focusing Mirror (FM) relative to normal (deg)
a.es=10;                % Exit slit vertical opening (um)
a.CMs=3.13;             % Collimating Mirror slope error in Sagittal direction (urad)
a.FMs=5;                % Focusing Mirror slope error in Sagittal direction (urad)
a.PMt=0.16;             % Plane Mirror slope error in Tangential direction (urad)
a.GRt=0.25;             % GRating slope error in Tangential direction (urad)
a.ZP_dr=21;             % Width of outermost zone on zone plate (nm)
a.ZP_p=7.7;               % Distance between exit slit and zone plate (m) 
a.ZP_eps=0.5;           % Central beam stop diameter (or OSA diameter) divided by ZP_D 

% Create nanoARPES class from parameters defined above
tic
c=nanoARPES(a)
toc

Find minimum beam size at sample position
min(c.s_v,[],'all')

% Provide index of element in up to 22-dimensional array which gives minimum
% sample beamsize
[I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I13,I14,I15,I16,I17,I18,I19,I20,I21,I22]= ...
    ind2sub(size(c.s_v),find(c.s_v==min(c.s_v,[],'all')));
format shortG
[a.energy(I1);a.el_h(I2);a.el_v(I3);a.el_hdiv(I4);a.el_vdiv(I5);a.L(I6);a.order(I7);...
    a.N(I8);a.cff(I9);a.r1(I10);a.r2(I11);a.r3(I12);a.CMa(I13);a.FMa(I14);a.es(I15); ...
    a.CMs(I16);a.FMs(I17);a.PMt(I18);a.GRt(I19);a.ZP_dr(I20);a.ZP_p(I21);a.ZP_eps(I22)]'
format short
% Print sample beamsize
c.s_v(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I13,I14,I15,I16,I17,I18,I19,I20,I21,I22)