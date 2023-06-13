a.energy=3000;          % Photon energy (eV)
a.el_h=35.91;           % Horizontal electron source size [sigma] (um)
a.el_v=5.292;           % Vertical electron source size [sigma] (um)
a.el_hdiv=4.327;        % Horizontal electron source divergence [sigma] (urad)
a.el_vdiv=1.512;        % Vertical electron source divergence [sigma] (urad)
a.L=5;                  % Length of the photon source (m)
a.order=1;              % Diffraction order (normally 1)
a.N=1600;               % Line density of grating (lines/mm)
a.cff=2.4;              % cos(beta)/cos(alpha) [alpha/beta is Grating incidence/diffraction angle]
a.r1=25.500;            % Distance from the source to the collimating mirror (m)
a.r2=15.000;            % Distance from the collimating mirror to the exit slit (m)
a.r3=9.000;             % Distance from the VFM to the exit slit (m)
a.r4=6.000;             % Distance from the HFM to the exit slit (m)
a.CMa=89.4;             % Incidence angle for the Collimating Mirror (CM) relative to normal (deg)
a.VFMa=89.3;            % Incidence angle for the Focusing Mirror (FM) relative to normal (deg)
a.HFMa=89.4;            % Incidence angle for the Focusing Mirror (FM) relative to normal (deg)
a.es=10;                % Exit slit vertical opening (um)
a.CMs=5;            	% Collimating Mirror slope error in Sagittal direction (urad)
a.VFMs=5;               % Focusing Mirror slope error in Sagittal direction (urad)
a.HFMs=5;               % Focusing Mirror slope error in Sagittal direction (urad)
a.PMt=0.2;              % Plane Mirror slope error in Tangential direction (urad)
a.GRt=0.2;              % GRating slope error in Tangential direction (urad)
a.ZP_dr=20;             % Width of outermost zone on zone plate (nm)
a.ZP_p=4;               % Distance between exit slit and zone plate (m) 
a.ZP_eps=0.5;           % Central beam stop diameter (or OSA diameter) divided by ZP_D 

% Create nanoARPES class from parameters defined above
tic
c=CSXID(a)
toc

% Find minimum beam size at sample position
min(c.s_v,[],'all')

% Provide index of element in up to 22-dimensional array which gives minimum
% sample beamsize
[I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I13,I14,I15,I16,I17,I18,I19,I20,I21,I22,I23,I24,I25]= ...
    ind2sub(size(c.s_v),find(c.s_v==min(c.s_v,[],'all')));
format shortG
[a.energy(I1);a.el_h(I2);a.el_v(I3);a.el_hdiv(I4);a.el_vdiv(I5);a.L(I6); ...
    a.order(I7);a.N(I8);a.cff(I9);a.r1(I10);a.r2(I11);a.r3(I12);a.r4(I13); ...
    a.CMa(I14);a.VFMa(I15);a.HFMa(I16);a.es(I17);a.CMs(I18);a.VFMs(I19); ...
    a.HFMs(I20);a.PMt(I21);a.GRt(I22);a.ZP_dr(I23);a.ZP_p(I24);a.ZP_eps(I25)]'
format short
% Print sample beamsize
c.s_v(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12, ...
    I13,I14,I15,I16,I17,I18,I19,I20,I21,I22,I23,I24,I25)
ZP_eff_ratio=c.ZP_eff_ratio(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12, ...
    I13,I14,I15,I16,I17,I18,I19,I20,I21,I22,1,I24)