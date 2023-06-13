classdef RP_cPGM < photon & pgm & phsource
    properties
        r1      % Distance from the source to the collimating mirror (m)
        r2    	% Distance from the collimating mirror to the exit slit (m)
        r3  	% Distance from the focusing mirror to the exit slit (m)
        CMa 	% Incidence angle for the Collimating Mirror (CM) relative to normal (deg)
        FMa 	% Incidence angle for the Focusing Mirror (FM) relative to normal (deg)
        es    	% Exit slit vertical opening (um)
        CMs  	% Collimating Mirror slope error in Sagittal direction (urad)
        FMs     % Focusing Mirror slope error in Sagittal direction (urad)
        PMt     % Plane Mirror slope error in Tangential direction (urad)
        GRt     % GRating slope error in Tangential direction (urad)
    end
    properties (Dependent)
       Vss
       Vcm
       Vpm
       Vgr
       Vfm
       Ves
       Vac
%        Vab
       V
       Rss
       Rcm
       Rpm
       Rgr
       Rfm
       Res
       Rac
%        Rab
       R
       deltaE
    end
    properties (Dependent, Hidden)
        U           % parameter which enables easier 
                    % calculation of RPs
    end
    properties (Constant, Hidden)
        FWHM=2*sqrt(2*log(2));	% FWHM is the scale factor
        % between rms and FWHM (~2.355)
    end
    methods
        function value=get.Vss(RP_cPGM)
            value = RP_cPGM.r3.*RP_cPGM.ph_v.*1e-6./(RP_cPGM.r1.*RP_cPGM.cff);
        end
        function value=get.Vcm(RP_cPGM)
            value = 2.*RP_cPGM.r3.*RP_cPGM.CMs.*1e-6.*cosd(RP_cPGM.CMa)./RP_cPGM.cff;
        end
        function value=get.Vpm(RP_cPGM)
            value = 2.*RP_cPGM.r3.*RP_cPGM.PMt.*1e-6./RP_cPGM.cff;
        end
        function value=get.Vgr(RP_cPGM)
            value = RP_cPGM.r3.*RP_cPGM.GRt.*1e-6.*(1+(1./RP_cPGM.cff));
        end
        function value=get.Vfm(RP_cPGM)
            value = 2.*RP_cPGM.r3.*RP_cPGM.FMs.*1e-6.*cosd(RP_cPGM.FMa);
        end
        function value=get.Ves(RP_cPGM)
            value = 1e-6.*RP_cPGM.es./RP_cPGM.FWHM;
        end
        function value=get.Vac(RP_cPGM)
            value = 0.5.*RP_cPGM.cff.*RP_cPGM.ph_hdiv.*1e-6.*RP_cPGM.ph_vdiv.*1e-6.* ...
                (RP_cPGM.r1.*RP_cPGM.r2.*realmin./RP_cPGM.r3).*tand(RP_cPGM.CMa); % realmin makes this term very small for now
        end
        function value=get.V(RP_cPGM)
            value = 1e6*sqrt(RP_cPGM.Vss.^2+RP_cPGM.Vcm.^2+RP_cPGM.Vpm.^2+ ...
                RP_cPGM.Vgr.^2+RP_cPGM.Vfm.^2+RP_cPGM.Ves.^2+RP_cPGM.Vac.^2);
        end
        function value=get.U(RP_cPGM)
            value = RP_cPGM.N.*1e3.*RP_cPGM.order.*RP_cPGM.wavelength.*1e-9.*RP_cPGM.r3./ ...
                (RP_cPGM.FWHM.*cosd(RP_cPGM.beta));
        end
        function value=get.Rss(RP_cPGM)
            value = RP_cPGM.U./(RP_cPGM.Vss);
        end
        function value=get.Rcm(RP_cPGM)
            value = RP_cPGM.U./(RP_cPGM.Vcm);
        end
        function value=get.Rpm(RP_cPGM)
            value = RP_cPGM.U./(RP_cPGM.Vpm);
        end
        function value=get.Rgr(RP_cPGM)
            value = RP_cPGM.U./(RP_cPGM.Vgr);
        end
        function value=get.Rfm(RP_cPGM)
            value = RP_cPGM.U./(RP_cPGM.Vfm);
        end
        function value=get.Res(RP_cPGM)
            value = RP_cPGM.U./(RP_cPGM.Ves);
        end
        function value=get.Rac(RP_cPGM)
            value = RP_cPGM.U./(RP_cPGM.Vac);
        end
        function value=get.R(RP_cPGM)
            value = ((1./RP_cPGM.Rss).^2+(1./RP_cPGM.Rcm).^2+...
                (1./RP_cPGM.Rpm).^2+(1./RP_cPGM.Rgr).^2+...
                (1./RP_cPGM.Rfm).^2+(1./RP_cPGM.Res).^2+...
                (1./RP_cPGM.Rac).^2).^(-0.5);
        end
        function value=get.deltaE(RP_cPGM)
            value = 1e3.*RP_cPGM.energy./RP_cPGM.R;
        
        end
        function obj = RP_cPGM(val)
            obj@photon(val);
            obj@pgm(val);
            obj@phsource(val);
            if exist('val','var')
               for i=2:19
                   eval(['if isfield(val,''' obj.inputpropnames{i} ''') || isprop(val,''' obj.inputpropnames{i} ''')' newline ...
                       'obj.string=repmat(' '''1,''' ',1,i-1);' newline ...
                       'end'])
                   eval(['if isfield(val,''' obj.inputpropnames{i} ''') || isprop(val,''' obj.inputpropnames{i} ''')' newline ...
                       'obj.' obj.inputpropnames{i} '=zeros(' obj.string 'numel(val.' obj.inputpropnames{i} '));' newline ...
                       'obj.' obj.inputpropnames{i} '(' obj.string ':)=val.' obj.inputpropnames{i} ';' newline ...
                       'end'])
               end  
            end
        end
   end
end

% 12th August 2019
% EXAMPLE USAGE
% a.r1=26.162;a.r2=26.162;a.r3=11;a.CMa=87;a.FMa=87;a.es=5:5:100;
% a.CMs=5;a.FMs=5;a.PMt=[1e-10:0.01:1];a.GRt=[1e-10:0.01:0.5];
% a.order=1;a.N=400:100:2000;a.cff=2:0.25:5;a.energy=[20:5:200];
% a.el_h=178.4;a.el_v=12.6;a.el_hdiv=16.5;a.el_vdiv=2.2;a.L=5;
% b=RP_cPGM(a);

% [X,Y]=meshgrid(squeeze(b.GRt),squeeze(b.PMt));
% h=pcolor(X,Y,squeeze(b.V(1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1, 1,:,:)));
% set(h, 'EdgeColor', 'none');
% colorbar;

% [X,Y]=meshgrid(squeeze(b.r1),squeeze(b.ph_v(1,1,1,1,:)));
% h=pcolor(X,Y,squeeze(b.V(1,1,1,1,:,1,1,1,:)));
% set(h, 'EdgeColor', 'none');
% colorbar;


% 13th August 2019
% a.r1=26.162;a.r2=26.162;a.r3=11;a.CMa=87;a.FMa=87;a.es=5;
% a.CMs=5;a.FMs=5;a.PMt=0.16;a.GRt=[0.05:0.01:0.5];
% a.order=1;a.N=400:100:3000;a.cff=2:0.1:5;a.energy=[20:5:200];
% a.el_h=178.4;a.el_v=12.6;a.el_hdiv=16.5;a.el_vdiv=2.2;a.L=5;
% b=RP_cPGM(a);
