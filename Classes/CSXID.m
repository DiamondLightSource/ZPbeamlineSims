classdef CSXID < RP_cPGM_CSXID & photon_CSXID & pgm_CSXID & phsource_CSXID
    properties (Constant)
        ZP_q_min=0.003;             % Minimum image distance for zone plate (m)
        ZP_OSA_sample_min=0.0015;	% Minimum distance between OSA and sample (m)
        ZP_D_max=0.29;              % Maximum diameter of zone plate (mm)
        M1_M3_dist_min=4;           % Minimum distance between M1 toroid and M3 cylinder (m)
        BL_length_max=60;         % Maximum beamline length (m)
        ZP_eff_ratio_min=0.0001;     % Minimum value for (ZP area)/(beam size at zone plate)
        ZP_h_min=0.001;              % Minimum horizontal beam size at ZP to ensure ZP is fully-illuminated
        ZP_v_min=0.001;              % Minimum horizontal beam size at ZP to ensure ZP is fully-illuminated
        ZP_p_min=1.0;                 % Minimum p for possible steering mirror (m)
    end
    properties (Constant,Hidden)
        Rayl_percentage=85.171; % If there is no central beam stop, 85.171% of the total intensity 
                                % lies inside the diameter defined by the Rayleigh criterion 
                                % (s_Rayl = 1.21966989*ZP_dr)
    end
    properties
        ZP_dr           % Width of outermost zone on zone plate (nm)
        ZP_p            % Object distance for zone plate (m)
        ZP_eps          % Central beam stop diameter (or OSA diameter) divided by ZP_D
    end
    properties (Dependent)
        ZP_f            % focal length for zone plate (mm)
        ZP_q            % Image distance for zone plate (m)
        ZP_OSA_sample   % Distance between Order Sorting Aperture and sample (mm)
        ES_h            % Horizonal beam size at exit slit [FWHM] (mm)
        ZP_h            % Horizonal beam size at zone plate [FWHM] (mm)
        ZP_v            % Vertical beam size at zone plate [FWHM] (mm)
        ZP_D            % Diameter of zone plate (mm)
        ZP_N            % Number of zones in zone plate
        s_Rayl          % Rayleigh resolution limit for zone plate image (nm)
        s_source_h      % Horizontal source contribution for zone plate image (nm)
        s_source_v      % Vertical source contribution for zone plate image (nm)
        s_chrom     	% Chromatic aberration contribution for zone plate image (nm)
        s_h             % Horizontal size of zone plate image (nm)  
        s_v             % Vertical size of zone plate image (nm)
        s_A             % Cross-sectional area of zone plate image (nm^2)
        BL_length       % Total length of beamline (m)
        ZP_eff_ratio    % (ZP area)/(beam size at zone plate)
    end
    methods
     	function value=get.ZP_f(CSXID)
            value = (CSXID.ZP_D.*CSXID.ZP_dr-1e-6.*CSXID.ZP_dr.^2)./CSXID.wavelength;
        end
        function value=get.ZP_q(CSXID)
            value = 1./((1./(1e-3.*CSXID.ZP_f))-1./CSXID.ZP_p);
            value(value<=CSXID.ZP_q_min)=NaN;
        end
        function value=get.ZP_OSA_sample(CSXID)
            value = CSXID.ZP_eps.*CSXID.ZP_q;
        end
        function value=get.ES_h(CSXID)
            D = CSXID.r4./(CSXID.r1+CSXID.r2-CSXID.r4);
            value = 1e-3.*CSXID.FWHM.*D.*CSXID.ph_h;
        end
        function value=get.ZP_h(CSXID)
            % Calculates assuming exit slit is wide open in horizontal
%             value = 1e-3.*sqrt((CSXID.FWHM.*CSXID.ph_h.*CCSXID.r4./(CSXID.r1+CSXID.r2-CSXID.r4)).^2+(CSXID.FWHM.*CSXID.ZP_p.*CSXID.ph_hdiv.*(CSXID.r1+CSXID.r2-CSXID.r4)./CSXID.r4)).^2);
%              Calculates assuming exit slit is a square (es by es)
            value = 1e-3.*sqrt(CSXID.es.^2+(CSXID.FWHM.*CSXID.ZP_p.*CSXID.ph_hdiv.*(CSXID.r1+CSXID.r2-CSXID.r4)./CSXID.r4).^2);
        end
        function value=get.ZP_v(CSXID)
            value = 1e-3.*sqrt(CSXID.es.^2+(CSXID.FWHM.*CSXID.ZP_p.*CSXID.ph_vdiv.*CSXID.cff.*CSXID.r1./CSXID.r3).^2);
        end
        function value=get.ZP_D(CSXID)
            % Sets ZP diameter equal to TWICE vertical coherence length 
            % considering slit fully illuminated by fully incoherent light 
            % (see "Soft X-rays and Extreme Ultraviolet Radiation" by 
            % Attwood for more information)
%             value = CSXID.ZP_p.*CSXID.wavelength./(pi.*CSXID.es);
%             value = CSXID.ZP_p.*CSXID.wavelength./(pi.*(CSXID.es/CSXID.FWHM));
            % Sets ZP diameter equal to TWICE horizontal coherence length 
            % considering fully incoherent Gaussian source at exit slit 
            % (assumes exit slit is completely open horizontally)  
%             value = CSXID.FWHM.*CSXID.ZP_p.* ...
%             CSXID.wavelength./(2.*pi.*1e3.*CSXID.ES_h);
%             Sets ZP diameter equal to FWHM of degree of coherence 
%             considering slit fully illuminated by fully incoherent light
            value = 4.44.*CSXID.ZP_p.*CSXID.wavelength./(2.*pi.*CSXID.es);
%             % Sets ZP diameter equal to FWHM of vertical beam at ZP
%             value=CSXID.ZP_v;
%             % Sets ZP diameter equal to FWHM of horizontal beam at ZP
%             value=CSXID.ZP_h;
%             % Sets ZP diameter equal to maximum value
%             value=CSXID.ZP_D_max;
%             value(value>CSXID.ZP_D_max)=NaN;
        end
        function value=get.ZP_N(CSXID)
            value = 1e6*(sqrt(4*CSXID.ZP_f.^2+CSXID.ZP_D.^2)-2*CSXID.ZP_f)./CSXID.wavelength;
        end
        function value=get.s_Rayl(CSXID)
%             value = 1.21966989.*CSXID.wavelength.*CSXID.ZP_f./CSXID.ZP_D;
            dr=CSXID.ZP_dr; eps=CSXID.ZP_eps;
            r_2D=repmat((1e-10:0.1:20*max(dr))',size(dr));
            r=repmat(r_2D,size(eps));
            x=pi.*r./dr;
            term1=2.*besselj(1,x)./(r+realmin);
            term2=eps.^2*2.*besselj(1,eps.*x)./(eps.*r+realmin);
            I=(term1-term2).^2;
            run_int=cumsum(I,1)./sum(I,1);
            value=0.*dr.*eps;
            for i=1:length(dr)
                for j=1:length(eps)
                    value(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,i,:,j) = ...
                        2.*r(logical(diff(run_int(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,i,:,j)>0.01*CSXID.Rayl_percentage)));
                end
            end
        end
        function value=get.s_source_h(CSXID)
              % Calculates assuming exit slit is wide open in horizontal
%             value = 1e6.*CSXID.ES_h.*CSXID.ZP_q./CSXID.ZP_p;
             % Calculates assuming exit slit is a square (es by es)
            value = 1e3.*CSXID.es.*CSXID.ZP_q./CSXID.ZP_p;
        end
        function value=get.s_source_v(CSXID)
            value = 1e3.*CSXID.es.*CSXID.ZP_q./CSXID.ZP_p;
        end
        function value=get.s_chrom(CSXID)
            value = 1e6.*CSXID.ZP_D./CSXID.R;
        end
        function value=get.s_h(CSXID)
            cond1=CSXID.BL_length<=CSXID.BL_length_max;
         	cond2=(CSXID.r2-CSXID.r3)>=CSXID.M1_M3_dist_min;
         	cond3=CSXID.ZP_eff_ratio>=CSXID.ZP_eff_ratio_min;
            cond4=CSXID.ZP_h>=CSXID.ZP_h_min;
            cond5=CSXID.ZP_v>=CSXID.ZP_v_min;
            cond6=CSXID.ZP_OSA_sample>=CSXID.ZP_OSA_sample_min;
            cond7=CSXID.ZP_p>=CSXID.ZP_p_min;
            value = sqrt(CSXID.s_Rayl.^2+ ...
                        CSXID.s_source_h.^2+ ...
                        CSXID.s_chrom.^2)./ ...
                        (cond1.*cond2.*cond3.*cond4.*cond5.*cond6.*cond7);
        end
        function value=get.s_v(CSXID)
            cond1=CSXID.BL_length<=CSXID.BL_length_max;
         	cond2=(CSXID.r2-CSXID.r3)>=CSXID.M1_M3_dist_min;
         	cond3=CSXID.ZP_eff_ratio>=CSXID.ZP_eff_ratio_min;
            cond4=CSXID.ZP_h>=CSXID.ZP_h_min;
            cond5=CSXID.ZP_v>=CSXID.ZP_v_min;
            cond6=CSXID.ZP_OSA_sample>=CSXID.ZP_OSA_sample_min;
            cond7=CSXID.ZP_p>=CSXID.ZP_p_min;
            value = sqrt(CSXID.s_Rayl.^2+ ...
                        CSXID.s_source_v.^2+ ...
                        CSXID.s_chrom.^2)./ ...
                        (cond1.*cond2.*cond3.*cond4.*cond5.*cond6.*cond7);
        end
        function value=get.s_A(CSXID)
            value = pi*CSXID.s_h.*CSXID.s_v./4;
        end
        function value=get.BL_length(CSXID)
            value = CSXID.r1+CSXID.r2+...
                    CSXID.ZP_p+CSXID.ZP_q;
        end
        function value=get.ZP_eff_ratio(CSXID)
%             value = CSXID.ZP_D.^2./(CSXID.ZP_h.*CSXID.ZP_v);
            cond1=CSXID.BL_length<=CSXID.BL_length_max;
         	cond2=(CSXID.r2-CSXID.r3)>=CSXID.M1_M3_dist_min;
            cond4=CSXID.ZP_h>=CSXID.ZP_h_min;
            cond5=CSXID.ZP_v>=CSXID.ZP_v_min;
            cond6=CSXID.ZP_OSA_sample>=CSXID.ZP_OSA_sample_min;
            cond7=CSXID.ZP_p>=CSXID.ZP_p_min;
            value = CSXID.ZP_D.^2./(CSXID.ZP_h.*CSXID.ZP_v)./ ...
                        (cond1.*cond2.*cond4.*cond5.*cond6.*cond7);

        
        end
        function obj = CSXID(val)
            obj@RP_cPGM_CSXID(val);
            obj@photon_CSXID(val);
            obj@pgm_CSXID(val);
            obj@phsource_CSXID(val);
            if exist('val','var')
               for i=2:numel(obj.inputpropnames)
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
