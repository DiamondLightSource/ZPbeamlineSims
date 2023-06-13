classdef photon_CSXID
	properties
        energy          % Photon energy (eV)
    end
    properties (Dependent)
        wavelength	% Photon wavelength (nm)
    end
    properties (Constant, Hidden)
        factor=1239.841984;     % factor=h*c/e, where
                                % h=6.62607015e-34, c=299792458 and
                                % e=1.602176634e-19
        inputpropnames={'energy','el_h','el_v','el_hdiv','el_vdiv','L','order', ...
            'N','cff','r1','r2','r3','r4','CMa','VFMa','HFMa', ...
            'es','CMs','VFMs','HFMs','PMt','GRt', ...
            'ZP_dr','ZP_p','ZP_eps'};  % Names of the input properties
    end
    properties (Hidden)
      	string    	% temporary string used to create multi-dimensional arrays
    end
    methods
        function value=get.wavelength(photon)
            value=photon.factor./photon.energy;
        end
        function obj = photon_CSXID(val)
            if exist('val','var')
                if isfield(val,'energy') || isprop(val,'energy')
                    obj.energy=zeros(numel(val.energy),1);
                    obj.energy(:,1)=val.energy;
                end
            end
        end
    end
end