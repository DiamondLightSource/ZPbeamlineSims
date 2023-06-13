classdef phsource < photon
    properties
        el_h            % Horizontal electron source size [sigma] (um)
        el_v            % Vertical electron source size [sigma] (um)
        el_hdiv         % Horizontal electron source divergence [sigma] (urad)
        el_vdiv         % Vertical electron source divergence [sigma] (urad)
        L               % Length of the photon source (m)
    end
    properties (Dependent)
        ph_h % Horizontal photon source size [sigma] (um)
        ph_v % Vertical photon source size [sigma] (um)
        ph_hdiv % Horizontal photon source divergence [sigma] (urad)
        ph_vdiv % Vertical photon source divergence [sigma] (urad)
    end
   methods
       function value=get.ph_h(phsource)
            value = 1e6.*sqrt((2.*phsource.wavelength.*1e-9.*phsource.L./((2.*pi).^2))+(phsource.el_h.*1e-6).^2);
       end
       function value=get.ph_v(phsource)
            value = 1e6.*sqrt((2.*phsource.wavelength.*1e-9.*phsource.L./((2.*pi).^2))+(phsource.el_v.*1e-6).^2);
       end
       function value=get.ph_hdiv(phsource)
            value = 1e6.*sqrt((phsource.wavelength.*1e-9./(2*phsource.L))+(phsource.el_hdiv*1e-6)^2);
       end
       function value=get.ph_vdiv(phsource)
            value = 1e6.*sqrt((phsource.wavelength.*1e-9./(2*phsource.L))+(phsource.el_vdiv*1e-6)^2);
       end
       function obj = phsource(val)
           obj@photon(val);
            if exist('val','var')
               for i=2:6
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

% jim=pgm
% jim.energy=[60:0.1:64]'
% jim.order=1
% jim.N=[700:50:1000]
% jim.cff(1,1,:)=[2:0.1:3]