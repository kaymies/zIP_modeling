%% CLASSLINK Helper class for generating lumped parameters for
%  the double-inverted-pendulum model of human standing
%
% Author: Jongwoo Lee, PhD.
% Modified by: Rika Sugimoto-Dimitrova (rikasd@mit.edu)
% 2025-07-01

classdef ClassLink

    properties
        L_ratio_F % link length as a fraction of total body height for female model
        L_ratio_M % link length as a fraction of total body height for male model
        m_ratio_F % link mass as a fraction of total body mass for female model
        m_ratio_M % link mass as a fraction of total body mass for male model
        c_ratio_F % link center-of-mass position as a fraction of total body height for female model
        c_ratio_M % link center-of-mass position as a fraction of total body height for male model
        r_ratio_frt_F % link frontal-plane radius of gyration normalized by total body height for female model
        r_ratio_frt_M % link frontal-plane radius of gyration normalized by total body height for male model
        r_ratio_sgt_F % link sagittal-plane radius of gyration normalized by total body height for female model
        r_ratio_sgt_M % link sagittal-plane radius of gyration normalized by total body height for male model
        
        m % link mass (kg)
        L % link length (m)
        c % link center-of-mass position from joint (m)
        j % link moment of inertia about center-of-mass (kgm^2)
    end

    methods
        function obj = ClassLink(L_F, L_M, m_ratio_F, m_ratio_M, c_ratio_F, c_ratio_M, r_ratio_frt_F, r_ratio_frt_M, r_ratio_sgt_F, r_ratio_sgt_M)
            HeightM_mm = 1741; % total body height for average young male (de Leva 1996)
            HeightF_mm = 1735; % total body height for average young female (de Leva 1996)
            if nargin > 1
                obj.L_ratio_F = L_F/HeightF_mm;
                obj.L_ratio_M = L_M/HeightM_mm;
                obj.m_ratio_F = m_ratio_F/100;
                obj.m_ratio_M = m_ratio_M/100;
                obj.c_ratio_F = c_ratio_F/100;
                obj.c_ratio_M = c_ratio_M/100;
                obj.r_ratio_frt_F = r_ratio_frt_F/100;
                obj.r_ratio_frt_M = r_ratio_frt_M/100;
                obj.r_ratio_sgt_F = r_ratio_sgt_F/100;
                obj.r_ratio_sgt_M = r_ratio_sgt_M/100;
            end
        end

        function obj = assignLink(obj, totalMass, totalHeight, gender, plane)
            switch gender
                case 'F'
                    obj.m = obj.m_ratio_F*totalMass;
                    obj.L = obj.L_ratio_F*totalHeight;
                    obj.c = obj.c_ratio_F*obj.L;
                    switch plane
                        case 'sgt'
                            obj.j = obj.m*(obj.L*obj.r_ratio_sgt_F)^2;
                        case 'frt'
                            obj.j = obj.m*(obj.L*obj.r_ratio_frt_F)^2;
                    end
                case 'M'
                    obj.m = obj.m_ratio_M*totalMass;
                    obj.L = obj.L_ratio_M*totalHeight;
                    obj.c = obj.c_ratio_M*obj.L;
                    switch plane
                        case 'sgt'
                            obj.j = obj.m*(obj.L*obj.r_ratio_sgt_M)^2;
                        case 'frt'
                            obj.j = obj.m*(obj.L*obj.r_ratio_frt_M)^2;
                    end
            end
        end

    end
end
