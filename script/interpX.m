function H = interpX( t0, H0, t, method )
%INTERPX implements interpolations of SE(3) and PCG(3): SO(3)xR^3
% Inputs:
%     t0     ----- original time steps
%     H0     ----- original frames
%     t      ----- new time steps
%     method ----- 'SE' or 'PCG' interpolations
%     
% Output:
%     H      ----- interpolated frames
%
% Author: Sipu Ruan, ruansp@jhu.edu, September 2017

H = nan(size(H0,1), size(H0,2), length(t));
for i = 1:length(t)
    % interpolate between two adjacent time steps
    dt = t(i)-t0;
    [~,idx] = min(dt(dt>=0));
    if isempty(idx)
        continue;
    end
    if idx ~= length(t0)
        Ho_oa = H0(:,:,idx);
        Ho_ob = H0(:,:,idx+1);
        Ho_ab = Ho_ob / Ho_oa;
        dt0 = t0(idx+1)-t0(idx);
        
        if strcmp(method, 'SE')
            % Interpolate on SE(3)
            Xi = (dt(idx)/dt0) * logm(Ho_ab);
            
            H(:,:,i) = expm(Xi) * Ho_oa;

        elseif strcmp(method, 'PCG')
            % Extract rotation and translation parts
            Ro_ab = Ho_ab(1:3,1:3); do_ab = Ho_ab(1:3,4);
            Ro_oa = Ho_oa(1:3,1:3); to_oa = Ho_oa(1:3,4); 
            to_ab = do_ab - (eye(3)-Ro_ab)*to_oa;
            
            % Interpolate rotation and translation parts separately
            Rpcg(:,:,i) = expm(dt(idx)/dt0 * logm(Ro_ab)) * Ro_oa;
            tpcg(:,:,i) = dt(idx)/dt0 * to_ab + to_oa;
            
            H(:,:,i) = [Rpcg(:,:,i) tpcg(:,:,i); zeros(1,3) 1];
        end
    else
        H(:,:,i) = H0(:,:,idx);
    end
end
end

