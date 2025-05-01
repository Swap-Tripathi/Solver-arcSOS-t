function [status,GS0Wout, GS12Wout, GVout] = solveSDPAlternatived3(dimsGV,r,d,L)
%for d=4 use the below statement rather than the above one
%function [status,GS0Wout, GS12Wout, GS23Wout, GS13Wout, GVout] = solveSDPAlternative(dimsGV,r,d,L)

    dimsGS0W = dimsGV + L*ones ;           % n_v+L.1
    dimsGSLW = dimsGV + (L-r)*ones ;       % n_v+L.1-r.1
    
    cvx_begin sdp quiet
    cvx_precision low

    variable GS0W(prod(dimsGS0W+ones), prod(dimsGS0W+ones)) hermitian
    variable GS12W(prod(dimsGSLW+ones), prod(dimsGSLW+ones)) hermitian
    %variable GS23W(prod(dimsGSLW+ones), prod(dimsGSLW+ones)) hermitian
    %variable GS13W(prod(dimsGSLW+ones), prod(dimsGSLW+ones)) hermitian
    variable GV(prod(dimsGV+ones), prod(dimsGV+ones)) hermitian
    
    bb = max(dimsGS0W);                             
    B = cell(1,d-1);
    [B{:}] = ndgrid(-bb:bb); 
    B = cellfun(@(M) M(:), B, 'uniform', 0);
    S = [B{end:-1:1}]';                             % All the steps in total produce all possible index vectors for the trig poly corresponding to S_l^w as columns
    
    Iep=eye(prod(dimsGV+ones));
    Iep(1,1)=0;

    %% Define the SDP for Arcs as existence problem.

    
    minimize(1); % +10*norm(GS12W(:),1) - 0.1*log_det(GS12W + epsilon*eye)
    subject to
        
        epsilon = 0.0; 
        (GS0W)>=0;
        (GS12W)>=0;
        %(GS23W)>=0;
        %(GS13W)>=0;
        (GV)>=epsilon*Iep;
        
        %% Define terms appearing in the main equality corresponding to Lyapunov-like condition eq (21)

        for k = S(:,1:(end+1)/2)
           
            %calculate LHS of theorem: Local stability on Arcs

            lhs = 0;
           
           lhs = lhs + (TrFind(GS0W, dimsGS0W, k))+ 0.5 * (TrFind(GS12W, dimsGSLW, k+[-r,r])+TrFind(GS12W, dimsGSLW, k+[r,-r]));
                                                          %(TrFind(GS12W, dimsGSLW, k+[-r,r,0])+TrFind(GS12W, dimsGSLW, k+[r,-r,0]) ...
                                                          %+ TrFind(GS23W, dimsGSLW, k+[0,-r,r])+TrFind(GS23W, dimsGSLW, k+[0,r,-r]) ...
                                                          %+ TrFind(GS13W, dimsGSLW, k+[-r,0,r])+TrFind(GS13W, dimsGSLW, k+[r,0,-r]));
           
           %calculate RHS of theorem: Local stability on Arcs
           
           rhs = 0;
           
           for c = 1:d-1
               for j = S(:,:)
                   % Compute the RHS of eq 26 as per the summation
                   rhs = rhs - 1i * j(c)* F(c,k-j) * ( TrFind(GV, dimsGV, j) );
               end
           end
       
           % Now add the constraint to the CVX problem
           lhs == rhs;
        end
       
        % Add other constraint to matrices so that the solution does not go
        % close to zero matrix
        trace(GS0W) >= 64;
        trace(GS12W) >= 36;
        %trace(GS23W) >= 12;
        %trace(GS13W) >= 12;
        sum(GS0W(:)) == 0;  %S0W(theta)=0 at origin
        sum(GS12W,2) == 0;
        %sum(GS12W(:)) == 0; %S12W(theta)=0 at origin
        %sum(GS23W(:)) == 0; %S23W(theta)=0 at origin
        %sum(GS13W(:)) == 0; %S13W(theta)=0 at origin
        sum(GV(:))==0;

    cvx_end
    %% SDP for Arcs ends.
    %[P12,D12]   = eig(GS12W);
    %Check_exactly_single_zero= (max(abs(P12))-min(abs(P12)))/max(abs(P12));
    %[P0,D0]     = jordan(GS0W);
    %opt_value   = cvx_optval;
    status      = cvx_status;
    GS0Wout     = GS0W;
    GS12Wout    = GS12W;
    %GS23Wout    = GS23W;
    %GS13Wout    = GS13W;
    GVout       = GV;
end