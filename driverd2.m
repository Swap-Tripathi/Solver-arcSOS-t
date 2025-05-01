clc; clear; clear cvx; close; close all;


%% example for d  = 3

r=3;                                    % to check stability in arc(pi/2r)
d=3;                                    % order of original system
L=3;                                    % number of Harmonics of coupling function
harm=4;                                 % >=r-L; will be later used to define nv=harm*[1,1]
global alphaVec;                        % global variable storing alpha_c for c=1,...,d 
global betaVec;                         % global variable storing beta_c for c=1,..,d
alphaVec=[4 2 4];                 % length equal to L. All entries nonnegative
betaVec=[pi/8 -pi/4 -pi/3];               % all betas have absolute values less than pi/2

if ~isequal(size(alphaVec),size(zeros(1,L))) || ~isequal(size(betaVec),size(zeros(1,L)))
    error ('check sizes of alphaVec and betaVec')
else
    for c=1:L
        if alphaVec(c)<0 || abs(betaVec(c))>=pi/2
            error ('check bounds for alphas and betas')
        end
    end
end
    dimsGV = harm * [1,1];               % dimsGV=n_v. will be used to define matrix size
    [val_four,GS0W, GS12W, GV]  = solveSDPAlternatived3(dimsGV, r, d, L);

    % for d=4, use the following instead of the above
    % [status,GS0W, GS12W, GS23W, GS13W, GV]  = solveSDPAlternative(dimsGV, r, d, L);
    disp(harm)
    disp(val_four)
    disp('===========')