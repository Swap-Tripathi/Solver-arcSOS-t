function y = F(l,i)

   
global alphaVec;                        % global variable storing alpha_c for c=1,...,L 
global betaVec;                         % global variable storing beta_c for c=1,..,L


    %for one oscillator, use the following
    %FF   = [-1      1       -1j*alphaVec(1) * exp(1j * betaVec(1))/2        1j*alphaVec(1)*exp(-1j*betaVec(1))/2;
    %         1      0       1j*alphaVec(1)*cos(betaVec(1))                  1j*alphaVec(1)*exp(1j*betaVec(1))/2;
    %         0      1       1j*alphaVec(1)*exp(1j*betaVec(1))/2             1j*alphaVec(1)*cos(betaVec(1))];

    
    %for three oscillators, use the following
    FF   = [-1      1       -1j*alphaVec(1) * exp(1j * betaVec(1))/2        1j*alphaVec(1)*exp(-1j*betaVec(1))/2;
            -2      2       -1j*alphaVec(2) * exp(1j * betaVec(2))/2        1j*alphaVec(2)*exp(-1j*betaVec(2))/2;
            -3      3       -1j*alphaVec(3) * exp(1j * betaVec(3))/2        1j*alphaVec(3)*exp(-1j*betaVec(3))/2;
             1      0       1j*alphaVec(1)*cos(betaVec(1))                  1j*alphaVec(1)*exp(1j*betaVec(1))/2;
             2      0       1j*alphaVec(2)*cos(betaVec(2))                  1j*alphaVec(2)*exp(1j*betaVec(2))/2;
             3      0       1j*alphaVec(3)*cos(betaVec(3))                  1j*alphaVec(3)*exp(1j*betaVec(3))/2;
             0      1       1j*alphaVec(1)*exp(1j*betaVec(1))/2             1j*alphaVec(1)*cos(betaVec(1));
             0      2       1j*alphaVec(2)*exp(1j*betaVec(2))/2             1j*alphaVec(2)*cos(betaVec(2));
             0      3       1j*alphaVec(3)*exp(1j*betaVec(3))/2             1j*alphaVec(3)*cos(betaVec(3))];



    % In FF, columns from 1 to d-1 (dimension of phase-difference system) holds the non-zero
    % coefficient indexes,
    % columns form (d) to (2d-2) represents the non-zero coefficients


    
    % we indexed the non-zero coefficienttics in FF.
    % if " i th " coefficient is non-zero, return it in y.
    % return zero for zero coefficients.
    if (isempty(find(ismember(FF(:,1:2), i', 'rows'))) && isempty(find(ismember(FF(:,1:2), -i', 'rows'))))
        y = 0; % zero coefficients,
    elseif ~isempty(find(ismember(FF(:,1:2), i', 'rows')))
        % return non zero coeffs
        y = FF(find(ismember(FF(:,1:2), i', 'rows')),2+l);
    else
        % return conjugates of the coefficients for non-existing indexes.
        y = conj(FF(find(ismember(FF(:,1:2), -i', 'rows')),2+l));

    end
end