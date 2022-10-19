function [ Mat_A,Data_ant] = Generate_Matrix( receivedFields_mag, ...
    receivedFields_pha, greensFunctions_mag, greensFunctions_pha, measMag, ...
  measPha,  numFreqs, numTR, omegas, eps0, PairTR, tauP, bbSize)
% this function is used for generating the matrix of linear equation.
% multiple frequency approach

for nf = 1:numFreqs
    Mat_A = zeros(numTR*numFreqs*2,bbSize*3);
    Data_ant = zeros(numTR*numFreqs*2,1);
    % output of function project_FDTD_2D
    calcMag = squeeze(receivedFields_mag(:,:,nf));
    calcPha = squeeze(receivedFields_pha(:,:,nf));
    % get vector of complex fields on E-mesh inside bounding
    % box at frequency point for this source antenna
    internalFields = squeeze(greensFunctions_mag(:,:,nf).* ...
        exp(1i.*greensFunctions_pha(:,:,nf)));
    %multiplier to convert field to Green's Function
    multiplier = 1i.*omegas(nf).*eps0;
    
    mmMag = zeros(numTR,1);
    mmPha = zeros(numTR,1);
    ccMag = zeros(numTR,1);
    ccPha = zeros(numTR,1);
    % for all source-receiver pairs
    % multiplier*internalFields is the real green function
    A = multiplier*internalFields(PairTR(:,1),:).*internalFields(PairTR(:,2),:);
    for i = 1:numTR
        % store meas and calc antenna observations for TR pair
        mmMag(i) = measMag(PairTR(i,1),PairTR(i,2),nf);
        mmPha(i) = measPha(PairTR(i,1),PairTR(i,2),nf);
        ccMag(i) = calcMag(PairTR(i,1),PairTR(i,2));
        ccPha(i) = calcPha(PairTR(i,1),PairTR(i,2));
    end
    % calculate complex data at freq nf (diff between meas and calc)
    data = ccMag.*exp(1i.*ccPha) - mmMag.*exp(1i.*mmPha); % see (1) in "Application of the DBIM-TwIST Algorithm toExperimental Microwave Imaging Data"
    %============================================================================
    % real part of the multi-frequency matrix
    numer1 = omegas(nf).*tauP;
    %denom1 = 1;
    denom1 = 2*pi*1e9./omegas(nf);
    denom2 = 1+numer1.^2;
    %------------------------------
    Mat_A(1+(nf-1)*numTR*2:(nf-1)*numTR*2+numTR,1:end/3)= real(A);
    Mat_A(1+(nf-1)*numTR*2:(nf-1)*numTR*2+numTR,end/3+1:2*end/3) ...
        = real(A).*(1./denom2)+imag(A).*(numer1./denom2);
    Mat_A(1+(nf-1)*numTR*2:(nf-1)*numTR*2+numTR,2*end/3+1:end) ...
        = imag(A).*denom1;
    % imaginary part of the multi-frequency matrix
    Mat_A(1+(nf-1)*numTR*2+numTR:nf*numTR*2,1:end/3)= imag(A);
    Mat_A(1+(nf-1)*numTR*2+numTR:nf*numTR*2,end/3+1:2*end/3) ...
        = imag(A).*(1./denom2)+real(A).*(-numer1./denom2);
    Mat_A(1+(nf-1)*numTR*2+numTR:nf*numTR*2,2*end/3+1:end) ...
        = -real(A).*denom1;
    % place calculations into multifrequency data
    Data_ant(1+(nf-1)*numTR*2:nf*numTR*2) = [real(data);imag(data)];
end % loop over freqs

end

