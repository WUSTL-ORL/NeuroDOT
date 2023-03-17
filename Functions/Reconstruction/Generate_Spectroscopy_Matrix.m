function Eout=Generate_Spectroscopy_Matrix(lambdas,spectraIN,ExtCoeffs)
%
% Make new spectroscopy matrix
%
% Inputs:
% lambdas   the list of system wavelengths for which we will make the
%           matrix. Assumed to be in nm.
% spectra   power distribution for spectra of source for a given system, 
%           if available. if not spectrum is available,then 1 is passed
%           in to simply read off the extiction coeff at that wavelength.
%           To allow for differing lengths, this is a cell array with 2
%           columns each: wavelengths and spectra.
% ExtCoeffs     Input structure containing:
%               lambda_nm - wavelenths
%               HbO2_loge - HbO2 log-e extinction coeffs
%               HbR_loge - HbR log-e extinction coeffs
%           *** An example can be loaded 'Extinction_Coefficients.m' ***
%           If so, ExtCoeffs.Prahl or ExtCoeffs.Wray gets passed in.

% Outputs:
% Eout      Extinction coefficient matrix (wavelengths X chromophores)
%
% EG Usage:
% load('Extinction_Coefficients.m')
% lambdas=[750,850]; 
% spectra{1}=1;spectra{2}=1;
% Eout=Generate_Spectroscopy_Matrix([750,850],spectra,ExtCoeffs.Prahl)


%% Parameters and Initialization
numled=numel(lambdas);
Eout=zeros(numled,2);
lamda_min=min(ExtCoeffs.lambda_nm);
lamda_max=max(ExtCoeffs.lambda_nm);


%% Generate extinction coefficient matrix
for n=1:numled      
    
    % Check for spectra content
    if spectraIN{n}==1
        spectra{n}(:,1)=[floor(lamda_min):ceil(lamda_max)];
        spectra{n}(:,2)=zeros(size(spectra{n}(:,1)));
        spectra{n}(spectra{n}(:,1)==round(lambdas(n)),2)=1;
    else
        spectra{n}=spectraIN{n};
    end
    
    % Interpolate from Spectrometer Wavelengths to Reference Wavelengths
    ledpower=interp1(spectra{n}(:,1),spectra{n}(:,2),...
        ExtCoeffs.lambda_nm,'pchip');
    
    if ~any(ledpower)
        [~,idx]=min(abs(ExtCoeffs.lambda_nm-lambdas(n)));
        ledpower(idx)=1;
    end
    
    % Normalize
    ledpower=ledpower./max(ledpower);
    
    % Zero Out Noise
    ledpower(ledpower<0.01)=0;

    % Normalize
    ledpower=ledpower./sum(ledpower);
    
    % Spectroscopy Matrix
    Eout(n,1)=sum(ExtCoeffs.HbO2_loge.*ledpower);
    Eout(n,2)=sum(ExtCoeffs.HbR_loge.*ledpower);

end
