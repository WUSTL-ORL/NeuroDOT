function Properties=OptProp(lambdas)

% the point of this function is to interpolate and extrapolate absorption
% and scattering coefficients at various wavelengths based on published data 
%
% all values are in inverse mm.

% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht
% Other Contributor(s): Zachary E. Markow
% Eggebrecht et al., 2014, Nature Photonics; Zeff et al., 2007, PNAS.
%
% Washington University hereby grants to you a non-transferable, 
% non-exclusive, royalty-free, non-commercial, research license to use 
% and copy the computer code that is provided here (the Software).  
% You agree to include this license and the above copyright notice in 
% all copies of the Software.  The Software may not be distributed, 
% shared, or transferred to any third party.  This license does not 
% grant any rights or licenses to any other patents, copyrights, or 
% other forms of intellectual property owned or controlled by Washington 
% University.
% 
% YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS 
% PROVIDED AS IS, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR 
% IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY 
% OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY 
% THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  
% IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON 
% UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR 
% CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH 
% THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER 
% IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

if ~exist('lambdas','var'),lambdas=[750,850];end

Range_mu=600:1:1000;


%% Published optical data

% Okada et al, 1997
Mua_scalp_skull_Okada=0.04;
Mua_csf_Okada=0.001;
Mua_gm_Okada=0.025;
Mua_wm_Okada=0.005;

Mus_scalp_skull_Okada=2.0;
Mus_csf_Okada=0.01;
Mus_gm_Okada=2.5;
Mus_wm_Okada=6.0;

% Bevilacqua et al, 1999 via FDPM (freq domain photon mig)
Wavelengths_Bevilacqua_1=[674 811 849 956];
Mua_Bev_gm_F=[0.0173,0.0182,0.0185,0.0206];
Mua_Bev_gm_T=[0.0179,0.0190,0.0179,0.0218];

Mus_Bev_gm_F=[1.12,0.74,0.74,0.8];
Mus_Bev_gm_T=[0.99,0.48,0.45,0.42];

Wavelengths_Bevilacqua_2=[674 849 956];
Mua_Bev_skull=[0.0208,0.0215,0.0355];
Mua_Bev_wm=[0.0165,0.0132,0.0299];

Mus_Bev_skull=[1.19,0.91,0.77];
Mus_Bev_wm=[1.34,0.98,0.84];

% Primary data: Strangman et al, 2002 ('Factors affecting accuracy of NIRS
% concentration calculations for focal changes in oxygenation parameters').
Wavelengths=[690 760 780 830]; % in nm
Mua_Scalp=[0.0159, 0.0177, 0.0164, 0.0191];
Mua_Skull=[0.0101, 0.0125, 0.0115, 0.0136];
Mua_Brain=[0.0178, 0.0195, 0.0170, 0.0186];

Mus_Scalp=[0.8, 0.73, 0.71, 0.66];
Mus_Skull=[1.0, 0.93, 0.91, 0.86];
Mus_Brain=[1.25, 1.18, 1.16, 1.11];


% PLOT raw data sets for comparison over 600:900 nm (in 10 nm steps)
figure;
subplot(2,2,1)
plot(Range_mu,repmat(Mua_scalp_skull_Okada,1,length(Range_mu)),'-m',...
    Range_mu,repmat(Mua_gm_Okada,1,length(Range_mu)),'-b',...
    Range_mu,repmat(Mua_wm_Okada,1,length(Range_mu)),'--b',...
    Wavelengths_Bevilacqua_1,Mua_Bev_gm_F,'-b+',...
    Wavelengths_Bevilacqua_1,Mua_Bev_gm_T,'-.b+',...
    Wavelengths_Bevilacqua_2,Mua_Bev_skull,'-g+',...
    Wavelengths_Bevilacqua_2,Mua_Bev_wm,':b+',...
    Wavelengths,Mua_Scalp,'-ro',...
    Wavelengths,Mua_Skull,'-go',...
    Wavelengths,Mua_Brain,'-bo');
title('Absorption coefficients of scalp (red) skull (green) and brain (blue).')
xlabel('Wavelength [nm]')
ylabel('\mu_a [mm^-^1]')
ylim([0 0.05]);

subplot(2,2,2)
plot(Range_mu,repmat(Mus_scalp_skull_Okada,1,length(Range_mu)),'-m',...
    Range_mu,repmat(Mus_gm_Okada,1,length(Range_mu)),'-b',...
    Range_mu,repmat(Mus_wm_Okada,1,length(Range_mu)),'--b',...
    Wavelengths_Bevilacqua_1,Mus_Bev_gm_F,'-b+',...
    Wavelengths_Bevilacqua_1,Mus_Bev_gm_T,'-.b+',...
    Wavelengths_Bevilacqua_2,Mus_Bev_skull,'-g+',...
    Wavelengths_Bevilacqua_2,Mus_Bev_wm,':b+',...
    Wavelengths,Mus_Scalp,'-ro',...
    Wavelengths,Mus_Skull,'-go',...
    Wavelengths,Mus_Brain,'-bo');
title('Reduced scattering coefficients of scalp (red) skull (green) and brain (blue).')
xlabel('Wavelength [nm]')
ylabel('\mu_s_'' [mm^-^1]')
legend('Okada scalp-skull','Okada gm','Okada wm','Bevilacqua gm frontal',...
    'Bevilacqua gm temporal','Bevilaqua skull','Bevilacqua wm',...
    'Strangman scalp','Strangman skull','Strangman brain',...
    'Location','NorthEast');
ylim([0 6.05]);

%% Fit the primary data with linear functions

% Absorption
p_mua_scalp=polyfit(Wavelengths,Mua_Scalp,1);
p_mua_skull=polyfit(Wavelengths,Mua_Skull,1);
p_mua_brain=polyfit(Wavelengths,Mua_Brain,1);
p_mua_gm=polyfit(Wavelengths_Bevilacqua_1,(Mua_Bev_gm_F+Mua_Bev_gm_T)./2,1);
p_mua_wm=polyfit(Wavelengths_Bevilacqua_2,Mua_Bev_wm,1);

Mua_Scalp_fit=polyval(p_mua_scalp,Range_mu);
Mua_Skull_fit=polyval(p_mua_skull,Range_mu);
Mua_Brain_fit=polyval(p_mua_brain,Range_mu);
Mua_GM_fit=polyval(p_mua_gm,Range_mu);
Mua_WM_fit=polyval(p_mua_wm,Range_mu);

subplot(2,2,3)
plot(Wavelengths,Mua_Scalp,'ro',Range_mu,Mua_Scalp_fit,'--r',...
    Wavelengths,Mua_Skull,'go',Range_mu,Mua_Skull_fit,'--g',...
    Wavelengths,Mua_Brain,'bo',Range_mu,Mua_Brain_fit,'--b',...
    Wavelengths_Bevilacqua_1,(Mua_Bev_gm_F+Mua_Bev_gm_T)./2,'ks',Range_mu,Mua_GM_fit,'--k',...
    Wavelengths_Bevilacqua_2,Mua_Bev_wm,'ms',Range_mu,Mua_WM_fit,':m');
hold on;
title([{'Linear fits to Strangman et al and Bevilacqua et al'};...
    {'absorption coefficients.'}])
xlabel('Wavelength [nm]')
ylabel('\mu_a [mm^-^1]')
ylim([0 0.05]);
xlim([600 1000])

% Scattering
p_mus_scalp=polyfit(Wavelengths,Mus_Scalp,1);
p_mus_skull=polyfit(Wavelengths,Mus_Skull,1);
p_mus_brain=polyfit(Wavelengths,Mus_Brain,1);
p_mus_gm=polyfit(Wavelengths_Bevilacqua_1,(Mus_Bev_gm_F+Mus_Bev_gm_T)./2,1);
p_mus_wm=polyfit(Wavelengths_Bevilacqua_2,Mus_Bev_wm,1);

Mus_Scalp_fit=polyval(p_mus_scalp,Range_mu);
Mus_Skull_fit=polyval(p_mus_skull,Range_mu);
Mus_Brain_fit=polyval(p_mus_brain,Range_mu);
Mus_GM_fit=polyval(p_mus_gm,Range_mu);
Mus_WM_fit=polyval(p_mus_wm,Range_mu);

subplot(2,2,4)
plot(Wavelengths,Mus_Scalp,'ro',Range_mu,Mus_Scalp_fit,'--r',...
    Wavelengths,Mus_Skull,'go',Range_mu,Mus_Skull_fit,'--g',...
    Wavelengths,Mus_Brain,'bo',Range_mu,Mus_Brain_fit,'--b',...
    Wavelengths_Bevilacqua_1,(Mus_Bev_gm_F+Mus_Bev_gm_T)./2,'ks',Range_mu,Mus_GM_fit,'--k',...
    Wavelengths_Bevilacqua_2,Mus_Bev_wm,'ms',Range_mu,Mus_WM_fit,':m');
hold on;
title([{'Linear fits to Strangman et al and Bevilacqua et al'};...
    {'reduced scattering coefficients.'}])
xlabel('Wavelength [nm]')
ylabel('\mu_s_'' [mm^-^1]')
legend('Strangman scalp data','Strangman scalp fit','Strangman skull data',...
    'Strangman skull fit','Strangman brain data','Strangman brain fit',...
    'Bevilacqua gm data','Bevilacqua gm fit','Bevilacqua wm data','Bevilaqua wm fit',...
    'Location','NorthEast');    
ylim([0 6.05]);
xlim([600 1000])

% Populate relavant fitted points
Properties.Scalp_mua(1,:)=lambdas;
Properties.Skull_mua(1,:)=lambdas;
Properties.Brain_mua(1,:)=lambdas;
Properties.GM_mua(1,:)=lambdas;
Properties.WM_mua(1,:)=lambdas;
Properties.Scalp_mus(1,:)=lambdas;
Properties.Skull_mus(1,:)=lambdas;
Properties.Brain_mus(1,:)=lambdas;
Properties.GM_mus(1,:)=lambdas;
Properties.WM_mus(1,:)=lambdas;


for i=1:length(lambdas)
    Properties.Scalp_mua(2,i)=Mua_Scalp_fit(Range_mu==lambdas(i));
    Properties.Skull_mua(2,i)=Mua_Skull_fit(Range_mu==lambdas(i));
    Properties.Brain_mua(2,i)=Mua_Brain_fit(Range_mu==lambdas(i));
    Properties.GM_mua(2,i)=Mua_GM_fit(Range_mu==lambdas(i));
    Properties.WM_mua(2,i)=Mua_WM_fit(Range_mu==lambdas(i));
    
    Properties.Scalp_mus(2,i)=Mus_Scalp_fit(Range_mu==lambdas(i));
    Properties.Skull_mus(2,i)=Mus_Skull_fit(Range_mu==lambdas(i));
    Properties.Brain_mus(2,i)=Mus_Brain_fit(Range_mu==lambdas(i));
    Properties.GM_mus(2,i)=Mus_GM_fit(Range_mu==lambdas(i));
    Properties.WM_mus(2,i)=Mus_WM_fit(Range_mu==lambdas(i));    
end

subplot(2,2,3);
plot(Properties.Scalp_mua(1,:),Properties.Scalp_mua(2,:),'kx','MarkerSize',10)
plot(Properties.Skull_mua(1,:),Properties.Skull_mua(2,:),'kx','MarkerSize',10)
plot(Properties.Brain_mua(1,:),Properties.Brain_mua(2,:),'kx','MarkerSize',10)
plot(Properties.GM_mua(1,:),Properties.GM_mua(2,:),'kx','MarkerSize',10)
plot(Properties.WM_mua(1,:),Properties.WM_mua(2,:),'kx','MarkerSize',10)
ylim([0 0.05]);
xlim([600 1000])

subplot(2,2,4);
plot(Properties.Scalp_mus(1,:),Properties.Scalp_mus(2,:),'kx','MarkerSize',10)
plot(Properties.Skull_mus(1,:),Properties.Skull_mus(2,:),'kx','MarkerSize',10)
plot(Properties.Brain_mus(1,:),Properties.Brain_mus(2,:),'kx','MarkerSize',10)
plot(Properties.GM_mus(1,:),Properties.GM_mus(2,:),'kx','MarkerSize',10)
plot(Properties.WM_mus(1,:),Properties.WM_mus(2,:),'kx','MarkerSize',10)
ylim([0 6.05]);
xlim([600 1000])
