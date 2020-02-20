%%
x = importdata('Protocol1Edit.csv');
A = x.data(:,2);
n = 3 %Number of decomposition levels
% Here A is the vector,n is the decompression levels you want, the third is the wave family. Mine one is db10.  
[C1,L1] = wavedec(A,n,'sym2');
%Then wenergy gives the percentage value which corresponds to the 
%decompression level plus the energy values which corresponds to the decomposition levels.  
[Ea1,Ed1] = wenergy(C1,L1)
% Ea1 is the percentage of energy of the approximation, and Ed1 is the
% percentage energy vector of details
%%%%%%%%%%%%Reconstruction
A3 = wrcoef('a',C1,L1,'sym2',n);
%%
figure(5)
[C1,L1] = wavedec(A,n,'sym2'); %Decompose to n levels using the db10 wavelet
[Ea1,Ed1] = wenergy(C1,L1);% Looks at the energy proportions represented by each decomposition level
for i = 1:n
subplot(1,n,i)
sth = wrcoef('a',C1,L1,'sym2',i); %Reconstruction of the i-th level
%%'a' for approximation and 'd' for detail reconstruction
plot(A,'k-','linewidth',3);hold on;plot(sth,'r-','linewidth',3)
xlabel('Sample','fontsize',16)
ylabel('Voltage from Myoware Sensor (V)','fontsize',16)
set(gca,'fontsize',16)
end

figure(6);%Same in 3D
CWTcoeffs = cwt(A,1:20:500,'sym2','plot'); colormap jet; colorbar;
xlabel('Sample (time)','fontsize',16);
ylabel('Scale level (Inverse of frequency)','fontsize',16)

figure(7);%Same in 3D
surf(CWTcoeffs); colormap jet;colorbar;
shading('interp'); view(-30,10);
xlabel('Sample (time)','fontsize',16);
ylabel('Scale level (Inverse of frequency)','fontsize',16)
zlabel('Intensity (coefficients)')