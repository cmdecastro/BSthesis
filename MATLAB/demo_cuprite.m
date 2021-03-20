%  demo_cuprite
%
%  This demo illustrates the MVSA unmixing a subscene of the cuprite 
%  data set.
%
%  ------------------------------------------------------------------------
%   Observation model: linear mixing model
%     
%     Y = MS + N    with   S>= 0,   and   sum(S) = ones(1,np)
%
% -----------------------------------------------------------------------
%  Optimization problem (MVSA criterion): 
%
%   1-Project Y onto a p-dimensional subspace containing the data set y
%
%            Yp = Up'*Y;      Up is an isometric matrix (Up'*Up=Ip)
%
%   2- solve the   optimization problem
%
%     Q^* = arg min_Q  -\log abs(det(Q))
%
%      subject to: Q*Yp >= 0 and ones(1,p)*Q=mq,
%
%     where mq = ones(1,N)*Yp'inv(Yp*Yp)
%
%   3- Compute
%
%      M = Up*inv(Q^*);
%
%
% Author: Jun Li, Jose Bioucas Dias, January 2015
%

clear all
close all


%--------------------------------------------------------------------------
% Load cuprite data set 
%--------------------------------------------------------------------------
load cuprite_ref.mat
Y = x;
%  Remove the first 5 bands (too noisy)
band_initial = 5;
Y=Y(band_initial:end-1,:);
BANDS = BANDS(band_initial:end-1);
clear x;

% 
% Notes about the data set:
%
%  -  the spectral reflectances are in vector Y
%
%  -  the original data set constains 224 bands. 
%     36+(band_initial-1) were removed due to low SNR. The retained bands 
%     are saved in the vector BANDS
%
%  -  The vector wavelen  contains the map between bands and wavelenghts; 
%     this is the same we have for the library USGS
%  

[L,np] = size(Y);   

%--------------------------------------------------------------------------
% Estimate the signal subspace with Hysime
%--------------------------------------------------------------------------
[w Rw] = estNoise(Y);
[kf,Ek,E,delta_p]=hysime(Y,w,Rw);
clear w;
% project in the subspace
p = kf;
Ek = Ek(:,1:p);
Yp = Ek*Ek'*Y;

% dpft (perpective)  projection. Works better than 'otrh' in tthis data set
[Yp,Up,my,angles,scales] = affineProj(Yp,p,'proj_type','dpft');

% note: Up sapans the same subspace as Ek



%% scattergram centered at the mass center
figure(10)
ym = mean(Yp,2);
Yc = Yp-repmat(ym,1,np);
plot(Yc(1,:), Yc(2,:),'.');
title('Projection of Y onto the firt two PCA components')


%--------------------------------------------------------------------------
% Estimate the mixing matrix with VCA
%--------------------------------------------------------------------------
vol_ant = -inf;
VCA_runs = 30;
small = 1e-4;
for i=1:VCA_runs
    Aux = VCA(Yp,'Endmembers',p,'SNR',1,'verbose','yes');
    vol = sum(log(svd(Aux) + small));
    if vol > vol_ant
        Avca = Aux;
        vol_ant = vol;
    end
end

% estimate abundances 
% This is a modified version of sunsal whic allows the abundances to be
% negative.  IN THIS WAY, AND IF THERE ARE SAMPLES OUSIDE THE SIMPLEX
% COMPUTED BY THE PURE PIXEL ALGORITHM, WE MAY PROJECT THEN  ON AN 
% INFLATED SIMPLEX AND THEN RUN MVSA, WHITCH WILL WORK WITHOUT PROBLEMS
% BECAUSE THE FACETS ARE REGULARIZED


% PROJECT ON THE INFLATED SIMPEX  (ABUNDANCES LARGER THAN -0.01) 
slack = -0.01;
Xvca = sunsal_mod(Avca,Yp,'POSITIVITY',slack,'VERBOSE','yes','ADDONE','yes', ...
    'lambda', 0,'AL_ITERS',1000, 'TOL', 1e-8);
% 


%--------------------------------------------------------------------------
% Estimate the mixing matrix with MVSA using the inflacted  simplex
%--------------------------------------------------------------------------
%%
%Ahat = mvsa(Ypp,p,'spherize','no', 'M0',Avca ,'verbose',1);
Amvsa = mvsa(Avca*Xvca,p,'spherize','yes' ,'verbose',1);

figure(30);
plot(Up*Amvsa)
title('Estimated  spectral signatures (MVSA)')

figure(35);
plot(Up*Avca)
title('Estimated spectral signatures (VCA)')

figure(40)
plot(Yp(1,:), Yp(2,:),'.')
hold on
plot(Avca(1,:),Avca(2,:),'ro',Amvsa(1,:),Amvsa(2,:),'gs','LineWidth',3)
hold off
legend('samples','VCA', 'MVSA')

%  abundances wrt to the Amvsa mixing matrix
Xmvsa = sunsal_mod(Amvsa,Yp,'POSITIVITY',0,'VERBOSE','yes','ADDONE','yes', ...
    'lambda', 0,'AL_ITERS',1000, 'TOL', 1e-8);
% 




figure(45)
imagesc(Xmvsa)
title('Estimated abundances (MVSA)')
max(Xmvsa')


%%
%--------------------------------------------------------------------------
% Compare with library signatures
%--------------------------------------------------------------------------
%load USGS_1995_Library
%aux_lib = datalib(BANDS,4:end);

load USGU_sublibrary.mat
aux_lib = Areal(band_initial:end-1,:);

       % project onto the subspace span{Upp}
        aux_lib = Up'*aux_lib;    
        % perpective projection  (to be comparable with dataset) 
        u = ones(L,1)/(ones(L,1)'*my);
        aux_lib = aux_lib./repmat((Up'*u)'*aux_lib, p,1);
        %lift
        aux_lib = Up*aux_lib;
    
    [~,n_sig] = size(aux_lib);
    %names = names(4:end,:);
    
    %used to avoid choosing the same spactral signatute twice
    indicator_set = zeros(1,n_sig);

    big = 1e3;

    for i=1:p
        end_sig = Up*Amvsa(:,i);
        [~,index] = min( sum((end_sig*ones(1,n_sig)-aux_lib).^2)+ indicator_set);
        indicator_set(i) = 0;
       
        indexes_lib(i) = index; 
        st(i,:)=sprintf('%s',char(names(index,:)));        
     
    end
    

%--------------------------------------------------------------------------
% plot abundances
%--------------------------------------------------------------------------

for i=1:p
   figure(100+i);
   imagesc(reshape(Xmvsa(i,:),Lines,Columns))
   title(sprintf('abundance map %d',i))
   title(st(i,:))
       set(gcf,'Position',[100 100 300 260])
    set(gca,'units','pixels','Visible','off'); 
    q=get(gca,'position'); 
    q(1)=0; q(2)=0;  
    set(gcf,'position',q-5)
    set(gca,'position',q-5)
    
        colorbar
        colormap(gray)
end


%unscale abundances 
load USGU_sublibrary.mat
aux_lib = Areal(band_initial:end-1,:);
sig_lib = aux_lib(:,indexes_lib);
figure(300)
plot(sig_lib)

% compute innner products/|| Amvsa ||^2
scales = sum((Up*Amvsa).*sig_lib)./sum(Amvsa.^2);
Amvsa_sca = Amvsa.*repmat(scales,p,1); 

load AVIRIS_wavelenght

 x1 = nan*ones(224,1); 
Bands = [7:104 116:149 171:221];


nband = 1:224;

for i=1:p
    
    xr=x1;        
    xmvsa = x1;
    

          if i<= 5
            figure(50)
            subplot(strcat('5,1,',num2str(i)));
            
        elseif i>5 & i<=10
            figure(60)
            subplot(strcat('5,1,',num2str(i-5)));
        else
                        figure(70)
            subplot(strcat('5,1,',num2str(i-10)));
          end
        
              
    xmvsa(Bands) =Up*Amvsa_sca(:,i);    xmvsa([2 3 104 116 149]) = xmvsa(1);
  
    xr(Bands) =sig_lib(:,i); 

    
        plot(wavebands,xr,'b-.','Linewidth',1.5)
        hold on
        plot(wavebands,xmvsa,'r','Linewidth',1.5)
        hold off
        if i ==3
ylabel('Reflectance (%)','Fontsize',14)
        end
        
        if i ==5
xlabel('Bands','Fontsize',14)
        end
        
        if i ==8
ylabel('Reflectance (%)','Fontsize',14)
        end
        
        if i ==10
xlabel('Bands','Fontsize',14)
        end
        if i ==13
ylabel('Reflectance (%)','Fontsize',14)
        end
        
             if i ==10
xlabel('Wavelength(\mum)','Fontsize',14)
        end
        if i ==13
ylabel('Reflectance (%)','Fontsize',14)
        end
                
                if i ==p
  xlabel('Wavelength(\mum)','Fontsize',14)
                end
                  if i ==1
legend('Real','MVSA','Fontsize',14)
                  end
          
        axis([0.4 2.5 0 1])
        st(i,:)=sprintf('%s',char(names(index,:)));
        st(i,:)
        
end


norm(Up*Amvsa_sca-sig_lib, 'fro')




