%% mvsa_demo_easy
%
% This demo runs the MVSA algorithm [1] (hard constraint version)  and plots
% results in an easy problem. 
% 
% For a soft constraint version see SISAL algorithm [3].
%
% MVSA: Minimum volume simplex analysis
%
% [3] J. Bioucas-Dias, "A variable splitting augmented Lagrangian approach
%     to linear spectral unmixing", in  First IEEE GRSS Workshop on 
%     Hyperspectral Image and Signal Processing-WHISPERS'2009, Grenoble, 
%     France,  2009. Available at http://arxiv.org/abs/0904.4635v
%
% 
% NOTE:  VCA (Vertex Component Analysis) is used to initialize MVSA. However,
%        VCA is a pure-pixel based algorithm and thus it is not suited to 
%        the data sets herein considered.  Nevertheless, we plot VCA results, 
%        to highlight the advantage of non-pure-pixel based algorithms over the 
%        the pure-pixel based ones.
%
%
% VCA: Vertex component analysis
%
% [4] J. Nascimento and J. Bioucas-Dias, "Vertex componet analysis",
%     IEEE Transactions on Geoscience and Remote Sensing, vol. 43, no. 4,
%     pp. 898-910, 2005.
%
%
% Authors: Jose M. Bioucas-Dias (May, 2014)
%          Jun Li
%          Alexander Agathos
%
%
%% -------------- PROBLEM DESCRIPTION ------------------------------------
%
%   MVSA estimates the vertices  M={m_1,...m_p} of the
%  (p-1)-dimensional   simplex of minimum volume containing the vectors
%  [y_1,...y_N]. 
%
%   MVSA assumes that y belongs to an affine space. It may happen,
%   however, that the data supplied by the user is not in an affine
%   set. For this reason, the first step this code implements
%   is the estimation of the affine set the best represent
%   (in the l2 sense) the data.
%
%   Any vector y_i   belongs  thus to the   convex hull of  the 
%   columns of M; i.e.,
%
%                   y_i = M*x_i
%
%  where x_i belongs to the probability (p-1)-simplex.
%

%%
%--------------------------------------------------------------------------
%       Beginning of the demo
%-------------------------------------------------------------------------

clear all,
clc
close all
verbose= 1;


%%
%--------------------------------------------------------------------------
%        Simulation parameters
%-------------------------------------------------------------------------
% p                         -> number of endmembers
% N                         -> number of pixels
% SNR                       -> signal-to-noise ratio (E ||y||^2/E ||n||^2) in dBs
% SIGNATURES_TYPE           -> see below
% L                         -> number of bands (only valid for SIGNATURES_TYPE = 5,6)
% COND_NUMBER               -> conditioning number of the mixing matrix (only for SIGNATURES_TYPE = 5,6)
% DECAY                     -> singular value decay rate 
% SHAPE_PARAMETER           -> determines the distribution of spectral
%                              over the simplex (see below)
% MAX_PURIRY                -> determines the maximum purity of  the
%                              mixtures.  If  MAX_PURIRY < 1, ther will be
%                              no pure pixels is the data


% SIGNATURES_TYPE;
%  1 - sampled from the USGS Library
%  2 - not available
%  3 - random (uniform)
%  4 - random (Gaussian)
%  5 - diagonal with conditioning number COND_NUMBER and DECAY exponent
%  6 - fully populated with conditioning number COND_NUMBER and DECAY exponent

% NOTE: For SIGNATURES_TYPE 5 or 6, the difficulty of the problem is
% determined by the parameters
% COND_NUMBER              
% DECAY                     


% Souces are Dirichlet distibuted
% SHAPE_PARAMETER ;           
%   = 1   -  uniform over the simplex
%   > 1   -  samples moves towards the center
%            of the simplex, corresponding to higly mixed materials and thus 
%            the unmixing is  more difficult
%   ]0,1[ -  samples moves towards the facets
%            of the simplex. Thus the unmixing is easier.


%
%--------------------------------------------------------------------------
%        SELECTED PARAMETERS FOR AN EASY PROBLEM 
%-------------------------------------------------------------------------

SIGNATURES_TYPE = 3;        % Uniform in [0,1]
p = 10;                      % number of endmembers
N = 2000;                   % number of pixels
SNR = 50;                   % signal-to-noise ratio (E ||y||^2/E ||n||^2) in dBs
L = 200;                    % number of bands (only valid for SIGNATURES_TYPE = 2,3)
% COND_NUMBER  = 1;           % conditioning number (only for SIGNATURES_TYPE = 5,6)
% DECAY = 1;                  % singular values decay rate  (only for SIGNATURES_TYPE = 5,6)
SHAPE_PARAMETER  = 1   ;    % uniform over the simplex
MAX_PURIRY = 0.8;           % there are pure pixels in the data set
OUTLIERS  = 0;              % Number of outliers in the data set


%%
%--------------------------------------------------------------------------
%        Begin the simulation
%-------------------------------------------------------------------------
switch SIGNATURES_TYPE
    case 1
        %rand('seed',5);
        load('USGS_1995_Library')
        wavlen=datalib(:,1);    % Wavelengths in microns
        [L n_materiais]=size(datalib);
        % select randomly
        sel_mat = 4+randperm(n_materiais-4);
        sel_mat = sel_mat(1:p);
        M = datalib(:,sel_mat);
        % print selected endmembers
%         fprintf('endmembers:\n')
%         for i=1:p
%             aux = names(sel_mat(i),:);
%             fprintf('%c',aux);
%             st(i,:) = aux;
%         end
        clear datalib wavelen names aux st;
    case 2
        error('type not available')
    case 3
        M = rand(L,p);
    case 4
        M = randn(L,p);
    case 5
        L=p;
        M = diag(linspace(1,(1/(COND_NUMBER)^(1/DECAY)),p).^DECAY);
    case 6
        L=p;
        M = diag(linspace(1,(1/(COND_NUMBER)^(1/DECAY)),p).^DECAY);
        A = randn(p);
        [U,D,V] = svd(A);
        M = U*M*V';
        clear A U D V;
    otherwise
        error('wrong signatute type')
end


%--------------------------------------------------------------------------
%        Set noise parameters (to be used in spectMixGen function)
%-------------------------------------------------------------------------
% white_noise = [0 1 1];     % white noise
% % non-white noise parameters
% eta   = 10;                % spread of the noise shape
% level = 10;                % floor lever
% gauss_noise = [level L/2 eta]; % Gaussian shaped noise centered at L/2 with spread eta
% % and floor given by level
% rect_noise  = [level L/2 eta]; % Rectangular shaped noise centered at L/2 with spread eta
% % and floor given by level

%%
%--------------------------------------------------------------------------
%        Generate the data set
%-------------------------------------------------------------------------
%
%   Sources are Diriclet distributed (shape is controled by 'Source_pdf' and
%   'pdf_pars': 'Source_pdf' = 'Diri_id' and 'pdf_pars' = 1 means a uniform
%   density over the simplex).  The user may change the parameter to
%   generate other shapes. Mixtures are aldo possible.
%
%   'max_purity' < 1 means that there are no pure pixels
%
[Y,x,noise] = spectMixGen(M,N,'Source_pdf', 'Diri_id','pdf_pars',SHAPE_PARAMETER,...
    'max_purity',MAX_PURIRY*ones(1,p),'no_outliers',OUTLIERS, ...
    'violation_extremes',[1,1.2],'snr', SNR, ...
    'noise_shape','uniform');

%%
%--------------------------------------------------------------------------
%        Remove noise  (optional)
%-------------------------------------------------------------------------
%   noise_hat = estNoise(Y);
%   Y = Y-noise_hat;
%   clear noise_hat




%%
%--------------------------------------------------------------------------
%       Project  on the  affine set defined by the data in the sense L2
%-------------------------------------------------------------------------
%
%   The application of this projection ensures that the data is in
%   an affine set.
%
%   Up is an isometric matrix that spans the subspace where Y lives
[Y,Up,my,sing_val] = dataProj(Y,p,'proj_type','affine');



%%
%--------------------------------------------------------------------------
%        Degree of Difficulty of the problem
%-------------------------------------------------------------------------
%% compute original subspace
sing_vects = svds(M,p);

% Condition number gives an idea of the difficulty in inferring
% the subspace
fprintf('Conditioning number of M = %2f \n', sing_vects(1)/sing_vects(end))
% fprintf('\n Hit any key: \n ');
% pause;

Cx = Up'*(M*x)*(M*x)'*Up/N;
Cn = Up'*noise*noise'*Up/N;
[U,D]=svd(Cx);

% compute the SNR along the direction corresponding the smaller eigenvalue 

LOWER_SNR= D(p,p)/(U(:,p)'*Cn*U(:,p));
fprintf('\nSNR along the signal smaller eigenvalue = %f \n', LOWER_SNR);
if LOWER_SNR < 20
   fprintf('\nWARNING: This problem is too hard and the results may be inaccurate \n') 
end

clear noise x;
    
    



%%
%--------------------------------------------------------------------------
%         ALGORITHMS
%-------------------------------------------------------------------------


%%
%--------------------------------------------------------------------------
%         MVSA - Minimum volume simple analysis
%-------------------------------------------------------------------------

%start timer
tic
A_est = mvsa(Y,p,'spherize','yes');
Mvsa = Up'*A_est;
%stop timer
t(1) = toc;


%%
%--------------------------------------------------------------------------
%         VCA  - Vertex component analysis
%-------------------------------------------------------------------------
%
% start timer
tic
[Ae, indice,ys]= VCA(Y,'Endmembers',p);
Mvca = Up'*Ae;
% stop timer
t(2) = toc;


%%
%--------------------------------------------------------------------------
%         Project the original mixing matrix and the data set on the
%         identified affine set.
%-------------------------------------------------------------------------
Mtrue = Up'*M;
Y=Up'*Y;


%%
%--------------------------------------------------------------------------
%        Display the results
%-------------------------------------------------------------------------

% selects axes  to display
%
%

I = 1;
J = 2;
K = 3;

% canonical orthogonal directions
E_I = eye(p);

v1 = E_I(:,I);
v2 = E_I(:,J);
v3 = E_I(:,K);

% original axes

Q = inv(Mtrue);
% v1 = Q(I,:)';
% v2 = Q(J,:)';
% v3 = Q(K,:)';



Y = [v1 v2 v3]'*Y;
m_true = [v1 v2 v3]'*Mtrue;



% legend
leg_cell = cell(1);
leg_cell{end} = 'data points';
H_1=figure(1);
plot(Y(1,:),Y(2,:),'k.','Color',[ 0 0 1])

hold on;
plot(m_true(1,[1:p 1]), m_true(2,[1:p 1]),'*', 'Color',[0 0 0])
leg_cell{end +1} = 'true';



m_vsa  = [v1 v2 v3]'*Mvsa;
plot(m_vsa(1,[1:p 1]), m_vsa(2,[1:p 1]),'O','Color',[1 0 0])
leg_cell{end +1} = 'MVSA';

m_vca  = [v1 v2 v3]'*Mvca;
plot(m_vca(1,[1:p 1]), m_vca(2,[1:p 1]),'p', 'Color',[0  1 0])
leg_cell{end +1} = 'VCA';



xlabel('v1''*Y'),ylabel('v2''*Y');
legend(leg_cell)
title('Endmembers and data points (2D projection)')


%%
%--------------------------------------------------------------------------
%       3D plot, if p>= 4
%-------------------------------------------------------------------------

% if p>= 4
%     % legend
%     leg_cell = cell(1);
%     leg_cell{end} = 'data points';
%     figure;
%     plot3(Y(1,:),Y(2,:),Y(3,:),'k.','Color',[ 0 0 1])
%
%     hold on;
%     plot3(m_true(1,[1:p 1]), m_true(2,[1:p 1]),m_true(3,[1:p 1]),'*', 'Color',[0 0 0])
%     leg_cell{end +1} = 'true';
%
%
%     plot3(m_vsa(1,[1:p 1]), m_vsa(2,[1:p 1]),m_vsa(3,[1:p 1]),'O','Color',[0 1 0])
%     leg_cell{end +1} = 'MVSA';
%
%        
%     plot3(m_vca(1,[1:p 1]), m_vca(2,[1:p 1]),m_vca(3,[1:p 1]),'p', 'Color',[0  0 0])
%     leg_cell{end +1} = 'VCA';

%
%
%     xlabel('v1''*Y'),ylabel('v2''*Y'),zlabel('v3''*Y');
%     legend(leg_cell)
%     title('Mixing matrices ')
%
% end

fprintf('\nTIMES (sec):\n    MVSA = %3.2f, \n    VCA = %3.2f\n', t(1), t(2))

%--------------------------------------------------------------------------
%        Display errors
%-------------------------------------------------------------------------

%% compute correspondences

% mvca
angles = Mtrue'*Mvca./(repmat(sqrt(sum(Mtrue.^2)),p,1)'.*(repmat(sqrt(sum(Mvca.^2)),p,1)));
P = zeros(p);
for i=1:p
    [dummy,j] = max(angles(i,:));
    P(j,i) = 1;
    angles(:,j) = -inf;
end
% permute colums
Mvca = Mvca*P;

% mvsa
angles = Mtrue'*Mvsa./(repmat(sqrt(sum(Mtrue.^2)),p,1)'.*(repmat(sqrt(sum(Mvsa.^2)),p,1)));
P = zeros(p);
for i=1:p
    [dummy,j] = max(angles(i,:));
    P(j,i) = 1;
    angles(:,j) = -inf;
end
% permute colums
Mvsa = Mvsa*P;

MVSA_ERR =norm(Mtrue-Mvsa,'fro')/norm(Mtrue,'fro');
VCA_ERR =norm(Mtrue-Mvca,'fro')/norm(Mtrue,'fro');

fprintf('\nERROR(mse relative):\n   MVSA = %f  \n   VCA  = %f\n', MVSA_ERR, VCA_ERR);


%--------------------------------------------------------------------------
%        Plot signatures
%-------------------------------------------------------------------------
% Choose signatures

leg_cell = cell(1);
H_2=figure(2);
hold on
clear p_H;

% plot signatures
p_H(1) = plot(1:L,(Up*Mtrue(:,1))','b');
leg_cell{end} = 'Mtrue';


p_H(end+1)= plot(1:L,(Up*Mvsa(:,1))','r');
leg_cell{end+1} = 'Mvsa';




for i=2:p
    plot(1:L,(Up*Mtrue(:,i))','b');
    plot(1:L,(Up*Mvsa(:,i))','r');
end

legend(leg_cell)
title('First endmember')
xlabel('spectral band')

pos1 = get(H_1,'Position');
pos1(1)=pos1(1)-300;
set(H_1,'Position', pos1)

pos2 = get(H_2,'Position');
pos2(1)=pos2(1)+300;
set(H_2,'Position', pos2)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %