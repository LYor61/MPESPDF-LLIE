function C_out = MPESPD(imgSeqColor,r1,noise,Cl,n,PGway,Cfun)
[D1,i_mean1,aa1,N1,Im,Var] = scale_fine(imgSeqColor,r1,Cl,PGway,Cfun);
if noise==2 % BM3D
y = D1;
noise_type =  'gw';
noise_var = 0.008; % Noise variance
[~, PSD, ~] = getExperimentNoise(noise_type, noise_var, 0, size(y));
y_est = CBM3D(y, PSD);
D1 = y_est;
elseif noise==1 %Weight Guided Filter
G = D1(:,:,1)*0.299+0.587*D1(:,:,2)+0.114*D1(:,:,3);
Smean=mean2(imgSeqColor(:));
RR=abs(log(Smean));
% n=my_contrast(D1);
C=0.01+Var;
for ii=1:3
    D1(:, :, ii) = Weight_Guided_Filter( G,D1(:, :, ii),ceil(RR), 10e-4,C); 
end 
end
%% the intermediate  scale
[w,h,~,~]=size(imgSeqColor);
nlev = floor(log(min(w,h)) / log(2))-n;
D2 = cell(nlev,1);
aa2= cell(nlev,1);
N2= cell(nlev,1);
Im2 = cell(nlev,1);
for ii=1:nlev
    r2=4;
    [ D2{ii},i_mean2,aa2{ii},N2{ii},Im2{ii}] = scale_interm(i_mean1,r2,ii,Cl,PGway,Cfun);
    i_mean1=i_mean2;
end
%% the coarsest  scale
r3=4;
[fI3,~,~,~,~] = scale_coarse(i_mean2,r3,Cl,PGway,Cfun);
%% reconstruct
%% Intermediate layers
for ii=nlev:-1:1
    temp=aa2{ii};
    fI=zeros(size(temp));
    fI(1:2:size(temp,1),1:2:size(temp,2))=fI3;
    B2=boxfilter(fI.*Im2{ii}, r2)./ N2{ii}+D2{ii};
    fI3=B2;
end
%% finest layers
fI=zeros(size(aa1));
fI(1:2:size(aa1,1),1:2:size(aa1,2))=B2;
B1=boxfilter(fI.*Im, r1)./ N1;
a=0.75;b=1.2;
C_out=D1*a+repmat(B1*b,[1 1 3]);
end

function [fI3, i_meant,aa,N1,Im,Var] = scale_fine(imgSeqColor,r,Cl,PGway,Cfun)
[h,w,c]=size(imgSeqColor);
N = boxfilter(ones(h, w), r);
tem=ones(h, w);
tem(:,2:2:w)=0;
tem(2:2:h,:)=0;
N1= boxfilter(tem, r);
p=4;
[Im,imgSeqColor] = Enhancement(imgSeqColor,r,1,Cl,PGway);
[WD, Cmax,i_mean2,Var]= weight_cal_detail(imgSeqColor,p,r,Cfun); %w -> β,β=c^(p-1)/(c^p+eps)
WD=WD.*Cmax; % β*c' if c'=c1, WD=1
WD = WD.*Im; % 
% WD = WD; %
F_temp2_detail=zeros(h,w,c);
i_meant=zeros(ceil(h/2),ceil(w/2));
aa=i_mean2(:,:).*tem;
i_meant(:,:)=aa(1:2:h,1:2:w);
W_D1=boxfilter(i_mean2(:,:).*WD(:,:), r)./ N;
% W_D1=i_mean2(:,:).*WD(:,:);
W_D2=boxfilter(WD(:,:), r)./ N;
F_temp2_detail(:,:,:)=repmat(W_D2,[1 1 3]).*imgSeqColor(:,:,:)-repmat(W_D1,[1 1 3]); %WD*(X-U)   
% S=F_temp2_detail./Cmax; % s*Im
%     figure,imshow(S,[]),title(['S',num2str(0)]);
%     figure,imshow(Cmax,[]),title(['C',num2str(0)]);
% Cmax = mat2gray(Cmax);
% for i=1:3
% % SS(:,:,i) = mat2gray(S(:,:,i));
% end
% imwrite(S,'Se.png'),imwrite(Cmax,'Ce.png')
fI3=F_temp2_detail;
end

function [fI3, i_meant,aa,N1,Im] = scale_interm(imgSeqColor,r,ii,Cl,PGway,Cfun)
[h,w]=size(imgSeqColor);
N = boxfilter(ones(h, w), r);
tem=ones(h, w);
tem(:,2:2:w)=0;
tem(2:2:h,:)=0;
N1= boxfilter(tem, r);
p=4;
[Im,imgSeqColor] = Enhancement(imgSeqColor,r,ii,Cl,PGway);
[WD, Cmax,i_mean2]= weight_cal_detail(imgSeqColor,p,r,Cfun);
WD=WD.*Cmax;
WD=WD.*Im;
F_temp2_detail=zeros(h,w);
i_meant=zeros(ceil(h/2),ceil(w/2));
aa=i_mean2(:,:).*tem;
i_meant(:,:)=aa(1:2:h,1:2:w);
tmp=i_mean2(:,:);
W_D1=boxfilter(tmp.*WD(:,:), r)./ N;
W_D2=boxfilter(WD(:,:), r)./ N;
% W_D1=tmp.*WD(:,:);
% W_D2=WD(:,:);
F_temp2_detail(:,:)=W_D2.*(imgSeqColor(:,:))-W_D1;
fI3=sum(F_temp2_detail,3);
end

function [fI3, i_meant,aa,N1,Im] = scale_coarse(imgSeqColor,r,Cl,PGway,Cfun)
[h,w]=size(imgSeqColor);
N = boxfilter(ones(h, w), r);
tem=ones(h, w);
tem(:,2:2:w)=0;
tem(2:2:h,:)=0;
N1= boxfilter(tem, r);
p=4;
[Im,imgSeqColor] = Enhancement(imgSeqColor,r,1,Cl,PGway);
[WD, Cmax,i_mean2]= weight_cal_detail(imgSeqColor,p,r,Cfun);
WD=WD.*repmat(Cmax,[1 1 ]);
WD=WD.*repmat(Im,[1 1 ]);
F_temp2_detail=zeros(h,w);
F_temp2_base=zeros(h,w);
F_temp2=zeros(h,w);
i_meant=zeros(ceil(h/2),ceil(w/2));
aa=i_mean2(:,:).*tem;
i_meant(:,:)=aa(1:2:h,1:2:w);
templ=exp((i_mean2(:,:))-0.2);
% templ=i_mean2(:,:).^0.1;
%     templ=i_mean2(:,:).*(well_exposedness(i_mean2(:,:))+0.3);
W_D1=boxfilter(i_mean2(:,:).*WD(:,:), r)./ N;
W_D2=boxfilter(WD(:,:), r)./ N;
W_B2=boxfilter(templ.*Im, r)./ N;
% W_D1=i_mean2(:,:).*WD(:,:);
% W_D2=WD(:,:);
% W_B2=templ.*Im;
F_temp2_detail(:,:)=W_D2.*(imgSeqColor(:,:))-W_D1;
F_temp2_base(:,:)= W_B2;
F_temp2(:,:)=F_temp2_detail(:,:)+F_temp2_base(:,:);
fI3=F_temp2;
end

function [Im,Iout1] = Enhancement(I,r,ii,Cl,PGway)
%% Ci is the minimum threshold for image processsing, which is designed to suppress noise
[h,w,~,~] = size(I);
N1 = boxfilter(ones(h, w), r);
Lamda = 10;
switch PGway
    case 'dehaze'
        if size(I,3)==1 % grayscale graph
            img = I;
            img(img<=Cl) = Cl;
            i_mean2= boxfilter(img, r)./ N1 + eps;
             input = imcomplement(img);
            output = imreducehaze(input);
            Iout = imcomplement(output);
        else % color graph
            img = I;
            img(img<=Cl) = Cl;
            i_mean2=boxfilter(img(:,:,1), r)./ N1 + boxfilter(img(:,:,2), r)./ N1 + boxfilter(img(:,:,3), r)./ N1;
            i_mean2= i_mean2./3 + eps;
            input = imcomplement(img);
            output = imreducehaze(input);
            Iout = imcomplement(output);
        end
    case 'lime'
        if size(I,3)==1 % grayscale graph
            img = I;
            img(img<=Cl) = Cl;
            i_mean2= boxfilter(img, r)./ N1 + eps;
            tmp =lime(img);
            Iout(:,:) = img./tmp;
        else % color graph
            img = I;
            img(img<=Cl) = Cl;
            i_mean2=boxfilter(img(:,:,1), r)./ N1 + boxfilter(img(:,:,2), r)./ N1 + boxfilter(img(:,:,3), r)./ N1;
            i_mean2= i_mean2./3 + eps;
            tmp =lime(img);
            Iout = img./repmat(tmp,[1 1 3]);
        end
    case 'pespd'
        if size(I,3)==1 % grayscale graph
            img = I;
            img(img<=Cl) = Cl;
            i_mean2= boxfilter(img, r)./ N1 + eps;
            tt = i_mean2(:,:);
            tmp = i_mean2(:,:);
            tmp(tmp > 0.5) = 0.5;
            tmp = 1./(log(Lamda*0.5./tmp)/log(Lamda))*max(tt(:));
            Iout(:,:) = img./tmp;
        else % color graph
            img = I;
            img(img<=Cl) = Cl;
            i_mean2=boxfilter(img(:,:,1), r)./ N1 + boxfilter(img(:,:,2), r)./ N1 + boxfilter(img(:,:,3), r)./ N1;
            i_mean2= i_mean2./3 + eps;
            tt = i_mean2(:,:);
            tmp = i_mean2(:,:);
            tmp(tmp > 0.5) = 0.5;
            tmp = 1./(log(Lamda*0.5./tmp)/log(Lamda))*max(tt(:));
            Iout(:,:,:) = img./repmat(tmp,[1 1 3]);
        end
    case 'mpespdf'
        if size(I,3)==1 % grayscale graph
        img = I;
        img(img<=Cl) = Cl;
        i_mean2= boxfilter(img, r)./ N1 + eps;
        Iout =My_enhance(img,Cl*0.1);
    else % color graph
        img = I;
        img(img<=Cl) = Cl;
        i_mean2=boxfilter(img(:,:,1), r)./ N1 + boxfilter(img(:,:,2), r)./ N1 + boxfilter(img(:,:,3), r)./ N1;
        i_mean2= i_mean2./3 + eps;
        Iout = My_enhance(img,Cl*0.1);
        end
end
Iout1=nol(Iout,1);% standardization
tmp2 = mean(i_mean2,3);
tt2 = max(tmp2(:));
Im=tt2;
% Im = sqrt(Im);
Im = min([Im,0.85]);
Im = max([Im,0.55]);
end

function img_out = My_enhance(I0,Cli)
I00 = 1- I0; 
I = removeHaze(I00,Cli); 
I2 = 1 - I;
img_out = I2;
end

function  I_out = removeHaze( I, Cli )
c=size(I,3);
t0=Cli;
al =mean2(I);
% al =min(I(:))
% al=0.5;
A=ones(size(I)); 
if c==3
J =min(I,[],3);
T_est = 1 - al*J;
T_est=t0+(1-t0)*T_est;
dehazed = zeros(size(I));
for c = 1:3
    dehazed(:,:,c) =(I(:,:,c) - A(:,:,c))./T_est+ A(:,:,c);
end
elseif c==1
J =0.001;
T_est = 1 - al*J;
T_est=t0+(1-t0)*T_est;
dehazed = zeros(size(I));
dehazed(:,:) =(I(:,:) - A(:,:))./T_est+ A(:,:);
end
I_out = dehazed;
end

function out = nol(in,m)
[h,w,c]=size(in);
% m=min(in(:))+1;
out=zeros([h,w,c]);
for p=1:c
    tmp=in(:,:,p);
    mi=max(min(tmp(:)),0);
    out(:,:,p)=(tmp-mi)./(max(tmp(:))-mi)*m;
end
end

function q = Weight_Guided_Filter(I, p, r, eps, C)
%% parameter description
% I is guided image
% p is pending image
% r is the size of the filter
% eps&C determines the weight of this filter

    [hei, wid]  =  size(I);
    N = boxfilter(ones(hei, wid), r); % the size of each local patch; N=(2r+1)^2 except for boundary pixels.
    mean_I   =  boxfilter(I, r) ./ N;
    mean_p   =  boxfilter(p, r) ./ N;
    mean_Ip  =  boxfilter(I.*p, r) ./ N;
    cov_Ip   =  mean_Ip - mean_I .* mean_p; % this is the covariance of (I, p) in each local patch.

    mean_II  =  boxfilter(I.*I, r) ./ N;
    var_I    =  mean_II - mean_I .* mean_I; %variance

    a = cov_Ip./(var_I+eps.*C);
    b = mean_p-a.*mean_I;
    mean_a = boxfilter(a, r) ./ N;
    mean_b = boxfilter(b, r) ./ N;
    q = mean_a .* I + mean_b; 
end

function [sMap,Cmax,i_mean2,Variance] = weight_cal_detail(I,p,r,Cfun)
R=2*r+1;

[h,w,~,~]=size(I);
N1 = boxfilter(ones(h, w), r);
   
if size(I,3)==1
    img = I;
    i_mean2=boxfilter(img, r)./ N1;
    i_var2= boxfilter(img.*img, r) ./ N1- i_mean2.* i_mean2;
    i_var2=sqrt(max(i_var2,0));
    C = i_var2 * sqrt( R^2  )+ 1e-12; % signal strengh
else
    img = I;
    i_mean2=boxfilter(img(:,:,1), r)./ N1+boxfilter(img(:,:,2), r)./ N1+boxfilter(img(:,:,3), r)./ N1;
    i_mean2= i_mean2./3;
    i_var2= boxfilter(img(:,:,1).*img(:,:,1), r)./ N1+ boxfilter(img(:,:,2).*img(:,:,2), r)./ N1+...
        boxfilter(img(:,:,3).*img(:,:,3), r)./ N1- i_mean2.* i_mean2*3;

    i_var2 = i_var2./3;
    i_var2=sqrt(max(i_var2,0));
    C = i_var2 * sqrt( 3*R^2  )+ 1e-12; % signal strengh
end
Variance=i_var2;
% the C of Cmax and sMap must be sure to match
% Cm=histeq(C);
% core= fspecial('gaussian',[7,7],100);
% C=imfilter(C,core,"replicate","same");
k=ceil(max(C(:)));
switch Cfun
    case 'threesigma'
        C1=fsigma(C,3,k);
    case 'log'
        C1=log(C+k);
    case 'gamma'
        gamma = 0.5;
        C1 = (C/k).^gamma;
    case 'multik'
        C1 = C*k;
    case 'exp'
        C1 = exp(C+k);
    case 'constk'
        C1=ones(size(C))*k;
end
% C=10*C;
% C(C<0.1) = 0.1;
% C2=nol(C,1);
% b=1.2-C2;
% b=mean(C(:));
b=sqrt(mean(C(:)));
% b=mean2(C);
% b=mean(C1(:));
% b=1.5;
% Cmax=C1;
Cmax=C1+b;
sMap1 = C1.^p; % signal structure weighting map
sMap2 = C1.^(p-1)+eps;
% sMap2(sMap2==0) = eps;
sMap1 = sMap1+eps ;
% normalizer = sum(sMap1,3)+eps; % here sMap only has one channel
% sMap = sMap2 ./ normalizer ;
sMap = sMap2 ./ sMap1 ;
% sMaptmp = sMap2 ./ sMap1 ;
end




