function C_out = fusion(imgSeqColor,r1,noise,Cl,n,PGway,Cfun)
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
% m=0.5*mean(imgSeqColor(:))
% a=1-m;b=1+m;
% a=0.75;b=1.2;
a=0.75;b=1.2;
C_out=D1*a+repmat(B1*b,[1 1 3]);
% imwrite(B1,"L1.png")
% imwrite(D1,"H2.png")
% imwrite(B2,"L2.png")
% imwrite(D2{1},"H2.png")
end

