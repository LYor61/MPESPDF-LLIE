% clear;
% addpath('bm3d');
% imgSeqColor=im2double(imread("D:\myDatasets\lowlight\loldata\our485\low\2.png"));
% Img=im2double(imread("D:\myDatasets\lowlight\loldata\our485\high\2.png"));
% % r1=4; %patch size p=r1*2+1
% % noise=2;%denoising ways
% % Cl=0.005;%low threshold
% % n=4;%Z-scale
% % PGway='mpespdf';%'dehaze','lime','pespd','mpespdf' in perception gain
% % Cfun='log';%'threesigma','log','gamma','multik','exp','constk'
% 
% tic
% I_f= MPESPDF_LLIE(imgSeqColor,r1,noise,Cl,n,PGway,Cfun);
% toc
% AB=mean2(I_f);
% figure,imshow(I_f,[])
% PSNR=psnr(I_f,Img);
% SSIM=ssim(I_f,Img);
% title(['PSNR: ' num2str(PSNR) ',   SSIM: ' num2str(SSIM) ])


% % % % % bacth processing
clear
tic
addpath('bm3d');
Input_path ="C:\Users\28371\Desktop\OLT修改\SICE100\";
namelist = dir(strcat(Input_path,'*.jpg'));  
len = length(namelist);
DirCell = struct2cell(namelist);
Dir = sort_nat(DirCell(1,:)); 
r1=4; %patch size p=r1*2+1
noise=1;%denoising ways
Cl=0.005;%low threshold
n=4;%Z-scale
PGway='mpespdf';%'dehaze','lime','pespd','mpespdf' in perception gain
Cfun='log';%'threesigma','log','gamma','multik','exp','constk'
for i = 1:len
    Dir1 = Dir{i};
    image =  imread(strcat(Input_path,Dir1));% input image
    fprintf('NO.%d processing...  %s\n',i,strcat(Input_path,Dir1));
    img_in=im2double(image); 
    K2=MPESPD(img_in,r1,noise,Cl,n,PGway,Cfun);
    imwrite(K2,['SICE\',char(Dir1)]); % output path
end
toc
             