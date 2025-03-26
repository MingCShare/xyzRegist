% Batch processing xyz align zstack images to a reference zstack
% XY align were using turboreg by Miji package/ matlab-fiji,  
% instruction for Miji installation pls goto https://imagej.net/plugins/miji 
% Z align  by maximize correlation coeffcient for corresponding Z frames
% To increase speed, image were binned to calculate the CC value
%
%  mingchen@lglab.ac.cn,  @202406


% edit parameter here 

% fdir ='G:\Xian data backup\2P\20240606_gfap_mouse1\tiff\green_channel\';   % file directory

fdir ='G:\Xian data backup\2P\20241125GFAP_KGC_mCherry\Tiff\laser_ablation';

refimgname ='2_GFAP_KGC_before_z12to72_2umstep.tif';


dualchannel = 1;  % 0 for single channel and;  1 for dual channel data, 
refchn = 1;    % reference channel for turboreg

algmethod = 'translation';    % 1, 'translation' ; 2, 'rigid'


xybin = 6;        % xy binning factor for CC calculation

hfstep = 10;  % maximum shift pixel in binned image, ie, 15 in xybin 4 equal to 60 pixel in original images

shiftrange = 15;  % maximum shift range in Z axis;
upscalefactor = 10;  % up scale step to calculate shift , step equal to 1/upscalefactor

overwirteMode = 0;

%===

if ~endsWith(fdir,'\')
    fdir = strcat(fdir,'\');
end


resfolder = strcat(fdir,'analyseV2\');  % creat result folder
if ~exist(resfolder,'dir')
   mkdir(resfolder);
end

resTifffolder = strcat(fdir,'cropTiff\');  % creat result folder
if ~exist(resTifffolder,'dir')
   mkdir(resTifffolder);
end

tempfilename = struct2cell(dir([fdir,'*.tif']));

f1 = tempfilename{1,1};

fnum = size(tempfilename,2);

fznum = zeros(fnum,1);

t0=tic;

zshiftn = zeros(fnum,3);
fullstep = hfstep*2; 

for fi = 1:fnum

%           f1name = tempfilename{1,1};                          % first file name
            f1name = refimgname;
            
            f2name = tempfilename{1,fi};                         % second file name
                
               
            fprintf('Searching for best Z shift for %s... %d/%d\r',f2name,fi,fnum);
              
            ldim1 = loadtif(strcat(fdir,f1name),xybin);
            ldim2 = loadtif(strcat(fdir,f2name),xybin);            
            
            if dualchannel
                
                im1 = squeeze(ldim1(:,:,refchn,:));
                im2 = squeeze(ldim2(:,:,refchn,:));  
            else

                im1 = ldim1;
                im2 = ldim2; 
                
            end
            
            
            sz1 = size(im1);
            sz2 = size(im2);
            
            
            fznum (fi)= sz2(3);
            
            ccmat = zeros(sz1(3),sz2(3),fullstep,fullstep);
            
            scaleF = sz1(1)/256;
            xycenter = round([sz1(2)/2,sz1(1)/2]);  % x, y in image dimension 
            rectxy = round(100*scaleF);

            CenterRect = zeros(1,4);  % [y0, y1, x0, x1]
            CenterRect(1) = max(hfstep,xycenter(2)-rectxy) + 1;
            CenterRect(2) = min(sz1(2)-hfstep,xycenter(2)+rectxy);
            CenterRect(3) = max(hfstep,xycenter(1)-rectxy) + 1;
            CenterRect(4) = min(sz1(1)-hfstep,xycenter(1)+rectxy);

            if sum(hfstep>xycenter-rectxy)>0
                fprintf('pls use smaller shiftstep - %d in current xybin - %d\r',hfstep,xybin);
            end

           
            f1pre = regexp(f1name,'\.','split');
            f2pre = regexp(f2name,'\.','split');
            savematname= strcat(resfolder,f1pre{1},'-',f2pre{1},'-xybin',num2str(xybin),'_CC.mat');
            
            
            if exist(savematname,'file')==2 && overwirteMode==0
            
               loadcc = load(savematname);
               ccmat = loadcc.ccmat;
               fprintf('load existed CC mat, skip calculation\r');
            
            else
            
            
                    for ii = 1:sz1(3)
                       refimg = im1(CenterRect(1):CenterRect(2),CenterRect(3):CenterRect(4),ii);
                    
                       for jj = 1:sz2(3)
                    
                           for xi = 1:fullstep
                               for yi = 1:fullstep
                                   MovRect = CenterRect-hfstep + [xi,xi,yi,yi];
                                   movimg = im2(MovRect(1):MovRect(2),MovRect(3):MovRect(4),jj);
                                   ccmat(ii,jj,xi,yi) = corr2(refimg,movimg);
                               end
                           end 
                       end
                    
                       if ii==1 || mod(ii,10)==0
                            disp(['Calculate CC, frame....',num2str(ii),' / ',num2str(sz1(3))]);
                       end
                    
                    end
                    
                    
                    save(savematname,'ccmat');
            
            end
            
            mcc = squeeze(mean(mean(ccmat,2),1));
            
            [r,s]=find(mcc==max(mcc(:)));
            
            
            f1=figure;
            set(gcf,'units','points','position',[100,100,1500,600])
            subplot(1,2,1);
            hold on;
            contour(imresize(mcc,3));
            colormap winter
           
            
            plot([1,fullstep*3],[hfstep,hfstep]*3,'k-');            
            plot([hfstep,hfstep]*3,[1,fullstep*3],'k-');
            plot(hfstep*3,hfstep*3,'ro');
            text(hfstep*3+1,hfstep*3+3,"0,0",'fontsize',14)
            
            set(gca,'XTick',[1,hfstep*3,hfstep*6],'XTickLabel',{strcat('-',num2str(hfstep*xybin)),'0',num2str(hfstep*xybin)});
            set(gca,'YTick',[1,hfstep*3,hfstep*6],'YTickLabel',{strcat('-',num2str(hfstep*xybin)),'0',num2str(hfstep*xybin)});
            title('average xyshift');
            xlabel('X-shiftPix');
            ylabel('Y-shiftPix');
            axis equal
            set(gca,'fontsize',14)
            hold off
            
            
                       
            maxmat =squeeze(ccmat(:,:,r,s));
            
%             if maxmat(1,sz2(3))>= maxmat(sz1(3),1)
%                  [~,sftmat]= max(maxmat,[],1);
%                  fidx = [1:sz2(3);sftmat]';
%             elseif  maxmat(1,sz2(3))< maxmat(sz1(3),1)
%                  [~,sftmat]= max(maxmat,[],2);
%                  fidx = [sftmat';1:sz1(3)]';
%             end
%                 
% %             fidx = [1:sz1(3);sftmat]';           
%                         
%             if min(sftmat)==1
%                  frm = find(sftmat>1,1,'first');
%             
%                  fidx(1:frm-2,:)=[];
%             
%             else
%                  frm = find(sftmat==sz1(3),1,'first');
%             
%                  fidx(frm+1:end,:)=[];
%             
%             end
%             
%             
%             ft = fittype( 'poly1' );
%             [fitresult, gof] = fit( fidx(:,1), fidx(:,2), ft );           
%             k0 = abs(fitresult.p1);            
%             sft0 = fitresult.p2;            
%             zsft = round(abs(sft0));
%             
%             linexx = 1:0.2:sz1(3);
%             lineyy = fitresult(linexx);
            
            
            
            
             rszmat = imresize(maxmat,upscalefactor);
             nstep = 1/upscalefactor;
             sftzi = -shiftrange:nstep:shiftrange;

             slen = length(sftzi);

             mmval = zeros(slen,1);

             xx = 1:size(rszmat,2);
            for si = 1:slen

               yy = round(xx + sftzi(si)*upscalefactor);

               xx(yy<=0)=[];
               yy(yy<=0)=[];

               xx(yy>size(rszmat,1))=[];
               yy(yy>size(rszmat,1))=[];
                temp = zeros(length(yy),1);

                for ti = 1:length(yy)
                   temp(ti) = rszmat(yy(ti),xx(ti));
                end
               mmval(si) = mean(temp);

            end
         [mval,midx]=max(mmval);

%             ft = fittype( 'poly1' );
% %             [fitresult, gof] = fit( fidx(:,1), fidx(:,2), ft );  
%             [fitresult, gof] = fit( mx, my, ft );   
%             k0 = abs(fitresult.p1);            
%             sft0 = fitresult.p2;  
            k0=1;
            sft0 =sftzi(midx);
            zsft = round(abs(sft0));
            
            linexx = 1:0.2:sz2(3);
            lineyy = linexx + sft0;
            
            
            
            
            
            subplot(1,2,2);
            hold on
            imagesc(imresize(maxmat,3));
            colormap default
            colorbar; 
            caxis([0.4 1]);
                        
%             plot((1:sz1(3))*3,sftmat*3,'w.');           
            linekp = ((lineyy>0).*(lineyy<sz1(3)))>0;            
            plot(linexx(linekp)*3,lineyy(linekp)*3,'b-','linewidth',1);
            
%             text(5,20*3,strcat("y = ",num2str(fitresult.p1)," * X + (",num2str(fitresult.p2),")"),'fontsize',14,'color','red');
            text(5,20*3,strcat("y = X + (",num2str(sft0),")"),'fontsize',14,'color','red');        
            set(gca,'xlim',[0,sz1(3)*3+1],'ylim',[0,sz1(3)*3+1]);            
            set(gca, 'YDir','reverse');            
            set(gca,'XTick',[1,sz1(3)*3],'XTickLabel',{'1',strcat(num2str(sz1(3)))});
            set(gca,'YTick',[1,sz1(3)*3],'YTickLabel',{'1',strcat(num2str(sz1(3)))});            
                       
            title('fitted max CC');
            set(gca,'fontsize',14)
            
             xlabel('T2 Frame');
             ylabel('T1 Frame');
             axis equal
             hold off
            

           sgtitle(strcat('zAlign',f1pre{1},'-',f2pre{1}),'fontsize',16);

           savename = strcat(resfolder,'zAlign',f1pre{1},'-',f2pre{1},'.Tiff');
           f=getframe(gcf);
           imwrite(f.cdata,savename,'tif');
           close(f1);       
            
           zshiftn(fi,:)=[zsft,k0,sft0]; 

           fprintf('Frame 1 in %s equal to Frame %4.1f in %s\r\n', f1name, -sft0 + 1, f2name);

%             end

           te = toc(t0);
          
           fprintf('Searching Cost %6.1f min\r\n', te/60);


end

saveshift = strcat(resfolder,'xybin',num2str(xybin),'_shift.mat');

% if exist(saveshift,'file')==2 && overwirteMode==0        
%            loadcc = load(saveshift);
%            zshiftn = loadcc.zshiftn;
%            fprintf('load existed zshiftn mat, skip calculation\r'); 
% 
%            figure;plot(zshiftn(:,3));
%            title('zshift for all images');
%            set(gca,'fontsize',14);          
%            xlabel('image file');
%            ylabel('z shift index');
% else
           save(saveshift,'zshiftn');            
           f2=figure;plot(zshiftn(:,3));
           title('zshift for all images');
           set(gca,'fontsize',14);          
           xlabel('image file');
           ylabel('z shift index');            
           savename = strcat(resfolder,'zShift_xybin',num2str(xybin),'.Tiff');
           f=getframe(gcf);
           imwrite(f.cdata,savename,'tif');
           close(f2);  
% end



frmend = min(fznum - zshiftn(:,1).*(zshiftn(:,3)<0));

rmRef =  round(max(zshiftn(:,3)));  % delete initial unmatched frames in ref stack

CommFrame = frmend-rmRef;

% rmAlg = round(-min(zshiftn(:,3)));    % delete initial unmatched frames in align stack
% rmRef =  round(max(zshiftn(:,3)));  % delete initial unmatched frames in ref stack
% 
% CommFrame = sz1(3)-rmAlg-rmRef;


Miji
inter = ij.macro.Interpreter;
inter.batchMode = true;    

fprintf('Image XY Alignment... \r\n');


for fi = 1:fnum

         Refname = refimgname;                          % first file name            

         Aligname = tempfilename{1,fi};                         % second file name
         
         Aligpre = regexp(Aligname,'\.','split');

         fprintf('Align image %s  ...   %d/%d \r', Aligname, fi,fnum);

         zsft = zshiftn(fi,1); 
         k0 = zshiftn(fi,2); 
         sft0 = zshiftn(fi,3); 

        if abs(k0-1)<0.3

                imf1 = loadtif(strcat(fdir,Refname),1);
                imf2 = loadtif(strcat(fdir,Aligname),1);

               

                sz = size(imf1); 
                
                if length(sz)==3
                    
                   imf1 = reshape(imf1,size(imf1,1),size(imf1,2),1,size(imf1,3));
                   imf2 = reshape(imf2,size(imf2,1),size(imf2,2),1,size(imf2,3));
                   refchn = 1;
                    
                end

                regall = imf1*0;
                
                switch algmethod
                    case 'translation'
                       regxy = zeros(CommFrame,4);
                    case 'rigid'
                       regxy = zeros(CommFrame,3,4);
                end
                
                
                
                
                for ii = 1:CommFrame
                                           
                       refimg = imf1(:,:,refchn,ii + rmRef);
                       shiftmove = zsft* ((sft0<0)-0.5)*2;
                       movimg = imf2(:,:,refchn,ii + rmRef + shiftmove); 
                    
                       [~,regres] = turboreg(refimg,movimg,sz(2),sz(1),algmethod);

%                        regall(:,:,ii) = regframe;  
                     if size(regres,1)>1
                          regxy(ii,:,:) = regres;
                     else
                          regxy(ii,:) = regres;
                     end

                   if ii==1 || mod(ii,10)==0
                       disp(['Align frame....',num2str(ii),' / ',num2str(sz1(3))]);
                   end
                end
                
                regcordinate = squeeze(mean(regxy(1:CommFrame ,:,:),1));

                for ii = 1: CommFrame                    
                       for ci = 1:size(imf1,3)
                           shiftmove = zsft* ((sft0<0)-0.5)*2;
                           movimg = imf2(:,:,ci,ii + rmRef + shiftmove);
                        
                           [regframe] = turbotransform(movimg,sz(2),sz(1),regcordinate,algmethod);
                           regall(:,:,ci,ii) = regframe;       
                       end
                end
              
                regall(:,:,:,ii+1:end)=[];
                
                regall = squeeze(regall);               
                                          
                WriteTiff(regall,strcat(resTifffolder,Aligpre{1},'_Reg_Rm',num2str(rmRef),'_',num2str(zsft),'.tif'),16);              

        end

end
    inter.batchMode = false;
    MIJ.exit 





function img=loadtif(imgname,binnum)
   
    finfo =imfinfo(imgname); 
    t = finfo(1);
    
    xpix = t.Width;
    ypix = t.Height;
    
    imgdesc = t.ImageDescription;
    
    imgpa = regexp(imgdesc,'\W','split');
    
    timgdesc = imgpa{1,find(strcmp('images',imgpa)==1)+1};
    if ~isempty(timgdesc)      
         numimg = str2double(imgpa{1,find(strcmp('images',imgpa)==1)+1});
    else
         numimg = length(finfo);
        
    end

    if sum(strcmp('channels',imgpa)==1)>0
          numchannel = str2double(imgpa{1,find(strcmp('channels',imgpa)==1)+1});
    else
        numchannel = 1;
    end
    if sum(strcmp('slices',imgpa)==1)>0
         numzstep = str2double(imgpa{1,find(strcmp('slices',imgpa)==1)+1});
    else
         numzstep = 1;

    end
    if sum(strcmp('frames',imgpa)==1)>0
            numTframe = str2double(imgpa{1,find(strcmp('frames',imgpa)==1)+1});
    else
            numTframe = 1;
    end
    
    
    img = zeros(ypix,xpix,numimg);
    for mi = 1:numimg
    
         img(:,:,mi) = imread(imgname,mi);
    
    end


    
     if binnum>1
       ylen = floor(ypix/binnum)*binnum;
       xlen = floor(xpix/binnum)*binnum;
       timg = img(1:ylen,1:xlen,1:numimg); 
       txpix = size(timg,2);
       typix = size(timg,1);

       cimg = reshape(timg,binnum,typix/binnum,txpix,numimg);
       cimg = squeeze(mean(cimg,1));

       cimg = reshape(cimg,typix/binnum,binnum,txpix/binnum,numimg);
       cimg = squeeze(mean(cimg,2));


       img = cimg;
       imy = floor(ypix/binnum);
       imx = floor(xpix/binnum);
     else
         imy = ypix;
         imx = xpix;
         cimg= img;
     end
     img=squeeze(reshape(cimg, imy, imx,numchannel,numzstep,numTframe));
end


function [regimg,res] = turboreg(ref,move,alignx,aligny,method)

       MIJ.createImage('source',uint16(move),true);                           

       MIJ.createImage('target',uint16(ref),true);          

       alignSource = strcat(' -window source',{' 1 1 '},num2str(alignx),{' '},num2str(aligny));
       alignTarget = strcat(' -window target',{' 1 1 '},num2str(alignx),{' '},num2str(aligny));

       switch method

           case 'translation'

                   alignMethod = strcat(' -translation ',{' '},num2str(alignx/2),{' '},num2str(aligny/2),{' '},num2str(alignx/2),{' '},num2str(aligny/2),' -showOutput');

           case 'rigid'

               
                   alignMethod = strcat(' -rigidBody ',{' '},num2str(alignx/2),{' '},num2str(aligny/2),{' '},num2str(alignx/2),{' '},num2str(aligny/2),...
                       {' '},num2str(alignx/2),{' '},num2str(aligny/8),{' '},num2str(alignx/2),{' '},num2str(aligny/8),...
                       {' '},num2str(alignx/2),{' '},num2str(aligny*7/8),{' '},num2str(alignx/2),{' '},num2str(aligny*7/8),...
                       ' -showOutput');
       end
            
                   MIJ.run('TurboReg ',strcat(' -align ',alignSource,alignTarget,alignMethod)); 
                 
                   temp = MIJ.getCurrentImage;   
                   
                   res = MIJ.getResultsTable;  % sX ,sY, tX,tY   



       regimg = uint16(temp(:,:,1));            
             
       MIJ.run('Close All')      

end


function [timg] = turbotransform(move,alignx,aligny,regcordinate,method)

 MIJ.createImage('source',uint16(move),true);
 tfmSource = strcat(' -window source',{' '},num2str(alignx),{' '},num2str(aligny));

  switch method

           case 'translation'
               tfmMethod = strcat(' -translation  ',{' '},num2str(regcordinate(1)), {' '},num2str(regcordinate(2)), ...
                             {' '},num2str(regcordinate(3)),{' '},num2str(regcordinate(4)), ' -showOutput');

          case 'rigid'

               tfmMethod = strcat(' -rigidBody ',...
                   {' '},num2str(regcordinate(1,1)), {' '},num2str(regcordinate(1,2)),  {' '},num2str(regcordinate(1,3)),{' '},num2str(regcordinate(1,4)),...
                   {' '},num2str(regcordinate(2,1)), {' '},num2str(regcordinate(2,2)),  {' '},num2str(regcordinate(2,3)),{' '},num2str(regcordinate(2,4)),... 
                   {' '},num2str(regcordinate(3,1)), {' '},num2str(regcordinate(3,2)),  {' '},num2str(regcordinate(3,3)),{' '},num2str(regcordinate(3,4)),...  
                             ' -showOutput');


  end

 MIJ.run('TurboReg ', strcat(' -transform ', tfmSource,tfmMethod ));
 
 temp = MIJ.getCurrentImage; 
 timg = uint16(temp(:,:,1));
 
 MIJ.run('Close All')

end



function WriteTiff(data,filename,bitdepth)




sz = size(data);
nframe = prod(sz)/sz(1)/sz(2);

arg1 = [' ImageJ=1.52p' ... 
            'images=' num2str(nframe)... 
            ];

if length(sz) == 2
    arg2 = [];
elseif length(sz) == 3

    arg2 = [' frames=' num2str(sz(3)) ... 
            ];

elseif length(sz) == 4

    arg2 = [' channels=' num2str(sz(3)) ...
            ' frames=' num2str(sz(4)) ... 
            ];

elseif length(sz) == 5
    arg2 = [' channels=' num2str(sz(3)) ...
            ' slices=' num2str(sz(4)) ...
            ' frames=' num2str(sz(5)) ... 
            ];

end
           
arg3 = [    ' hyperstack=true' ...
            ' mode=grayscale' ...  
            ' loop=false' ...  
            ' min=0.0' ...      
            ' max=65535.0'];  % change this to 256 if you use an 8bit image

fiji_descr = strcat(arg1,arg2,arg3);
            
% t = Tiff('test.tif','w')
tagstruct.ImageLength = sz(1);
tagstruct.ImageWidth = sz(2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.Compression = Tiff.Compression.LZW;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
tagstruct.ImageDescription = fiji_descr;


tf = Tiff(filename,'w');

if length(sz)==2
      tf.setTag(tagstruct);
       switch bitdepth
        case 8
             tf.write(uint8(data)); 
        case 16
             tf.write(uint16(data));
        otherwise
            fprintf('please specify tiff file bitdepth is 8 or 16');
       end
        tf.writeDirectory();   

elseif length(sz)==3
        
    for frameid=1:sz(3)

         tf.setTag(tagstruct);

         Fn = data(:,:,frameid);
            switch bitdepth
                case 8
                     tf.write(uint8(Fn)); 
                case 16
                     tf.write(uint16(Fn));
                otherwise
                    fprintf('please specify tiff file bitdepth is 8 or 16');
            end

         tf.writeDirectory();   

    end

elseif length(sz)==4
        
    for dim3=1:sz(3)
         for dim4 = 1:sz(4) 
             tf.setTag(tagstruct);
    
             Fn = data(:,:,dim3,dim4);
                switch bitdepth
                    case 8
                         tf.write(uint8(Fn)); 
                    case 16
                         tf.write(uint16(Fn));
                    otherwise
                        fprintf('please specify tiff file bitdepth is 8 or 16');
                end
    
              tf.writeDirectory();   

        end
    end

end

 tf.close;
end