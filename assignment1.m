lena = imread('lena.bmp');

% factor = 2;
% newsubsample = subsample(lena, factor);
% imshow(newsubsample);

% newshrink = shrink(lena);
% imshow(newshrink);

% newzoom = zoom(lena);
% imshow(newzoom);

% newrotate = myrotate(lena);
% imshow(newrotate);

% newreflect = reflect(lena);
% imshow(newreflect);

% fraction = 1/5;
% newdim = dim(lena,fraction);
% imshow(newdim);

newcontrastcompress = contrast_compress(lena);
imshow(newcontrastcompress);

function newimage = subsample(image, factor)  
    newimage = image(1:factor:end, 1:factor:end, :);
end

function newimage = shrink(image)  
    newimage = uint8(blockproc(image,[2 2],@(A) mean(mean(A.data), 2)));
end

function newimage = zoom(image)
      newimage = image(:,:,:);
%     [sizex, sizey, sizez]= size(image);
%     newsize = [2*sizex 2*sizey sizez];
%     newblank = zeros(newsize);
%     newimage = newblank(:,:,:);
end

function newimage = myrotate(image)
   newimage = rot90(image,-1);
end

function newimage = reflect(image)
    newimage = fliplr(image);
end

function newimage = dim(image,fraction)
    newimage = fraction*image;
end

function newimage = contrast_compress(image)
    rchan = sqrt(double(image(:,:,1)));
    gchan = sqrt(double(image(:,:,2)));
    bchan = sqrt(double(image(:,:,3)));
    
    rimg = rchan - min(rchan(:)) ;
    rnorm = rimg / max(rimg(:)) ;
    
    gimg = gchan - min(gchan(:)) ;
    gnorm = gimg / max(gimg(:)) ;
    
    bimg = bchan - min(bchan(:)) ;
    bnorm = bimg / max(bimg(:)) ;
    
    newimage = rnorm + gnorm + bnorm;
    
%     img = uint8((rchan + gchan + bchan));
%     newimage = uint8((rchan + gchan + bchan) / (3*sqrt(255)));
end