function c=convb(g,r)
%This function conv g and out of focus blurring kernel with reflection boundary.
%function c=convb_outfocus(g,r);
%r is the radius of the point spread function.

PSF = fspecial('disk', r);
p = (size(PSF, 1) - 1) / 2;
g = padarray(g,[p p], 'symmetric');
c = conv2(g, PSF, 'valid');