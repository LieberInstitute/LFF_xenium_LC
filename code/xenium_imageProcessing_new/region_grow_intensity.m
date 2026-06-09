function BWrg = region_grow_intensity(I, seedMask, lo, hi)
% I: grayscale (double 0..1). seedMask: logical seeds (e.g., your BW).
% lo/hi: intensity band to grow through (e.g., 0.1, 0.9)

I = im2double(I);
pass = I >= lo & I <= hi;          % where growth is allowed
marker = seedMask & pass;          % seeds that satisfy the band
BWrg = imreconstruct(marker, pass);
BWrg = imfill(BWrg, 'holes');
end