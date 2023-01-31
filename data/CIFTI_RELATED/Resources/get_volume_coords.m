function coords = get_volume_coords(volumefile)

data = load_untouch_nii(volumefile,1);
data.img(:) = 1;
[xs,ys,zs] = ind2sub(data.hdr.dime.dim(2:4),find(data.img));
xs = (xs - 1) .* data.hdr.hist.srow_x(1) + data.hdr.hist.srow_x(4);
ys = (ys - 1) .* data.hdr.hist.srow_y(2) + data.hdr.hist.srow_y(4);
zs = (zs - 1) .* data.hdr.hist.srow_z(3) + data.hdr.hist.srow_z(4);

coords = [xs ys zs];