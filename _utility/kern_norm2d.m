function pfield_smooth = kern_norm2d(pfield,kern)

[fieldx, fieldy] = size(pfield);
pfield_smooth = NaN(fieldx, fieldy);

[kernx, kerny] = size(kern);
for r = 1:fieldx
    for c = 1:fieldy;
        target = pfield(r,c);
        if ~isnan(target)
            kx1 = (kernx+1)/2 - 1;
            kx2 = (kernx+1)/2 - 1;
            ky1 = (kerny+1)/2 - 1;
            ky2 = (kerny+1)/2 - 1;
            if r-kx1<=0
                kx1 = r-1;
            elseif kx2+r>fieldx
                kx2 = fieldx-r;
            end
            if c-ky1<=0
                ky1 = c-1;
            elseif ky2+c>fieldy
                ky2 = fieldy-c;
            end
            neigh = pfield([r-kx1:r+kx2],[c-ky1:c+ky2]);
            ktemp = kern([(kernx+1)/2 - kx1: (kernx+1)/2 + kx2], [(kerny+1)/2 - ky1: (kerny+1)/2 + ky2]);
            pfield_smooth(r,c) = nansum(nansum(((neigh.*ktemp)/sum(ktemp(:)))));            
        end
    end
end