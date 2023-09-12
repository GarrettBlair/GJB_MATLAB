function cmap_out = shift_colormap(cmap, kappa)
if kappa>0
cmap_out = cmap + (1-cmap)/kappa;
elseif kappa<0
cmap_out = cmap - (cmap)/abs(kappa);
cmap_out = cmap_out - min(cmap_out);
end