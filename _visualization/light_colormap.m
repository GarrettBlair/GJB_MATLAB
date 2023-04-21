function cmap_out = light_colormap(cmap_in, kappa)
if kappa>0
cmap_out = cmap_in + (1-cmap_in)/kappa;
elseif kappa<0
cmap_out = cmap_in - (cmap_in)/abs(kappa);
end