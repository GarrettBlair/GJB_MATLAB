function cmap = shift_colormap(cmap, kappa)
if kappa>0
cmap = cmap + (1-cmap)/kappa;
elseif kappa<0
cmap = cmap - (cmap)/abs(kappa);
cmap = cmap - min(cmap);
end