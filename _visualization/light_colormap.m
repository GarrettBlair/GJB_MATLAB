function cmap = light_colormap(cmap, kappa)
if kappa>0
cmap = cmap + (1-cmap)/kappa;
elseif kappa<0
cmap = cmap - (1-cmap)/abs(kappa);
end