import matplotlib.pyplot as plt
import numpy as np
from readData.DataLoader import DataLoader
import analysis.analyze as analyze
import util.utilities as util
import ctypes
import visualization.colors as viscolors

def hist(path, snap_num, part_types, side_length=1, thickness=1, center=[0,0,0]):

    keys = ['Coordinates', 'Masses']
    cat = DataLoader(path, snap_num, part_types, keys)

    side_length *= cat.boxsize
    thickness *= cat.boxsize

    if type(center)==type([]):
        center = np.array(center)

    offset_pos = np.zeros((0,3))
    mass = np.zeros(0)
    for i in range(6):
        key = f'PartType{i}/Masses'
        if key in cat.data.keys():
            mass = np.concatenate((mass, cat[key]))
            offset_pos = np.concatenate((offset_pos, cat[f'PartType{i}/Coordinates']))

    is_in_box = analyze.get_box_cut(offset_pos, np.array([0,0,0]), side_length)
    in_box = offset_pos[is_in_box]
    mass_in_box = mass[is_in_box]

    is_in_thickness_1 = analyze.get_box_cut(offset_pos[is_in_box], np.array([0,0,0]), [side_length, side_length, thickness])
    is_in_thickness_2 = analyze.get_box_cut(offset_pos[is_in_box], np.array([0,0,0]), [side_length, thickness, side_length])
    is_in_thickness_3 = analyze.get_box_cut(offset_pos[is_in_box], np.array([0,0,0]), [thickness, side_length, side_length])

    xy_data = in_box[is_in_thickness_1]
    xz_data = in_box[is_in_thickness_2]
    yz_data = in_box[is_in_thickness_3]

    xy_weights = mass_in_box[is_in_thickness_1]
    xz_weights = mass_in_box[is_in_thickness_2]
    yz_weights = mass_in_box[is_in_thickness_3]

    return xy_data, xz_data, yz_data, xy_weights, xz_weights, yz_weights

## 
## routine to use 'raytrace_projection_compute' to make mock stellar images, 
##   treating the starlight as sources and accounting for gas extinction via ray-tracing
## 
## important to set KAPPA_UNITS appropriately: code loads opacity (kappa) for 
##   the bands of interest in cgs (cm^2/g), must be converted to match units of input 
##   mass and size. the default it to assume gadget units (M=10^10 M_sun, l=kpc)
##
## from torreylabtools/Python/visualization/raytrace_projection.py
def stellar_raytrace( BAND_IDS, \
        stellar_x, stellar_y, stellar_z, \
        stellar_mass, stellar_age, stellar_metallicity, stellar_hsml, \
        gas_x, gas_y, gas_z, gas_mass, gas_metallicity, gas_hsml, \
        xrange=0, yrange=0, zrange=0, pixels=720, 
        KAPPA_UNITS=2.08854068444, \
        IMF_CHABRIER=1, IMF_SALPETER=0 , \
        ADD_BASE_METALLICITY=0.0, ADD_BASE_AGE=0.0 ):
        
    Nbands=len(np.array(BAND_IDS)); Nstars=len(np.array(stellar_mass)); Ngas=len(np.array(gas_mass));
    if (Nbands != 3): print("stellar_raytrace needs 3 bands, you gave",Nbands); return -1,-1,-1,-1;
    ## check if stellar metallicity is a matrix
    if (len(stellar_metallicity.shape)>1): stellar_metallicity=stellar_metallicity[:,0];
    if (len(gas_metallicity.shape)>1): gas_metallicity=gas_metallicity[:,0];

    ## get opacities and luminosities at frequencies we need:
    stellar_metallicity[stellar_metallicity>0] += ADD_BASE_METALLICITY;
    gas_metallicity[gas_metallicity>0] += ADD_BASE_METALLICITY;
    stellar_age += ADD_BASE_AGE;
    kappa=np.zeros([Nbands]); lums=np.zeros([Nbands,Nstars]);
    for i_band in range(Nbands):
        nu_eff = util.colors_table(np.array([1.0]),np.array([1.0]), \
            BAND_ID=BAND_IDS[i_band],RETURN_NU_EFF=1);
        kappa[i_band] = util.opacity_per_solar_metallicity(nu_eff);
        l_m_ssp = util.colors_table( stellar_age, stellar_metallicity/0.02, \
            BAND_ID=BAND_IDS[i_band], CHABRIER_IMF=IMF_CHABRIER, SALPETER_IMF=IMF_SALPETER, CRUDE=1, \
            UNITS_SOLAR_IN_BAND=1); ## this is such that solar-type colors appear white
        l_m_ssp[l_m_ssp >= 300.] = 300. ## just to prevent crazy values here 
        l_m_ssp[l_m_ssp <= 0.] = 0. ## just to prevent crazy values here 
        lums[i_band,:] = stellar_mass * l_m_ssp
    gas_lum=np.zeros(Ngas); ## gas has no 'source term' for this calculation
    stellar_mass_attenuation = np.zeros(Nstars); ## stars have no 'attenuation term'
    gas_mass_metal = gas_mass * (gas_metallicity/0.02);
    kappa *= KAPPA_UNITS;
    
    ## combine the relevant arrays so it can all be fed into the ray-tracing
    x=np.concatenate([stellar_x,gas_x]); y=np.concatenate([stellar_y,gas_y]); z=np.concatenate([stellar_z,gas_z]);
    mass=np.concatenate([stellar_mass_attenuation,gas_mass_metal]);
    hsml=np.concatenate([stellar_hsml,gas_hsml]);
    wt1=np.concatenate([lums[0,:],gas_lum]); wt2=np.concatenate([lums[1,:],gas_lum]); wt3=np.concatenate([lums[2,:],gas_lum]);
    k1=kappa[0]; k2=kappa[1]; k3=kappa[2];
       
    return raytrace_projection_compute(x,y,z,hsml,mass,wt1,wt2,wt3,k1,k2,k3,\
        xrange=xrange,yrange=yrange,zrange=zrange,pixels=pixels,TRIM_PARTICLES=1);

def checklen(x):
    return len(np.array(x,ndmin=1));
def int_round(x):
    return np.int(np.round(x));
def ok_scan(input,xmax=1.0e30,pos=0):
    if (pos==0):
        return (np.isnan(input)==False) & (np.isfinite(input)) & (np.fabs(input)<=xmax);
    if (pos==1):
        return (np.isnan(input)==False) & (np.isfinite(input)) & (np.fabs(input)<=xmax) & (input > 0.);
    if (pos==2):
        return (np.isnan(input)==False) & (np.isfinite(input)) & (np.fabs(input)<=xmax) & (input >= 0.);
def fcor(x):
    return np.array(x,dtype='f',ndmin=1)
def vfloat(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));
def single_vec_sorted(x,reverse=False):
    return sorted(np.reshape(x,x.size),reverse=reverse);
def clip_256(x,max=255,min=2):
    x = x*(max-min) + min;
    x[x >= max]=max;
    x[x <= min]=min;
    x[np.isnan(x)]=min;
    x /= 256.;
    return x;


##
##  Wrapper for raytrace_rgb, program which does a simply line-of-sight projection 
##    with multi-color source and self-extinction along the sightline: here called 
##    from python in its most general form: from the c-code itself:
##
##  int raytrace_rgb(
##    int N_xy, // number of input particles/positions
##    float *x, float *y, // positions (assumed already sorted in z)
##    float *hsml, // smoothing lengths for each
##    float *Mass, // total weight for 'extinction' part of calculation
##    float *wt1, float *wt2, float *wt3, // weights for 'luminosities'
##    float KAPPA1, float KAPPA2, float KAPPA3, // opacities for each channel
##    float Xmin, float Xmax, float Ymin, float Ymax, // boundaries of output grid
##    int Xpixels, int Ypixels, // dimensions of grid
##    float *OUT0, float *OUT1, float *OUT2, float*OUT3 ) // output vectors with final weights
##
## from torreylabtools/Python/visualization/raytrace_projection.py
def raytrace_projection_compute( x, y, z, hsml, mass, wt1, wt2, wt3, \
    kappa_1, kappa_2, kappa_3, xrange=0, yrange=0, zrange=0, pixels=720, \
    TRIM_PARTICLES=1 ):

    #print zrange

    ## define bounaries
    if(checklen(xrange)<=1): xrange=[np.min(x),np.max(x)];
    if(checklen(yrange)<=1): yrange=[np.min(y),np.max(y)];
    if(checklen(zrange)<=1): zrange=[np.min(z),np.max(z)];
    xr=xrange; yr=yrange; zr=zrange;
    x00=0.5*(xr[1]+xr[0]); y00=0.5*(yr[1]+yr[0]); z00=0.5*(zr[1]+zr[0]); 
    tolfac = 1.0e10;
    if (TRIM_PARTICLES==1): tolfac = 0.05; 

    ## clip to particles inside those
    xlen=0.5*(xr[1]-xr[0]); ylen=0.5*(yr[1]-yr[0]); zlen=0.5*(zr[1]-zr[0]);
    x-=x00; y-=y00; z-=z00; dx=xlen*(1.+tolfac*2.); dy=ylen*(1.+tolfac*2.); dz=zlen*(1.+tolfac*2.);
    ok=ok_scan(x,xmax=dx) & ok_scan(y,xmax=dy) & ok_scan(z,xmax=dz) & \
        ok_scan(hsml,pos=1) & ok_scan(mass+wt1+wt2+wt3,pos=1)
        #& ok_scan(mass) & ok_scan(wt1) & ok_scan(wt2) & ok_scan(wt3);
    x=x[ok]; y=y[ok]; z=z[ok]; hsml=hsml[ok]; mass=mass[ok]; wt1=wt1[ok]; wt2=wt2[ok]; wt3=wt3[ok];
    N_p=checklen(x); xmin=-xlen; xmax=xlen; ymin=-ylen; ymax=ylen;
    if(N_p<=1): 
        print(' UH-OH: EXPECT ERROR NOW (raytrace_projection), there are no valid source/gas particles to send!'); return -1,-1,-1,-1;

    ## now sort these in z (this is critical!)
    s=np.argsort(z);
    x=x[s]; y=y[s]; z=z[s]; hsml=hsml[s]; mass=mass[s]; wt1=wt1[s]; wt2=wt2[s]; wt3=wt3[s];
    ## cast new copies to ensure the correct formatting when fed to the c-routine:
    x=fcor(x); y=fcor(y); z=fcor(z); hsml=fcor(hsml); mass=fcor(mass); wt1=fcor(wt1); wt2=fcor(wt2); wt3=fcor(wt3);

    ## load the routine we need
    exec_call=util.return_python_routines_cdir()+'/RayTrace_RGB/raytrace_rgb.so'
    routine=ctypes.cdll[exec_call];
    
    ## cast the variables to store the results
    aspect_ratio=ylen/xlen; Xpixels=int_round(pixels); Ypixels=int_round(aspect_ratio*np.float(Xpixels));

    print("some diags in raytrace projection ")
    print("aspect_ratio = {:f}".format(aspect_ratio))
    print("Xpixels = {:d}".format(Xpixels))
    print("Ypixels = {:d}".format(Ypixels))
    N_pixels=Xpixels*Ypixels; out_cast=ctypes.c_float*N_pixels; 
    out_0=out_cast(); out_1=out_cast(); out_2=out_cast(); out_3=out_cast();

    ## main call to the calculation routine
    routine.raytrace_rgb( ctypes.c_int(N_p), \
        vfloat(x), vfloat(y), vfloat(hsml), vfloat(mass), \
        vfloat(wt1), vfloat(wt2), vfloat(wt3), \
        ctypes.c_float(kappa_1), ctypes.c_float(kappa_2), ctypes.c_float(kappa_3), \
        ctypes.c_float(xmin), ctypes.c_float(xmax), ctypes.c_float(ymin), ctypes.c_float(ymax), \
        ctypes.c_int(Xpixels), ctypes.c_int(Ypixels), \
        ctypes.byref(out_0), ctypes.byref(out_1), ctypes.byref(out_2), ctypes.byref(out_3) );

    ## now put the output arrays into a useful format 
    out_0 = np.copy(np.ctypeslib.as_array(out_0));
    out_1 = np.copy(np.ctypeslib.as_array(out_1));
    out_2 = np.copy(np.ctypeslib.as_array(out_2));
    out_3 = np.copy(np.ctypeslib.as_array(out_3));
    out_0 = out_0.reshape([Xpixels,Ypixels]);
    out_1 = out_1.reshape([Xpixels,Ypixels]);
    out_2 = out_2.reshape([Xpixels,Ypixels]);
    out_3 = out_3.reshape([Xpixels,Ypixels]);

    return out_0, out_1, out_2, out_3;





def make_threeband_image_process_bandmaps(r,g,b, \
    dont_make_image=0, maxden=0, dynrange=0, pixels=720, \
    color_scheme_nasa=1, color_scheme_sdss=0 , \
    filterset = ['r','g','b'], **kwargs ):

    ## now clip the maps and determine saturation levels
    cmap_m=np.zeros((checklen(r[:,0]),checklen(r[0,:]),3),dtype='f');
    cmap_m[:,:,0]=r; cmap_m[:,:,1]=g; cmap_m[:,:,2]=b; 
    if (dont_make_image==1): return cmap_m;

    if (maxden<=0):
        f_saturated=0.005 ## fraction of pixels that should be saturated 
        x0=int_round( f_saturated * (np.float(checklen(r)) - 1.) );
        for rgb_v in [r,g,b]: 
            rgbm=single_vec_sorted(rgb_v,reverse=True);
            if(rgbm[x0]>maxden): maxden=rgbm[x0]

    if (dynrange<=0):
        f_zeroed=0.1 ## fraction of pixels that should be black 		
        x0=int_round( f_zeroed * (np.float(checklen(r)) - 1.) ); minden=np.max(r); 
        rgbm=single_vec_sorted(r+g+b,reverse=False); 
        if(rgbm[x0]<minden): minden=rgbm[x0];
        #for rgb_v in [r,g,b]: 
        #    rgbm=single_vec_sorted(rgb_v,reverse=False); 
        #    if(rgbm[x0]<minden): minden=rgbm[x0];
        if (minden<=0):
            minden = np.min(np.concatenate((r[r>0.],g[g>0.],b[b>0.])));
        dynrange = maxden/minden;

    ## now do the color processing on the maps
    maxnorm=maxden; minnorm=(maxden/dynrange);
    #print 'maxnorm == ',maxnorm,' dynrange == ',dynrange,' minnorm == ',minnorm;

    i = (r+g+b)/3.
 
    #print "min/max/mean/median of i"
    #print np.min(i), np.max(i), np.mean(i), np.median(i)
    #print "min/max/mean/median of r"
    #print np.min(r), np.max(r), np.mean(r), np.median(r)
    #print "min/max/mean/median of g"
    #print np.min(g), np.max(g), np.mean(g), np.median(g)
    #print "min/max/mean/median of b"
    #print np.min(b), np.max(b), np.mean(b), np.median(b)

    if np.isfinite( np.log10(maxnorm/minnorm) ) and np.log10(maxnorm/minnorm) != 0:
        f_i = np.log10(i/minnorm) / np.log10(maxnorm/minnorm);	# nominally goes from 0 to 1
    else:
        f_i = np.log10(i/minnorm)

    #print "number of saturated pixels"
    #print np.sum( i >=maxnorm )

    f_i[i>=maxnorm]=1.; f_i[i<=minnorm]=0.; f_i[i==0]=0.;
    bad=(i<=0)|(np.isnan(i))|(np.isnan(f_i))
    f_i[bad]=0.;

    #print "min/max/mean/median of f_i"
    #print np.min(f_i), np.max(f_i), np.mean(f_i), np.median(f_i)

    if (color_scheme_sdss==1):
        q=9.; alpha=0.3;
        f_i = np.arcsinh( alpha * q * (i/minnorm) ) / q; 
        wt=f_i/i; r*=wt; g*=wt; b*=wt;
    if (color_scheme_nasa==1):
        r[r>0] = np.log(r[r>0]/minnorm) / np.log(maxnorm/minnorm);
        g[g>0] = np.log(g[g>0]/minnorm) / np.log(maxnorm/minnorm);
        b[b>0] = np.log(b[b>0]/minnorm) / np.log(maxnorm/minnorm);

    #print "min/max/mean/median of r"
    #print np.min(r), np.max(r), np.mean(r), np.median(r)


    ## rescale to saturation limit
    if False:
        bad=(i<=0.); maxrgb=0.;
        if (checklen(i[bad])>0): r[bad]=0.; g[bad]=0.; b[bad]=0.;
        f_saturated=0.0004  ## fraction of pixels that should be saturated
        f_saturated=0.0001  ## fraction of pixels that should be saturated
        x0=int_round( f_saturated * (np.float(checklen(r)) - 1.) ); 

        for rgb_v in [r,g,b]: 
            rgbm=single_vec_sorted(rgb_v,reverse=True); 
            if(rgbm[x0]>maxrgb): maxrgb=rgbm[x0]

    #this is an arbitrary scaling to saturate the image...
#    if (maxrgb > 1.): r/=maxrgb; g/=maxrgb; b/=maxrgb;

    ## rescale to 256-colors to clip the extremes (rescales back to 0-1):
    max_c=255; min_c=2;


    r=clip_256(r,max=max_c,min=min_c);
    g=clip_256(g,max=max_c,min=min_c);
    b=clip_256(b,max=max_c,min=min_c);
    
    #print np.min(r), np.max(r), np.mean(r), np.median(r)


    image24=np.zeros((checklen(r[:,0]),checklen(r[0,:]),3),dtype='f');
    image24[:,:,0]=r; image24[:,:,1]=g; image24[:,:,2]=b; 

    ## ok have r, g, b -- really just three re-scaled maps. no reason they 
    ##   have to map to r, g, b colors: use the filter set given to map them ::
    if True:

       image24_new = 0.*image24
       viscolors.load_my_custom_color_tables();
       for i in [0,1,2]:
         im=image24[:,:,i]
         if filterset[i]=='r': image24_new[:,:,0] = im
         if filterset[i]=='g': image24_new[:,:,1] = im
         if filterset[i]=='b': image24_new[:,:,2] = im
         if (filterset[i] != 'r') & (filterset[i] != 'g') & (filterset[i] != 'b'):
            my_cmap = matplotlib.cm.get_cmap(filterset[i])
            rgb_im = my_cmap(im)
            image24_new += rgb_im[:,:,0:3] ## dropping the alpha channel here!
       image24 = image24_new

    return image24, cmap_m; ## return both processed image and massmap
