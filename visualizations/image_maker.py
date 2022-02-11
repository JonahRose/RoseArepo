import numpy as np
import ctypes




def fcor(x):
    return np.array(x,dtype='f',ndmin=1)
def vfloat(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));
def cfloat(x):
    return ctypes.c_float(x);
def checklen(x):
    return len(np.array(x,ndmin=1));



def make_surface_density_image(parttype=0, xpixels=256, fov=10.0, axisratio=1.0, ma=None, mi=None ):

    


    ### To be re-written by Jonah using RoseArepo ###
    data = rs.snapshot_data( snapshot, parttype, **kwargs ) # possible tags: cosmo, fof_num, sub_num
    data.snapdir = snapdir
    data.snapnum = snapnum

    print "... centering data ..."
    data = idm.center_data( data, **kwargs )


    phi, theta = util.determine_rotation_angle( data, **kwargs )
    print "... rotating data ..."
    data = idm.rotate_data( data, phi, theta )

    print "... image bounds ..."
    data = ip.determine_image_bounds( data, **kwargs )

    print "... image stretch ..."
    data = ip.determine_image_stretch( data, **kwargs )

    print "... clipping ..."
    data = idm.clip_particles( data,  data.xr, data.yr, data.zr,  **kwargs )

    print "... rescaling ..."
    data = idm.rescale_hsml(   data, **kwargs)

    ### At this point, have data object with required partice data.




#    massmap,image = \
#                cmakepic.simple_makepic(data.pos[:,0],data.pos[:,1],weights=getattr(data,weight_type),hsml=data.h,\
#                xrange=data.xr,yrange=data.yr,
#                set_maxden=data.set_maxden , set_dynrng=data.set_dynrng,
#                pixels=pixels )


  ## load the routine we need
    exec_call=util.return_python_routines_cdir()+'/SmoothedProjPFH/allnsmooth.so'
    smooth_routine=ctypes.cdll[exec_call];
    ## make sure values to be passed are in the right format

   ### Jonah would need to edit, to ref data fields loaded above (e.g., weight = mass)

    N=checklen(x); x=fcor(x); y=fcor(y); M1=fcor(weight); M2=fcor(weight2); M3=fcor(weight3); H=fcor(hsml)

    # could be packaged into a sep class dealing with image properties
    xpixels=np.int(xpixels); ypixels=np.int(xpixels * axisratio)
    xmin = -fov; xmax = fov;
    ymin = -fov; ymax = fov

    ## check for whether the optional extra weights are set
    NM=1;
    if(checklen(M2)==checklen(M1)):
        NM=2;
        if(checklen(M3)==checklen(M1)):
            NM=3;
        else:
            M3=np.copy(M1);
    else:
        M2=np.copy(M1);
        M3=np.copy(M1);
    ## initialize the output vector to recieve the results
    XYpix=xpixels*ypixels; MAP=ctypes.c_float*XYpix; MAP1=MAP(); MAP2=MAP(); MAP3=MAP();
    ## main call to the imaging routine
    smooth_routine.project_and_smooth( \
        ctypes.c_int(N), \			# number of particles
        vfloat(x), vfloat(y), 			# x/y pos of particles
        vfloat(H), \				# hsml of particle 
        ctypes.c_int(NM), \			# Number of map dimensions
        vfloat(M1), vfloat(M2), vfloat(M3), \	# The 3 possible map weights
        cfloat(xmin), cfloat(xmax), 		# x min/max for the image
	cfloat(ymin), cfloat(ymax), \		# y min/max for the image
        ctypes.c_int(xpixels), ctypes.c_int(ypixels), \ # just the number of pixels in each direction
        ctypes.byref(MAP1), ctypes.byref(MAP2), ctypes.byref(MAP3) );	# "empty" arrays to hold the output

    ## now put the output arrays into a useful format
    MassMap1=np.ctypeslib.as_array(MAP1).reshape([xpixels,ypixels]);
    MassMap2=np.ctypeslib.as_array(MAP2).reshape([xpixels,ypixels]);
    MassMap3=np.ctypeslib.as_array(MAP3).reshape([xpixels,ypixels]);

    # Here, the MassMaps should contain actual mass maps with units of [M1] / ( [dimensions of pos]^2 )
    # For example, if mass is in units of 10^{10} M_solar and pos is in kpc:  10^{10} M_solar / kpc ^2


    # set boundaries and do some clipping
    MassMap = np.copy(MassMap1);
    print "MassMap : max: ", np.max(MassMap), "   min: ", np.min(MassMap)

####### could put similar functionality back in later.
#    if (set_percent_maxden !=0) or (set_percent_minden !=0):
#        print 'percent max/min = ',set_percent_maxden,set_percent_minden
#        Msort=np.sort(np.ndarray.flatten(MassMap));
#        if (set_percent_maxden != 0): ma=Msort[set_percent_maxden*float(checklen(MassMap)-1)];
#        mi=ma/set_dynrng;
#        if (set_percent_minden != 0): mi=Msort[set_percent_minden*float(checklen(MassMap)-1)];
#        if (set_percent_maxden == 0): ma=mi*set_dynrng;
#        ok=(Msort > 0.) & (np.isnan(Msort)==False)
#        print 'min == ',mi
#        print 'max == ',ma
#        print 'element = ',set_percent_maxden*float(checklen(MassMap)-1)
#        print Msort
#        print Msort[set_percent_minden*float(checklen(MassMap)-1)]
#        if (mi <= 0) or (np.isnan(mi)): mi=np.min(Msort[ok]);
#        if (ma <= 0) or (np.isnan(ma)): ma=np.max(Msort[ok]);
#############

    if(ma is None):  ma = np.max( MassMap )
    if(mi is None):  mi = np.min( MassMap )

    print "Clipping the weighted maps at   ma= ", ma, " mi= ", mi
    MassMap[MassMap < mi]=mi; MassMap[MassMap > ma]=ma;


    # now set colors
    cols=255. # number of colors
    Pic = (np.log(MassMap/mi)/np.log(ma/mi)) * cols 	###  (cols-3.) + 2.  # old code; double check this works
    # Pic should manifestly go from 0 to 255

    #  This should be unnecessary b/c of clipping + image stretching
    ###	Pic[Pic > 255.]=255.; Pic[Pic < 1]=1.;


####    backgrd = np.where((Pic<=2) | (np.isnan(Pic)))
#   Setting background color
####    Pic[backgrd] = 0; # black

    if (NM>1):
        return MassMap1,MassMap2,MassMap3, Pic
    else:
        return MassMap, Pic



