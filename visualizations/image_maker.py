import numpy as np
import ctypes

from readData.DataLoader import DataLoader
from analysis import analyze

class ImageData():

    def __init__(self, cat):

        #load in values from DataLoader
        self.pos = cat['Coordinates'] / cat.h * cat.time
        self.vels = cat['Velocities'] * np.sqrt(cat.time) 
        self.masses = cat['Masses'] * 1e10 / cat.h

        #check if hsml is stored in the snapshot.
        #if not, calculate it ourselves
        if 'SubfindHsml' in cat:
            self.hsml = cat['SubfindHsml'] / cat.h * cat.time
        else:
            self.hsml = analyze.get_particle_hsml(self.pos[:,0], self.pos[:,1], self.pos[:,2]) 

        #check if this in an image of a single galaxy/group
        if cat.sub_idx != -1:
            self.gal_pos = cat['SubhaloPos'] / cat.h * cat.time
        elif cat.fof_idx != -1:
            self.gal_pos = cat['GroupPos'] / cat.h * cat.time
        else:
            self.gal_pos = None

        #initialize fields to be filled in later
        self.phi = None
        self.theta = None

        self.xr = None
        self.yr = None
        self.zr = None

        self.dynrange = None
        self.maxden = None
    
        self.is_centered = False

        return

    #center data so that either the center of the halo is [0,0,0]
    #or the center of the particles that have been loaded 
    def center_data(self):
        if self.gal_pos is not None:
            self.pos -= self.gal_pos
        else:
            for i in range(3):
                self.pos[:,i] -= (np.max(self.pos[:,i]) + np.min(self.pos[:,i])) / 2
        self.is_centered = True
        return

    #not necessary so far, but if you want to obtain these you can 
    def get_phi_theta(self, is_face=True, is_edge=False):
        phi, theta = analyze.get_rotate_data(self.pos, self.vels, self.masses, edge_on=is_edge, face_on=is_face, get_pt=True)
        self.phi = phi
        self.theta = theta
        return phi, theta

    #This function will calculate phi, theta automatically if you want a face- or edge-on image
    #otherwise you must supply the phi and theta you want 
    def rotate_data(self, phi=0, theta=0, is_face=True, is_edge=False):
        self.pos, self.vels = analyze.get_rotate_data(self.pos, self.vels, self.masses, phi, theta, is_edge, is_face)
        return

    #if not set explicetly, set implicetly as maximum 1d length in x, y, or z
    def set_image_bounds(self, xrange=None, yrange=None, zrange=None):
        if xrange is not None:
            if yrange is None:
                yrange = xrange
            if zrange is None:
                zrange= xrange
        else:
            x_stretch = np.max(self.pos[:,0]) - np.min(self.pos[:,0])
            y_stretch = np.max(self.pos[:,1]) - np.min(self.pos[:,1])
            x_mid     = 0.5*( np.max( self.pos[:,0]) + np.min(self.pos[:,0]) )
            y_mid     = 0.5*( np.max( self.pos[:,1]) + np.min(self.pos[:,1]) )
            z_mid     = 0.5*( np.max( self.pos[:,2]) + np.min(self.pos[:,2]) )

            stretch = np.max( [x_stretch, y_stretch] )

            xrange = [ x_mid - stretch/2.0 , x_mid + stretch/2.0 ]
            yrange = [ y_mid - stretch/2.0 , y_mid + stretch/2.0 ]
            zrange = [ z_mid - stretch/2.0 , z_mid + stretch/2.0 ]       

        self.xr = xrange
        self.yr = yrange
        self.zr = zrange
            
        return

    #If not set explicitly, assumed as min max values
    def set_image_stretch(self, dynrange=None, maxden=None):
        self.dynrange = dynrange
        self.maxden = maxden
        return

    #cut out a box based on set_image_bounds
    def clip_particles(self):
        cut = ( self.pos[:,0] > self.xr[0]) & \
              ( self.pos[:,0] < self.xr[1]) & \
              ( self.pos[:,1] > self.yr[0]) & \
              ( self.pos[:,1] < self.yr[1]) & \
              ( self.pos[:,2] > self.zr[0]) & \
              ( self.pos[:,2] < self.zr[1]) 

        self.pos = self.pos[cut]
        self.vels = self.vels[cut]
        self.masses = self.masses[cut]
        self.hsml = self.hsml[cut]

        return 

    def rescale_hsml(self, h_rescale_factor=1):
        self.hsml *= h_rescale_factor
        return

    def fcor(self, x):
        return np.array(x,dtype='f',ndmin=1)

    def vfloat(self, x):
        return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    def cfloat(self, x):
        return ctypes.c_float(x)

    def checklen(self, x):
        return len(np.array(x,ndmin=1))

def make_surface_density_image(path, snap_num, parttype=0, xpixels=256, fov=10.0, axisratio=1.0, ma=None, mi=None, sub_idx=-1, fof_idx=-1):

    print("... loading data ...")

    keys = ['Coordinates', 'Velocities', 'Masses', 'SubhaloPos', 'GroupPos', 'SubfindHsml']
    cat = DataLoader(path, snap_num, part_types=parttype, keys=keys, sub_idx=sub_idx, fof_idx=fof_idx)

    data = ImageData(cat)

    ### To be re-written by Jonah using RoseArepo ###
    #data = rs.snapshot_data( snapshot, parttype, **kwargs ) # possible tags: cosmo, fof_num, sub_num
    #data.snapdir = snapdir
    #data.snapnum = snapnum

    print("... centering data ...")
    #data = idm.center_data( data, **kwargs )
    data.center_data()

    #phi, theta = util.determine_rotation_angle( data, **kwargs )
    print("... rotating data ...")
    #data = idm.rotate_data( data, phi, theta )
    data.rotate_data()

    print("... image bounds ...")
    #data = ip.determine_image_bounds( data, **kwargs )
    data.set_image_bounds()

    print("... image stretch ...")
    #data = ip.determine_image_stretch( data, **kwargs )
    data.set_image_stretch()

    print("... clipping ...")
    #data = idm.clip_particles( data,  data.xr, data.yr, data.zr,  **kwargs )
    data.clip_particles()

    print("... rescaling ...")
    #data = idm.rescale_hsml(   data, **kwargs)
    data.rescale_hsml()

    ### At this point, have data object with required partice data.


#    massmap,image = \
#                cmakepic.simple_makepic(data.pos[:,0],data.pos[:,1],weights=getattr(data,weight_type),hsml=data.h,\
#                xrange=data.xr,yrange=data.yr,
#                set_maxden=data.set_maxden , set_dynrng=data.set_dynrng,
#                pixels=pixels )


  ## load the routine we need
    #exec_call=util.return_python_routines_cdir()+'/SmoothedProjPFH/allnsmooth.so'
    exec_call = '/home/j.rose/torreylabtools_git/C/c_libraries/SmoothedProjPFH/allnsmooth.so'
    smooth_routine=ctypes.cdll[exec_call];
    ## make sure values to be passed are in the right format

    ### Jonah would need to edit, to ref data fields loaded above (e.g., weight = mass)

    #N=checklen(x); x=fcor(x); y=fcor(y); M1=fcor(weight); M2=fcor(weight2); M3=fcor(weight3); H=fcor(hsml)

    N=data.checklen(data.pos[:,0]) #not used
    
    x=data.fcor(data.pos[:,0])
    y=data.fcor(data.pos[:,1])
    
    weight = data.masses
    weight2 = weight3 = None

    M1=data.fcor(weight)
    M2=data.fcor(weight2)
    M3=data.fcor(weight3)
    
    H=data.fcor(data.hsml)

    # could be packaged into a sep class dealing with image properties
    xpixels=np.int(xpixels)
    ypixels=np.int(xpixels * axisratio)
    
    xmin = -fov
    xmax = fov
    ymin = -fov
    ymax = fov

    ## check for whether the optional extra weights are set
    NM=1
    if(data.checklen(M2)==data.checklen(M1)):
        NM=2
        if(data.checklen(M3)==data.checklen(M1)):
            NM=3
        else:
            M3=np.copy(M1)
    else:
        M2=np.copy(M1)
        M3=np.copy(M1)

    ## initialize the output vector to recieve the results
    XYpix=(xpixels*ypixels)
    MAP=ctypes.c_float*XYpix
    MAP1=MAP()
    MAP2=MAP()
    MAP3=MAP()

    ## main call to the imaging routine
    smooth_routine.project_and_smooth( \
        ctypes.c_int(N), 			# number of particles
        data.vfloat(x), data.vfloat(y), 			# x/y pos of particles
        data.vfloat(H), 				# hsml of particle 
        ctypes.c_int(NM), 			# Number of map dimensions
        data.vfloat(M1), data.vfloat(M2), data.vfloat(M3), 	# The 3 possible map weights
        data.cfloat(xmin), data.cfloat(xmax), 		# x min/max for the image
        data.cfloat(ymin), data.cfloat(ymax), 		# y min/max for the image
        ctypes.c_int(xpixels), ctypes.c_int(ypixels),  # just the number of pixels in each direction
        ctypes.byref(MAP1), ctypes.byref(MAP2), ctypes.byref(MAP3) )	# "empty" arrays to hold the output

    ## now put the output arrays into a useful format
    MassMap1=np.ctypeslib.as_array(MAP1).reshape([xpixels,ypixels])
    MassMap2=np.ctypeslib.as_array(MAP2).reshape([xpixels,ypixels])
    MassMap3=np.ctypeslib.as_array(MAP3).reshape([xpixels,ypixels])

    # Here, the MassMaps should contain actual mass maps with units of [M1] / ( [dimensions of pos]^2 )
    # For example, if mass is in units of 10^{10} M_solar and pos is in kpc:  10^{10} M_solar / kpc ^2

    # set boundaries and do some clipping
    MassMap = np.copy(MassMap1)
    print("MassMap : max: ", np.max(MassMap), "   min: ", np.min(MassMap))

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

    if ma is None:  
        ma = np.max( MassMap )
    if mi is None:  
        mi = np.min( MassMap )

    print("Clipping the weighted maps at   ma= ", ma, " mi= ", mi)
    MassMap[MassMap < mi]=mi
    MassMap[MassMap > ma]=ma

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



