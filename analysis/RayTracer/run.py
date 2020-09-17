# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 18:26:43 2020

@author: Jonah Rose
"""

#import sys
#sys.path.insert(0, "C:\\Users\\astro\\OneDrive\\Desktop\\RoseArepo")

from RoseArepo.analysis.RayTracer.Sphere import Sphere
from RoseArepo.analysis.RayTracer.Screen import Screen
from RoseArepo.analysis.RayTracer.Light import Light
from RoseArepo.analysis.RayTracer.Camera import Camera
from RoseArepo.analysis.RayTracer.Ray import Ray

import matplotlib.pyplot as plt
import numpy as np

def main():
    
    pixels = (10, 10)
    
    objs = []
    objs.append(Sphere([50,50,50], 5))
    
    cam = Camera()
    
    lights = []
    lights.append(Light([0,0,0], [1,1,1]))
    
    screen = Screen(pixels, objs, lights, cam)
    
    for i in range(pixels[0]):
        for j in range(pixels[1]):
            ray = Ray(screen.cam_pos, screen.obj_pos)
            if ray.hit:
                sray = Ray(ray.inter_pos, screen.lig_pos)
                if sray.hit:
                    sray.shadow = True
                    screen.pixels[i, j] = np.zeros(3)
                else:
                    obj_color = screen.objs[ray.obj_idx].color
                    lig_color = screen.lights[sray.lig_idx].color
                    screen.pixels[i,j] = obj_color * lig_color
    
    plt.figure()
    plt.imshow(screen.pixels)
    plt.savefig("test.pdf")
    return

if __name__=="__main__":
    main()