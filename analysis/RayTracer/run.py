# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 18:26:43 2020

@author: Jonah Rose
"""

#import sys
#sys.path.insert(0, "C:\\Users\\astro\\OneDrive\\Desktop\\")
#print(sys.path)

from RoseArepo.analysis.RayTracer.Sphere import Sphere
from RoseArepo.analysis.RayTracer.Screen import Screen
from RoseArepo.analysis.RayTracer.Light import Light
from RoseArepo.analysis.RayTracer.Camera import Camera
from RoseArepo.analysis.RayTracer.Ray import Ray
from RoseArepo.analysis.RayTracer.Plane import Plane

import matplotlib.pyplot as plt
import numpy as np

def main():
    
    pixels = (2000, 2000)
    
    objs = []
    objs.append(Sphere([500,500,500], 100, color=[255,0,0]))
    #objs.append(Plane([100,0,0], 100, 100, [255,255,255]))
    
    cam = Camera()
    
    lights = []
    lights.append(Light([0,0,0], [1,1,1]))
    l = Light([0,25,50], [1,1,1])
    
    screen = Screen(pixels, objs, lights, cam) #on y-z plane
    
    for i in range(pixels[0]):
        if i % 100 == 0:
            print(i)
        for j in range(pixels[1]):
            ray = Ray(np.array([0,j,i]), np.array([1000,j, i])) #parallel ray
            hit_result = screen.hit(ray)
            if hit_result is not None: #hit_result is a coordinate
                sray = Ray(hit_result, l.pos, shadow=True)
                sray_hit_result = screen.hit(sray)
                if sray_hit_result is not None:
                    screen.pixels[i, j] = np.zeros(3)
                else:
                    obj_color = screen.objs[ray.obj_idx].color
                    lig_color = l.color
                    screen.pixels[i,j] = obj_color * lig_color
            else:
                screen.pixels[i,j] = np.zeros(3)

    fig, ax = plt.subplots()
    ax.imshow(screen.pixels)
    fig.savefig("C:\\Users\\astro\\OneDrive\\Desktop\\test.pdf")
    return

if __name__=="__main__":
    main()