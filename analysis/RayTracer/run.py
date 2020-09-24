# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 18:26:43 2020

@author: Jonah Rose
"""

#import sys
#sys.path.insert(0, "C:\\Users\\astro\\OneDrive\\Desktop\\")
#print(sys.path)
#stop

from RoseArepo.analysis.RayTracer.Sphere import Sphere
from RoseArepo.analysis.RayTracer.Screen import Screen
from RoseArepo.analysis.RayTracer.Light import Light
from RoseArepo.analysis.RayTracer.Camera import Camera
from RoseArepo.analysis.RayTracer.Ray import Ray
from RoseArepo.analysis.RayTracer.Box import Box

import matplotlib.pyplot as plt
import numpy as np

def main():
    
    pixels = (200, 200)
    
    objs = []
    #objs.append(Sphere([100,100,100], 5, color=[255,255,255]))
    objs.append(Box([0,0,0], 100, 100, 10, [255,255,255]))
    
    cam = Camera()
    
    lights = []
    lights.append(Light([50,50,50], [1,1,1]))
    #lights.append(Light([200,200,200], [0,0,1]))
    #lights.append(Light([0,200,0], [0,1,0]))
    #lights.append(Light([0,0,200], [1,1,0]))
    
    screen = Screen(pixels, objs, lights, cam) #on y-z plane
    
    for i in range(pixels[0]):
        if i % 100 == 0:
            print(i)
        for j in range(pixels[1]):
            ray = Ray(np.array([0,j,i]), np.array([1000,j, i])) #parallel ray
            hit_result = screen.hit(ray)
            if hit_result is not None: #hit_result is a coordinate
                obj_color = []
                lig_color = []
                for l in screen.lights:
                    sray = Ray(hit_result, l.pos, shadow=True)
                    sray_hit_result = screen.hit(sray)
                    if sray_hit_result is None:
                        obj_color.append(screen.objs[ray.obj_idx].color)
                        lig_color.append(l.color)
                ocolor = np.average(obj_color, axis=0)
                lcolor = np.average(lig_color, axis=0)
                if type(ocolor) != type(np.zeros(0)):
                    screen.pixels[i,j] = np.zeros(3, dtype=int)
                else:
                    screen.pixels[i,j] = 1.0 * ocolor * lcolor

    fig, ax = plt.subplots()
    ax.imshow(screen.pixels)
    fig.savefig("C:\\Users\\astro\\OneDrive\\Desktop\\test.pdf")
    return

if __name__=="__main__":
    main()