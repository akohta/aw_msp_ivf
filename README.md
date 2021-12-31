# aw_msp_ivf
This is the analysis program of sound wave scattering by a spherical object in a inviscid fluid. 
This program can analyze multiple scattering between spherical objects, and can analyze radiation force acting on each spherical object. 
The sound pressure analysis program "multi_aw" is used to analyze incident field.
Gnu Scientific Library and libpng are required.

## Usage of example code  
1. type 'make' command to compile.  
   The executable aw_msp_solver, example1.out, example2.out, example3.out, example4.out are created.
   The aw_msp_solver is the solver of field coefficients. 
   The example1.out is the executable of source code example1.c, it shows a simplest example using "aw_msp_ivf". 
   The example2.out is the executable of source code example2.c, it shows a example of sound pressure intensity analysis. 
   The example3.out is the executable of source code example3.c, it shows a example of outputting the instantaneous value of the sound pressure as an image.
   The example4.out is the executable of source code example4.c, it shows a example of far-field intensity analysis.  

2. type './aw_msp_solver' with an argument of output datafile name.  
   For './aw_msp_solver ex.dat'. 
   The beam datafile "mfb.txt" (focused beam is defined) and the sphere datafile "msphr.txt" (radius 5mm silicon oil droplet is defined) are used. 
   This executable calcluates field coefficients, outputs them to binary file with the specified file name. 
   This program searches for a sphere datafile in current directory using the default datafile name "msphr.txt". 
   As a simple representation of the analysis model, the nodes used for the surface integral are output as point cloud data. 
   In this case, the file "ex.particles" is output, and the visualization result is "particle.png" ( using gnuplot script gscript_particles.plt ).

3. type './example1.out' with an argument of datafile name.  
   For example, './example1.out ex.dat'. 
   This executable calculates sound pressure, particle velocity, radiation force and torque. 
   It also outputs the verification results using boundary conditions and so on.
   This is the simplest example using this code.
   
4. type './example2.out' with an argument of datafile name.  
   For example, './example2.out ex.dat'. 
   This executable calculates sound pressure intensity distributions, outputs them to text files. 
   The I_example2.png is the visualization result of sound pressure intensity distributions, created by using Gnuplot script gscript_example2.plt
   (converted eps to png by using ImageMagick).
   
5. type './example3.out' with an argument of datafile name.  
   For example, './example3.out ex.dat'. 
   This executable calculates instantaneous value of the sound pressure, outputs them to png image files. 
   The image files are output to the folder which has a name adding "images" to the datafile name specified in the argument (file-extension is excluded). 
   Each image file has a name that indicates the cross section, field component, and number of time steps (ex. xz_Ex_014.png). 
   The color bar is output as color_bar.png in the same folder.
   The range of color bar in each cross section is output to the info.txt file (ex. xy_info.txt for z=0 plane). 
   The xz_p.gif and xy_p.gif are animated gifs that concatenate the png files created by using the shell script gif_animation.sh.  
   
6. type './example4.out' with an argument of datafile name.  
   For example, './example4.out ex.dat'. 
   This executable calculates scattered field (sound puressure) intensity distributions in far-field and outputs them to text files. 
   The I_example4.png is the visualization result of the intensity distributions, created by gnuplot script gscript_example4.plt.
   The I_example4_3d.png is the visualization result of the intensity distributions on a spherical surface, created by gnuplot script gscript_example4_3d.plt.  

Please see msp_src/aw_msp_ivf.h for detail of functions, maw_src/multi_aw.h for detail of incident field functions. 
The aw_msp_solver, example2.out, example3.out, example4.out are parallelized using OpenMP. 
The number of threads is controlled by the environment variable OMP_NUM_THREADS. 
The analysis example of multiple spheres is in the folder analysis_sample1.
The analysis example of multiple focused beams is in the folder analysis_sample2.

![point cloud data](particles.png "nodes for surface integral (particles.png)")  
![intensity distributions](I_example2.png "intensity distributions (I_example2.png)")  
![xz_p.gif](xz_p.gif "instantaneous value of the p on y=0 plane (xz_p.gif)")![xy_p.gif](xy_p.gif "instantaneous value of the p on z=0 plane (xy_p.gif)")  
![far-field 3d](I_example4_3d.png "far-field intensity distribution (I_example4_3d.png)")  

## Analysis sample of radiation force  

The codes in the folder radiation_force1 ~ radiation_force3 are the analysis samples for radiation force. 
To use it, type 'make' and './radiation_force_force*.out' in that directory. 
The radiation_force1 analyzes the force acting on the sphere irradiated by a plane wave, 
the radiation_force2 analyzes the force acting on the sphere irradiated by a Bessel beam. 
The poisition of the sphere is on the acoustic axis.
The radiation_force3 analyzes the radiation force vector on each axis cross section. 
In this example, the incident field is the two focused beams crossing at origin (the same as analysis_sample2).
The figure below was created using the gnuplot script in that folder. 

* irradiated by the plane wave
![radiation force 1](radiation_force1/radiation_force1.png "result of radiation_force1 (radiation_force1/radiation_foece1.png)")  
  
* irradiated by the Bessel beam
![radiation force 2](radiation_force2/radiation_force2.png "result of radiation_force2 (radiation_force2/radiation_foece2.png)")   
  
* irradiated by the two focused beams crossing at origin
![radiation force 3](radiation_force3/radiation_force3.png "result of radiation_force3 (radiation_force3/radiation_foece3.png)")  


## References

1. GNU Scientific Library [GSL](https://www.gnu.org/software/gsl/)  
2. The command-line driven graphing utility [gnuplot](http://www.gnuplot.info/)  
3. The utilities for manipulating images [ImageMagick](https://imagemagick.org/)  
4. The official PNG reference library [libpng](http://www.libpng.org/pub/png/libpng.html)  
5. The soud pressure analysis program [multi_aw](https://github.com/akohta/multi_aw/)  
6. Andrade, José Henrique Araújo Lopes de. "Acoustic radiation force and torque on suspended objects in an inviscid fluid." (2014).
   (Please note that there are some misprints in the formula.)
