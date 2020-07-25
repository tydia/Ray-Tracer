# Ray Tracer

## What is this project?
A ray tracer that generate images from .scene files

## Usage
Make sure you have the `pic` library compiled. Then in RayTracer folder, type `make` to compile. To run the program, type `./RayTracer <scenefile> [jpeg_name] [AA_ON/AA_OFF]`. Then last argument is for turn on/off antialiasing. 

## Details about Functionalities:
- Ray tracing triangles
- Ray tracing spheres
- Triangle Phong Shading
- Sphere Phong Shading
- Shadow Rays
- Anti-aliasing: done by supersampling rays with a factor of three. The final color of a output pixel is taken as the average of the 9 of its coresponding pixels in the supersampled ray-traced image. To use anti-aliasing, you should have 4th input argument typed as "AA_ON". i.e. `./RayTracer test.scene test.jpg AA_ON`

## Demos
See `demos` directory.
