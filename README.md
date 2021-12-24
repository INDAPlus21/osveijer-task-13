# osveijer Ray Tracer

Can render spheres and planes. Spheres can have metal and lambertian materials. Planes can have both of these but can also have a texture material which behaves like a lambertian material but also has a texture to control the albedo value.

Renders into ```.png``` format. Texture files are 50x50  ```.png```s.

```troll_render.png``` contains the rendered result of the current code setup. Shapes can be added, removed or moved by modifying the ```spheres``` and ```planes``` variables in ```main()```.