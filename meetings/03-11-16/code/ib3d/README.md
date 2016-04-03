#IBM3D - A 3D Immersed Boundary Method implementation in Matlab


## Build instructions
--
Run the following command on the Matlab GUI to generate vtk data in the folder `./op_sphere` (or any other folder with appropriate read/write/execute permissions):

```
> ib3d
```

After the simulation completes, we suggest to render the .vtk files via POV-Ray by running the python script in terminal:

```
> python render.py
```
which will create png files in the `./rendered`  folder.
