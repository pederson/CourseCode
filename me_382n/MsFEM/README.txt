***********README for multiscale FEM code*****************

To run the code:

1) Navigate to the 'MsFEM/Code' folder in matlab
2) The file 'multiscale.m' contains all the user-defined values
3) Set user-defined variables as you wish (see below)
4) Run multiscale.m by clicking "Run" in matlab, or by typing "multiscale"
   on the matlab command prompt
5) various plots will be displayed, along with some relevant data to the 
   matlab terminal


******* User-defined variables guide ********
The user-defined variables are mostly found at the top of the script
All variables are in SI units, unless noted otherwise
They are defined as follows:

    conductivity_*      = conductivity of material
    density_*           = density of material
    specific_heat_*     = specific heat of material
    x_min               = left-most x value on the domain
    x_max               = right-most x value on the domain
    spacing_coarse      = node spacing for the coarse grid
    spacing_fine        = node_spacing for the fine grid
    T_left              = Temperature at left boundary
    T_right             = Temperature at right boundary
    dt                  = time step
    t_end               = terminating time
    T_0                 = vector of initial temperatures on the coarse grid


*********** Other goodies ************
throughout the script, you can change certain things by uncommenting
and commenting out the appropriate lines. Most notably, in the last section
titled "Unsteady Problem" you can change the time stepping method by 
uncommenting certain lines.
