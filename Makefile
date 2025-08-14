##############################################################################
#    Makefile for CitcomSVE-3.0 with all *.c and *.h in one directory
#    To compile to generate executable CitcomSFull, type
#    $make CitcomSFull
##############################################################################
##############################################################################
#	Individual Machine Variations
##############################################################################

COMPRESS= /usr/bin/compress
LIB_PATH=
LIB_LIST=
LIB= $(LIB_PATH) $(LIB_LIST) -lm


###################################################################
#	Operating System Variables
###################################################################

CC=mpicc -Wno-deprecated-non-prototype -Wno-incompatible-pointer-types
#CC=mpicc
CPP=

CEXT=c
FEXT=F   # which implies further action of cpp
OBJEXT=o
FOBJEXT=o
OBJFLAG=-c
FOBJFLAG=-c


###################################
# Choose your machine from here.
###################################

LinuxFLAGS=
LinuxLDFLAGS=
LinuxOPTIM=-O3

####################################

FLAGS= $(LinuxFLAGS) -DCOMPRESS_BINARY=\"$(COMPRESS)\"
LDFLAGS= $(LinuxLDFLAGS)
OPTIM= $(LinuxOPTIM)

####################################################################
#	C files .....
####################################################################

CFILES= Advection_diffusion.c\
	BC_util.c\
	Checkpoints.c\
	Citcom.c\
	Citcom_init.c\
	CitcomSFull.c\
	Composition_related.c\
	Construct_arrays.c\
	Convection.c\
	Determine_net_rotation.c\
	Drive_solvers.c\
	Element_calculations.c\
	Free_surfaces.c\
	Full_boundary_conditions.c\
	Full_geometry_cartesian.c\
	Full_lith_age_read_files.c\
	Full_parallel_related.c\
	Full_read_input_from_files.c\
	Full_solver.c\
	Full_sphere_related.c\
	Full_tracer_advection.c\
	Full_version_dependent.c\
	General_matrix_functions.c\
	Ggrd_handling.c\
	Global_operations.c\
	Initial_temperature.c\
	Instructions.c\
	Interuption.c\
	Lith_age.c\
	Material_properties.c\
	Mineral_physics_models.c\
	Nodal_mesh.c\
	Output.c\
	Output_gzdir.c\
	Output_h5.c\
	Output_vtk.c\
	Pan_problem_misc_functions.c\
	Parallel_util.c\
	Parsing.c\
	Phase_change.c\
	Problem_related.c\
	Process_buoyancy.c\
	Regional_boundary_conditions.c\
	Regional_geometry_cartesian.c\
	Regional_lith_age_read_files.c\
	Regional_parallel_related.c\
	Regional_read_input_from_files.c\
	Regional_solver.c\
	Regional_sphere_related.c\
	Regional_tracer_advection.c\
	Regional_version_dependent.c\
	Sea_Level_Equation.c\
	Shape_functions.c\
	Size_does_matter.c\
	Solve_equations_motions.c\
	Solver_conj_grad.c\
	Solver_multigrid.c\
	Sphere_harmonics.c\
	Sphere_util.c\
	Topo_gravity.c\
	Tracer_setup.c\
	Visco_elastic.c\
	Viscosity_structures.c

HEADER = advection_diffusion.h\
	 advection.h\
	 checkpoints.h\
	 citcom_init.h\
	 composition_related.h\
	 convection_variables.h\
	 drive_solvers.h\
	 element_definitions.h\
	 global_defs.h\
	 hdf5_related.h\
	 initial_temperature.h\
	 interuption.h\
	 lith_age.h\
	 material_properties.h\
	 output.h\
	 output_h5.h\
	 parallel_related.h\
	 parsing.h\
	 phase_change.h\
	 solver.h\
	 sphere_communication.h\
	 tracer_defs.h\
	 ve_loading.h\
	 viscosity_descriptions.h


OBJFILES=$(CFILES:.c=.o)
#
#global_defs.h: viscosity_descriptions.h advection.h\
#		Convection_variables.h
		
# The following entries can probably be automated from $CFILES etc

Advection_diffusion.o: $(HEADER) Advection_diffusion.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Advection_diffusion.c
#
BC_util.o: $(HEADER) BC_util.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  BC_util.c
#	
Checkpoints.o: $(HEADER) Checkpoints.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Checkpoints.c
#	
Citcom.o: $(HEADER) Citcom.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Citcom.c
#	
Citcom_init.o: $(HEADER) Citcom_init.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Citcom_init.c
#	
CitcomSFull.o: $(HEADER) CitcomSFull.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  CitcomSFull.c
#	
Composition_related.o: $(HEADER) Composition_related.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Composition_related.c
#	
Construct_arrays.o: $(HEADER) Construct_arrays.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Construct_arrays.c
#	
Convection.o: $(HEADER) Convection.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Convection.c
#	
Determine_net_rotation.o: $(HEADER) Determine_net_rotation.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Determine_net_rotation.c
#	
Drive_solvers.o: $(HEADER) Drive_solvers.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Drive_solvers.c
#	
Element_calculations.o: $(HEADER) Element_calculations.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Element_calculations.c
	
Free_surfaces.o: $(HEADER) Free_surfaces.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Free_surfaces.c
	
Full_boundary_conditions.o: $(HEADER) Full_boundary_conditions.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Full_boundary_conditions.c
	
Full_geometry_cartesian.o: $(HEADER) Full_geometry_cartesian.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Full_geometry_cartesian.c
	
Full_lith_age_read_files.o: $(HEADER) Full_lith_age_read_files.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Full_lith_age_read_files.c
	
Full_parallel_related.o: $(HEADER) Full_parallel_related.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Full_parallel_related.c
	
Full_read_input_from_files.o: $(HEADER) Full_read_input_from_files.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Full_read_input_from_files.c
	
Full_solver.o: $(HEADER) Full_solver.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Full_solver.c
	
Full_sphere_related.o: $(HEADER) Full_sphere_related.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Full_sphere_related.c
	
Full_tracer_advection.o: $(HEADER) Full_tracer_advection.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Full_tracer_advection.c
	
Full_version_dependent.o: $(HEADER) Full_version_dependent.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Full_version_dependent.c
	
General_matrix_functions.o: $(HEADER) General_matrix_functions.c
	$(CC) $(OPTIM) $(FLAGS) $(OBJFLAG)   General_matrix_functions.c
		
Ggrd_handling.o: $(HEADER) Ggrd_handling.c
	$(CC) $(OPTIM) $(FLAGS) $(OBJFLAG)   Ggrd_handling.c
		
Global_operations.o: $(HEADER) Global_operations.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Global_operations.c
	
Initial_temperature.o: $(HEADER) Initial_temperature.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Initial_temperature.c
	
Instructions.o: $(HEADER) Instructions.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Instructions.c
	
Interuption.o: $(HEADER) Interuption.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Interuption.c
	
Lith_age.o: $(HEADER) Lith_age.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Lith_age.c
	
Material_properties.o: $(HEADER) Material_properties.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Material_properties.c
	
Mineral_physics_models.o: $(HEADER) Mineral_physics_models.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG) Mineral_physics_models.c
	
Nodal_mesh.o: $(HEADER) Nodal_mesh.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Nodal_mesh.c

Output.o: $(HEADER) Output.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Output.c
	
Output_gzdir.o: $(HEADER) Output_gzdir.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Output_gzdir.c
	
Output_h5.o: $(HEADER) Output_h5.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Output_h5.c
	
Output_vtk.o: $(HEADER) Output_vtk.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Output_vtk.c
	
Pan_problem_misc_functions.o: $(HEADER)  Pan_problem_misc_functions.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Pan_problem_misc_functions.c
		
Parallel_util.o: $(HEADER) Parallel_util.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Parallel_util.c
	
Parsing.o: $(HEADER) Parsing.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Parsing.c

Phase_change.o: $(HEADER) Phase_change.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Phase_change.c

Problem_related.o: $(HEADER) Problem_related.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Problem_related.c

Process_buoyancy.o: $(HEADER) Process_buoyancy.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Process_buoyancy.c

Regional_boundary_conditions.o: $(HEADER) Regional_boundary_conditions.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Regional_boundary_conditions.c

Regional_geometry_cartesian.o: $(HEADER) Regional_geometry_cartesian.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Regional_geometry_cartesian.c

Regional_lith_age_read_files.o: $(HEADER) Regional_lith_age_read_files.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Regional_lith_age_read_files.c

Regional_parallel_related.o: $(HEADER) Regional_parallel_related.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Regional_parallel_related.c

Regional_read_input_from_files.o: $(HEADER) Regional_read_input_from_files.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Regional_read_input_from_files.c

Regional_solver.o: $(HEADER) Regional_solver.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Regional_solver.c

Regional_sphere_related.o: $(HEADER) Regional_sphere_related.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Regional_sphere_related.c

Regional_tracer_advection.o: $(HEADER) Regional_tracer_advection.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Regional_tracer_advection.c

Regional_version_dependent.o: $(HEADER) Regional_version_dependent.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Regional_version_dependent.c

Sea_Level_Equation.o: $(HEADER) Sea_Level_Equation.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Sea_Level_Equation.c
	
Shape_functions.o: $(HEADER) Shape_functions.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Shape_functions.c
	
Size_does_matter.o: $(HEADER) Size_does_matter.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Size_does_matter.c
	
Solve_equations_motions.o: $(HEADER) Solve_equations_motions.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG) Solve_equations_motions.c
	
Solver_conj_grad.o: $(HEADER) Solver_conj_grad.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Solver_conj_grad.c
	
Solver_multigrid.o: $(HEADER) Solver_multigrid.c
	$(CC) $(OPTIM) $(FLAGS) $(OBJFLAG)   Solver_multigrid.c

Sphere_harmonics.o: $(HEADER) Sphere_harmonics.c
	$(CC) $(OPTIM) $(FLAGS) $(OBJFLAG)   Sphere_harmonics.c

Sphere_util.o: $(HEADER) Sphere_util.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Sphere_util.c
	
Topo_gravity.o: $(HEADER) Topo_gravity.c 
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Topo_gravity.c
	
Tracer_setup.o: $(HEADER) Tracer_setup.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Tracer_setup.c
	
Visco_elastic.o: $(HEADER) Visco_elastic.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Visco_elastic.c
	
Viscosity_structures.o: $(HEADER) Viscosity_structures.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Viscosity_structures.c
	

smaller: 
	compress $(CFILES)

larger:
	uncompress $(CFILES)


CitcomSFull: $(OBJFILES)  $(HEADER) Makefile
	$(CC) $(OPTIM) $(FLAGS) $(LDFLAGS)  -o CitcomSFull $(OBJFILES)  $(FFTLIB)  $(LIB)
