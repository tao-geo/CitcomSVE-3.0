/******************************************************************************
 *
 * CitcomSVE-3.0 is a finite element package that solves for dynamic
 * deformation of planetary mantle and crust in a spherical shell under either
 * surface or tidal loads. The mantle and crust may assume a fully compressible
 * or incompressible medium with 3-D viscoelastic structures. The package works
 * with parallel computers with the number of cores or CPUs greater than 12
 * (test runs with >6,000 cores). This version was built on CitcomS-3.1.1,
 * which is a code for planetary mantle convection of a purely viscous flow.
 * This heading was added to the files in CitcomSVE-3.0 that has been modified
 * extensively from or new to CitcomS. Following CitcomS, CitcomSVE has been
 * publicly available as an open source package at github since 2022.
 *
 *               Shijie Zhong, University of Colorado
 ******************************************************************************/



#include "global_defs.h"
#include "element_definitions.h"
#include "parsing.h"
#include <math.h>
#include <sys/types.h>

#include <stdio.h>

/* ========================================================================== 
 * set_free_surfaces(E)
 * Called only once in read_instructions (in Instructions.c)
 * ========================================================================== */
void set_free_surfaces(E)
    struct All_variables *E;
{

    void parallel_process_termination();
    void  get_surface_loads();
    int lev,i,node1,node2,m;

    /* mark the vertical d.o.f. for nodes on density boundaries */

    if (E->parallel.me_loc[3] == 0)                      { 
        for (lev=E->mesh.gridmax;lev>=E->mesh.gridmin;lev--)
        for (m=1;m<=E->sphere.caps_per_proc;m++)  
        for (i=1;i<=E->lmesh.NEL[lev];i++)    
        if ((i-1)%E->lmesh.ELZ[lev]==0)   {
            E->ELEMENT[lev][m][i] = E->ELEMENT[lev][m][i] | BOUNDARIES;
            E->ELEMENT[lev][m][i] = E->ELEMENT[lev][m][i] | BOUNDARY2;    /*CMB*/
            E->NODE[lev][m][E->IEN[lev][m][i].node[1]] = E->NODE[lev][m][E->IEN[lev][m][i].node[1]] | RESTORE;    /*CMB*/
            E->NODE[lev][m][E->IEN[lev][m][i].node[2]] = E->NODE[lev][m][E->IEN[lev][m][i].node[2]] | RESTORE;    /*CMB*/
            E->NODE[lev][m][E->IEN[lev][m][i].node[3]] = E->NODE[lev][m][E->IEN[lev][m][i].node[3]] | RESTORE;    /*CMB*/
            E->NODE[lev][m][E->IEN[lev][m][i].node[4]] = E->NODE[lev][m][E->IEN[lev][m][i].node[4]] | RESTORE;    /*CMB*/
        }
    }
    if (E->parallel.me_loc[3] == E->parallel.nprocz-1)   {
        for (lev=E->mesh.gridmax;lev>=E->mesh.gridmin;lev--)
        for (m=1;m<=E->sphere.caps_per_proc;m++)  
        for (i=1;i<=E->lmesh.NEL[lev];i++)    
        if (i%E->lmesh.ELZ[lev]==0)   {
            E->ELEMENT[lev][m][i] = E->ELEMENT[lev][m][i] | BOUNDARIES;
            E->ELEMENT[lev][m][i] = E->ELEMENT[lev][m][i] | BOUNDARY1;    /*Surface*/
            E->NODE[lev][m][E->IEN[lev][m][i].node[5]] = E->NODE[lev][m][E->IEN[lev][m][i].node[5]] | RESTORE;
            E->NODE[lev][m][E->IEN[lev][m][i].node[6]] = E->NODE[lev][m][E->IEN[lev][m][i].node[6]] | RESTORE;
            E->NODE[lev][m][E->IEN[lev][m][i].node[7]] = E->NODE[lev][m][E->IEN[lev][m][i].node[7]] | RESTORE;
            E->NODE[lev][m][E->IEN[lev][m][i].node[8]] = E->NODE[lev][m][E->IEN[lev][m][i].node[8]] | RESTORE;
        }
    }

	// get the load (single harmonic, iceModel, or ice caps):
    /* get the initial rho*g*h for nodes on density boundaries */
    get_surface_loads(E);

    return;
}


/* ======================================================= 
  before assembling a force vector, determine topography
  at each density interface and convert it into stress 
  drho*g*h
 =======================================================  */

void get_boundary_forces(struct All_variables *E, int count) 
//    struct All_variables *E;
//    int count;
{
    int m,ll,mm,i,j,p;
    double density_cmb,density_surf,total, length;
    double modified_plgndr_a(),con,t1;
    void remove_average();
    void parallel_process_termination();
    double total_surface_integral();
    double *TG[4],temp1,temp2;

  int e,n, node;
  double density_jump, stress_scaling;
  const int onedp = onedvpoints[E->mesh.nsd];
  const int vpts=vpoints[E->mesh.nsd];
  const int lev = E->mesh.levmax;

    void load_to_CM();
    void load_to_CM_grav();
    void load_to_CM_grav_comp();

    density_surf = 1.0;
    density_cmb = (E->data.density_below-E->data.density)/E->data.density;

        /* If this is the second (or later) self-gravity iteration for this
         * timestep, put the deformation from the previous iteration into
         * slice.surf/botm[2], then zero the CM for these loads (these loads
         * will be used in calculate_potential).
         * Then add them to the initial slice.load[0/2] to get the loads for 
         * this timestep, slice.load[1/3] (although I believe these are just
         * over-written before they're used).      */

    if (count==0)   {   
        /* slice.load[0] and slice.load[2] are the surf/cmb loads
         * from the end of the previous timestep, so return. */
      return;
      }
    else if (count)   {  //>>>> modify by tao:
		
	if(E->ve_data_cont.compressible == 0){ /* comment by tao: use surf/botm when incompress */
        // surface:
	        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)   {
	            for (m=1;m<=E->sphere.caps_per_proc;m++)
	            for (j=1;j<=E->lmesh.nsf;j++)    {
	                i = j*E->lmesh.noz;
	                E->slice_ve.surf[2][m][j] = E->U[m][E->id[m][i].doff[3]]
	                                               * E->ve_data_cont.surf_scaling
	                                        + E->slice_ve.dynamic_oceanload[m][j];
	            }
	            remove_average(E,E->slice_ve.surf[2],1);
	        }
	        // cmb:
	        if (E->parallel.me_loc[3]==0)   {
	            for (m=1;m<=E->sphere.caps_per_proc;m++)
	            for (j=1;j<=E->lmesh.nsf;j++)    {
	                i = (j-1)*E->lmesh.noz + 1;
	                E->slice_ve.botm[2][m][j] = E->U[m][E->id[m][i].doff[3]];
	                E->slice_ve.botm[2][m][j] *= E->ve_data_cont.botm_scaling;
	            }
	            remove_average(E,E->slice_ve.botm[2],0);
	        }
		}
	else {  // for compressible
		stress_scaling =  E->data.density*E->data.grav_acc*E->sphere.dradius/E->ve_data_cont.shear_mod;
			
		for (m=1;m<=E->sphere.caps_per_proc;m++){
		  for (node=1;node<=E->lmesh.nno;node++) {
                       j = (node-1)%E->lmesh.noz + E->lmesh.nzs;
		       E->slice_ve.all_stress[m][node] = E->rho_jump1D[E->mesh.levmax][j]*E->U[m][E->id[m][node].doff[3]]; /* i do not multiply this by grav acc because I only need mass */
		  }
			
		  /* add dynamic_oceanload to the Earth's surface */
	       if (E->ve_data_cont.SLE){
		  if ( E->parallel.me_loc[3]==E->parallel.nprocz-1 ) {
			for (j=1;j<=E->lmesh.nsf;j++)	 {
			  node = j*E->lmesh.noz;
			  E->slice_ve.all_stress[m][node] += E->slice_ve.dynamic_oceanload[m][j]/stress_scaling; 
 /* note that the rho used in dynamic_oceanload is actually rho_ice */
 /* note that dynamic_oceanload is non-dim stress, so convert it to non-dim mass - rho_water * u */			
		}
	  }
	  }   // end for SLE 
	}  // end for m
    }  // end for compressible

  if(E->ve_data_cont.compressible){
	  load_to_CM_grav_comp(E, 0, 0);
	}
  else{
	  load_to_CM_grav( E, E->slice_ve.surf[2], E->slice_ve.botm[2], 0 ); 
	}

        // todo_tao: 

    }// end if count
    	//<<<< end tao

    return;
}



/* ======================================================= 
   for each element el on boundaries, determine contributions
   from topographic deflection to elt_f
 =======================================================  */

void load_boundary_forces(struct All_variables *E, int el, double elt_f[24],int m)
//    struct All_variables *E;
 //   int el,m;
//    double elt_f[24];
{

    void get_global_1d_shape_fn_2();
//    void construct_c3x3matrix_el();
    void get_rtf_at_vpts();

    double radius2,temp,x[3][3],force[9],force_at_gs[9];
    int ic,nn[5],e, i, k, p, p2, lnode[5],node[9];
    double sinaa[9],ct[9],rtf[4][9];
    double div_w[4][9][9];
    double *gnxx,*gnda;

    const int lev = E->mesh.levmax;
    const int onedp = onedvpoints[E->mesh.nsd];
    const int vpts=vpoints[E->mesh.nsd];
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int ends=enodes[dims];

    struct Shape_function1 GM;
    struct Shape_function1_dA dGammax;

    static struct CC Cc;
    static struct CCX Ccx;

	//>>>> add by tao:
    double con;
    int a,j;
    void construct_c3x3matrix_el();
	
	//<<<< end tao

     // for incompressible
  if (!E->ve_data_cont.compressible)  {

    if (E->parallel.me_loc[3]==0 || E->parallel.me_loc[3]==E->parallel.nprocz-1)  {

        if (el%E->lmesh.elz==0 && E->parallel.me_loc[3]==E->parallel.nprocz-1)   {
            e = el/E->lmesh.elz;
            for(k=1;k<=onedvpoints[E->mesh.nsd];k++)   {
                lnode[k] = E->sien[m][e].node[k];
                force[k] = E->slice_ve.load[1][m][ lnode[k] ];
            }
            ic = 1;     /*      top  */
        }
        else if ((el-1)%E->lmesh.elz==0 && E->parallel.me_loc[3]==0)  {
            e = (el-1)/E->lmesh.elz+1;
            for(k=1;k<=onedvpoints[E->mesh.nsd];k++)   {
                lnode[k] = E->sien[m][e].node[k];
                force[k] = E->slice_ve.load[3][m][ lnode[k] ];
            }
            ic = 0;     /*      bottom */
        }
        else {
            return;
        }

        for(k=1;k<=onedvpoints[E->mesh.nsd];k++)   
            nn[k] = k+ic*4;

        for (i=1;i<=onedvpoints[E->mesh.nsd];i++)   {
            force_at_gs[i] = 0.0;
            for(k=1;k<=onedvpoints[E->mesh.nsd];k++)
                force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];
        }

        for(k=1;k<=onedp;k++) {
            for (i=1;i<=onedp;i++)  
                elt_f[3*nn[k]-1] += E->B_R[lev][m][nn[k]][(e-1)*onedp+i] * force_at_gs[i];
            }
    }
  }
  else if (E->ve_data_cont.compressible) {
    con = 0.0;
    if (E->ve_data_cont.SELFG)
       con = E->ve_data_cont.Rsg;

    for(k=1;k<=onedvpoints[E->mesh.nsd];k++)   {
      lnode[k] = E->ien[m][el].node[k+onedp];
      force[k] = E->slice_ve.all_load[m][lnode[k]];
    }
   
    radius2 = E->sx[m][3][lnode[1]]*E->sx[m][3][lnode[1]];

    for (i=1;i<=onedvpoints[E->mesh.nsd];i++)   {
      force_at_gs[i] = 0.0;
      for(k=1;k<=onedvpoints[E->mesh.nsd];k++)
        force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];
    }

    /* elt_f has been set to zero before this subroutine is called, grep assemble_forces */
    /* and only the r components of elt_f are changed */
    for(k=1;k<=onedp;k++)
      for (i=1;i<=onedp;i++)
        elt_f[3*(k+onedp)-1] += /*radius2* */E->B_RR[lev][m][k+onedp][(el-1)*onedp+i] * force_at_gs[i]; 
// modified by tao: before shijie used radius2*E->B_R, which will be same as E->B_RR.

    /* also pick the bottom surface if it is at CMB */
    if ((el-1)%E->lmesh.elz==0 && E->parallel.me_loc[3]==0)  {
      for(k=1;k<=onedvpoints[E->mesh.nsd];k++)   {
        lnode[k] = E->ien[m][el].node[k];
        force[k] = E->slice_ve.all_load[m][lnode[k]];
      }
      for (i=1;i<=onedvpoints[E->mesh.nsd];i++)   {
        force_at_gs[i] = 0.0;
        for(k=1;k<=onedvpoints[E->mesh.nsd];k++)
          force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];

      }
      for(k=1;k<=onedp;k++)
        for (i=1;i<=onedp;i++)
          elt_f[3*k-1] += E->B_RR[lev][m][k][(el-1)*onedp+i] * force_at_gs[i];
    }   // no need for multiplying radius^2 here for the CMB, as it was already, 

	// also add in the div_w*rho*phi*dV term 
	//  compute the rho*phi at gaussian points 

    for (a=1;a<=ends;a++) {
      node[a] = E->ien[m][el].node[a];
      force[a] = con*E->slice_ve.total_potential[m][node[a]];
    }
    for (i=1;i<=vpts;i++) {
      force_at_gs[i] = 0.0;
      for (k=1;k<=ends;k++)
        force_at_gs[i] += E->erho[lev][m][el]*force[k]*E->N.vpt[GNVINDEX(k,i)];
    }
	// compute div_w at each dimension 
    gnxx = E->gNX[m][el].vpt;
    gnda = E->gDA[m][el].vpt;
    get_rtf_at_vpts(E,m,E->mesh.levmax,el,rtf);

    if ((el-1)%E->lmesh.elz==0)
      construct_c3x3matrix_el(E,el,&Cc,&Ccx,E->mesh.levmax,m,0);

    for (j=1;j<=vpts;j++) {
	      sinaa[j] = sin(rtf[1][j]);
	      ct[j] = cos(rtf[1][j])/sinaa[j];
	      for (i=1;i<=ends;i++) {
	        for (k=1;k<=dims;k++){
	          div_w[k][i][j] = rtf[3][j]*
	                 (gnxx[GNVXINDEX(0,i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]
	                 +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(1,k,1,i,j)]
	                 +E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(3,k,i,j)])
	                          +rtf[3][j]*
	                (E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]*ct[j]
	                +E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(3,k,i,j)]
	                +(gnxx[GNVXINDEX(1,i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
	                 +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(2,k,2,i,j)])/sinaa[j])
	                          +gnxx[GNVXINDEX(2,i,j)]*Cc.vpt[BVINDEX(3,k,i,j)];
	        }
	      }
    }
	//      add to elf_f 
    for(a=1;a<=ends;a++) {
	      p = dims*(a-1);
	      for(j=1;j<=vpts;j++) {
	        elt_f[p  ] -= div_w[1][a][j]*force_at_gs[j]*gnda[j]*g_point[j].weight[dims-1];
	        elt_f[p+1] -= div_w[2][a][j]*force_at_gs[j]*gnda[j]*g_point[j].weight[dims-1];
	        elt_f[p+2] -= div_w[3][a][j]*force_at_gs[j]*gnda[j]*g_point[j].weight[dims-1];
	      }
	    }

  }


    return;
}


/* ========================================================================== 
 * get_surface_loads(E)
 * Called only once in set_free_surfaces (above).
 * Gets the ice load (single harmonic, iceModel, or ice caps).
 * ========================================================================== */
void get_surface_loads(struct All_variables *E) 
 //   struct All_variables *E;
{
    FILE *fp;
    int oldstep,m,ll,mm,p,i,j,n,k;
    double tau,total, length, stress_scaling;
    double dist,distance(),modified_plgndr_a(),con,t1,f1,r1;
    double x1,y1,z1,density_ice,r,t,f,ri,ro,icet,icef,icer,oldampt;
    int m1 = E->parallel.me;
    int node;

    double itemp[100], temp[100];

    char outfile[255];
    void sphere_expansion_output();
    void remove_average();
    void get_iceModel();  
    void parallel_process_termination();
    void apply_new_loads(); 

    void get_static_oceanload();
    void load_to_CM() ;
    void load_to_CM_grav() ;
    void load_to_CM_grav_comp();
    void shift_V_to_CM();
    void calculate_potential_incomp();
    double * temp_load_surf[NCS]; 
    double * temp_load_botm[NCS]; 
    double * temp_potential_surf[NCS];
    double * temp_potential_botm[NCS];

    E->ve_data_cont.potential_scaling = 4.0*M_PI*E->data.grav_const*E->data.density
        *E->sphere.dradius*E->sphere.dradius;
    E->ve_data_cont.rotation_rate = E->ve_data_cont.rotation_rate*2.0*M_PI/(24.0*3600.0);

        E->ve_data_cont.Rsg = 4.0*M_PI*E->data.grav_const*E->data.density
            *E->data.density*E->sphere.dradius*E->sphere.dradius
            /E->ve_data_cont.shear_mod;
        // the following are used to convert nondim height -> nondim stress
        E->ve_data_cont.botm_scaling = (E->data.density_below-E->data.density)
                              *E->data.grav_acc
                              *E->sphere.dradius/E->ve_data_cont.shear_mod;
        E->ve_data_cont.surf_scaling = E->data.density*E->data.grav_acc
                              *E->sphere.dradius /E->ve_data_cont.shear_mod;

    density_ice = 917.4;
    E->ve_data_cont.ice_stress_scale = E->ve_data_cont.surf_scaling*density_ice/E->data.density;
    stress_scaling = E->data.density*E->data.grav_acc*E->sphere.dradius/E->ve_data_cont.shear_mod;

    if (E->parallel.me==0)  {
        fprintf(E->fp,"scale %.5e %.5e %.5e %.5e\n",
                      E->ve_data_cont.botm_scaling,
                      E->ve_data_cont.surf_scaling,
                      E->ve_data_cont.ice_stress_scale,E->ve_data_cont.tau_in_years);
        fflush(E->fp); 
    }

// the following three lines are moved to Instructions.c
//    input_int("stages",&(E->ve_data_cont.stages),"1",m1);
//    input_int_vector("step",E->ve_data_cont.stages,(E->ve_data_cont.stages_step),m1);
//    input_double_vector("timestep",E->ve_data_cont.stages,(E->ve_data_cont.stages_time),m1);

    E->ve_data_cont.DIRECT=1;
    if (E->ve_data_cont.Heaviside==1)   { 
        E->monitor.elapsed_time=0.0;
        }
    else  {
        // set the current timstep to 1, and the current elapsed_time to dt:
        E->monitor.elapsed_time=E->ve_data_cont.stages_timestep[0];
        E->advection.timestep=E->ve_data_cont.stages_timestep[0];
        E->monitor.solution_cycles = 1;
    }

    // Set the initial loads:
    
    E->ve_data_cont.change_of_load = 1;  // used by apply_new_loads

    if (E->ve_data_cont.Heaviside==1)    {
        /* loading at a single harmonic */
        mm = E->convection.perturb_mm[0];
        ll = E->convection.perturb_ll[0];
        con = E->convection.perturb_mag[0];

        if (E->ve_data_cont.apply_potential==0)  {   // surface load
          for (m=1;m<=E->sphere.caps_per_proc;m++)
          for (j=1;j<=E->lmesh.nsf;j++)    {
            i = j*E->lmesh.noz;
            t1 = E->SX[E->mesh.levmax][m][1][i];
            f1 = E->SX[E->mesh.levmax][m][2][i];
            E->Xsurf[3][m][j] = con*modified_plgndr_a(ll,mm,t1)*cos(mm*f1);
            }
        } // end for surface load
        else if (E->ve_data_cont.apply_potential==1) {   // potential load
          if (E->ve_data_cont.compressible==0)  {  // incompressible
                  for (m=1;m<=E->sphere.caps_per_proc;m++)
                  for (j=1;j<=E->lmesh.nsf;j++)    {
                    i = j*E->lmesh.noz;
                    t1 = E->SX[E->mesh.levmax][m][1][i];
                    f1 = E->SX[E->mesh.levmax][m][2][i];
                    E->Xsurf[3][m][j] = 0.0;
                    E->init_potential[0][m][j] = con*(pow(E->sphere.ro,ll))*modified_plgndr_a(ll,mm,t1)*cos(mm*f1);
                    E->init_potential[1][m][j] = con*(pow(E->sphere.ri,ll))*modified_plgndr_a(ll,mm,t1)*cos(mm*f1);

              }
            }  // end for incompressible
          else if (E->ve_data_cont.compressible==1)  {  // compressible
                  for (m=1;m<=E->sphere.caps_per_proc;m++)
                  for (j=1;j<=E->lmesh.nsf;j++)    {
                    E->Xsurf[3][m][j] = 0.0;
                    for (k=1;k<=E->lmesh.noz;k++)    {
                      i = (j-1)*E->lmesh.noz + k;
                      t1 = E->SX[E->mesh.levmax][m][1][i];
                      f1 = E->SX[E->mesh.levmax][m][2][i];
                      r1  = E->SX[E->mesh.levmax][m][3][i];
                      E->slice_ve.init_total_potential[m][i] = con*(pow(r1,ll))*modified_plgndr_a(ll,mm,t1)*cos(mm*f1);
                }
            }
          }   // end for compressible
        }  // end for apply potential

        remove_average(E,E->Xsurf[3],1);

        sphere_expansion_output(E,1,E->Xsurf[3],
                                E->sphere.sphc[0],E->sphere.sphs[0],
                                E->monitor.solution_cycles,"init_surf");

        /* Get the potential of surface load, this is inside of Heaviside==1; Note this function only gets called once */
        if (E->ve_data_cont.apply_potential==0)  {
          
          for (m=1;m<=E->sphere.caps_per_proc;m++){
            temp_load_surf[m] = (double *)malloc((E->lmesh.nsf+1)*sizeof(double));
            temp_load_botm[m] = (double *)malloc((E->lmesh.nsf+1)*sizeof(double));
            temp_potential_surf[m] = (double *)malloc((E->lmesh.nsf+1)*sizeof(double));
            temp_potential_botm[m] = (double *)malloc((E->lmesh.nsf+1)*sizeof(double));
          }

          for (m=1;m<=E->sphere.caps_per_proc;m++){
            for (j=1;j<=E->lmesh.nsf;j++) {
              temp_load_surf[m][j] = E->Xsurf[3][m][j]*E->ve_data_cont.surf_scaling;  // m is always 1
              temp_load_botm[m][j] = 0.0;
            }
          }

          // here the load is Xsurf[3] * surf_scaling, non-dim stress (although in the calculate_potential, it is converted to non-dim mass)
          calculate_potential_incomp(E, temp_load_surf, temp_load_botm,
                               temp_potential_surf, temp_potential_botm, 0); // icon=0 so not considering the centrifugal potential

          // get sph expansion of potential
          sphere_expansion_output(E,1,temp_potential_surf,
            E->sphere.sphc_init_load,E->sphere.sphs_init_load,
            E->monitor.solution_cycles,"init_load_pttl");
          
          // done, free the memory
          for (m=1;m<=E->sphere.caps_per_proc;m++){
            free(temp_load_surf[m]);
            free(temp_load_botm[m]);
            free(temp_potential_surf[m]);
            free(temp_potential_botm[m]);
          }
        }    // end for apply surface load
        else if (E->ve_data_cont.apply_potential==1)  {
          for (m=1;m<=E->sphere.caps_per_proc;m++){
            temp_potential_surf[m] = (double *)malloc((E->lmesh.nsf+1)*sizeof(double));
          }
          if (E->ve_data_cont.compressible==0){
            for (m=1;m<=E->sphere.caps_per_proc;m++){
              for (j=1;j<=E->lmesh.nsf;j++){
                temp_potential_surf[m][j] = E->init_potential[0][m][j];
              }
            }
          }
          else if (E->ve_data_cont.compressible==1) {
            for (m=1;m<=E->sphere.caps_per_proc;m++)
              for (j=1;j<=E->lmesh.nsf;j++)    {
                node = j*E->lmesh.noz;
                temp_potential_surf[m][j] = E->slice_ve.init_total_potential[m][node];
              }
          }
          // get sph expansion of potential
          sphere_expansion_output(E,1,temp_potential_surf,
            E->sphere.sphc_init_load,E->sphere.sphs_init_load,
            E->monitor.solution_cycles,"init_tidal_load_pttl");
          
          // done, free the memory
          for (m=1;m<=E->sphere.caps_per_proc;m++){
            free(temp_potential_surf[m]);
          }
        }    //end for apply potential

        if(E->ve_data_cont.compressible) {  // add by tao: compressible version
            
          for (m=1;m<=E->sphere.caps_per_proc;m++)
            for (i=1;i<=E->lmesh.nno;i++)
              E->slice_ve.all_stress[m][i] = 0.0;

          if(E->parallel.me_loc[3]==E->parallel.nprocz-1)
            for (m=1;m<=E->sphere.caps_per_proc;m++)
              for (j=1;j<=E->lmesh.nsf;j++)    {
                node = j*E->lmesh.noz;
                E->slice_ve.all_stress[m][node] = E->Xsurf[3][m][j]; // comment by tao: here I follow geruo, all_stress is rho*u, or dynamic_oceanload / stress_scaling.
                    // so do not *E->ve_data_cont.surf_scaling.
              }

          for(i=1; i<=1; i++){
            load_to_CM_grav_comp(E, 1, 1);
            if(E->parallel.me==0)  {
                printf("debug deg-1 motion: now with shift_V_to_CM\n");
                printf("CM_incr: %e %e %e\n",E->ve_data_cont.CM_incr[0],E->ve_data_cont.CM_incr[1],E->ve_data_cont.CM_incr[2]);
            }
            // add by tao: follow geruo: note for ice loading case, shift_V_to_CM is not called.
            if(E->parallel.me==0)
              printf("shift_V_to_CM called\n");

            shift_V_to_CM(E);

          }

          for (m=1;m<=E->sphere.caps_per_proc;m++)
            for (i=1;i<=E->lmesh.nno;i++)
              E->slice_ve.total_load[m][i] = 0.0;
              
//	  fprintf(E->fp_out,"oooka \n");
          for (m=1;m<=E->sphere.caps_per_proc;m++)
            for (j=1;j<=E->lmesh.nno;j++)	 {
              node = j;
              E->slice_ve.total_load[m][node] = E->slice_ve.all_stress[m][node] * E->grav[E->mesh.levmax][m][node]*stress_scaling;
//	      fprintf(E->fp_out,"%d %g\n",j,E->slice_ve.total_load[m][node]);
            }
        }
        // incompressible
        else {
            for (m=1;m<=E->sphere.caps_per_proc;m++)
              for (j=1;j<=E->lmesh.nsf;j++)    {
                  E->slice_ve.surf[2][m][j] = E->Xsurf[3][m][j]*E->ve_data_cont.surf_scaling;
                  E->slice_ve.botm[2][m][j] = 0.0 ;
              }

            load_to_CM_grav( E, E->slice_ve.surf[2], E->slice_ve.botm[2], 1 );
            for (m=1;m<=E->sphere.caps_per_proc;m++)
            for (j=1;j<=E->lmesh.nsf;j++)	 {
              E->slice_ve.load[0][m][j] = E->slice_ve.surf[2][m][j] ;
              E->slice_ve.load[2][m][j] = E->slice_ve.botm[2][m][j] ;
            }
        }
          
    }  // end for Heaviside==1
      
    else if (E->ve_data_cont.Heaviside==2)  {   /* ice-model */
        // read 1st epoch of iceModel and (if appropriate) the oceanload:
        get_iceModel(E,E->slice_ve.iceload[0]);  
        if (E->ve_data_cont.SLE) get_static_oceanload(E);

    } // end ice loading

    // Put the initial loads into slice.load[0/2], and set init_potential:
    apply_new_loads(E);

    if (E->ve_data_cont.Heaviside==1)  
        E->ve_data_cont.change_of_load = 0; // no more load changes in this case
					    
    return;
}

/* =============================================================================
 * double distance(x0,y0,z0,t1,f1)
 *    dist = 2.0*arcsin(dist/diameter) = 2.0*arcsin(0.5*dist/radius) 
 * ===========================================================================*/
double distance(x0,y0,z0,t1,f1)
    double x0,y0,z0,t1,f1;
{

    double dist,x2,y2,z2;

    x2 = sin(t1)*cos(f1); 
    y2 = sin(t1)*sin(f1); 
    z2 = cos(t1);
    dist = sqrt( (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0) );
    dist = 0.5*dist;
    dist = 2.0*asin(dist);
    /* dist = 2.0*arcsin(dist/diameter) = 2.0*arcsin(0.5*dist/radius) */
    return (dist);
}


/* ======================================================= 
  for self-G cases, the relation between stress srr and boundary
  topography needs to add the effects of potential for CMB.
 =======================================================  */
void  append_boundary_forces(E)
    struct All_variables *E;
{
    int m,is,ib,j;
    double con,con1;
    int debug_write_load = 0;
    void sphere_expansion_output();

    if (E->ve_data_cont.SELFG)  {// for conversion of potential -> geoid 
        con1 = E->ve_data_cont.Rsg; // surface
        con =  E->ve_data_cont.Rsg  // CMB
               * (E->data.density_below-E->data.density) / E->data.density;
    }
    else  con = con1 = 0.0 ;

    if (E->parallel.me_loc[3]==E->parallel.nprocz-1)          // surface
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)    {
            E->slice_ve.load[1][m][j] = -(  E->slice_ve.load[0][m][j] 
                                       + E->slice_ve.dynamic_oceanload[m][j] )
                                     + con1*E->potential[0][m][j];
            /* E->slice_ve.load[1][m][j] = -E->slice_ve.load[1][m][j] 
                                     + con1*E->potential[m][is];*/
        }
    if (E->parallel.me_loc[3]==0)                              // CMB
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)    {
            /* the negative sign out of ( ... ) is from the surface
               normal at the bottom is negative or points to the
               -r direction */
            E->slice_ve.load[3][m][j] = - ( E->slice_ve.load[2][m][j] 
                                         - con*E->potential[1][m][j] ) ; 
            /* E->slice_ve.load[3][m][j] = - ( E->slice_ve.load[3][m][j] 
                                        - con*E->potential[m][ib]) ; */
        }

    if (debug_write_load) 
        sphere_expansion_output(E,1,E->slice_ve.load[1],
                                E->sphere.sphc[0],E->sphere.sphs[0],
                                E->monitor.solution_cycles,"surfaceload");

    return;
}

// only useful for compressible models
void  append_multi_boundary_forces(struct All_variables * E)
//    struct All_variables *E;
{
    int m,is,ib,j,i,e;
    double con,con1;
    int debug_write_load = 0;
    void sphere_expansion_output();

    int n,node;
    double density_jump, con2;
    const int onedp = onedvpoints[E->mesh.nsd];
    const int vpts=vpoints[E->mesh.nsd];
    const int lev = E->mesh.levmax;

    if (E->ve_data_cont.SELFG)  { 
        con = E->ve_data_cont.Rsg; 
    }
    else  con = 0.0 ;

    for (m=1;m<=E->sphere.caps_per_proc;m++) {
      for (node=1;node<=E->lmesh.nno;node++) {
           j = (node-1)%E->lmesh.noz + E->lmesh.nzs;
           E->slice_ve.all_load[m][node] = -E->slice_ve.total_load[m][node] + con*E->rho_jump1D[lev][j]*E->slice_ve.total_potential[m][node];
      } // end node loop

      /* add dynamic_oceanload to the Earth's surface */
      if ( E->parallel.me_loc[3]==E->parallel.nprocz-1 ) {
        for (j=1;j<=E->lmesh.nsf;j++)    {
          node = j*E->lmesh.noz;
          E->slice_ve.all_load[m][node] -= E->slice_ve.dynamic_oceanload[m][j];
          // ideally, should mupltiply by non-dim surface g, but it is 1.0
        }
      }
    }  

//    fprintf(E->fp_out,"in append_multi %d %d %d\n",E->monitor.count,E->monitor.solution_cycles,E->lmesh.nno);
//    for (m=1;m<=E->sphere.caps_per_proc;m++)
//          for (node=1;node<=E->lmesh.nno;node++)    {
//                  fprintf(E->fp_out,"%d %g %g\n",node,E->slice_ve.all_load[m][node],E->slice_ve.total_potential[m][node]);
//          }
//     fflush(E->fp_out);

    return;
}


/*==============================================================================
 * calculate_potential (E,X_surf,X_cmb,potential_surf,potential_cmb)
 * =============================================================================
 * Takes the surface masses given in X_surf/cmb, computes their gravitational
 * potential at the surface and CMB, and puts it in potential_surf/cmb. 
 * If E->control.polar_wander is true, it further changes the potential arrays
 * (in spherical harmonic l,m=2,1) to via polar_wander_effects_incomp.
 * Notes:
 *     X_surf/cmb must be nondimensional stresses (such as E->iceload).
 *     This routine will not change the values of X_surf and X_cmb.
 *=========================================================================== */
void calculate_potential_incomp(E,X_surf,X_cmb,potential_surf,potential_cmb,icon)
    struct All_variables *E;
    double **X_surf;
    double **X_cmb;
    double **potential_surf;
    double **potential_cmb;
    int icon;
{
    int ib,is,m,ll,ll1,ll2,mm,p,i,j,k,n;
    int isource,ipolarw;
    double density_surf,density_cmb;
    static double *potential_temp[NCS];

    static int been=0;
    void  sphere_expansion_VE();
    void  polar_wander_effects_incomp();
    void  exchange_sphcs();
    void calculate_potential_final();

    density_surf = 1.0;
    density_cmb = (E->data.density_below-E->data.density)/E->data.density;

    if (been==0)  {
        been = 1;
      for (m=1;m<=E->sphere.caps_per_proc;m++)
         potential_temp[m]=(double *) malloc((E->lmesh.nno+2)*sizeof(double)); 
    }

   for (m=1;m<=E->sphere.caps_per_proc;m++)
   for (j=1;j<=E->lmesh.nsf;j++) {
        potential_surf[m][j] = 0.0; 
        potential_cmb[m][j] = 0.0; 
        E->Xsurf[1][m][j] = 0.0;
        E->Xsurf[2][m][j] = 0.0;
        }

   for (k=1;k<=E->mesh.noz;k++)   {
   for (p=0;p<=E->sphere.hindice;p++)   {
            E->sphere.sphc_c[k][p] = 0.0;
            E->sphere.sphs_c[k][p] = 0.0;
	    }
	}

    // Ylm expansion of surface mass:
    if (E->parallel.me_loc[3]==E->parallel.nprocz-1)   {
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)   {
            // scale nondim stress -> nondim height of rock
            E->Xsurf[1][m][j] = X_surf[m][j]/E->ve_data_cont.surf_scaling;
        }
        sphere_expansion_VE(E,E->Xsurf[1],
                         E->sphere.sphc[0],E->sphere.sphs[0],E->output.llmax);
        for (p=0;p<=E->sphere.hindice;p++)   {
            E->sphere.sphc[0][p] *= density_surf;
            E->sphere.sphs[0][p] *= density_surf;
        }
    }

    // Ylm expansion of CMB mass:
    if (E->parallel.me_loc[3]==0)  {
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)   {
            // scale nondim stress -> nondim height of cmb boundary
            E->Xsurf[2][m][j] = X_cmb[m][j]/E->ve_data_cont.botm_scaling;
        }
        sphere_expansion_VE(E,E->Xsurf[2],
                          E->sphere.sphc[1],E->sphere.sphs[1],E->output.llmax);
        for (p=0;p<=E->sphere.hindice;p++)   {
            E->sphere.sphc[1][p] *= density_cmb;
            E->sphere.sphs[1][p] *= density_cmb;
        }
    }

    if (E->parallel.nprocz>1) { // if >1 cap in z-direction (probably not)
        exchange_sphcs(E,E->sphere.sphc[0],E->sphere.sphs[0],
                         E->sphere.sphc[1],E->sphere.sphs[1]);
                         }
    

    if (E->parallel.me_loc[3]==E->parallel.nprocz-1 || E->parallel.me_loc[3]==0)   {
        for (p=0;p<=E->sphere.hindice;p++)   {
            E->sphere.sphc_c[E->mesh.noz][p] = E->sphere.sphc[0][p];
            E->sphere.sphs_c[E->mesh.noz][p] = E->sphere.sphs[0][p];
            E->sphere.sphc_c[1][p] = E->sphere.sphc[1][p];
            E->sphere.sphs_c[1][p] = E->sphere.sphs[1][p];
        }
    }

    isource=0;
    calculate_potential_final(E,E->sphere.sphc_c,E->sphere.sphs_c,potential_temp,E->output.llmax,E->ve_data_cont.compressible,isource,icon);

    if (E->parallel.me_loc[3]==0)  {
      for (m=1;m<=E->sphere.caps_per_proc;m++)
      for (j=1;j<=E->lmesh.nsf;j++) {
        n = (j-1)*E->lmesh.noz + 1;
        potential_cmb[m][j] = potential_temp[m][n];
        }
      }
    if (E->parallel.me_loc[3]==E->parallel.nprocz-1)   {
      for (m=1;m<=E->sphere.caps_per_proc;m++)
      for (j=1;j<=E->lmesh.nsf;j++) {
        n = j*E->lmesh.noz;
        potential_surf[m][j] = potential_temp[m][n];
        }
      }

    // alter the Y_21 potential, if polar_wander is on:
    if (E->ve_data_cont.polar_wander && icon==1)  
        // Xsurf[1/2] have already been converted to nondim height of boundary
        polar_wander_effects_incomp( E, X_surf, X_cmb, 
                                 potential_surf, potential_cmb );

    return;
}




//>>>> add by tao: only calculate the phi, the CM effect should be corrected before separately
// todo_tao: not implemented.

void exchange_sphcs_to_all(struct All_variables *E, double **sphc, int i) // i means the sending layer, index of sphs, global z index.
{

 void parallel_process_termination();
 int p,target_proc,msginfo[8];


 static int been_here = 0;
 static int idb,sizeofk;
 static double *RV, *RV_positive, *RV_negative;

 MPI_Status status[100];
 MPI_Status status1;
 MPI_Request request[100];

 if (been_here ==0 )   {

   idb =  E->sphere.hindice;
   sizeofk = (idb+1)*sizeof(double);
	RV = (double *)malloc( sizeofk );
   been_here ++;
   }

    MPI_Allreduce(sphc[i], RV, idb, MPI_DOUBLE, MPI_SUM , E->parallel.vertical_comm);
    
    for (p=0;p<E->sphere.hindice;p++)  {
      sphc[i][p] = RV[p];
      }

	return;
}

void exchange_depth_to_all(struct All_variables *E, double *depth) 
{

 void parallel_process_termination();
 int p,target_proc,msginfo[8];


 static int been_here = 0;
 static int idb,sizeofk;
 static double *RV;

 MPI_Status status[100];
 MPI_Status status1;
 MPI_Request request[100];

 if (been_here ==0 )   {

   idb =  E->mesh.noz + 1;
   sizeofk = idb*sizeof(double);

   RV =(double *)malloc( sizeofk );
   been_here ++;
   }

    MPI_Allreduce(depth, RV, idb, MPI_DOUBLE, MPI_SUM , E->parallel.vertical_comm);
    
    for (p=1;p<=E->mesh.noz;p++)  {
		depth[p] = RV[p];
      }
	
	return;
}


void calculate_potential_comp(struct All_variables *E, double ** all_potential, int load_only)
{
	void calculate_potential_comp_OR_correct_CM();

	calculate_potential_comp_OR_correct_CM(E,all_potential,load_only,1,0,0);
	return;
}


/*
  This subroutine severs two distinct purposes:
  A) if Only_Correct_CM is false, it is used to calculate potential. (used in calculate_potential_comp)
  B) if Only_Correct_CM is true, it is used to correct CM motion, i.e., updating the all_stress variable. (used in load_to_CM_grav_comp)
      add_incr_CM only works with Only_Correct_CM
  
  icon is associated with centrifugal potential, that is, if icon = 0, this subroutine will ignore calculation about the centrifugal potential.
*/
void calculate_potential_comp_OR_correct_CM(struct All_variables *E, double ** all_potential, int load_only, int icon, int Only_Correct_CM, int add_incr_CM)
{
    int ib,is,m,ll,ll1,ll2,mm,p,i,j,k,n;
    double con,r,t,f,rir,ri,ro,modified_plgndr_a();
    double density_surf,density_cmb;

    static int been=0;
    void  sphere_expansion_VE();
    void  exchange_sphcs();
    void remove_average();
    void construct_c3x3matrix_el();
    void calculate_potential_final();
    void polar_wander_effects_comp();
    void center_mass_effects_comp();

    int global_z_id,noz,node,llmax;
    static double *global_r; 

    static double * rho_div_u[NCS];
    static double M_conv, C_conv, omega;

    double *gnxx, *gnda;
    double rtf[4][9],sinaa[9],ct[9],div_u_gp[9],div_u[4][9][9];
    static struct CC Cc;
    static struct CCX Ccx;
    double temp;
    int a, el;
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int vpts=vpoints[E->mesh.nsd];
    const int ends=enodes[E->mesh.nsd];
    
    void get_rtf_at_vpts();
    double phi_lm_cs, phi_lm_sn, tmp, dlayer;
    double M[3],centri_pot;

    static int been_here = 0;
    double  potential[3],prod[3];
    static double *oldpotl_surf[NCS], *oldpotl_cmb[NCS];
    int procz;
    double cm[3];
    static double Earth_mass;
    double xy_mag, r0, t0, f0, density_jump, t1, f1, rel_cos;


    ri = E->sphere.ri;
    ro = E->sphere.ro;
    density_surf = 1.0;
    density_cmb = (E->data.density_below-E->data.density)/E->data.density;

    if (been==0)  {
        been = 1;
        rir = 1.0;

        for(m=1; m<=E->sphere.caps_per_proc;m++){
          rho_div_u[m] = (double *) malloc((E->lmesh.nno+2)*sizeof(double));
        }

        omega = 2.0*M_PI/(24.0*3600.0);  // rotation_rate
        M_conv = (3.0*E->data.grav_const*E->data.density)/(omega*omega*E->ve_data_cont.kf)*5.0*sqrt(4.0*M_PI/15.0); // 3G*rho/(w^2*kf)*5sqrt(4pi/15)
          // convert non-dim phi into non-dim M
        C_conv = -omega*omega/(4.0*M_PI*E->data.grav_const*E->data.density);
          // convert non-dim r^2(m1*cosf+m2*sinf)*sintcost into non-dim centrifugal potential, which can be added into the potential terms;
          // to dimensionalize the potential, multiply by 4*pi*G*rho*R^2.
    }

    //comment by tao: should set to 0 every time called, since only layer in this processor is modified, others will be lefted.
    for (i = 1; i <= E->mesh.noz; ++i) { // commented by tao: only worked for compressible, otherwise out of range
      for (p = 0; p < E->sphere.hindice+2; ++p) {
        E->sphere.sphc_c[i][p] = 0.;
        E->sphere.sphs_c[i][p] = 0.;
      }
    }

    if(!Only_Correct_CM)
      llmax=E->output.llmax;
    else
      llmax=2; // although 1 is enough.

  /* 
    Part I.2. get sphc_v and sphs_v. calculate div_u 's sph expansion, 
              only needed for !load_only 
            // !!! compute rho*div_u*dV at each gaussian points, then spherical expansion
            // here will get sphere.sphc_v and sphs_v as geruo use sc_v ss_v.
  */

    for (m=1;m<=E->sphere.caps_per_proc;m++)
         for(i=1;i<=E->lmesh.nno;i++)
              rho_div_u[m][i] = 0.0;

    if (!load_only){
            for (m=1;m<=E->sphere.caps_per_proc;m++)
            for (el=1;el<=E->lmesh.nel;el++) {
              gnxx = E->gNX[m][el].vpt;
              gnda = E->gDA[m][el].vpt;
              get_rtf_at_vpts(E,m,E->mesh.levmax,el,rtf); // get r,t,f at gaussian points
              if ((el-1)%E->lmesh.elz==0) {
                construct_c3x3matrix_el(E,el,&Cc,&Ccx,E->mesh.levmax,m,0);
              }
              temp = 0.;
              for (j=1;j<=vpts;j++) {
                // compute div_u at gaussian points
                sinaa[j] = sin(rtf[1][j]);
                ct[j] = cos(rtf[1][j])/sinaa[j];
                div_u_gp[j] = 0.0;
                for (a=1;a<=ends;a++) {
                  node = E->ien[m][el].node[a];
                  for (k=1;k<=dims;k++) {
                    div_u[k][a][j] = rtf[3][j]*
                            (gnxx[GNVXINDEX(0,a,j)]*Cc.vpt[BVINDEX(1,k,a,j)]
                            +E->N.vpt[GNVINDEX(a,j)]*Ccx.vpt[BVXINDEX(1,k,1,a,j)]
                            +E->N.vpt[GNVINDEX(a,j)]*Cc.vpt[BVINDEX(3,k,a,j)])
                                    +rtf[3][j]*
                            (E->N.vpt[GNVINDEX(a,j)]*Cc.vpt[BVINDEX(1,k,a,j)]*ct[j]
                            +E->N.vpt[GNVINDEX(a,j)]*Cc.vpt[BVINDEX(3,k,a,j)]
                            +(gnxx[GNVXINDEX(1,a,j)]*Cc.vpt[BVINDEX(2,k,a,j)]
                            +E->N.vpt[GNVINDEX(a,j)]*Ccx.vpt[BVXINDEX(2,k,2,a,j)])/sinaa[j])
                                    +gnxx[GNVXINDEX(2,a,j)]*Cc.vpt[BVINDEX(3,k,a,j)];
                    div_u_gp[j] +=div_u[k][a][j]*E->U[m][E->id[m][node].doff[k]];
                  }  
                }

                temp += div_u_gp[j]*E->erho[E->mesh.levmax][m][el];

              } // end j loop
              temp = temp/vpts;
              for(j=1;j<=ends;j++)                {
                node = E->ien[m][el].node[j];
                rho_div_u[m][node] += E->TWW[E->mesh.levmax][m][el].node[j]*temp;
              }
            } // end el loop

            (E->exchange_node_d) (E,rho_div_u,E->mesh.levmax);

            for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(node=1;node<=E->lmesh.nno;node++)
               rho_div_u[m][node] *= E->MASS[E->mesh.levmax][m][node];

            for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(node=1;node<=E->lmesh.nno;node++)  {
	      i = (node-1)%E->lmesh.noz+E->lmesh.nzs;
              dlayer = (E->sphere.rr[i+1]-E->sphere.rr[i-1])*0.5;
              if (i==E->mesh.noz) dlayer = (E->sphere.rr[i]-E->sphere.rr[i-1])*0.5;
              else if (i==1) dlayer = (E->sphere.rr[i+1]-E->sphere.rr[i])*0.5;
              rho_div_u[m][node] = rho_div_u[m][node]*dlayer;
	       // rho_div_u is now the surface density like all_stress
	    }

	 }   // end for if (!load_only)
	    
        noz = E->lmesh.noz-1;
        if (E->parallel.me_loc[3]==E->parallel.nprocz-1) noz = E->lmesh.noz;
            // noz as defined avoids double counting in exchange_sphcs_to_all
        for (i=1;i<=noz;i++)  {
              // save the ith layer into E->Xsurf[1]
          if (!load_only)  {
            for (m=1;m<=E->sphere.caps_per_proc;m++) 
              for (j=1;j<=E->lmesh.nsf;j++) {
                  node = (j-1)*E->lmesh.noz+i;
                  E->Xsurf[1][m][j] = E->slice_ve.all_stress[m][node] 
			            - rho_div_u[m][node];
                  }
	      }
	  else {
            for (m=1;m<=E->sphere.caps_per_proc;m++) 
              for (j=1;j<=E->lmesh.nsf;j++) {
                  node = (j-1)*E->lmesh.noz+i;
                  E->Xsurf[1][m][j] = E->slice_ve.all_stress[m][node];
                  }
	      }
          global_z_id = E->lmesh.nzs + i -1;
          sphere_expansion_VE(E,E->Xsurf[1],
                  E->sphere.sphc_c[global_z_id],E->sphere.sphs_c[global_z_id],llmax);
         }  // end for loop i 

    if (E->parallel.nprocz>1)  // if >1 cap in z-direction (probably not)
      for(k=1; k<=E->mesh.noz; k++) {  //comment by tao: global z;
        exchange_sphcs_to_all(E,E->sphere.sphc_c, k);
        exchange_sphcs_to_all(E,E->sphere.sphs_c, k);
      }

  /*
    The following part of this subroutine is either calculating potential or correcting CM motion, depending on the value of Only_Correct_CM.
  */

  if  (!Only_Correct_CM)  {     // if Only_Correct_CM is false
      // Now calculate potential:
      calculate_potential_final(E,E->sphere.sphc_c,E->sphere.sphs_c,all_potential,llmax,1,1,icon);

       // note M[0] and M[1] were computed already in calculate_potential_final() above
       // for me_loc[3] at the surface
      if  (E->ve_data_cont.polar_wander && icon == 1)
          polar_wander_effects_comp(E,all_potential,icon);

     }    // end if !Only_Correct_CM

//    Part III. This part functions to correct CM motion. (for Only_Correct_CM)
//      This part corrects the all_stress to account for the effect of CM.
//      Having it here is not ideal, but is convienent for now. 
//      The reason is it requires the potential to be calculated.

  if (Only_Correct_CM)  { 
      center_mass_effects_comp(E,all_potential,add_incr_CM);
   }

  return;
}

void calculate_potential_deg1_2(E,X_surf,X_cmb)
    struct All_variables *E;
    double **X_surf;
    double **X_cmb;
{
    int llmax,ib,is,m,ll,ll1,ll2,mm,p,i,j,k,n;
    double con,r,t,f,rir,ri,ro,modified_plgndr_a();
    double density_surf,density_cmb;

    static int been=0;
    static double brll1[65],brll2[65];
    static double srll1[65],srll2[65];
    void  sphere_expansion_VE();
    void  exchange_sphcs();

    ri = E->sphere.ri;
    ro = E->sphere.ro;
    density_surf = 1.0;
    density_cmb = (E->data.density_below-E->data.density)/E->data.density;

    llmax = 2;
    if (been==0)  {
        been = 1;
        for (ll=1;ll<=llmax;ll++)    {
            srll1[ll] = pow(ro,(double)(ll))/(2.0*ll+1.0); 
            srll2[ll] = ri*pow(ri,(double)(ll+1))/(2.0*ll+1.0); 
        }
    }

    // Ylm expansion of surface mass:
    if (E->parallel.me_loc[3]==E->parallel.nprocz-1)   {
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)   {
            // scale nondim stress -> nondim height of rock
            E->Xsurf[1][m][j] = X_surf[m][j]/E->ve_data_cont.surf_scaling;
        }
        sphere_expansion_VE(E,E->Xsurf[1],
                             E->sphere.sphc[0],E->sphere.sphs[0],llmax);
        for (ll=0;ll<=llmax;ll++)
        for (mm=0;mm<=ll;mm++)   {
            p = E->sphere.hindex[ll][mm];
            E->sphere.sphc[0][p] *= density_surf;
            E->sphere.sphs[0][p] *= density_surf;
        }
    }

    // Ylm expansion of CMB mass:
    if (E->parallel.me_loc[3]==0)  {
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)   {
            // scale nondim stress -> nondim height of cmb boundary
            E->Xsurf[2][m][j] = X_cmb[m][j]/E->ve_data_cont.botm_scaling;
        }
        sphere_expansion_VE(E,E->Xsurf[2],
                             E->sphere.sphc[1],E->sphere.sphs[1],llmax);
        for (ll=0;ll<=llmax;ll++)
        for (mm=0;mm<=ll;mm++)   {
            p = E->sphere.hindex[ll][mm];
            E->sphere.sphc[1][p] *= density_cmb;
            E->sphere.sphs[1][p] *= density_cmb;
        }
    }

    if (E->parallel.nprocz>1)  // if >1 cap in z-direction (probably not)
        exchange_sphcs(E,E->sphere.sphc[0],E->sphere.sphs[0],
                         E->sphere.sphc[1],E->sphere.sphs[1]);

    // Now calculate potential:
	//comment by tao: seems here sphere.sphc[0][1] is only implemented in top and bottom processors? todo_tao
    ll=1; mm=0;  
    p = E->sphere.hindex[ll][mm]; 
    E->ve_data_cont.CM_pot[0]= srll1[ll]*E->sphere.sphc[0][p]+srll2[ll]*E->sphere.sphc[1][p];

    ll=1; mm=1;  
    p = E->sphere.hindex[ll][mm]; 
    E->ve_data_cont.CM_pot[1]= srll1[ll]*E->sphere.sphc[0][p]+srll2[ll]*E->sphere.sphc[1][p];
    E->ve_data_cont.CM_pot[2]= srll1[ll]*E->sphere.sphs[0][p]+srll2[ll]*E->sphere.sphs[1][p];

    ll=2; mm=1;  
    p = E->sphere.hindex[ll][mm]; 
    E->ve_data_cont.PW_pot[0]= srll1[ll]*E->sphere.sphc[0][p]+srll2[ll]*E->sphere.sphc[1][p];
    E->ve_data_cont.PW_pot[1]= srll1[ll]*E->sphere.sphs[0][p]+srll2[ll]*E->sphere.sphs[1][p];

  return;

  }

/*===========================================================================
 * get_potential(E,count)
 *===========================================================================
 * Sets E->potential[0/1], the nondim grav'l potential at the surface/cmb.
 *               |  pred      (if there is no self-grav iteration)
 *   potential = |  init      (for the first step in a self-grav iteration)
 *               |  init+incr (if in the middle of self-grav interation)
 *=========================================================================== */
void get_potential(struct All_variables *E, int count)
//    struct All_variables *E;
 //   int count;
{
    int ib,is,m,ll,ll1,ll2,mm,p,i,j,k,n;
    void calculate_potential_incomp();
    void calculate_potential_comp();
    void get_potential_from_3D_rho();
    void parallel_process_termination();
    double total_surface_integral();
    static int been_here=0;

/*    double temp1,*TG[4];

  for (m=1;m<=E->sphere.caps_per_proc;m++)
    TG[m] = (double *)malloc((E->lmesh.nsf+1)*sizeof(double));
*/

    if (E->control.CONVECTION) {  // compute potential from internal loads
				 
       if (count==0 && been_here==0) {  // only computed once
          get_potential_from_3D_rho(E,E->buoyancy,E->slice_ve.potential_int_loads,E->ve_data_cont.compressible);  
	   // potential_int_loads is the potential caused by mantle buoyancy.  
	   //
	  if (E->ve_data_cont.compressible)  {
	    for (m=1;m<=E->sphere.caps_per_proc;m++)
              for (n=1;n<=E->lmesh.nno;n++)   {
		E->slice_ve.init_total_potential[m][n] += E->slice_ve.potential_int_loads[m][n];
	      }
	  }
	  else {  // incompressible
	    for (m=1;m<=E->sphere.caps_per_proc;m++)
              for (j=1;j<=E->lmesh.nsf;j++)   {
		n = j*E->lmesh.noz;
		k = n-E->lmesh.noz+1;
               	E->init_potential[1][m][j] += E->slice_ve.potential_int_loads[m][k];
               	E->init_potential[0][m][j] += E->slice_ve.potential_int_loads[m][n];
            	}
	  }

	  been_here=1;
       }

    }

    if (count==0)   {  // => first (perhaps only) visit here for this timestep
            // We are performing self-grav iteration. So initialize E->potential
            // with the potential at the beginning of the timestep.

			//>>>> modify by tao:
            if (E->ve_data_cont.compressible)  {
		for(m=1;m<=E->sphere.caps_per_proc;m++)
			for (i=1; i<=E->lmesh.nno; ++i) {
				E->slice_ve.all_potential[m][i] = 0.0;
				E->slice_ve.total_potential[m][i] = E->slice_ve.init_total_potential[m][i];
				}
            	}
	    else  {
		for (m=1;m<=E->sphere.caps_per_proc;m++)
            		for (j=1;j<=E->lmesh.nsf;j++)   {
                	E->potential[1][m][j] = E->init_potential[1][m][j];
                	E->potential[0][m][j] = E->init_potential[0][m][j];
            		}
		}			
			//<<<< end tao
         }

    else  {   // second time or later in self-grav iteration
		//>>>> add by tao:
	if (E->ve_data_cont.compressible)  {
		calculate_potential_comp(E, E->slice_ve.all_potential,0); 

		for(m=1;m<=E->sphere.caps_per_proc;m++)
		  	for (i =1; i <=E->lmesh.nno; ++i) {
			  E->slice_ve.total_potential[m][i] = E->slice_ve.init_total_potential[m][i]
					+ E->slice_ve.all_potential[m][i];
				  }
			}		
		//<<<< end tao
		
       else  { // comment by tao: here is original incompressible version.
	 // calculate E->incr_potential[0/1] from E->slice.surf/botm[2]:
	      calculate_potential_incomp( E, E->slice_ve.surf[2], E->slice_ve.botm[2], 
	                       E->incr_potential[0], E->incr_potential[1], 1);

	      for (m=1;m<=E->sphere.caps_per_proc;m++)
	        for (j=1;j<=E->lmesh.nsf;j++)   {
		            E->potential[1][m][j] =  E->init_potential[1][m][j] 
		                                   + E->incr_potential[1][m][j];  // cmb
		            E->potential[0][m][j] =  E->init_potential[0][m][j] 
		                                   + E->incr_potential[0][m][j];  // surface
	        }
	}
    }

if(E->ve_data_cont.compressible){
		    // save the surface potential as in the original version; used in dynamic ocean calculation.
		    for (m=1;m<=E->sphere.caps_per_proc;m++)
		      for (j=1;j<=E->lmesh.nsf;j++)   {
		        E->potential[0][m][j] = E->slice_ve.total_potential[m][j*E->lmesh.noz];
		      }

		}

/*
for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++)  { 
         TG[m][j] = fabs(E->potential[0][m][j]);
            }
temp1=total_surface_integral(E,TG,1);
if (E->parallel.me==0) fprintf(stderr,"in get_potential %g count %d\n",temp1,count);
    parallel_process_termination();

  free ((void *)TG);
*/

    return;
}

/* ======================================================= 
   for each element el on boundaries, determine contributions
   from incremental displacement. Only the matrix is determined
   here.
 =======================================================  */
void  construct_B_R(E)
  struct All_variables *E;
 {

  void get_global_1d_shape_fn_2();
  void  exchange_id_d();

  const int onedp=onedvpoints[E->mesh.nsd];

  struct Shape_function1 GM;
  struct Shape_function1_dA dGammax;

  double con,slope;
  double temp,x[3][3],force[9],force_at_gs[9];
  int lev,m,ic,nn[5],e,el, i, k, p, p2, lnode[5];

 if (E->parallel.me_loc[3]==E->parallel.nprocz-1 || E->parallel.me_loc[3]==0)  {
 for (lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)     
 for (m=1;m<=E->sphere.caps_per_proc;m++)   {

   if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  
     for (e=1;e<=E->lmesh.SNEL[lev];e++)    {
       ic = 1;     /*top  */
       el = e*E->lmesh.ELZ[lev];

       get_global_1d_shape_fn_2(E,el,&GM,&dGammax,ic,lev,m);

       for(k=1;k<=onedp;k++)
          nn[k] = k+ic*onedp;

       for(k=1;k<=onedp;k++)
         for (i=1;i<=onedp;i++)
           E->B_R[lev][m][nn[k]][(e-1)*onedp+i] =
                     E->M.vpt[GMVINDEX(k,i)]
                     * dGammax.vpt[GMVGAMMA(0,i)];
       }
   if (E->parallel.me_loc[3]==0)  
     for (e=1;e<=E->lmesh.SNEL[lev];e++)     {
       ic = 0;     /*bottom  */
       el = (e-1)*E->lmesh.ELZ[lev]+1;

       get_global_1d_shape_fn_2(E,el,&GM,&dGammax,ic,lev,m);

       for(k=1;k<=onedp;k++)
          nn[k] = k+ic*onedp;

       for(k=1;k<=onedp;k++)
         for (i=1;i<=onedp;i++)
           E->B_R[lev][m][nn[k]][(e-1)*onedp+i] =
                     E->M.vpt[GMVINDEX(k,i)]
                     * dGammax.vpt[GMVGAMMA(0,i)];
       }
    }
  }

 return;
 }

void  construct_B_RR(struct All_variables *E)
 // struct All_variables *E;
 {

  void get_global_1d_shape_fn_2();
  void  exchange_id_d();

  const int onedp=onedvpoints[E->mesh.nsd];

  struct Shape_function1 GM;
  struct Shape_function1_dA dGammax;

  double con,slope;
  double temp,x[3][3],force[9],force_at_gs[9];
  int lev,m,ic,nn[5],e,el, i, k, p, p2, lnode[5], j;

  if (1)  {
   for (lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)
     for (m=1;m<=E->sphere.caps_per_proc;m++)   {
       for (e=1;e<=E->lmesh.SNEL[lev];e++)    {
         for (j=1;j<=E->lmesh.ELZ[lev];j++) {
            ic = 1;
            el = (e-1)*E->lmesh.ELZ[lev] + j;
            get_global_1d_shape_fn_2(E,el,&GM,&dGammax,ic,lev,m);
            
            for(k=1;k<=onedp;k++)
              nn[k] = k+ic*onedp;
            for(k=1;k<=onedp;k++)
              for (i=1;i<=onedp;i++)
                E->B_RR[lev][m][nn[k]][(el-1)*onedp+i] =
                     E->M.vpt[GMVINDEX(k,i)]
                     * dGammax.vpt[GMVGAMMA(0,i)];
         }
       }
       if (1 /*E->parallel.me_loc[3]==0*/) {
         for (e=1;e<=E->lmesh.SNEL[lev];e++)     {
           ic = 0;     /*bottom  */
           el = (e-1)*E->lmesh.ELZ[lev]+1;
           get_global_1d_shape_fn_2(E,el,&GM,&dGammax,ic,lev,m);

           for(k=1;k<=onedp;k++)
             nn[k] = k+ic*onedp;
           for(k=1;k<=onedp;k++)
             for (i=1;i<=onedp;i++)
               E->B_RR[lev][m][nn[k]][(el-1)*onedp+i] =
                     E->M.vpt[GMVINDEX(k,i)]
                     * dGammax.vpt[GMVGAMMA(0,i)];
         }
       }
     }
  }
  return;
}

/*===========================================================================
 * add_restoring(E,lev,m,Au,u)
 *===========================================================================*/
void add_restoring(E,lev,m,Au,u)
    struct All_variables *E;
    double **Au,**u;
    int lev,m;
{

    const int onedp=onedvpoints[E->mesh.nsd];

    double temp,x[3][3],force[9],force_at_gs[9];
    int ic,nn[5],el,e, i, k, p, p2, lnode[5];

    if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {

        for (e=1;e<=E->lmesh.SNEL[lev];e++)   {
            el=e*E->lmesh.ELZ[lev];

            for(k=1;k<=onedp;k++)   {
                nn[k] = E->ID[lev][m][E->IEN[lev][m][el].node[k+onedp]].doff[3];
                force[k] = E->ve_data_cont.surf_scaling*u[m][nn[k]];
            }

            for (i=1;i<=onedp;i++)   {
                force_at_gs[i] = 0.0;
                for(k=1;k<=onedp;k++)
                    force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];
            }

            for (k=1;k<=onedp;k++)
            for (i=1;i<=onedp;i++)
                Au[m][nn[k]] += E->B_R[lev][m][k+onedp][(e-1)*onedp+i] * force_at_gs[i];

        }
    }

    if (E->parallel.me_loc[3]==0)  {

        for (e=1;e<=E->lmesh.SNEL[lev];e++)   {
            el=(e-1)*E->lmesh.ELZ[lev]+1;

            for(k=1;k<=onedp;k++)   {
                nn[k] = E->ID[lev][m][E->IEN[lev][m][el].node[k]].doff[3];
                force[k] = E->ve_data_cont.botm_scaling*u[m][nn[k]];
            }

            for (i=1;i<=onedp;i++)   {
                force_at_gs[i] = 0.0;
                for(k=1;k<=onedp;k++)
                    force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];
            }

            for(k=1;k<=onedp;k++)
            for (i=1;i<=onedp;i++)
                Au[m][nn[k]] += E->B_R[lev][m][k][(e-1)*onedp+i] * force_at_gs[i];

        }
    }

    return;
}

//add by tao: from geruo
void add_restoring_comp(struct All_variables *E, int lev, int m, double **Au, double** u)
//    struct All_variables *E;
//    double **Au,**u;
//    int lev,m;
{

    void get_global_1d_shape_fn_2();
    const int onedp=onedvpoints[E->mesh.nsd];
    const int vpts=vpoints[E->mesh.nsd];

    double temp,x[3][3],force[9],force_at_gs[9];
    int ic,nn[5],el,e, i, k, p, p2, lnode[5];
    int j,n;
    double density_jump, grav, stress_scaling;

    struct Shape_function1 GM;
    struct Shape_function1_dA dGammax;
    static int been=0;

    stress_scaling = E->data.density*E->data.grav_acc*E->sphere.dradius/E->ve_data_cont.shear_mod;

    // fprintf(E->fp_out,"kkkk density\n");
        for (e=1;e<=E->lmesh.SNEL[lev];e++)  {
          for (j=1;j<=E->lmesh.ELZ[lev];j++) {
            el = (e-1)*E->lmesh.ELZ[lev] + j;
            i = E->lmesh.NZS[lev] + j -1;
	    density_jump = E->rho_jump1D[lev][i+1];  // at the top of each element

            /* find the forces acting on the top nodes of the specific element */
            for (k=1;k<=onedp;k++) {
              nn[k] = E->ID[lev][m][E->IEN[lev][m][el].node[k+onedp]].doff[3];
              grav = E->grav[lev][m][E->IEN[lev][m][el].node[k+onedp]];
            /* always pick the top surface of each element so: k+onedp; doff[3] means the r-direction component*/
              force[k] = density_jump*grav*u[m][nn[k]]*stress_scaling;
            }

//	               fprintf(E->fp_out,"%d %d %g\n",j,el,density_jump*E->data.density);

            /* find the forces at each gaussian point */
            for (i=1;i<=onedp;i++)   {
                force_at_gs[i] = 0.0;
                for(k=1;k<=onedp;k++)
                    force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];
            }

            /* add the forcing terms ( drho*g*u_r, here it only changes the nn(k) terms in Au) to the Au vector */
            for (k=1;k<=onedp;k++)
            for (i=1;i<=onedp;i++)
                  Au[m][nn[k]] += E->B_RR[lev][m][k+onedp][(el-1)*onedp+i] * force_at_gs[i];
          }  // end for j=lmesh.ELZ
        }   // end for e=SNEL

	// do the same for the bottom nodes of bottom elements of bottom CPU, CMB
    if (E->parallel.me_loc[3] ==0)  {
      for (e=1;e<=E->lmesh.SNEL[lev];e++)   {
	    j = 1;  //only for bottom element of each cpu
            el=(e-1)*E->lmesh.ELZ[lev]+j;  
            i = E->lmesh.NZS[lev] + j -1;
	    density_jump = E->rho_jump1D[lev][i];  // at the bottom of this element
            for(k=1;k<=onedp;k++)   {
                nn[k] = E->ID[lev][m][E->IEN[lev][m][el].node[k]].doff[3];
                grav = E->grav[lev][m][E->IEN[lev][m][el].node[k]];
                force[k] = density_jump*grav*u[m][nn[k]]*stress_scaling;
            }

//	                fprintf(E->fp_out,"lll %d %d %g\n",1,el,density_jump*E->data.density);
            for (i=1;i<=onedp;i++)   {
                force_at_gs[i] = 0.0;
                for(k=1;k<=onedp;k++)
                    force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];
            }

            for(k=1;k<=onedp;k++)
            for (i=1;i<=onedp;i++)
                Au[m][nn[k]] += E->B_RR[lev][m][k][(el-1)*onedp+i] * force_at_gs[i];
        }    // end for e=SNEL
      }    // end for bottom CPU

    return;
}


/* =============================================================================
 * apply_new_loads (E)
 * =============================================================================
 * At the beginning of each timestep, just after setting the current value for
 * iceload (in get_iceModel, for example), this routine will prepare the following
 * for the stokes solver:
 *     slice.load[0/2]
 *     init_potential[1/0]
 * Also, the degree-1 iceload is changed via load_to_CM.
 * About the flag E->control.change_of_load:
 *     0 => incr loads are same as last timestep (eg, within an epoch)
 *     1 => incr loads are changing from prev timestep
 *     2 => incr loads are all zero (eg, the last 4kyrs, post-glacial)
 * ===========================================================================*/

void apply_new_loads(struct All_variables * E)  
 //   struct All_variables *E;
{
    int m,i,node;
    void calculate_potential_incomp();
    void calculate_potential_comp();
    double total_surface_integral();
    double temp1,temp2;
    void parallel_process_termination();
    void load_to_CM();
    void load_to_CM_grav();
	void load_to_CM_grav_comp();
        double stress_scaling;

   stress_scaling = E->data.density*E->data.grav_acc*E->sphere.dradius/E->ve_data_cont.shear_mod;

    if (E->ve_data_cont.change_of_load==0) {  // no more ice or static ocean load
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (i=1;i<=E->lmesh.nsf;i++)   {
            E->incr_potential[2][m][i] = 0.0 ;
            E->incr_potential[3][m][i] = 0.0 ;
        }
		//>>>> add by tao:
		if (E->ve_data_cont.compressible){
			for(m=1;m<=E->sphere.caps_per_proc;m++)
			for (i=1; i<=E->lmesh.nno; ++i)
				{
					E->slice_ve.all_load_potential[m][i] = 0.0;
				}
		}
		//<<<< end tao
		

        return;
    }

    if (E->ve_data_cont.Heaviside==2) {  // only for ice-model type loading

        // Put loads into slice.surf/botm[2]:
        if(E->ve_data_cont.compressible){
	        // reset all_stress to zero
	        for (m=1;m<=E->sphere.caps_per_proc;m++)
	          for (i=1;i<=E->lmesh.nno;i++)
	            E->slice_ve.all_stress[m][i] = 0.0;

	        // Put loads into the surface layer of all_stress:
	        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)
	            for (m=1;m<=E->sphere.caps_per_proc;m++)
	            for (i=1;i<=E->lmesh.nsf;i++){
	              node = i*E->lmesh.noz;
	              E->slice_ve.all_stress[m][node] = (E->slice_ve.iceload[0][m][i]
	                                        + E->slice_ve.static_oceanload[m][i])/stress_scaling; /* convert to non-dim rho*u from non-dim stress */
          }
        }
        else{
              if (E->parallel.me_loc[3]==E->parallel.nprocz-1)
                  for (m=1;m<=E->sphere.caps_per_proc;m++)  
                  for (i=1;i<=E->lmesh.nsf;i++)
                      E->slice_ve.surf[2][m][i] =  E->slice_ve.iceload[0][m][i]
                                              + E->slice_ve.static_oceanload[m][i];
              if (E->parallel.me_loc[3]==0)
                  for (m=1;m<=E->sphere.caps_per_proc;m++)  
                  for (i=1;i<=E->lmesh.nsf;i++)
                      E->slice_ve.botm[2][m][i] = 0.0;
        }
		
        // Adjust loads to prevent CM motion (the 1 flag will overwrite the
        // current data.CM_incr, starting fresh for this timestep):
        if(E->ve_data_cont.compressible){
          load_to_CM_grav_comp(E, 1, 1); // load_only=1; add_incr_CM=1, which means overwrite data.CM_incr
          for(i=0;i<3;i++)
            E->ve_data_cont.CM_incr_ice_static_ocean[i] = E->ve_data_cont.CM_incr[i]; // record ice and static ocean induced CM motion.			
        }
        else{
          load_to_CM_grav( E, E->slice_ve.surf[2], E->slice_ve.botm[2], 1 );
        }

		    // comment by tao: update total_load, which is delta_rho * g * Ur + ice and ocean loads. the delta_rho*g*Ur is add in deform_grid
        if(E->ve_data_cont.compressible){
        
          if (E->ve_data_cont.compressible)
          for (m=1;m<=E->sphere.caps_per_proc;m++)
            for (i=1;i<=E->lmesh.nno;i++){
              E->slice_ve.total_load[m][i] +=  E->slice_ve.all_stress[m][i]*E->grav[E->mesh.levmax][m][i]*stress_scaling; // here we need a non-dim stress
            }

          }
        else{
              // Add these loads into slice.load[0/2] (which is used in
              // assemble_forces to create force vector):
              if (E->parallel.me_loc[3]==E->parallel.nprocz-1)
                  for (m=1;m<=E->sphere.caps_per_proc;m++)  
                  for (i=1;i<=E->lmesh.nsf;i++)
                      E->slice_ve.load[0][m][i] +=  E->slice_ve.surf[2][m][i];
              if (E->parallel.me_loc[3]==0)
                  for (m=1;m<=E->sphere.caps_per_proc;m++)  
                  for (i=1;i<=E->lmesh.nsf;i++)
                      E->slice_ve.load[2][m][i] += E->slice_ve.botm[2][m][i];
          }

    } // end if ice-model

    if (E->ve_data_cont.SELFG)   {

        // if the loads have changed, find their new incremental potential due
        // to incremental loads from ice and static ocean loads:
        if (E->ve_data_cont.change_of_load==1 || E->ve_data_cont.SLE)
        {
            // Note: if we're using the SLE, the static_oceanload does not
            // change in a uniform way (due to the time-dependant ocean
            // function).  So we have to recalculate the incr_potential every
            // timestep.

            if (E->ve_data_cont.compressible) {
              calculate_potential_comp(E, E->slice_ve.all_load_potential,
                                       1); /* */
            } else{                         /*comment by tao: incompressible*/
              calculate_potential_incomp(E, E->slice_ve.surf[2], E->slice_ve.botm[2],
                                  E->incr_potential[2], E->incr_potential[3],
                                  1);
            }
        }

        // Add load's incr_potenial to init_ and pred_:
	 /* comment by tao: update init_total_potential with all_potential*/
       if (E->ve_data_cont.compressible)  {
//	       fprintf(E->fp_out,"ooookkkb\n");
		for(m=1;m<=E->sphere.caps_per_proc;m++)
		for (i =1; i <=E->lmesh.nno; ++i) {
		   E->slice_ve.init_total_potential[m][i] += E->slice_ve.all_load_potential[m][i];
//		   fprintf(E->fp_out,"%d %g %g\n",i,E->slice_ve.init_total_potential[m][i],E->slice_ve.all_load_potential[m][i]);
		}
	  }
        else {  /* comment by tao: incompressible*/
		for (m=1;m<=E->sphere.caps_per_proc;m++)
		for (i=1;i<=E->lmesh.nsf;i++)	{
				E->init_potential[0][m][i] += E->incr_potential[2][m][i];
				E->init_potential[1][m][i] += E->incr_potential[3][m][i];
		}
	  }
    }     // end if self-gravitation

    return;
}


/* =============================================================================
 * get_iceModel(E,iceload)  
 * =============================================================================
 * Specifically for iceModel loading.
 * Puts new load into iceload[m][n] (n = surface node). Specifically, 
 * iceload is the incremental (ie, per timestep) nondimensional stress, 
 * after removing its average.
 * Only sets (passed in) array 'iceload'.
 * About the flag E->control.change_of_load:
 *     0 => incr loads are same as last timestep (eg, within an epoch)
 *     1 => incr loads are changing from prev timestep
 *     2 => incr loads are all zero (eg, the last 4kyrs, post-glacial)
 * ========================================================================== */

void get_iceModel(E,iceload)  
    struct All_variables *E;
    double **iceload;
{
    FILE *fp;
    int m,i,j,n;
    static int ifile;  // this is the iceModel 'epoch'
    static int been_here=0;
    static int step_prev=0;
    void sphere_expansion_output();
    void parallel_process_termination();
    void read_reg_grids();
    void read_sve_grids();
    void read_sph_harmonics();
    void remove_average();
    double total,temp1, temp2;
    char outfile[255],input_s[200];

	void get_iceModel_withInitialLoad(struct All_variables * ,double ** );
	
	if(E->ve_data_cont.NonZero_Init_ICELoad == 1){
		get_iceModel_withInitialLoad(E,iceload);
		return;
	}

    if (been_here==0)  {
        // Allocate and initialize some arrays:
        for (m=1;m<=E->sphere.caps_per_proc;m++)   {
            E->slice_ve.ice_height_prev[m] = (double *)
                              malloc((E->lmesh.nsf+2)*sizeof(double));
            E->slice_ve.ice_height_curr[m] = (double *)
                              malloc((E->lmesh.nsf+2)*sizeof(double));
            for (n=1;n<=E->lmesh.nsf;n++)  {  
                E->slice_ve.ice_height_prev[m][n] = 0.0;
                E->slice_ve.ice_height_curr[m][n] = 0.0;
            }
        } // end allocate & initialize

    }

   if (E->ve_data_cont.DIRECT == 0 && been_here!=0) 
        return;

    been_here ++;

        // Read in the iceload for each epoch. 

    ifile = E->ve_data_cont.stage + 1;

    // The following flag will force a re-calculation of load potential in
    // apply_new_loads:
    E->ve_data_cont.change_of_load = 1;

    if (ifile==1)  {       // for the first epoch 
//        sprintf(outfile,"%s%dx%d.%d",
//            E->ve_data_cont.ice_file, E->sphere.elx, E->sphere.ely, ifile-1);
        sprintf(outfile,"%s.%d",E->ve_data_cont.ice_file, ifile-1);
        if (E->ve_data_cont.read_ice_format==0) 			 
           read_reg_grids(E,outfile,E->slice_ve.ice_height_prev);
        else if (E->ve_data_cont.read_ice_format==1)   // on SVE grid 
           read_sve_grids(E,outfile,E->slice_ve.ice_height_prev);
        else if (E->ve_data_cont.read_ice_format==2)  // Jackson form
           read_sph_harmonics(E,outfile,E->slice_ve.ice_height_prev);

        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  {
           E->slice_ve.ice_height_prev[m][j] = E->slice_ve.ice_height_prev[m][j]/E->sphere.dradius; // non-dimensionalized by the Earth's radius
           E->slice_ve.ice_height_prev[m][j] = 0.0; // reset to be ice free
           }
        step_prev = 0;
        }
    else if (ifile>1)  {       // if we're beyond the first epoch, curr->prev
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  
            E->slice_ve.ice_height_prev[m][j] = E->slice_ve.ice_height_curr[m][j];   
        step_prev = E->ve_data_cont.stages_step[E->ve_data_cont.stage-1];
        }

    // read new epoch data into slice.ice_height_curr
    //
    sprintf(outfile,"%s.%d",E->ve_data_cont.ice_file, ifile);
    if (E->ve_data_cont.read_ice_format==0) 			 
           read_reg_grids(E,outfile,E->slice_ve.ice_height_curr);
    else if (E->ve_data_cont.read_ice_format==1)  // on SVE grid
           read_sve_grids(E,outfile,E->slice_ve.ice_height_curr);
    else if (E->ve_data_cont.read_ice_format==2)  //Jackson form
           read_sph_harmonics(E,outfile,E->slice_ve.ice_height_curr);

    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++)  {
        E->slice_ve.ice_height_curr[m][j] = E->slice_ve.ice_height_curr[m][j]/E->sphere.dradius; // non-dimensionalized by the Earth's radius
    }


    // output ice_height_curr for debugging
    if (E->ve_data_cont.change_of_load == 1){
      sprintf(outfile, "%s.ice_height.%d.%d.%d", E->control.data_file,E->parallel.me, ifile, E->monitor.solution_cycles);
      fp = fopen(outfile,"w");
      for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  {
          fprintf(fp,"%.4e\n", E->slice_ve.ice_height_curr[m][j]* E->sphere.dradius);
      }
      fclose(fp);
    }

    
    // temp1 is height->stress scaling / number of timesteps in this epoch
    temp1 = E->ve_data_cont.ice_stress_scale/(E->ve_data_cont.stages_step[E->ve_data_cont.stage] - step_prev);

    // iceload now becomes incremental (ie, per timestep) nondimensional stress
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++) 
        iceload[m][j] = ( E->slice_ve.ice_height_curr[m][j]
                         -E->slice_ve.ice_height_prev[m][j])*temp1;

    // remove average and save the spherical harmonic coefficients:
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (n=1;n<=E->lmesh.nsf;n++)    
        E->Xsurf[3][m][n] = iceload[m][n];
    remove_average(E,E->Xsurf[3],1);
    sprintf(outfile,"iceload_stage%d",ifile);

    sphere_expansion_output( E, 1, E->Xsurf[3], 
                             E->sphere.sphc[0], E->sphere.sphs[0],
                             E->monitor.solution_cycles,outfile );
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (n=1;n<=E->lmesh.nsf;n++)
        iceload[m][n] = E->Xsurf[3][m][n];

    return;
}

// version of get_iceModel: have initial non-zero loading.
void get_iceModel_withInitialLoad(struct All_variables * E,double ** iceload)  
//    struct All_variables *E;
//    double **iceload;
{
    FILE *fp;
    int m,i,j,n;
    static int ifile;  // this is the iceModel 'epoch'
    static int been_here=0;
    static int step_prev=0;
    void sphere_expansion_output();
    void parallel_process_termination();
    void read_reg_grids();
    void remove_average();
    double total,temp1, temp2;
    char outfile[255],input_s[200];

    if (been_here==0)  {
        // Allocate and initialize some arrays:
        for (m=1;m<=E->sphere.caps_per_proc;m++)   {
            E->slice_ve.ice_height_prev[m] = (double *)
                              malloc((E->lmesh.nsf+2)*sizeof(double));
            E->slice_ve.ice_height_curr[m] = (double *)
                              malloc((E->lmesh.nsf+2)*sizeof(double));
            for (n=1;n<=E->lmesh.nsf;n++)  {  
                E->slice_ve.ice_height_prev[m][n] = 0.0;
                E->slice_ve.ice_height_curr[m][n] = 0.0;
            }
        } // end allocate & initialize

        ifile = 0;
        sprintf(outfile,"%s.%d",
            E->ve_data_cont.ice_file, ifile);
        read_reg_grids(E,outfile,E->slice_ve.ice_height_prev);
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  {
           E->slice_ve.ice_height_curr[m][j] = E->slice_ve.ice_height_prev[m][j]/E->sphere.dradius; // non-dimensionalized by the Earth's radius
           E->slice_ve.ice_height_prev[m][j] = 0.0; 
           }
        //step_prev = 0;
        E->ve_data_cont.change_of_load = 1;

        // set elapsed time and solution_cycle back to step 0, since the 1st epoch has non-zero load
        E->monitor.elapsed_time=0;
        E->monitor.solution_cycles = 0;

        E->advection.timestep=E->ve_data_cont.stages_timestep[0];


        for (m=1;m<=E->sphere.caps_per_proc;m++)
          for (j=1;j<=E->lmesh.nsf;j++) 
            iceload[m][j] = E->slice_ve.ice_height_curr[m][j] * E->ve_data_cont.ice_stress_scale;

        // remove average and save the spherical harmonic coefficients:
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (n=1;n<=E->lmesh.nsf;n++)    
            E->Xsurf[3][m][n] = iceload[m][n];
        remove_average(E,E->Xsurf[3],1);
        sprintf(outfile,"iceload_stage%d",ifile);

        sphere_expansion_output( E, 1, E->Xsurf[3], 
                                 E->sphere.sphc[0], E->sphere.sphs[0],
                                 E->monitor.solution_cycles,outfile );
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (n=1;n<=E->lmesh.nsf;n++)
            iceload[m][n] = E->Xsurf[3][m][n];
		
        been_here ++;

        return;
    }

   if (E->ve_data_cont.DIRECT == 0 && been_here!=0) 
        return;

    

        // Read in the iceload for each epoch. 

    ifile = E->ve_data_cont.stage + 1;

    // The following flag will force a re-calculation of load potential in
    // apply_new_loads:
    E->ve_data_cont.change_of_load = 1;

    // Actually, in this case, don't need to do special treatment for ifile==1, since the initial load is read already.
    if (ifile==1)  {       // for the first epoch 
        sprintf(outfile,"%s.%d",
            E->ve_data_cont.ice_file, ifile-1);
        read_reg_grids(E,outfile,E->slice_ve.ice_height_prev);
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  {
           E->slice_ve.ice_height_prev[m][j] = E->slice_ve.ice_height_prev[m][j]/E->sphere.dradius; // non-dimensionalized by the Earth's radius
           // E->slice_ve.ice_height_prev[m][j] = 0.0; // reset to be ice free
           }
        step_prev = 0;
        }
    else if (ifile>1)  {       // if we're beyond the first epoch, curr->prev
        for (m=1;m<=E->sphere.caps_per_proc;m++)
          for (j=1;j<=E->lmesh.nsf;j++)  
            E->slice_ve.ice_height_prev[m][j] = E->slice_ve.ice_height_curr[m][j];   
        step_prev = E->ve_data_cont.stages_step[E->ve_data_cont.stage-1];
        }

    // read new epoch data into slice.ice_height_curr
    //
    sprintf(outfile,"%s.%d",
                E->ve_data_cont.ice_file, ifile);
    read_reg_grids(E,outfile,E->slice_ve.ice_height_curr);

    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++)  {
        E->slice_ve.ice_height_curr[m][j] = E->slice_ve.ice_height_curr[m][j]/E->sphere.dradius; // non-dimensionalized by the Earth's radius
    }
    
    // temp1 is height->stress scaling / number of timesteps in this epoch
    temp1 = E->ve_data_cont.ice_stress_scale/(E->ve_data_cont.stages_step[E->ve_data_cont.stage] - step_prev);

    // iceload now becomes incremental (ie, per timestep) nondimensional stress
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++) 
        iceload[m][j] = ( E->slice_ve.ice_height_curr[m][j]
                         -E->slice_ve.ice_height_prev[m][j])*temp1;

    // remove average and save the spherical harmonic coefficients:
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (n=1;n<=E->lmesh.nsf;n++)    
        E->Xsurf[3][m][n] = iceload[m][n];
    remove_average(E,E->Xsurf[3],1);
    sprintf(outfile,"iceload_stage%d",ifile);

    sphere_expansion_output( E, 1, E->Xsurf[3], 
                             E->sphere.sphc[0], E->sphere.sphs[0],
                             E->monitor.solution_cycles,outfile );
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (n=1;n<=E->lmesh.nsf;n++)
        iceload[m][n] = E->Xsurf[3][m][n];

    return;
}


/////////////////////////////////////////////////////
//
 void read_reg_grids(E,outfile,field)
   struct All_variables *E;
   char *outfile;
   double **field;   
  {

char input_s[1000];
FILE *fp2;

void  parallel_process_termination();
int numtheta,numphi,nregnodes;
int kk,j;
int ntheta,nphi,iregnode,iregnode2;
int node;
int c1,c2,c3,c4;
int isurface_element,lev;

double del_degree;
double del_theta;
double del_phi;
double lat_min,lat_max;
double theta_min,theta_max;
double phi_min,phi_max;
double vphi;
double rlat,rlong;
double theta,phi;
double phiprime,thetaprime;
double xi,eta;
double shape1,shape2,shape3,shape4;

static double *regular_f;
static int been_here=0;

/* read from file */

   if ( (fp2=fopen(outfile,"r"))==NULL)   {
      fprintf(stderr,"ERROR(read_reg_grids)- file %s not found\n",outfile);
      fflush(stderr);
      parallel_process_termination();
   }

   fgets(input_s,1000,fp2);
   sscanf(input_s,"%d %d",&numphi,&numtheta);

//   numtheta=E->sphere.nox-1;       /* number of theta increments */
//   numphi=2*numtheta;               /* number of phi increments */

/* INPUT FILE SPECIFIC STUFF */
/* format longitude, latitude, vphi, vtheta (vtheta is opposite in convention) */
/* file sequence given in specific order (see below):                          */
/* inputfile specific parameters:                                              */

   del_degree=180.0/numtheta; /* increments given in input file */
   lat_min=-(90.0-del_degree/2.0);         /* minimum lattitude */
   lat_max=+(90.0-del_degree/2.0);         /* maximum lattitude */
   phi_min=del_degree/2.0;             /* minimum longitude */
   phi_max=360.0-del_degree/2.0;           /* maximum longitude */
 
   del_theta=del_degree*M_PI/180.0;
   del_phi=del_degree*M_PI/180.0;
   theta_min=-1.0*(lat_max-90.0)*M_PI/180.0;
   theta_max=-1.0*(lat_min-90.0)*M_PI/180.0;
   phi_min=phi_min*M_PI/180.0;
   phi_max=phi_max*M_PI/180.0;

   nregnodes=numtheta*numphi;


/* adjust parameters for ghosting - extra regular elements are added on all sides */
/* This is required to interpolate nodes between thetamin and thetamax, and phimin and phimax */
/* Domain is padded by extra +2 numphi and +2 numtheta */

   numtheta=numtheta+2;
   numphi=numphi+2;
   theta_min=theta_min-del_theta;
   theta_max=theta_max+del_theta;
   phi_min=phi_min-del_phi;
   phi_max=phi_max+del_phi;

   nregnodes=numtheta*numphi;

 if (been_here==0)  {
   regular_f = (double *)malloc((nregnodes+1)*sizeof(double));
   been_here = 1;
   }


   for (kk=1;kk<=nregnodes;kk++)
      regular_f[kk]=-999999;

/* input file specific order of reading  right here */
/* (note - leaving room for ghosting) */

//   for (ntheta=(numtheta-1);ntheta>1;ntheta--)
   for (ntheta=2;ntheta<=(numtheta-1);ntheta++)
     for (nphi=2;nphi<=(numphi-1);nphi++)     {

      iregnode=ntheta+(nphi-1)*numtheta;

      if ((iregnode>nregnodes)||(iregnode<1))   {
        fprintf(stderr,"ERROR(plate velocities)-wrong iregnode %d %d\n",iregnode,nregnodes);
        fflush(stderr);
        exit(10);
      }

      fgets(input_s,1000,fp2);
      sscanf(input_s,"%lf %lf %le",&rlong,&rlat,&vphi);

     regular_f[iregnode]=vphi;

   }

   fclose(fp2);

/* ghost sides of regular domain */

   nphi=1;
   for (ntheta=2;ntheta<=(numtheta-1);ntheta++)  {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=ntheta+((numphi-1)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   nphi=numphi;
   for (ntheta=2;ntheta<=(numtheta-1);ntheta++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=ntheta+(2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   ntheta=1;
   for (nphi=1;nphi<=(numphi/2);nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta+1)+(nphi+numphi/2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }
   for (nphi=(numphi/2)+1;nphi<=numphi;nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta+1)+(nphi-(numphi/2)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   ntheta=(numtheta);
   for (nphi=1;nphi<=(numphi/2);nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta-1)+(nphi+numphi/2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }
   for (nphi=(numphi/2)+1;nphi<=numphi;nphi++)  {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta-1)+(nphi-(numphi/2)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

/* interpolate real mesh */

   isurface_element=0;
   for (j=1;j<=E->sphere.caps_per_proc;j++)
   for (node=1;node<=E->lmesh.nno;node++)   {
     theta=E->sx[j][1][node];
     phi  =E->sx[j][2][node];
     if (node%E->lmesh.noz==0)  {
        isurface_element++;
/* find position on regular mesh */

        ntheta=((theta-theta_min)/del_theta)+1;
        nphi=((phi-phi_min)/del_phi)+1;

/* iregnode is the lower left node in the regular element */

        iregnode=ntheta+(nphi-1)*numtheta;

/* find theta and phi distance from node */

        thetaprime=theta-(theta_min+(ntheta-1)*del_theta);
        phiprime=phi-(phi_min+(nphi-1)*del_phi);

        if (thetaprime<0.0 || thetaprime>del_theta)  {
           fprintf(stderr,"ERROR(get plate velocities)-wierd thetaprime %f %f\n",thetaprime,del_theta);
           fflush(stderr);
           exit(10);
        }
        if (phiprime<0.0 || phiprime>del_phi)   {
           fprintf(stderr,"ERROR(get plate velocities)-wierd phiprime %f %f\n",phiprime,del_phi);
           fflush(stderr);
           exit(10);
        }

/* determine shape functions */

       xi=2.0/del_theta*thetaprime-1.0;
       eta=2.0/del_phi*phiprime-1.0;

       if (xi<(-1-1e-6) || xi > (1+1e-6) || eta < (-1-1e-6) || eta > (1+1e-6))  {
          fprintf(stderr,"bad local variables %f %f (%f %f) (%f %f)\n",xi,eta,thetaprime,phiprime,del_theta,del_phi);
          fflush(stderr);
          exit(10);
       }
       shape1=0.25*(1.0-xi)*(1.0-eta);
       shape2=0.25*(1.0+xi)*(1.0-eta);
       shape3=0.25*(1.0+xi)*(1.0+eta);
       shape4=0.25*(1.0-xi)*(1.0+eta);

       c1=iregnode;
       c2=iregnode+1;
       c3=iregnode+1+numtheta;
       c4=iregnode+numtheta;

       if ( (c3>(numtheta*numphi)) || (c1<1) )  {
         fprintf(stderr,"ERROR(get plate...) - %d %d %d %d\n",c1,c2,c3,c4);
         fflush(stderr);
         exit(10);
       }

       vphi = regular_f[c1]*shape1+regular_f[c2]*shape2
            + regular_f[c3]*shape3+regular_f[c4]*shape4;

       field[j][isurface_element]=vphi;

     }
   }

return;
}

// read ice heights that were already assigned to SVE grid

 void read_sve_grids(E,outfile,field)
   struct All_variables *E;
   char *outfile;
   double **field;   
  {

char input_s[1000],filename[255];
FILE *fp2;

void  parallel_process_termination();
int j,j1;
double temp;

 if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {
   sprintf( filename,"%s.%d",outfile,E->parallel.me);
   if ( (fp2=fopen(filename,"r"))==NULL)   {
      fprintf(stderr,"ERROR(read_reg_grids)- file %s not found\n",outfile);
      fflush(stderr);
      parallel_process_termination();
   }

  for (j=1;j<=E->sphere.caps_per_proc;j++)
    for (j1=1;j1<=E->lmesh.nsf;j1++)   {
      fgets(input_s,1000,fp2);
      sscanf(input_s,"%lf",&temp);
      field[j][j1]=temp;
    }
  fclose(fp2);
 }

  return;
  }

//
// read ice model in terms of spherical harmonics and assign ice to FE grid.
//
 void read_sph_harmonics(E,outfile,field)
   struct All_variables *E;
   char *outfile;
   double **field;   
  {

char input_s[1000];
FILE *fp2;

void  parallel_process_termination();
int llmax,ll,mm,i,j,j1,node,ll1,ll2;
double sqrt2,tempc,temps,temp1;
static double *clm,*slm,*plm_ll,*plm_ll_1,*plm_ll_2,*cosmf,*sinmf;
double modified_plgndr_a(),con,t1,f1,x1,y1;
double twoll, fourl2m1, fact1, fact2, llm1_2;
static int been_here=0;
const double two = 2.0;
const double three = 3.0;

double time1,time2,CPU_time0();

time1=CPU_time0();

   if ( (fp2=fopen(outfile,"r"))==NULL)   {
      fprintf(stderr,"ERROR(read_reg_grids)- file %s not found\n",outfile);
      fflush(stderr);
      parallel_process_termination();
   }

   sqrt2 = sqrt(2.0);

   fgets(input_s,1000,fp2);
   sscanf(input_s,"%d",&llmax);
   
//   llmax=200;

   j = (llmax+1)*(llmax+2)/2;
   if (been_here==0)  {
      clm = (double*)malloc((j+2)*sizeof(double));
      slm = (double*)malloc((j+2)*sizeof(double));
      plm_ll = (double*)malloc((llmax+2)*sizeof(double));
      plm_ll_1 = (double*)malloc((llmax+2)*sizeof(double));
      plm_ll_2 = (double*)malloc((llmax+2)*sizeof(double));
      cosmf = (double*)malloc((llmax+2)*sizeof(double));
      sinmf = (double*)malloc((llmax+2)*sizeof(double));
      been_here = 1;
   }

  //read in clm and slm
  i = 0;
  for (ll=0;ll<=llmax;ll++)
    for (mm=0; mm<=ll; mm++)   {
      i = i+1;
      fgets(input_s,1000,fp2);
      sscanf(input_s,"%d %d %lf %lf",&j,&j,&tempc,&temps);
 	// to turn Jackson Clm and Slm to SVE's
      if (mm==0) {
         temps=-temps;
      }
      else {
         tempc= sqrt2*tempc;
         temps=-sqrt2*temps;
      }
      clm[i] = tempc;
      slm[i] = temps;
    }
   fclose(fp2);


  for (j=1;j<=E->sphere.caps_per_proc;j++)
    for (j1=1;j1<=E->lmesh.nsf;j1++)   
      field[j][j1]=0.0;

  for (j=1;j<=E->sphere.caps_per_proc;j++)
  for (j1=1;j1<=E->lmesh.nsf;j1++)   {
      node = j1*E->lmesh.noz;
      t1=E->sx[j][1][node];
      f1=E->sx[j][2][node];
      x1 = cos(t1);
      for (mm=0; mm<=llmax; mm++)   {
         cosmf[mm] = cos(mm*f1);
         sinmf[mm] = sin(mm*f1);
         }
      i = 0;
      for (ll=0;ll<=llmax;ll++) {

	twoll= two*ll;
	fourl2m1 = twoll*twoll-1.0;    // 4l^2-1
	fact1 = (twoll+1.0)/(twoll-three);  // (2l+1)/(2l-3)
        llm1_2 = ll*ll-twoll+1.0;   // (l-1)*(l-1)
        for (mm=0; mm<=ll; mm++)   {
           i = i+1;
	   fact2 = (double)(ll-mm)/(double)(ll+mm);
	   if (mm==ll) 
 	      plm_ll[mm]=modified_plgndr_a(ll,mm,t1);
	   else if (mm==ll-1)
 	      plm_ll[mm]=x1*sqrt(two*mm+three)*plm_ll_1[mm];
	   else if (mm<ll-1) {
              plm_ll[mm]=x1*sqrt(fourl2m1*fact2)*plm_ll_1[mm] 
		      - sqrt(fact1*fact2*(llm1_2-mm*mm))*plm_ll_2[mm];
	      plm_ll[mm]=plm_ll[mm]/(ll-mm);
	      }

           field[j][j1]+=(clm[i]*cosmf[mm]+slm[i]*sinmf[mm])*plm_ll[mm];
        } // end for mm
	 
	// swap and redefine plm_ll_2, and plm_ll_1 for next ll, 
        for (mm=0; mm<ll; mm++) 
 	   plm_ll_2[mm]=plm_ll_1[mm];
        for (mm=0; mm<=ll; mm++)  
 	   plm_ll_1[mm]=plm_ll[mm];

     }   // end for ll
   }    // end for j1 or surface node
	//
  fprintf(E->fp_out,"icestage %d\n",E->ve_data_cont.stage);
  for (j=1;j<=E->sphere.caps_per_proc;j++)
    for (j1=1;j1<=E->lmesh.nsf;j1++)  { 
      fprintf(E->fp_out,"%d %.5e\n",j1,field[j][j1]);
    }
  fflush(E->fp_out);

time2=CPU_time0();
if (E->parallel.me==0) {
	fprintf(stderr,"CPU_time for ice model %g seconds\n",time2-time1);
        fprintf(E->fp,"CPU_time for ice model %g seconds\n",time2-time1);
	fflush(E->fp);
}

  return;
  }



void interpolate_to_FE_grid(E,field,regular_f,numtheta,numphi)
  struct All_variables *E;
  double *field,*regular_f;
  int numtheta, numphi;
{

void  parallel_process_termination();
int nregnodes;
int kk,j;
int ntheta,nphi,iregnode,iregnode2;
int node;
int c1,c2,c3,c4;
int isurface_element,lev;

double del_degree;
double del_theta;
double del_phi;
double lat_min,lat_max;
double theta_min,theta_max;
double phi_min,phi_max;
double vphi;
double rlat,rlong;
double theta,phi;
double phiprime,thetaprime;
double xi,eta;
double shape1,shape2,shape3,shape4;

   del_degree=180.0/(numtheta-2); /* increments given in input file */
   lat_min=-(90.0-del_degree/2.0);         /* minimum lattitude */
   lat_max=+(90.0-del_degree/2.0);         /* maximum lattitude */
   phi_min=del_degree/2.0;             /* minimum longitude */
   phi_max=360.0-del_degree/2.0;           /* maximum longitude */
 
   del_theta=del_degree*M_PI/180.0;
   del_phi=del_degree*M_PI/180.0;
   theta_min=-1.0*(lat_max-90.0)*M_PI/180.0;
   theta_max=-1.0*(lat_min-90.0)*M_PI/180.0;
   phi_min=phi_min*M_PI/180.0;
   phi_max=phi_max*M_PI/180.0;

   nregnodes=numtheta*numphi;

   theta_min=theta_min-del_theta;
   theta_max=theta_max+del_theta;
   phi_min=phi_min-del_phi;
   phi_max=phi_max+del_phi;

/* ghost sides of regular domain */

   nphi=1;
   for (ntheta=2;ntheta<=(numtheta-1);ntheta++)  {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=ntheta+((numphi-1)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   nphi=numphi;
   for (ntheta=2;ntheta<=(numtheta-1);ntheta++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=ntheta+(2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   ntheta=1;
   for (nphi=1;nphi<=(numphi/2);nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta+1)+(nphi+numphi/2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }
   for (nphi=(numphi/2)+1;nphi<=numphi;nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta+1)+(nphi-(numphi/2)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   ntheta=(numtheta);
   for (nphi=1;nphi<=(numphi/2);nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta-1)+(nphi+numphi/2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }
   for (nphi=(numphi/2)+1;nphi<=numphi;nphi++)  {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta-1)+(nphi-(numphi/2)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

/* interpolate real mesh */

   isurface_element=0;
   for (j=1;j<=E->sphere.caps_per_proc;j++)
   for (node=1;node<=E->lmesh.nno;node++)   {
     theta=E->sx[j][1][node];
     phi  =E->sx[j][2][node];
     if (node%E->lmesh.noz==0)  {
        isurface_element++;
/* find position on regular mesh */

        ntheta=((theta-theta_min)/del_theta)+1;
        nphi=((phi-phi_min)/del_phi)+1;

/* iregnode is the lower left node in the regular element */

        iregnode=ntheta+(nphi-1)*numtheta;

/* find theta and phi distance from node */

        thetaprime=theta-(theta_min+(ntheta-1)*del_theta);
        phiprime=phi-(phi_min+(nphi-1)*del_phi);

        if (thetaprime<0.0 || thetaprime>del_theta)  {
           fprintf(stderr,"ERROR(get plate velocities)-wierd thetaprime %f %f\n",thetaprime,del_theta);
           fflush(stderr);
           exit(10);
        }
        if (phiprime<0.0 || phiprime>del_phi)   {
           fprintf(stderr,"ERROR(get plate velocities)-wierd phiprime %f %f\n",phiprime,del_phi);
           fflush(stderr);
           exit(10);
        }

/* determine shape functions */

       xi=2.0/del_theta*thetaprime-1.0;
       eta=2.0/del_phi*phiprime-1.0;

       if (xi<(-1-1e-6) || xi > (1+1e-6) || eta < (-1-1e-6) || eta > (1+1e-6))  {
          fprintf(stderr,"bad local variables %f %f (%f %f) (%f %f)\n",xi,eta,thetaprime,phiprime,del_theta,del_phi);
          fflush(stderr);
          exit(10);
       }
       shape1=0.25*(1.0-xi)*(1.0-eta);
       shape2=0.25*(1.0+xi)*(1.0-eta);
       shape3=0.25*(1.0+xi)*(1.0+eta);
       shape4=0.25*(1.0-xi)*(1.0+eta);

       c1=iregnode;
       c2=iregnode+1;
       c3=iregnode+1+numtheta;
       c4=iregnode+numtheta;

       if ( (c3>(numtheta*numphi)) || (c1<1) )  {
         fprintf(stderr,"ERROR(get plate...) - %d %d %d %d\n",c1,c2,c3,c4);
         fflush(stderr);
         exit(10);
       }

       vphi = regular_f[c1]*shape1+regular_f[c2]*shape2
            + regular_f[c3]*shape3+regular_f[c4]*shape4;

       field[isurface_element]=vphi;

     }
   }

return;
}

// rho_3d is defined as density anomalies at each node in 3-D volume.
// This subroutine returns the gravitational potential caused by rho_3d. 
//
void get_potential_from_3D_rho(struct All_variables *E, double ** rho_3d, double ** potential_intl_loads, int icompress)
{
    static int been=0;
    void calculate_potential_from_3D_sigma();

    static double *sigma_3d[NCS];

    int el,m,n,i,j,k;
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int vpts=vpoints[E->mesh.nsd];
    const int ends=enodes[E->mesh.nsd];
    const int levmax=E->mesh.levmax;
    
    double tmp, tmp1,dr;

    static int been_here = 0;

    if (been==0)  {
        been = 1;
        for(m=1; m<=E->sphere.caps_per_proc;m++){
          sigma_3d[m] = (double *) malloc((E->lmesh.nno+2)*sizeof(double));
        }
    }

    // turn rho_3d to sigma_3d or surface density also on the nodes.
 
    for (m=1;m<=E->sphere.caps_per_proc;m++) {
      for (n=1;n<=E->lmesh.nno;n++)  {
	      i = (n-1)%E->lmesh.noz+E->lmesh.nzs;
	      dr = (E->sphere.rr[i+1]-E->sphere.rr[i-1])*0.5;
	      if (i==E->mesh.noz) dr = (E->sphere.rr[i]-E->sphere.rr[i-1])*0.5;
	      else if (i==1) dr = (E->sphere.rr[i+1]-E->sphere.rr[i])*0.5;
	      sigma_3d[m][n] = rho_3d[m][n]*dr;
              }
      }
      
    calculate_potential_from_3D_sigma(E,sigma_3d,potential_intl_loads,icompress);

  return;
}
//
void calculate_potential_from_3D_sigma(struct All_variables *E, double ** sigma_3d, double ** potential_intl_loads, int icompress)
{

    static int been=0;
    void  sphere_expansion_VE();
    void  exchange_sphcs();
    void calculate_potential_final();

    int m,i,j,k,p,noz,node,global_z_id,ipolarw,isource,llmax;

    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int vpts=vpoints[E->mesh.nsd];
    const int ends=enodes[E->mesh.nsd];
    
    if (been==0)  {
        been = 1;
    }

    llmax = E->output.llmax;

    for(k=1; k<=E->mesh.noz; k++)    {  //global z;
      for (p = 0; p < E->sphere.hindice+2; ++p) {
        E->sphere.sphc_c[k][p] = 0.;
        E->sphere.sphs_c[k][p] = 0.;
      }
    }

    noz = E->lmesh.noz-1;
    if(E->parallel.me_loc[3]==E->parallel.nprocz-1) noz = E->lmesh.noz;
       // noz as defined avoids double counting in exchange_sphcs_to_all 
    for (i=1;i<=noz;i++) {
          for (m=1;m<=E->sphere.caps_per_proc;m++) {
            for (j=1;j<=E->lmesh.nsf;j++) {
              node = (j-1)*E->lmesh.noz+i;
              E->Xsurf[1][m][j] = sigma_3d[m][node];
            }
          }
          global_z_id = E->lmesh.nzs + i -1;
          sphere_expansion_VE(E,E->Xsurf[1],
             E->sphere.sphc_c[global_z_id],E->sphere.sphs_c[global_z_id],llmax);
    }

    if (E->parallel.nprocz>1) { // if >1 cap in z-direction (probably not)
       for(k=1; k<=E->mesh.noz; k++)    {  //comment by tao: global z;
          exchange_sphcs_to_all(E,E->sphere.sphc_c, k);
          exchange_sphcs_to_all(E,E->sphere.sphs_c, k);
          }
      } // end for if E->parallel.nprocz>1

    ipolarw = 1;   // determine polar motion
    isource = 1;   // volume source
    calculate_potential_final(E,E->sphere.sphc_c,E->sphere.sphs_c,potential_intl_loads,llmax,icompress,isource,ipolarw);

    return;
}


// compute potential from the sources in sphc_c and sphs_c 
//
void  calculate_potential_final(struct All_variables *E, double ** sphc_c, double ** sphs_c, double ** potential, int llmax, int icompress,int isource, int icon)
{
    int ib,is,m,ll,ll1,ll2,mm,p,i,j,k,n,ik,node;
    static int been=0;
    int global_z_id;
    int noz,nz1,nz2,nz3,mz1,mz2,mz3;

    static double M_conv, C_conv, omega;
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int vpts=vpoints[E->mesh.nsd];
    const int ends=enodes[E->mesh.nsd];
    
    double phi_lm_cs, phi_lm_sn, tmp, tmp1, dlayer;
    double M[3],centri_pot;

    double cm[3];
    static double Earth_mass;

    if (been==0)  {
        been = 1;
        omega = 2.0*M_PI/(24.0*3600.0);  // rotation_rate
        M_conv = (3.0*E->data.grav_const*E->data.density)/(omega*omega*E->ve_data_cont.kf)*5.0*sqrt(4.0*M_PI/15.0); // 3G*rho/(w^2*kf)*5sqrt(4pi/15)
          // convert non-dim phi into non-dim M
        C_conv = -omega*omega/(4.0*M_PI*E->data.grav_const*E->data.density);
          // convert non-dim r^2(m1*cosf+m2*sinf)*sintcost into non-dim centrifugal potential, which can be added into the potential terms;
          // to dimensionalize the potential, multiply by 4*pi*G*rho*R^2.
    }

    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (node=1;node<=E->lmesh.nno;node++) 
        potential[m][node] = 0.0;

    // determine the depth at which the field is computed
    nz1 = 1;
    nz2 = E->lmesh.noz;
    nz3 = 1;
    if (icompress==0)  { // incompressible, only need CMB and surface potential
	  if (E->parallel.nprocz==1) {  // only 1 cpu in z direction
              nz3 = E->lmesh.noz-1;
	  }
	  else {  // >1 cpu in z
            if (E->parallel.me_loc[3]==0) {
              nz3 = E->lmesh.noz+1;  // only do for the CMB
	    }
	    else if (E->parallel.me_loc[3]==E->parallel.nprocz-1) {
              nz1 = E->lmesh.noz;  // only do for the surface 
	    }
	    else {
              nz1 = E->lmesh.noz+1;  // no calculations
	    }
	  }
     }

    // determine the depth at which the source is needed
    mz1 = 1;
    mz2 = E->mesh.noz;
    mz3 = 1;
    if (isource==0) { // only need the source at the surface and CMB as in incompressible case, most of the time
       mz3 = E->mesh.noz-1;
    }

    for (ll=1;ll<=llmax;ll++) // degree 0 term is supposed to be zero at surface, but here we are excluding it for each layer
      for (mm=0;mm<=ll;mm++) {
        p = E->sphere.hindex[ll][mm];
            
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (i=nz1;i<=nz2;i=i+nz3) { // field point
            // compute phi_lm(r)
            phi_lm_cs = 0.0;
            phi_lm_sn = 0.0;

            for (k=mz1;k<=mz2;k=k+mz3) { // source point
              ik = k + (i-1)*E->mesh.noz;
              tmp = E->sphere.tables_rr[ll][ik];
              phi_lm_cs += tmp * sphc_c[k][p];
              phi_lm_sn += tmp * sphs_c[k][p];
            }

	    //compute the potential at grid point for this depth or i
            for (j=1;j<=E->lmesh.nsf;j++) {
                node = (j-1)*E->lmesh.noz+i;
                potential[m][node] += (phi_lm_cs*E->Tbl_cs[m][mm][j] + phi_lm_sn*E->Tbl_sn[m][mm][j])
                            *E->Tbl_lm[m][p][j];          
            }

            if (ll==2 && mm==1)
              if (E->ve_data_cont.polar_wander && icon==1 && icompress && E->parallel.me_loc[3] == E->parallel.nprocz-1 && i==E->lmesh.noz){
                  M[0] = M_conv * phi_lm_cs * pow(E->sphere.R[E->mesh.levmax][i],3.0);
                  M[1] = M_conv * phi_lm_sn * pow(E->sphere.R[E->mesh.levmax][i],3.0);
                    E->ve_data_cont.PW_incr[0] = M[0];
                    E->ve_data_cont.PW_incr[1] = M[1];
            }
        } // end of (i=1;i<=E->lmesh.noz;i++) and m loop
      }  // end of ll and mm
	 //
return;
}
