#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include "Grid.h"

//*********************************************************************************//
//                                                                                 //
// originally taken from Probe Particle Model:                                     //
// https://github.com/ProkopHapala/ProbeParticleModel/blob/master/pyProbeParticle  //
//                                                                                 //
//*********************************************************************************//

// ================= INNER FUNCTIONS

 inline double sqr(double x){
     return x*x;
 }

 void copy_outer( int n_add, int * dims_small, double * V){
     /*
     copy the inner parts of the original (inner) V to the outer parts.
     */
     int sx = dims_small[0];
     int sy = dims_small[1];
     int sz = dims_small[2];
     int lx = dims_small[0]+2*n_add;
     int ly = dims_small[1];
     int lz = dims_small[2]+2*n_add;
     int xs = sy*sz ; // size of each x dimension of the 4D xyzgrid
     int ys = sz ;
     int xl = ly*lz ; // size of each x dimention of the 3D grid
     int yl = lz ;
     int inds = 0; // index on the side of the V array
     int indi = 0; // the original, point from where the number will be calculated
     for (int i=n_add; i< sx+n_add; i++){
         for (int j=0; j< ly; j++){
             for (int k=0; k<n_add; k++){
                 inds= i*xl+j*yl+k;
                 indi = inds + sz;
                 V[inds]= V[indi];
             }
             for (int k=sz+n_add; k<lz; k++){
                 inds= i*xl+j*yl+k;
                 indi = inds - sz;
                 V[inds]= V[indi];
             }
         }
     }
     for (int i=0; i<n_add; i++){
         for (int j=0; j< ly; j++){
             for (int k=0; k<lz; k++){
                 inds= i*xl+j*yl+k;
                 indi = inds + sx*xl;
                 V[inds]= V[indi];
             }
         }
     }
     for (int i=n_add+sx; i<lx; i++){
         for (int j=0; j< ly; j++){
             for (int k=0; k<lz; k++){
                 inds= i*xl+j*yl+k;
                 indi = inds - sx*xl;
                 V[inds]= V[indi];
             }
         }
     }

 }
void copy_outer_fast( int n_add, int * dims_small, double * V){
    /*
     copy the inner parts of the original (inner) V to the outer parts.
     */
    int sx = dims_small[0];
    int sy = dims_small[1];
    int sz = dims_small[2];
    int lx = dims_small[0]+2*n_add;
    int ly = dims_small[1];
    int lz = dims_small[2]+2*n_add;
    int xs = sy*sz ; // size of each x dimension of the 4D xyzgrid
    int ys = sz ;
    int xl = ly*lz ; // size of each x dimention of the 3D grid
    int yl = lz ;
    int inds = 0; // index on the side of the V array
    int indi = 0; // the original, point from where the number will be calculated
    for (int i=n_add; i< sx+n_add; i++){
        for (int j=0; j< ly; j++){
            for (int k=0; k<n_add; k++){  // both - front and back sides copies are combined together for speeding up
                inds= i*xl+j*yl+k;
                indi = inds + sz;
                V[inds]= V[indi];
            //}for (int k=sz+n_add; k<lz; k++){
                inds= i*xl+j*yl+k+sz+n_add;
                indi = inds - sz;
                V[inds]= V[indi];
            }
        }
    }
    for (int i=0; i<n_add; i++){
    //for (int i=n_add+sx; i<lx; i++){
        for (int j=0; j< ly; j++){
            for (int k=0; k<lz; k++){ // both - front and back sides copies are combined together for speeding up
                inds= i*xl+j*yl+k;
                indi = inds + sx*xl;
                V[inds]= V[indi];
                inds= (i+sx+n_add)*xl+j*yl+k;
                indi = inds - sx*xl;
                V[inds]= V[indi];
            }
        }
    }
    
}

double one_step_w_control(double * Vout, double * V, int * optind, int * dims_small, const int n_add, const int aidx, const int st, const int yl, const int xl){
    double diff = -1;
    double tmpd ;
    for (int i=1; i<aidx; i++){
        Vout[optind[i]] = 0.166666666667*(V[optind[i]+st]     +   V[optind[i]-st] + //z
                                          V[optind[i]+st*yl]  +   V[optind[i]-st*yl] +  //y
                                          V[optind[i]+st*xl]  +   V[optind[i]-st*xl]); //x
        tmpd = abs(Vout[optind[i]]-V[optind[i]]);
        diff = (tmpd > diff) ? tmpd:diff;
        //printf("debug: optind[i] %d, Vout[] %4.5f, V[] %4.5f, tmpd %4.5f, diff %4.5f; \n", optind[i],Vout[optind[i]],V[optind[i]],tmpd,diff );
    }
    for (int i=1; i<aidx; i++){
        V[optind[i]] = Vout[optind[i]];
    }
    // printf ("; diff %4.5f",diff);
    copy_outer_fast( n_add, dims_small, V);
    return diff;
}

void m_step (double * Vout, double * V, int * optind, int * dims_small, const int n_add, const int aidx, const int st, const int yl, const int xl, const int m_step){
    for (int instep; instep<m_step; instep++ ){
        for (int i=1; i<aidx; i++){
            Vout[optind[i]] = 0.166666666667*(V[optind[i]+st]     +   V[optind[i]-st] + //z
                                              V[optind[i]+st*yl]  +   V[optind[i]-st*yl] +  //y
                                              V[optind[i]+st*xl]  +   V[optind[i]-st*xl]); //x
            //printf("debug: optind[i] %d, Vout[] %4.5f, V[] %4.5f, tmpd %4.5f, diff %4.5f; \n", optind[i],Vout[optind[i]],V[optind[i]],tmpd,diff );
        }
        for (int i=1; i<aidx; i++){
            V[optind[i]] = Vout[optind[i]];
        }
        copy_outer_fast( n_add, dims_small, V);
    }
}
// Outer functions


extern "C" {
    void create_xyzGrid(double *xyzgrid, int * dims, double *vox_vec, double *g_null){
        // will create a 4 D grid [ix,iy,iz,{x,y,z}] of x, y, z, positions of every voxel in the given dimensions (has to be controlled before!!!) )]
        // vox_vec is 3x3 matrix with (in lines) vectors of a single voxel; g_null is a 3x vector of the original point (normally 0,0,0)
        int nx=dims[0];
        int ny=dims[1];
        int nz=dims[2];
        int xs = dims[1]*dims[2]*3 ; // size of each x dimension
        int ys = dims[2]*3 ;
        int zs = 3 ;
        for  (int i=0; i<nx; i++){
            for (int j=0; j<ny; j++){   
                for (int k=0; k<nz; k++){
                    xyzgrid[i*xs+j*ys+k*zs+0]=g_null[0]+i*vox_vec[0]+j*vox_vec[3]+k*vox_vec[6]; //x coordinate
                    xyzgrid[i*xs+j*ys+k*zs+1]=g_null[1]+i*vox_vec[1]+j*vox_vec[4]+k*vox_vec[7]; //y coordinate
                    xyzgrid[i*xs+j*ys+k*zs+2]=g_null[2]+i*vox_vec[2]+j*vox_vec[5]+k*vox_vec[8]; //z coordinate
                }
            }
        }
        
    }
    
    void fix_all(int * fixtop, int * fixbot, int *fixnum, double * xyzgrid, int * dims, int * nat, double * pos, int * top_at, int * top_fix_at, int * bot_at, int *  bot_fix_at){
        // out: fixtop 1D integer array (full of zeros), for storring those indexes (of the 3D V array, which will be fixed.
        // out: fixbot 1D integer array (full of zeros), for storring those indexes (of the 3D V array, which will be fixed.
        // out: fixnum 1D integer of -- 1 length of fixtop, 1 length of fixbot, then for each x and y min /j/ (bottom) and max /j/ (top)for the voltage ramp preparations
        double max = 1E6; double min = -1E6;
        int nta  = nat[0];
        int nba  = nat[1];
        int ntfa = nat[2];
        int nbfa = nat[3];
        for(int i=0; i<ntfa; i++){ // figuring out the lowest fixed atom (beware z is y )
            //printf("debug: z_pos - tip fixed %4.5f %d \n", pos[top_fix_at[i]*3+1],  top_fix_at[i] );
            if (pos[top_fix_at[i]*3+1] < max ){max = pos[top_fix_at[i]*3+1];}
        }
        for(int i=0; i<nbfa; i++){ // figuring out the highest fixed atom (beware z is y )
            //printf("debug: z_pos - bot fixed %4.5f %d \n", pos[bot_fix_at[i]*3+1] , bot_fix_at[i]  );
            if (pos[bot_fix_at[i]*3+1] > min ){min = pos[bot_fix_at[i]*3+1];}
        }
        printf("\n");
        printf("max (top): %20.20lf \n",max);
        printf("min (bot): %20.20lf \n \n",min);
        printf("-- pretty long fixing loop -- \n");
        int nx=dims[0];
        int ny=dims[1];
        int nz=dims[2];        
        int xs = dims[1]*dims[2]*3 ; // size of each x dimension of the 4D xyzgrid
        int ys = dims[2]*3 ;
        int zs = 3 ;
        int xf = dims[1]*dims[2] ; // size of each x dimention of the 3D grid
        int yf = dims[2] ; 
        int ifixt = 0; int ifixb = 0;
        int ifixn = 2; // bacause ifixt and ifixn are stored in the first 2 numbers of the array
        //int yfixt = 10000; int yfixb =-1;
        int tmpl = 0 ; int tmps = 0; //tmp indexes in the 4D (l) and 3D (s) grid
        for  (int i=0; i<nx; i++){
            for (int j=0; j<ny; j++){   
                for (int k=0; k<nz; k++){
                    tmpl=i*xs+j*ys+k*zs;
                    tmps=i*xf+j*yf+k;
                    ifixn=i*nz*2+k*2+2; //tmpi = (ifixn > tmpi) ? ifixn:tmpi;
                    if (xyzgrid[tmpl+1] < min){ // lowest z (y in this way) of the fixed atom
                        fixbot[ifixb] = tmps ; ifixb += 1 ; fixnum[ifixn] = (j > fixnum[ifixn]) ? j:fixnum[ifixn] ; //yfixb = (j > yfixb) ? j: yfixb ;
                        //printf("debug: fixb; %d %d %d %d %4.5f %4.5f %d %d %d %d \n", i,j,k,tmpl+1,xyzgrid[tmpl+1],min,fixnum[ifixn],ifixn,fixbot[ifixb],ifixb);
                    }else{if (xyzgrid[tmpl+1] > max){ // lowest z (y in this way) of the fixed atom
                        fixtop[ifixt] = tmps ; ifixt += 1 ; fixnum[ifixn+1] = (j < fixnum[ifixn+1]) ? j:fixnum[ifixn+1] ; // yfixt = (j < yfixt) ? j: yfixt ;
                        //printf("debug: fixt; %d %d %d %d %4.5f %4.5f %d %d %d %d \n", i,j,k,tmpl+1,xyzgrid[tmpl+1],min,fixnum[ifixn],ifixn,fixtop[ifixt],ifixt);
                    }else{for (int ia=0; ia<nba ; ia++){ // atoms of the bottom
                        if ( ( sqr(xyzgrid[tmpl]-pos[bot_at[ia]*3+0]) + sqr(xyzgrid[tmpl+1]-pos[bot_at[ia]*3+1]) + sqr(xyzgrid[tmpl+2]-pos[bot_at[ia]*3+2]) ) < 1  ){
                        //        (       x      -  ix[ia]        )**2 +    (        y     -    iy[ia]      )**2 +  (        z        -    iz[ia]      )**2     < 1**2 // No SQRT needed ..//
                            fixbot[ifixb] = tmps ; ifixb += 1 ; fixnum[ifixn] = (j > fixnum[ifixn]) ? j:fixnum[ifixn] ; break; // break - not to have the same index twice
                        };
                    };    for (int ia=0; ia<nta ; ia++){ // atoms of the tip
                        if ( ( sqr(xyzgrid[tmpl]-pos[top_at[ia]*3+0]) + sqr(xyzgrid[tmpl+1]-pos[top_at[ia]*3+1]) + sqr(xyzgrid[tmpl+2]-pos[top_at[ia]*3+2]) ) < 1  ){
                        //        (       x     -  ix[ia]        )**2 +    (        y      -    iy[ia]      )**2 +  (        z        -    iz[ia]      )**2     < 1**2 // No SQRT needed ..//
                            fixtop[ifixt] = tmps ; ifixt += 1 ; fixnum[ifixn+1] = (j < fixnum[ifixn+1]) ? j:fixnum[ifixn+1] ; break;
                        }
                    }}}
                }
            }
        }
        /*for  (int i=0; i<nx; i++){
            for (int k=0; k<nz; k++){
                ifixn=i*nz*2+k*2+2; 
                printf("i, k, ifixn, fixnum[ifixn] %d %d %d %d \n",i,k,ifixn,fixnum[ifixn]);
            }
        }*/
        //printf("ifix (top): %d \n ",ifixt);
        //printf("ifix (bot): %d \n ",ifixb);
        fixnum[0]=ifixt; fixnum[1]=ifixb; //fixnum[2]=yfixt; fixnum[3]=yfixb;
        printf("ifixn %d \n",ifixn);
        //return ifix ;
    }

    void prepare_V( int n_add, int * dims_small, int * fixtop, int * fixbot, int *fixnum, int *optind, double * V, double V_tip){
        //note: fixnum has first 2 numbers lenght of fixtop and fixbot
        //printf("debug: inside prepare_V \n");
        //printf("debug: %d %d %d %d \n", n_add, dims_small[0], dims_small[1], dims_small[2] );
        int sx = dims_small[0];
        int sy = dims_small[1];
        int sz = dims_small[2];
        int lx = dims_small[0]+2*n_add;
        int ly = dims_small[1];
        int lz = dims_small[2]+2*n_add;
        int xs = sy*sz ; // size of each x dimension of the 4D xyzgrid
        int ys = sz ;
        int xl = ly*lz ; // size of each x dimention of the 3D grid
        int yl = lz ; 
        int ifixt = 0; int ifixb = 0;
        int ifixn = 2; // first 2 numbers are ifixt and ifixb;
        int oidx  = 1; //we left oidx 0 for the amount of numbers i optind
        int lidx  = 0; int sidx  = 0; //indexes in the enlarged (l) -- V --  and orignal-sized (s) field
        double Vtmp=V_tip;
        //printf ("debug: sx, sy, sz, lx, ly, lz: %d %d %d %d %d %d \n", sx, sy, sz, lx, ly, lz);
        //printf ("debug: xs, ys, xl, yl,  n_add: %d %d %d %d    %d \n",xs, ys, xl, yl, n_add );
        //printf ("debug: V_tip, Vtmp: %d %d \n", V_tip, Vtmp);
        for  (int i=n_add; i<sx+n_add; i++){ // x - only the inner cell ;
            for (int j=0; j<sy; j++){ // y -- because   no need for bigger size in y ;
                for (int k=n_add; k<sz+n_add; k++){ // z - only the inner cell ;
                    sidx = (i-n_add)*xs+j*ys+(k-n_add);
                    lidx = i*xl+j*yl+k;
                    ifixn= (i-n_add)*sz*2+(k-n_add)*2+2; //tmpi = (ifixn > tmpi) ? ifixn:tmpi;
                    Vtmp = (j < fixnum[ifixn]) ? -0.5*V_tip:-0.5*V_tip+V_tip*(j-fixnum[ifixn])/(fixnum[ifixn+1]-fixnum[ifixn]);
                    
                    if( sidx == fixbot[ifixb]){
                        V[ lidx ]=-0.5*V_tip; ifixb +=1;
                    } else { if(sidx == fixtop[ifixt]){
                        V[ lidx ]=+0.5*V_tip; ifixt +=1;
                    } else {
                        optind[oidx]=lidx; oidx +=1;
                        V[ lidx ]= (j > fixnum[ifixn+1]) ? +0.5*V_tip:Vtmp;
                    }}
                    //printf ("debug: i,j,k,sidx,lidx,ifixn,fixnum[ifixn], V[lidx]: %d %d %d %d %d %d %d %4.5f \n", i,j,k,sidx,lidx,ifixn, fixnum[ifixn], V[lidx]);
                }
            }
        }
        optind[0]=oidx;
        //printf("ifixt,ifixb,fixnum[0],fixnum[1]: %d %d %d %d \n",ifixt,ifixb,fixnum[0],fixnum[1]);
        //printf("optind[0..3], %d %d %d %d \n",optind[0],optind[1],optind[2],optind[3]);
        printf("end of prepare_V \n");
        //printf("fixtop[ifixt-1],fixbot[ifixb-1]: %d %d \n",fixtop[ifixt-1],fixbot[ifixb-1]);
        //printf("fixtop[ifixt],fixbot[ifixb],tmpmax: %d %d %d \n",fixtop[ifixt],fixbot[ifixb],tmpmax);
        //printf("fixtop[ifixt+1],fixbot[ifixb+1]: %d %d \n",fixtop[ifixt+1],fixbot[ifixb+1]);
        printf("going to copy outer \n");
        copy_outer_fast( n_add, dims_small, V);
        printf("end of c++ \n");
    }

    void opt_V( int n_add, int * dims_small, int *optind, double * V, double * Vout , double prec, int precond, int innerstep, int idx){ // prec - precission
        printf("inside opt_V \n");
        int sx = dims_small[0];
        int sy = dims_small[1];
        int sz = dims_small[2];
        int lx = dims_small[0]+2*n_add;
        int ly = dims_small[1];
        int lz = dims_small[2]+2*n_add;
        int xs = sy*sz ; // size of each x dimension of the 4D xyzgrid
        int ys = sz ;
        int xl = ly*lz ; // size of each x dimention of the 3D grid
        int yl = lz ;
        int aidx = optind[0]; // amount of indexes
        //printf("after making dims \n"); -- some seg-faults if field too small
        double diff = -1.0; // printf("after diff \n");
        int st = 1; //printf("after st \n");
        int maxstep = 100000; //printf("after maxstep \n"); // outer step
        //printf("debug: prec %4.10f prec \n");
        printf("Go into opt: %d \n",idx);
        if (precond == 1){
            printf("going to preconditioner");
            for (int ist=n_add; ist>0; ist--){
                m_step( Vout, V, optind, dims_small, n_add, aidx, ist, yl, xl, innerstep);
                m_step( Vout, V, optind, dims_small, n_add, aidx, st, yl, xl, innerstep);
            }
        }
        for (int istep=0; istep<maxstep; istep++){//while (diff > prec){
            printf("INDEX: %d ; OptStep: %d ",idx, istep);
            m_step( Vout, V, optind, dims_small, n_add, aidx, st, yl, xl, innerstep);
            diff = one_step_w_control( Vout, V, optind, dims_small, n_add, aidx, st, yl, xl);
            /*
            for (int i=1; i<aidx; i++){
                Vout[optind[i]] = 0.166666666667*(V[optind[i]+st]     +   V[optind[i]-st] + //z
                                                  V[optind[i]+st*yl]  +   V[optind[i]-st*yl] +  //y
                                                  V[optind[i]+st*xl]  +   V[optind[i]-st*xl]); //x
                tmpd = abs(Vout[optind[i]]-V[optind[i]]);
                diff = (tmpd > diff) ? tmpd:diff;
                //printf("debug: optind[i] %d, Vout[] %4.5f, V[] %4.5f, tmpd %4.5f, diff %4.5f; \n", optind[i],Vout[optind[i]],V[optind[i]],tmpd,diff );
            }
            for (int i=1; i<aidx; i++){
                V[optind[i]] = Vout[optind[i]];
            }
            copy_outer_fast( n_add, dims_small, V);
            */
            printf ("; diff %4.10f vs. prec %4.10f \n",diff, prec);
            if (diff < prec){break;}
            
        }
        //free(Vout);
        printf("\n");
    }

    void writeCube( int n_at, int * mol_Z, double * mol_xyz, double * grid_origin, double * grid_vec , int * dims, double * data, char *fname){
        // number of atoms ; Z - atoms n_at vector; x,y,z - atoms n_at x 3; 3x vector         ; 3x3 matrix       ;  3x vector ; nx*ny*nz - one dimesional data - already rescalled in python; every distances already in bohrs from python
        FILE *f;
        int nx = dims[0]; int ny = dims[1]; int nz = dims[2];
        int iv = 0 ; // voxel index
        f=fopen(fname,"w");
        printf("starting to write the cube file %s \n", fname);
        fprintf(f,"Electric field for the KPFM simulations;\n Comment here \n");
        fprintf(f,"%d %4.8f %4.8f %4.8f \n",n_at,grid_origin[0],grid_origin[1],grid_origin[2]);
        for (int i=0; i<3; i++ ){
            fprintf(f,"%d %4.8f %4.8f %4.8f \n",dims[i],grid_vec[0+3*i],grid_vec[1+3*i],grid_vec[2+3*i]);
            }
        for (int i=0; i<n_at; i++ ){
            fprintf(f,"%d 0.0 %4.8f %4.8f %4.8f \n", mol_Z[i],mol_xyz[0+3*i],mol_xyz[1+3*i],mol_xyz[2+3*i]);
            }
        for (int i=0; i<nx; i++ ){
            for (int j=0; j<ny; j++ ){
                for (int k=0; k<nz; k++ ){
                    fprintf(f,"%4.10f ",data[iv]);
                    iv += 1;
                    if ((iv % 6) == 0){fprintf(f,"\n");}
                    }
                // fprintf(f,"\n"); // !! This is not in CP2K //
                }
            }
        fclose(f);
        printf("debug: nx*ny*nz: %d; iv: %d. \n", nx*ny*nz, iv);
    }
    
}

/* // original code left here for inspiration
GridShape gridShape;

extern "C" {

    int ReadNumsUpTo_C (char *fname, double *numbers, int * dims, int noline) {
        FILE *f;
        char line[5000]; // define a length which is long enough to store a line
        char *waste;
        int waste2;
        long i=0, j=0, k=0, tot=0; 
        int nx=dims[0];
        int ny=dims[1];
        int nz=dims[2];
        printf ("FileRead program: reading %s file\n", fname);
        printf ("XYZ dimensions are %d %d %d\n", dims[0], dims[1], dims[2]);
        f=fopen(fname, "r");
        if (f==NULL)        {
            fprintf(stderr, "Can't open the file %s", fname);
            exit (1); 
        }
        for (i=0; i<noline; i++) {   
            waste=fgets(line,5000, f);
        }
//       printf ("Line: %s", line);
        for  (tot=0, k=0; k<dims[2]; k++){
            for (j=0; j<dims[1]; j++){   
                for (i=0; i<dims[0]; i++){
                    waste2=fscanf(f,"%lf",&numbers[tot]);
//                    printf ("%20.20lf ", numbers[tot]);
                    tot++;
//                    if (tot > 5 ) exit(1);
                }
            }
        }
//       printf ("%lf %lf %lf %lf %lf\n", numbers[tot-1], numbers[tot-2], numbers[tot-3], numbers[tot-4], numbers[tot-5]);
        printf("Reading DONE\n");
        fclose(f);
        return 0;
    }

	void interpolate_gridCoord( int n, Vec3d * pos_list, double * data, double * out ){
		for( int i=0; i<n; i++ ){ 
			out[i] = interpolate3DWrap( data, gridShape.n, pos_list[i] ); 
		}
	}

	void interpolateLine_gridCoord( int n, Vec3d * p0, Vec3d * p1, double * data, double * out ){
		Vec3d dp,p; 
		dp.set_sub( *p1, *p0 );
		dp.mul( 1.0d/n );
		p.set( *p0 );
		for( int i=0; i<n; i++ ){ 
			//printf( " i, n  %i %i  pi0  %f %f %f  \n", i, n,  p.x,p.y,p.z );
			out[i] = interpolate3DWrap( data, gridShape.n, p );
			p.add( dp );
		}
	}

	void interpolateQuad_gridCoord( int * nij, Vec3d * p00, Vec3d * p01, Vec3d * p10, Vec3d * p11, double * data, double * out ){
		int ni = nij[0];
		int nj = nij[1];
		Vec3d dpi0,dpi1,pi0,pi1;
		dpi0.set_sub( *p10, *p00 ); dpi0.mul( 1.0d/ni );
		dpi1.set_sub( *p11, *p01 ); dpi1.mul( 1.0d/nj );
		pi0.set( *p00 );
		pi1.set( *p01 );
		for( int i=0; i<ni; i++ ){ 
			//printf( " i, ni, nj   %i %i %i  pi0  %f %f %f   pi1  %f %f %f \n", i, ni, nj,     pi0.x,pi0.y,pi0.z,     pi1.x,pi1.y,pi1.z );
			interpolateLine_gridCoord( nj, &pi0, &pi1, data, out + ( i * nj ) );
			pi0.add( dpi0 );
			pi1.add( dpi1 );
		}
	}

	void interpolate_cartesian( int n, Vec3d * pos_list, double * data, double * out ){
		for( int i=0; i<n; i++ ){ 
			Vec3d gpos;
			gridShape.cartesian2grid( pos_list[i], gpos );
			out[i] = interpolate3DWrap( data, gridShape.n, gpos ); 
		}
	}

	void setGridN( int * n ){
		//gridShape.n.set( *(Vec3i*)n );
		gridShape.n.set( n[2], n[1], n[0] );
		printf( " nxyz  %i %i %i \n", gridShape.n.x, gridShape.n.y, gridShape.n.z );
	}

	void setGridCell( double * cell ){
		gridShape.setCell( *(Mat3d*)cell );
		printf( " a     %f %f %f \n", gridShape. dCell.a.x,  gridShape. dCell.a.y,  gridShape. dCell.a.z );
		printf( " b     %f %f %f \n", gridShape. dCell.b.x,  gridShape. dCell.b.y,  gridShape. dCell.b.z );
		printf( " c     %f %f %f \n", gridShape. dCell.c.x,  gridShape. dCell.c.y,  gridShape. dCell.c.z );
		printf( " inv_a %f %f %f \n", gridShape.diCell.a.x,  gridShape. diCell.a.y, gridShape.diCell.a.z );
		printf( " inv_b %f %f %f \n", gridShape.diCell.b.x,  gridShape.diCell.b.y,  gridShape.diCell.b.z );
		printf( " inv_c %f %f %f \n", gridShape.diCell.c.x,  gridShape.diCell.c.y,  gridShape.diCell.c.z );
	}

}
*/
