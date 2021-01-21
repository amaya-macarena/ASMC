/*---------------------------------------------------------------------------*/
/* TIME_2D:                                                                  */
/*        FINITE-DIFFERENCES COMPUTATION OF FIRST TRAVEL TIMES IN 2D.        */
/*        P.Podvin <Pascal.Podvin@ensmp.fr>                                  */
/*        Centre de Recherche en Geophysique, Ecole des Mines de Paris,      */
/*        Fontainebleau, France.                                             */
/*                                                                           */
/*        The first version was released in October 1990.                    */
/*        A slightly better version was released 11 March 1999               */
/*        A "robustified and ANSIfied" version was dated 11 July 2003        */
/*        An "un-Sun-ified" version was dated 10 November 2003               */
/*        Corrections in initialization procedures were introduced in a      */
/*        version dated 16 October 2006.                                     */
/*        This version is dated 9 January 2007. An inessential upgrade       */
/*        was introduced to correct very rare errors that might occur at     */
/*        isolated nodes in models with very strong 1-cell-wide slow         */
/*        anomalies. The upgrade involves a negligible overhead and will     */
/*        change nothing in the results in almost all cases...               */
/*                                                                           */
/* PROTOTYPE : (see ../include/fdtimes.h)                                    */
/*                                                                           */
/* int time_2d(double *hs, double *t, int nx, int ny,                          */
/*             double xs, double ys, double eps_init, int messages);            */
/*                                                                           */
/* ARGUMENTS (C description; all FORTRAN arguments are pointers)             */
/*                                                                           */
/*       (int)     nx,ny           : dimensions of the time field (number of */
/*                                   grid points). Cells are square. No      */
/*                                   dimension may be lower than 2.          */
/*                                                                           */
/*       (double *) hs,t            : 1D arrays of nx*ny elements.            */
/*                                   t will contain the computed time field  */
/*                                   arranged as a succession of columns     */
/*                                   (x=ct,y=0,ny-1). hs contains model      */
/*                                   data (slowness*mesh spacing) organized  */
/*                                   in the same manner. Within the routine, */
/*                                   values stored in hs are implicitly      */
/*                                   ascribed to cell centers, whereas times */
/*                                   computed and stored in t refer to cell  */
/*                                   corners (grid points). Cells located at */
/*                                   coordinates x=nx-1 or y=ny-1 are thus   */
/*                                   dummy cells (out of the real model).    */
/*                                   The corresponding values will be        */
/*                                   ignored and treated as "infinite".      */
/*                              Note:                                        */
/*                                   Values stored in hs may not be negative */
/*                                   and must not exceed 0.499e+19 (see      */
/*                                   macro FD_HUGE defined below).           */
/*                                                                           */
/*       (double)   xs,ys         :   Point source coordinates referred to    */
/*                                   "upper left" corner of model (grid      */
/*                                   point with lowest indices), expressed   */
/*                                   in mesh spacing (h) unit.               */
/*                                   Licit range [0.0,nx-1.0][0.0,ny-1.0].   */
/*                                   If source point is found to be out      */
/*                                   of the licit range, the timefield is    */
/*                                   treated as already initialized in the   */
/*                                   calling program (extended source).      */
/*                              Note: although this is not required when     */
/*                                    source is licit, you should always     */
/*                                    initialize array t with whatever       */
/*                                    constant value (why not 0?), so that   */
/*                                    you'll be warned if a wrong (illicit)  */
/*                                    source location is entered.            */
/*                                                                           */
/*       (double)   eps_init      :   tolerance on relative inhomogeneity     */
/*                                   used (only) during initialization.      */
/*                                   Although licit values are in [0,1],     */
/*                                   relevant values should be <0.01.        */
/*                                   0.001 is a reasonable choice.           */
/*                                   (arg. ignored if source is illicit)     */
/*                                                                           */
/*       (int)     messages        : 0 makes the routine silent (except on   */
/*                                   diagnosed error); 1: small talk mode;   */
/*                                   >1: debug mode (verbose).               */
/*                                   A negative value is interpreted as 1.   */
/*                                                                           */
/* VALUE :                                                                   */
/*                                                                           */
/*       time_2d() returns a nonzero value if an error was detected.         */
/*       An explicit error message is printed on 'stderr'.                   */
/*                                                                           */
/* CALLING TIME_2D FROM A PROGRAM WRITTEN IN FORTRAN :                       */
/*                                                                           */
/*       The C-function time_2d_() is provided as the interface to Fortran.  */
/*       In this routine, all arguments are pointers as required by Fortran  */
/*       (where routine arguments are passed by address, not by value), and  */
/*       dimensions 'x' and 'y' are swapped to mimic standard Fortran memory */
/*       mapping (hs[i][j] in C "means" HS(J,I) in Fortran).                 */
/*       Normally, calling TIME_2D (without the trailing underscore) in      */
/*       Fortran will automatically result in a call to time_2d_() in C.     */
/*    Compiling :                                                            */
/*       With Sun compilers, nothing particular is required to obtain this   */
/*       automatic behaviour.                                                */
/*       With the GNU compilers (gcc, g77), you will need to compile your    */
/*       Fortran program with option -fno-second-underscore to obtain this.  */ 
/*       Because the C/FORTRAN interface is not standardized, you might      */
/*       have to browse your documentation on other platforms.               */
/*       If you get into trouble :                                           */
/*       Compiling this program with the option -DDEBUG_ARGS allows to check */
/*       what C-function is actually called (time_2d or time_2d_) and how    */
/*       its arguments are interpreted.                                      */
/*     Declarations in the calling program :                                 */
/*       As seen from Fortran, TIME_2D is an INTEGER FUNCTION.               */
/*       It should thus be declared as such in the calling program.          */
/*       Please note that arguments passed to TIME_2D MUST be declared with  */
/*       their correct types (e.g., XS,YS MUST NOT be declared INTEGER,      */
/*       while NX,NY MUST NOT be declared REAL !)                            */
/*       Not conforming to this may result into incomprehensible situations  */
/*       (once again, compile option -DDEBUG_ARGS may help...)               */
/*       All scalar arguments are read-only and may thus be declared         */
/*       constant (using PARAMETER statements).                              */
/*     Program skeleton :                                                    */
/*       INTEGER NX,NY                                                       */
/*       REAL HS(NX,NY),T(NX,NY)                                             */
/* C or  REAL HS(NX*NY),T(NX*NY)                                             */
/*       REAL XS,YS,EPS_INIT                                                 */
/*       INTEGER MESSAGES,TIME_2D,STATUS                                     */
/*       .....                                                               */
/*       STATUS=TIME_2D(HS,T,NX,NY,XS,YS,EPS_INIT,MESSAGES)                  */
/*       IF(STATUS.NE.0)                                                     */
/*         STOP "time_2d diagnosed a (presumably fatal) error"               */
/*       ENDIF                                                               */
/*       .....                                                               */
/*       STOP                                                                */
/*       END                                                                 */
/*                                                                           */
/* RECENT UPDATES :                                                          */
/*                                                                           */
/*        Although time_2d has been used by dozens of people worldwide       */
/*        for many years, some changes were recently introduced (2003).      */
/*        These changes do not affect the routine's usage (except for        */
/*        some values returned on error, and the treatment of "infinite"     */
/*        slowness values, now forbidden in hs). Their only justification    */
/*        is improved portability.                                           */
/*        The changes are the following :                                    */
/*        - The program now conforms to ANSI-C (you must include fdtimes.h   */
/*          in your calling program, if it is written in C).                 */
/*        - I decided to drop Sun-specific calls to routines implementing    */
/*          the IEEE standards for infinity and related tests. This is non   */
/*          portable, and would even create problems on Sun platforms when   */
/*          the calling program is in Fortran !                              */
/*          As a consequence, a finite threshold is now used to test whether */
/*          a double is treated as infinite. No value in array hs is allowed  */
/*          to exceed this threshold (see macro FD_HUGE below).              */
/*        - Unpredictible crashes were seldom encountered on Intel based PCs */
/*          due to the fact that the routine's behaviour strongly depended   */
/*          on tests comparing doubles in an "exact" manner (e.g., replacing  */
/*          a test x<y by x<=y altered significantly the code behaviour).    */
/*          These tests are now "fuzzified" with no significant consequence  */
/*          on precision (using macro EPS_FUZZY below).                      */
/*                                                                           */
/*        More recently (December 2005), A. Tryggvason & B. Bergman from     */
/*        Uppsala Univ. reported an unexpected dissymmetry in the results    */
/*        in 3D when the source is located at the vicinity of a severe       */
/*        velocity contrast. Revised versions correcting initialization      */
/*        procedures in both the 2D and 3D codes were released 2 Feb 2006.   */
/*        The revisions induce very modest changes in the results when the   */
/*        source-point is located in a heterogeneous region of the model.    */
/*                                                                           */
/*        Finally, as this code is being intensively used again in our lab   */
/*        with very contrasted models (air layer), minor bugs were revealed  */
/*        and corrected at the end of 2006, with no impact whatsoever on     */
/*        results in all "reasonable" cases... Thanks to Cedric Taillandier  */
/*        for very detailed examination of results in the weirdest cases !   */
/*                                                                           */
/* REFERENCE : Podvin & Lecomte, Geophys.J.Intern. 105, 271-284, 1991.       */
/*---------------------------------------------------------------------------*/

#include"time_2d.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef    M_SQRT2
#define    M_SQRT2     1.41421356237309504880
#endif
#define    min(x,y)    ( ((x)<(y))? (x):(y) )
#define    max(x,y)    ( ((x)>(y))? (x):(y) )
#define    NINT(x)     (int)floor((x)+0.5)
/* NINT is "lazy nearest integer", used here only with positive values */

#define    INFINITY    0.500e+19
#define    FD_HUGE     0.499e+19
#define    ISINF(x)    ((x)>FD_HUGE)
/* Note: INFINITY must be lower than the SQUARE ROOT of the highest */
/* acceptable double value (machine dependent) to prevent overflow.  */
/* FD_HUGE must be chosen such that INFINITY>FD_HUGE is true.       */
/* This last condition is actually tested by the program !          */

#define EPS_FUZZY       1.2e-07
/* this is a little higher than 2^-23, the smallest relative */
/* difference between two "normal" doubles encoded on 32 bits */
/* with standard encoding such as IEEE754 (23-bit mantissa). */
/* Please inform me if you find a reason to change this !    */

/*------------------------------------------------Static Functions----------*/
static void
    error(int flag),
    init_nearest(void),
    init_cell(double dx, double dy, int x, int y),
    propagate_point(void),
    send_x_headwave(int from, int to, int y, int y_s_now, int y_s_future),
    send_y_headwave(int from, int to, int x, int x_s_now, int x_s_future);

static int
    pre_init(void),
    init_point(void),
    recursive_init(void),
    x_side(int y,int future,int x_begin,int x_end),
    y_side(int x,int future,int y_begin,int y_end);

static double
    init_cellh(double vh, double vv, double hsc, double hsn);

/*------------------------------------------------Static Declarations-------*/

/* MODEL DESCRIPTION */

static    int       nmesh_x,nmesh_y;        /* Model dimensions            */
static    double    *hs_buf,**hs;            /* 1D and 2D slowness arrays   */
static    double    *hs_keep;                /* to save boundary values     */
static    double     eps_init;

/* TIMEFIELD DESCRIPTION */

static    int       nx,ny;                  /* Timefield dimensions        */
static    double    *t_buf,**t;              /* 1D and 2D time arrays       */
static    double     timeshift;              /* see comments: init_point()  */

/* SOURCE */

static    double     fxs,fys;                /* Point source coordinates    */
static    int       xs,ys;                  /* Nearest node                */

/* STATIC PARAMETERS */

#ifndef INIT_MIN
#define INIT_MIN        10              /* Initialize by re-discretized FD  */
                                        /* if radius of uniform region near */
                                        /* the source < INIT_MIN meshes.    */ 
                                        /* (ADJUSTABLE).                    */
#endif /* !INIT_MIN */
#define N_INIT_X    (4*INIT_MIN+3)
#define N_INIT      N_INIT_X*N_INIT_X

#ifndef INIT_RECURS_LIMIT
#define INIT_RECURS_LIMIT    1          /* Maximum level of recursivity     */
                                        /* during initialization.           */
                                        /* (ADJUSTABLE).                    */
#endif /* !INIT_RECURS_LIMIT */

static int  init_stage=0,               /* recursivity of initialization    */
            source_at_node=0,
            reverse_order=0,            /* recursivity of side computation  */
            current_side_limit,         /* current boundary of computations */
            X0,X1,Y0,Y1,                /* incl. boundaries of computed box */
            messages;                   /* message flag: 0: silent
                                                         1: smalltalk
                                                        >1: verbose.        */


/*----------------------------------------------- ERROR FLAGS ------------*/

#define NO_ERROR      0
#define ERR_INFBUG  (-1)
#define ERR_MULTUNI (-2)
#define ERR_MALLOC  (-3)
#define ERR_RECURS  (-4)
#define ERR_EPS     (-5)
#define ERR_RANGE   (-6)
#define ERR_PHYS    (-7)
#define ERR_DIM     (-8)

static char    *err_msg[]=
            {
            "\ntime_2d: Computations terminated normally.\n",
            "\ntime_2d: [Bug] macros INFINITY, FD_HUGE not properly set.\n",
            "\ntime_2d: Multiple source but input time map is uniform.\n",
            "\ntime_2d: Not enough virtual memory (malloc failed).\n",
            "\ntime_2d: Fatal error during recursive init.\n",
            "\ntime_2d: [Init] Illegal tolerance on inhomogeneity.\n",
            "\ntime_2d: Illegal ('infinite') value(s) in array hs.\n",
            "\ntime_2d: Non-physical negative value(s) in array hs.\n",
            "\ntime_2d: ciao dimension (nx,ny) is too small or negative.\n",
            };

/*---------------------------------------------------Error()---------------*/

static void error(int flag)
{    
    if(messages || flag) {
      fflush(stdout);
      fprintf(stderr,"%s",err_msg[-flag]);
      fflush(stderr);
    }
}

#define VERBOSE messages>1


/*-------------------------------------------------------TIME_2D()---------*/

int time_2d(double *HS, double *T, int NX, int NY,
            double XS, double YS, double EPS_INIT, int MESSAGES)

{  
  /*    printf("NX = %d NY= %d\n", NX, NY); */
  /*    printf("sourcex= %2f sourcez = %2f\n", XS, YS);  */
    
    int    i,j,signal;
    double  *pf;

#ifdef DEBUG_ARGS
    fprintf(stderr,"\n******** time_2d: Option DEBUG_ARGS is on.\n");
    if(init_stage) fprintf(stderr,"Recursively entering ");
    else fprintf(stderr,"Initially entering ");
    fprintf(stderr,"time_2d() in C-style, using `time_2d'.\n");
    fprintf(stderr,"Arrays dimensions: nx=%d ny=%d\n",NX,NY);
    fprintf(stderr,"Args HS, T are arrays[nx][ny], i.e. arrays[%d][%d]\n",
                    NX,NY);
    fprintf(stderr,"Licit src coordinate ranges: xs in [0.,%d] ys in [0.,%d]\n",
                    NX-1,NY-1);
    fprintf(stderr,"Args XS, YS (actual src coordinates) : xs=%g ys=%g\n",
                    XS,YS);
    fprintf(stderr,"Other Args: EPS_INIT=%g MESSAGES=%d\n",
                    EPS_INIT,MESSAGES);
    fprintf(stderr,"First elements of input arrays: HS[0][0]=%g, T[0][0]=%g\n",
                    *HS,*T);
    fprintf(stderr,"******** time_2d: Option DEBUG_ARGS done.\n");
    fflush(stderr);
#endif

/* This section merely copies arguments to internal variables. */
/* This is where you might change things in order to build     */
/* a customized calling interface (with an alternate arglist). */
/* If you choose to do so, you must keep time_2d() unchanged,  */
/* as it is needed internally (recursive init). Design your    */
/* customized interface as another function (having another    */
/* name...) and prototype it in fdtimes.h                      */
    hs_buf=HS;
    t_buf=T;
    nx=NX;
    ny=NY;
    fxs=XS;
    fys=YS;
    eps_init=EPS_INIT;
    if(MESSAGES<0) messages=1;
    else messages=MESSAGES;
/* You should change nothing below this  */

    if((signal=pre_init())!=NO_ERROR){
        error(signal);
        return signal;
    }
    if((signal=init_point())==NO_ERROR) propagate_point();

    if(!init_stage){
        int x,y;
/* FUZZIFIED COMPARISONS: see comments below, in init_point() */
        if(timeshift<0.0)
          for(x=0;x<nx;x++)
            for(y=0;y<ny;y++)
              t[x][y]+=timeshift;
        for(i=0,pf=hs_keep;i<nx;i++,pf++) hs[i][nmesh_y]= *pf;
        for(j=0;j<nmesh_y;j++,pf++) hs[nmesh_x][j]= *pf;
        free((char *)hs_keep);
    }
    free((char *)hs);
    free((char *)t);
    error(signal);
    return signal;
}

/*-------------------------------------------------------TIME_2D_()---------*/
/*---------------------------------------------FORTRAN INTERFACE------------*/
 
// int time_2d_(double *HS, double *T, int *NX, int *NY,
//              double *XS, double *YS, double *EPS_INIT, int *MESSAGES)
// {   
//     int    i,j,signal;
//     double  *pf;
// 
// #ifdef DEBUG_ARGS
//     fprintf(stderr,"\n******** time_2d: Option DEBUG_ARGS is on.\n");
//     if(init_stage) fprintf(stderr,"Recursively entering ");
//     else fprintf(stderr,"Initially entering ");
//     fprintf(stderr,"time_2d() in FORTRAN-style, using `time_2d_'.\n");
//     fprintf(stderr,"Arrays dimensions: nx=%d ny=%d\n",*NX,*NY);
//     fprintf(stderr,"Args HS, T are Fortran-arrays(nx,ny), i.e. arrays(%d,%d)\n",
//                     *NX,*NY);
//     fprintf(stderr,"Licit src coordinate ranges: xs in [0.,%d] ys in [0,%d]\n",
//                     *NX-1,*NY-1);
//     fprintf(stderr,"Args XS, YS (actual src coordinates) : xs=%g ys=%g\n",
//                     *XS,*YS);
//     fprintf(stderr,"Other Args: EPS_INIT=%g MESSAGES=%d\n",
//                     *EPS_INIT,*MESSAGES);
//     fprintf(stderr,"First elements of input arrays: HS(1,1)=%g, T(1,1)=%g\n",
//                     *HS,*T);
//     fprintf(stderr,"******** time_2d: Option DEBUG_ARGS done.\n");
//     fflush(stderr);
// #endif
//  
// /* This section merely copies values pointed to by scalar args */
// /* and the addresses of arrays HS, T to internal variables.    */
// /* This is where you might change things in order to build     */
// /* a customized calling interface (with an alternate arglist). */
// /* As this routine time_2d_ is not needed internally, you may  */
// /* straighforwardly edit it.                                   */
//     hs_buf=HS;
//     t_buf=T;
//     nx= *NY;
//     ny= *NX;
//     fys= *XS;
//     fxs= *YS;
//     eps_init= *EPS_INIT;
//     if(*MESSAGES<0) messages=1;
//     else messages= *MESSAGES;
// /* You should change nothing below this  */
//  
//     if((signal=pre_init())!=NO_ERROR){
//         error(signal);
//         return signal;
//     }
//     if((signal=init_point())==NO_ERROR) propagate_point();
// 
//     if(!init_stage){
//         int x,y;
// /* FUZZIFIED COMPARISONS: see comments below, in init_point() */
//         if(timeshift<0.0)
//           for(x=0;x<nx;x++)
//             for(y=0;y<ny;y++)
//               t[x][y]+=timeshift;
//         for(i=0,pf=hs_keep;i<nx;i++,pf++) hs[i][nmesh_y]= *pf;
//         for(j=0;j<nmesh_y;j++,pf++) hs[nmesh_x][j]= *pf;
//         free((char *)hs_keep);
//     }
//     free((char *)hs);
//     free((char *)t);
//     error(signal);
//     return signal;
// }


/*---------------------------------------------------Pre_init()------------*/

static int pre_init(void)
{    
    int      i,j,errtest;
    double    *pf;

    if(nx<2 || ny<2) return ERR_DIM;
    if(!(ISINF(INFINITY))) return ERR_INFBUG;
/* if you encounter this error, it probably means you played */
/* around with the values of macros FD_HUGE and INFINITY !   */

    nmesh_x=nx-1;
    nmesh_y=ny-1;

/* allocate 2D arrays */
    if(!(hs=(double **)malloc((unsigned)nx*sizeof(double *))))
        return ERR_MALLOC;
    if(!(t=(double **)malloc((unsigned)nx*sizeof(double *)))){
        free((char *)hs);
        return ERR_MALLOC;
    }
    for(i=j=0;i<nx;i++,j+=ny) {
        hs[i]=hs_buf+j;
        t[i]=t_buf+j;
    }

/* nothing more to do if this is a recursive call */
    if(init_stage) return NO_ERROR;

/* check that non-masked hs values are licit */
    errtest=NO_ERROR;
    for(i=0;i<nmesh_x;i++) {
        for(j=0;j<nmesh_y && errtest==NO_ERROR;j++){
          if(ISINF(hs[i][j])) errtest=ERR_RANGE;
          if(hs[i][j]<0.0) errtest=ERR_PHYS;
        }
    }
    if(errtest!=NO_ERROR){
        free((char *)hs);
        free((char *)t);
        return errtest;
    }

/* assign zero velocity borders to the dummy cells */
/* while keeping these masked values in hs_keep[]. */
    if(!(hs_keep=(double *)malloc((unsigned)(nx+nmesh_y)*sizeof(double )))){
        free((char *)t);
        free((char *)hs);
        return ERR_MALLOC;
    }
    for(i=0,pf=hs_keep;i<nx;i++,pf++){
        *pf=hs[i][nmesh_y];
        hs[i][nmesh_y]=INFINITY;
    }
    for(j=0;j<nmesh_y;j++,pf++){
        *pf=hs[nmesh_x][j];
        hs[nmesh_x][j]=INFINITY;
    }

    return NO_ERROR;
}


/*---------------------------------------------------Init_point()----------*/

static int init_point(void)
{
    int     x,y,xsc,ysc,signal,ntime,
            test,test_X0,test_X1,test_Y0,test_Y1;
    double   min_t,max_t,hs0,sq_dist,allowed_delta_hs,*pf;
 
/* test relevance of source position or locate minimum time source point */
    if(fxs>=0.0 && fxs<=nmesh_x && fys>=0.0 && fys<=nmesh_y){
/* valid single source */
/* if first pass, assign infinity to all times */
        ntime=nx*ny;
        if(init_stage==0){
            timeshift=0.0;
            if(eps_init<0.0 || eps_init>1.0){
                error(ERR_EPS);
                return ERR_EPS;
            }
            for(x=0,pf=t_buf;x<ntime;x++) *pf++=INFINITY;
        }
 
        xs=NINT(fxs);
        ys=NINT(fys);
        if(xs==fxs && ys==fys) source_at_node=1;
        if(messages)
          printf("\ntime_2d[Init%d]: Point source[%g,%g]. Nearest node[%d,%d].",
                 init_stage,fxs,fys,xs,ys);
    }
    else {
/* multiple source */
/* patch 16102006[2] : bug: xs,ys remained uninitialized in this loop */
        for(x=0,min_t=max_t=t[0][0],xs=ys=0;x<nx;x++)
            for(y=0;y<ny;y++){
                if(t[x][y]<min_t){
                    min_t=t[x][y];
                    xs=x;
                    ys=y;
                }
                if(t[x][y]>max_t) max_t=t[x][y];
            }
        if(max_t==min_t) return ERR_MULTUNI;
        source_at_node=1;
        X0=X1=xs;
        Y0=Y1=ys;
/* FUZZIFIED COMPARISONS: fuzzy tests take it for granted that times are */
/* non-negative. This is why global variable timeshift was introduced... */
/* If needed, time-shifting must be done only once at init_stage zero,   */
/* and finally undone at return time (see time_2d() and time_2d_()).     */
        if(init_stage==0){
          if(min_t<0.0){
            timeshift=min_t;
            for(x=0;x<nx;x++)
              for(y=0;y<ny;y++)
                t[x][y]-=timeshift;
          }/* shift all times to work with non-negative values */
          else timeshift=0.0;
        }
        if(messages)
            printf("\ntime_2d[Init%d]: Multiple source starting at node[%d,%d],\
 at time %g.",init_stage,xs,ys,min_t);
        return NO_ERROR;
    }

/* model properties at source */
    if(source_at_node) {
	xsc=(xs==nmesh_x)? xs-1:xs; 
	ysc=(ys==nmesh_y)? ys-1:ys; 
    }
    else{
        xsc=(fxs<xs) ? xs-1:xs;
        ysc=(fys<ys) ? ys-1:ys;
    }
    hs0=hs[xsc][ysc];

#ifdef SIMPLISTIC_INIT
    init_nearest();
    X0=X1=xs;
    Y0=Y1=ys;
    return NO_ERROR;
#else

/* search for the largest square box with constant slowness */
/* patch 281205[2] : initial boundaries set incorrectly
 *                   many changes in the do-loop !
 *  X0=max(xs-1,0);
 *  X1=min(xs+1,nmesh_x-1);
 *  Y0=max(ys-1,0);
 *  Y1=min(ys+1,nmesh_y-1);
 *
 *                   Initially, X0,X1,Y0,Y1 are inclusive boundaries of the 
 *                   indices of cells with quasi-identical slowness.    
 */
    allowed_delta_hs=hs0*eps_init;
    test_X0=test_X1=test_Y0=test_Y1=0;
    X0=X1=xsc;
    Y0=Y1=ysc;
    do{
        test=0;
        if(X0 && !test_X0){
            test++;
            x=--X0;
            for(y=Y0;y<=Y1 && !test_X0;y++)
                if(fabs(hs[x][y]-hs0)>allowed_delta_hs) test_X0=1;
            if(test_X0) X0++;
        }
        if(Y0 && !test_Y0){
            test++;
            y=--Y0;
            for(x=X0;x<=X1 && !test_Y0;x++)
                if(fabs(hs[x][y]-hs0)>allowed_delta_hs) test_Y0=1;
            if(test_Y0) Y0++;
        }
        if(X1<nmesh_x-1 && !test_X1){
            test++;
            x=++X1;
            for(y=Y0;y<=Y1 && !test_X1;y++)
                if(fabs(hs[x][y]-hs0)>allowed_delta_hs) test_X1=1;
            if(test_X1) X1--;
        }
        if(Y1<nmesh_y-1 && !test_Y1){
            test++;
            y=++Y1;
            for(x=X0;x<=X1 && !test_Y1;x++)
                if(fabs(hs[x][Y1]-hs0)>allowed_delta_hs) test_Y1=1;
            if(test_Y1) Y1--;
        }
/* patch 16102006[1] : force loop-break as soon as a heterogeneity is diagnosed
 * on any side of the box (so that initialization is performed in a quasi-
 * square box. When a very elongated initialization box is generated, it is
 * highly probable that critical conditions will be encountered along one of
 * the longest box edges and most of its successors during propagation of
 * computations (because this edge is generally 
 * orthogonal to mean velocity gradient in the model - e.g., think of a water
 * layer over a vertical constant-gradient medium -) This may lead to very
 * degraded performances with no significant improvement in precision,
 * simply because exactly initialized times are not at all first arrival
 * times and might be recomputed dozens of times according to the generation
 * of head waves... */
/* Note: could be implemented much more efficiently but I prefer this to be
         easily cancellable, and yet with negligible impact on performance */
        if(test) test=(test_X0+test_X1+test_Y0+test_Y1)? 0:1;
/* end of patch 16102006[1] */
    } while(test);
    X1++;
    Y1++;
/* X0,X1,Y0,Y1 are now the inclusive boundaries of the indices
 * of the nodes that will be initialized "exactly".
 * end 281205[2] */

/* decrement boundaries of homogeneous region (if not at model boundaries) */
/* so that heterogeneous interfaces are dealt with by the FD scheme.       */
/* patch 281205[3] : homogeneous zone shrinkage only at heterogeneous
 *                   boundaries (more straightforward...)
 *  if(X0) X0++;
 *  if(Y0) Y0++;
 *  if(X1<nmesh_x) X1--;
 *  if(Y1<nmesh_y) Y1--;
 */
    if(test_X0) X0++;
    if(test_Y0) Y0++;
    if(test_X1) X1--;
    if(test_Y1) Y1--;
/* end 281205[3] */

/* patch 281205[4] : once shrinked, the homogeneous region may endup not
 *                   containing the source point anymore (because it is
 *                   located at the immediate vicinity of a velocity
 *                   heterogeneity). In such case, only minimal initialization
 *                   will be performed (via init_nearest()).
 * The following 6 lines were added :
 */
    if(X0>fxs || X1<fxs || Y0>fys || Y1<fys){
        X0=xsc;
        Y0=ysc;
        X1=xsc+1;
        Y1=ysc+1;
    }
/* end 281205[4] */

/* initialize the time-field according to situation... */
    if(init_stage>=INIT_RECURS_LIMIT ||
        (   (X0==0 || (xs-X0)>=INIT_MIN) &&
            (Y0==0 || (ys-Y0)>=INIT_MIN) &&
            (X1==nmesh_x || (X1-xs)>=INIT_MIN) &&
            (Y1==nmesh_y || (Y1-ys)>=INIT_MIN)
        ) ) {
/* patch 281205[5] : this was simply wrong !
 *      if((X1-X0+1)*(Y1-Y0+1)==1) init_nearest();
 */
        if((X1-X0)*(Y1-Y0)==1) init_nearest();
/* end 281205[5] */

        else {
            for(x=X0;x<=X1;x++)
                for(y=Y0;y<=Y1;y++) {
                    sq_dist=(x-fxs)*(x-fxs)+(y-fys)*(y-fys);
                    t[x][y]=hs0*sqrt(sq_dist);
                }
            if(VERBOSE) 
                printf("\nHomogeneous region: x(%d->%d),y(%d->%d)",X0,X1,Y0,Y1);
        }
        signal=NO_ERROR;
    }
    else{
        if((signal=recursive_init())!=NO_ERROR) return signal;
        X0=max(xs-INIT_MIN,0);
        Y0=max(ys-INIT_MIN,0);
        X1=min(xs+INIT_MIN,nmesh_x);
        Y1=min(ys+INIT_MIN,nmesh_y);
    }

    return signal;

#endif /* !SIMPLISTIC_INIT */

}


/*----------------------------------------------------init_nearest()---------*/
/* last changes in this function 311205: manage headwaves properly and  */
/* account for all cells surrounding src when src coincides with node   */
/* changed init_nearest(), added functions init_cell() and init_cellh() */

static void init_nearest(void)

/* initialize the 1|4|6 nearest nodes when model is immediately inhomogeneous */

{
    int x,y;
    double dx,dy;

    if(VERBOSE) printf("\nInitializing closest node(s).");
    if(source_at_node){
/* formerly: t[xs][ys]=0.0; now also initializing the 8 first neighbours */
        if(xs<nmesh_x && ys<nmesh_y) init_cell(0.,0.,xs,ys);
        if(xs && ys<nmesh_y)         init_cell(1.,0.,xs-1,ys);
        if(xs<nmesh_x && ys)         init_cell(0.,1.,xs,ys-1);
        if(xs && ys)                 init_cell(1.,1.,xs-1,ys-1);
        return;
    }
    x=(fxs<xs) ? xs-1:xs;
    y=(fys<ys) ? ys-1:ys;
    dx=fxs-x;
    dy=fys-y;
    if(xs==fxs){
        init_cell(0.,dy,x,y);
        if(x) init_cell(1.,dy,x-1,y);
    }/* source is located on a x=Ct interface */
    else if(ys==fys){
        init_cell(dx,0.,x,y);
        if(y) init_cell(dx,1.,x,y-1);
    }/* source is located on a y=Ct interface */
    else {
        init_cell(dx,dy,x,y);
    }/* source is located within a cell */
}

static void init_cell(double dx, double dy, int x, int y)

/* x,y : cell indices ; dx,dy : offsets from src to corner at indices x,y */

{
    double hs0,est,hs1;
    hs0=hs[x][y];
    if((est=hs0*sqrt(dx*dx+dy*dy))<t[x][y])
        t[x][y]=est;
    if((est=hs0*sqrt((1.0-dx)*(1.0-dx)+dy*dy))<t[x+1][y])
        t[x+1][y]=est;
    if((est=hs0*sqrt(dx*dx+(1.0-dy)*(1.0-dy)))<t[x][y+1])
        t[x][y+1]=est;
    if((est=hs0*sqrt((1.0-dx)*(1.0-dx)+(1.0-dy)*(1.0-dy)))<t[x+1][y+1])
        t[x+1][y+1]=est;
    if(x && (hs1=hs[x-1][y])<hs0){
        if((est=init_cellh(dx,dy,hs0,hs1))<t[x][y]) t[x][y]=est;
        if((est=init_cellh(dx,1.0-dy,hs0,hs1))<t[x][y+1]) t[x][y+1]=est;
    }
    if(y && (hs1=hs[x][y-1])<hs0){
        if((est=init_cellh(dy,dx,hs0,hs1))<t[x][y]) t[x][y]=est;
        if((est=init_cellh(dy,1.0-dx,hs0,hs1))<t[x+1][y]) t[x+1][y]=est;
    }
    if(x<nmesh_y-1 && (hs1=hs[x+1][y])<hs0){
        if((est=init_cellh(1.0-dx,dy,hs0,hs1))<t[x+1][y]) t[x+1][y]=est;
        if((est=init_cellh(1.0-dx,1.0-dy,hs0,hs1))<t[x+1][y+1]) t[x+1][y+1]=est;
    }
    if(y<nmesh_y-1 && (hs1=hs[x][y+1])<hs0){
        if((est=init_cellh(1.0-dy,dx,hs0,hs1))<t[x][y+1]) t[x][y+1]=est;
        if((est=init_cellh(1.0-dy,1.0-dx,hs0,hs1))<t[x+1][y+1]) t[x+1][y+1]=est;
    }
}

static double init_cellh(double vh, double vv, double hsc, double hsn)

/* headwave initialization (returns headwave arrival time or INFINITY) */
/* vv: distance from src to interface                                  */
/* vh: distance from src to target grid point projected onto interface */
/* hsc,hsn: hs values in current,neighbour cell                        */
{
    double hsd;
    hsd=sqrt(hsc*hsc-hsn*hsn);
    if(vh*hsd>vv*hsn) return vh*hsn+vv*hsd; /* critical condition reached  */
    else return INFINITY;                   /* direct arrival not critical */
}


/*----------------------------------------------------recursive_init()-------*/

static int recursive_init(void)

/* time_2d() is used recursively to initialize the timefield on a */
/* small region around the source with a smaller mesh spacing.    */

{
    int     signal,
            nx_,ny_,
            xs_,ys_,
            messages_,
            n,d,
            i,ii,ihs,i0,
            j,jj,jhs,j0;
    double   
            fxs_,fys_,
            *hs_buf_,*t_buf_,
            HS[N_INIT],T[N_INIT];

/* save static parameters at this stage */
    nx_=nx;
    ny_=ny;
    hs_buf_=hs_buf;
    t_buf_=t_buf;
    xs_=xs;
    ys_=ys;
    fxs_=fxs;
    fys_=fys;
    messages_=messages;

/* increment count of recursivity level */
    init_stage++;

/* free (double **) pointers */
    free((char *)hs);
    free((char *)t);

/* build the re-discretized local model */
    for(i=0;i<N_INIT;i++) HS[i]=T[i]=INFINITY;
    nx=ny=N_INIT_X;
    xs=ys=2*INIT_MIN+1;
    i0=j0=1;
    ihs=xs_-INIT_MIN-1;
    if((d=INIT_MIN-xs_)>=0){
        ihs+=d+1;
        d=1+2*d;
        nx-=d;
        xs-=d;
        i0=0;
    }
    if((d=xs_+INIT_MIN-nx_+1)>=0) nx-=1+2*d;
    jhs=ys_-INIT_MIN-1;
    if((d=INIT_MIN-ys_)>=0){
        jhs+=d+1;
        d=1+2*d;
        ny-=d;
        ys-=d;
        j0=0;
    }
    if((d=ys_+INIT_MIN-ny_+1)>=0) ny-=1+2*d;
    for(i=ihs,n=ii=0;ii<nx;ii++){
        for(j=jhs,jj=0;jj<ny;jj++,n++){
            HS[n]=0.5*hs_buf_[i*ny_+j];
            if(jj%2!=j0) j++;
        }
        if(ii%2!=i0) i++;
    }/* No smoothing is associated with this re-discretization */

/* compute new source coordinates */
    fxs=xs+2.0*(fxs_-xs_);
    fys=ys+2.0*(fys_-ys_);

/* recursively compute times on this model */
    if(VERBOSE){
      printf("\nRecursive initialization: level %d",init_stage);
      signal=time_2d(HS,T,nx,ny,fxs,fys,eps_init,messages);
      printf("\nRecursive initialization: level %d : done.",init_stage);
    }
    else
      signal=time_2d(HS,T,nx,ny,fxs,fys,eps_init,0);

/* assign relevant times to final timefield */
    if(signal==NO_ERROR){
        for(i=ihs+i0,ii=i0;ii<nx;ii+=2,i++)
            for(j=jhs+j0,jj=j0;jj<ny;jj+=2,j++) t_buf_[i*ny_+j]=T[ii*ny+jj];
    }
    else {
      error(signal);
      signal=ERR_RECURS;
    }

/* retrieve old static parameters */
    nx=nx_;
    ny=ny_;
    hs_buf=hs_buf_;
    t_buf=t_buf_;
    xs=xs_;
    ys=ys_;
    fxs=fxs_;
    fys=fys_;
    messages=messages_;

/* reallocate (double **) pointers but do not re-initialize ! */
    if((i=pre_init())!=NO_ERROR){
      error(i);
      signal=ERR_RECURS;
    }

/* decrement count of recursivity level */
    init_stage--;

    return signal;
}


/*---------------------------------------------------Propagate_point()-----*/

static void propagate_point(void)

{
    int test;

    do {
        test=0;
        if(X0>0){
            X0--;
            if(VERBOSE) printf("\n<init#%d> side x=%d.",init_stage,X0);
            y_side(X0,-1,Y0,Y1);
            test++;
        }
        if(Y0>0){
            Y0--;
            if(VERBOSE) printf("\n<init#%d> side y=%d.",init_stage,Y0);
            x_side(Y0,-1,X0,X1);
            test++;
        }
        if(X1<nmesh_x){
            X1++;
            if(VERBOSE) printf("\n<init#%d> side x=%d.",init_stage,X1);
            y_side(X1,1,Y0,Y1);
            test++;
        }
        if(Y1<nmesh_y){
            Y1++;
            if(VERBOSE) printf("\n<init#%d> side y=%d.",init_stage,Y1);
            x_side(Y1,1,X0,X1);
            test++;
        }

    } while(test);

}


/*---------------------------------------------------y_side()---------------*/

static int y_side(int x,int future,int y_begin,int y_end)

/* propagates computations from row x-future to row x */
/* between inclusive y_begin and y_end coordinates.   */
/* y_side() returns the number of stencils adopted.   */

{
    int
        x0,          /* past side coordinate          */
        x_s0,        /* current mesh coordinate       */
        y,           /* current point coordinate      */
        y_old,       /* last local minimum coordinate */
        past,        /* ...opposite to future !       */
        updated,     /* counts accepted stencils      */
        longhead,    /* counts longitudinal headwaves */
        alert;       /* longhead is alert ? (0/1)     */
    double
        hs0,hs1,hs2, /* current slownesses            */
        t_est,       /* current time estimate         */
        dt;          /* time difference               */

    if(reverse_order==0) current_side_limit=x+future;

    updated=longhead=0;
    x0=x-future;
    if(future==1) x_s0=x0;
    else x_s0=x;

    for(y=y_begin;y<=y_end;){

/* Search the next local time minimum on row x-future, save its position */
        while(y<y_end && t[x0][y+1]<t[x0][y]) y++;
        y_old=y;

/* Compute time in front of this local minimum */
        hs1=hs[x_s0][y];
        if(y==0) hs0=INFINITY;
        else hs0=hs[x_s0][y-1];
        if((t_est=t[x0][y]+min(hs0,hs1))<t[x][y]) {
            t[x][y]=t_est;
            updated++;
        }/* timed by 1D transmission */

/* Proceed backwards until the last local time maximum (or y_begin) */
        y--;
        alert=0;
        while(y>=y_begin && (dt=t[x0][y]-t[x0][y+1])>=0.0) {
            hs0=hs[x_s0][y];
            if(dt<hs0/M_SQRT2 
                && (t_est=t[x0][y]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through vertical interface */
            dt=t[x][y+1]-t[x0][y+1];
            if(dt>=0.0 && dt<hs0/M_SQRT2
                && (t_est=t[x][y+1]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through horizontal interface */
            if(y!=0) {
                hs1=hs[x_s0][y-1];
                if((t_est=t[x0][y]+hs1)<t[x][y]){
                    t[x][y]=t_est;
                    updated++;
                }/* 1D transmission towards future */
            }
            if((t_est=t[x0][y+1]+hs0*M_SQRT2)<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* diffraction */
            if(x_s0+future>=0) {
                hs2=hs[x_s0+future][y];
/* FUZZIFIED COMPARISON (11 July 2003) */
                if(hs2<hs0){
                    t_est=t[x][y+1]+hs2;
                    dt=t[x][y]-t_est;
                    if(dt>EPS_FUZZY*t[x][y]){
                        t[x][y]=t_est;
                        updated++;
                        if(!alert){
                            longhead++;
                            send_y_headwave(y,y_begin,
                                x,x_s0,x_s0+future);
                            alert=1;
                        }
                    }/* 1D transmission along the current side. */
                    /* if adopted, implies reverse propagation  */
/* added 090107: must try headwave in the reverse direction too ! */
                    else { 
                        alert=0;
                        t_est=t[x][y]+hs2;
                        dt=t[x][y+1]-t_est;
                        if(dt>EPS_FUZZY*t[x][y+1]){
                            t[x][y+1]=t_est;
                            updated++;
                            send_y_headwave(y+1,y_end,x,x_s0,x_s0+future);
                            longhead++;
                        }/* this happens only VERY rarely... but it DOES ! */
                    }/* 1D forward transmission along the current side. */
                }
            }
            y--;
        }

/* Proceed forwards until the next local time maximum (or y_end) */
        if((y=y_old)==y_end) break;
        y++;
        alert=0;
        while(y<=y_end && (dt=t[x0][y]-t[x0][y-1])>=0.0) {
            hs0=hs[x_s0][y-1];
            if(dt<hs0/M_SQRT2
                && (t_est=t[x0][y]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through vertical interface */
            dt=t[x][y-1]-t[x0][y-1];
            if(dt>=0.0 && dt<hs0/M_SQRT2
                && (t_est=t[x][y-1]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through horizontal interface */
            hs1=hs[x_s0][y];
            if((t_est=t[x0][y]+hs1)<t[x][y]){
                t[x][y]=t_est;
                updated++;
            }/* 1D transmission towards future */
            if((t_est=t[x0][y-1]+hs0*M_SQRT2)<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* diffraction */
            if(x_s0+future>=0) {
                hs2=hs[x_s0+future][y-1];
/* FUZZIFIED COMPARISON (11 July 2003) */
                if(hs2<hs0){
                    t_est=t[x][y-1]+hs2;
                    dt=t[x][y]-t_est;
                    if(dt>EPS_FUZZY*t[x][y]){
                        t[x][y]=t_est;
                        updated++;
                        if(!alert){
                            longhead++;
                            send_y_headwave(y,y_end,
                                x,x_s0,x_s0+future);
                            alert=1;
                        }
                    }/* 1D transmission along the current side. */
                    /* if adopted, implies reverse propagation  */
/* added 090107: must try headwave in the reverse direction too ! */
                    else { 
                        alert=0;
                        t_est=t[x][y]+hs2;
                        dt=t[x][y-1]-t_est;
                        if(dt>EPS_FUZZY*t[x][y-1]){
                            t[x][y-1]=t_est;
                            updated++;
                            send_y_headwave(y-1,y_begin,x,x_s0,x_s0+future);
                            longhead++;
                        }/* this happens only VERY rarely... but it DOES ! */
                    }/* 1D backward transmission along the current side. */
                }
            }
            y++;
        }
    }

/* times are now computed at every point on the current side.  */
/* Reverse propagation may be necessary if a headwave has been */
/* generated along this side...                                */

    if(longhead) {

        reverse_order++;
        if(VERBOSE) 
            printf("\nReverse-propagation(order %d) from row %d:",
                reverse_order,x);
        past= -future;
        for(x=x0;x!=current_side_limit;x+=past) {
            if(x<0 || x>nmesh_x) break; /* The Horrific Case ! (see Note) */
            if(VERBOSE) printf("\nrow #%d: ",x);
            if(y_side(x,past,y_begin,y_end)==0) break;
            if(VERBOSE) printf(" (%d):updated.",x);
        }
        reverse_order--;

    }

    return updated;

}
/* Note(May 1993); New simple tests were added to prevent The Horrific Case. */
/* The Horrific Case may arise in very awkward (non-connex) models when      */
/* several sources are simulated in areas separated by zero-velocity regions */ 
/* and if one of these regions of propagation reaches model boundaries with  */
/* a particular geometry... Whaow !...                                       */
/* This case WAS encountered in a simulation of rupture along a fault-plane. */


/*---------------------------------------------------send_y_headwave()------*/

static void send_y_headwave(int from,int to,int x,int x_s_now,int x_s_future)

/* enforces propagation of a headwave to the end of the current row.   */
/* (this is only needed in very severe models where such headwave may  */
/* provide first arrival at points located FAR from the region where   */
/* critical conditions were encountered...e.g., the slow disk example).*/
/* These headwaves may "jump" over local maxima of the past side.      */
/* In such cases, this information would be lost.                      */

{
    int      y;
    double    hsnow,hsmin,t_est;

    if(from<to) for(y=from;y<to;y++){
        hsnow=hs[x_s_now][y];
        if(x_s_future<0) hsmin=hsnow;
        else hsmin=min(hsnow,hs[x_s_future][y]);
        if((t_est=t[x][y]+hsmin)<t[x][y+1]) t[x][y+1]=t_est;
    }
    else for(y=from;y>to;y--){
        hsnow=hs[x_s_now][y-1];
        if(x_s_future<0) hsmin=hsnow;
        else hsmin=min(hsnow,hs[x_s_future][y-1]);
        if((t_est=t[x][y]+hsmin)<t[x][y-1]) t[x][y-1]=t_est;
    }
}


/*---------------------------------------------------x_side()---------------*/

static int x_side(int y,int future,int x_begin,int x_end)

/* propagates computations from row y-future to row y */
/* between inclusive x_begin and x_end coordinates.   */
/* x_side() returns the number of stencils adopted.   */

{
    int
        y0,          /* past side coordinate          */
        y_s0,        /* current mesh coordinate       */
        x,           /* current point coordinate      */
        x_old,       /* last local minimum coordinate */
        past,        /* ...opposite to future !       */
        updated,     /* counts accepted stencils      */
        longhead,    /* counts longitudinal headwaves */
        alert;       /* longhead is alert ? (0/1)     */
    double
        hs0,hs1,hs2, /* current slownesses            */
        t_est,       /* current time estimate         */
        dt;          /* time difference               */

    if(reverse_order==0) current_side_limit=y+future;

    updated=longhead=0;
    y0=y-future;
    if(future==1) y_s0=y0;
    else y_s0=y;

    for(x=x_begin;x<=x_end;){

/* Search for the next local time minimum on row y-future, save its position */
        while(x<x_end && t[x+1][y0]<t[x][y0]) x++;
        x_old=x;

/* Compute time in front of this local minimum */
        hs1=hs[x][y_s0];
        if(x==0) hs0=hs1;
        else hs0=min(hs1,hs[x-1][y_s0]);
        if((t_est=t[x][y0]+hs0)<t[x][y]) {
            t[x][y]=t_est;
            updated++;
        }/* timed by 1D transmission */

/* Proceed backwards until the last local time maximum (or x_begin) */
        x--;
        alert=0;
        while(x>=x_begin && (dt=t[x][y0]-t[x+1][y0])>=0.0) {
            hs0=hs[x][y_s0];
            if(dt<hs0/M_SQRT2 
                && (t_est=t[x][y0]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through horizontal interface */
            dt=t[x+1][y]-t[x+1][y0];
            if(dt>=0.0 && dt<hs0/M_SQRT2
                && (t_est=t[x+1][y]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through vertical interface */
            if(x!=0) {
                hs1=hs[x-1][y_s0];
                if((t_est=t[x][y0]+hs1)<t[x][y]){
                    t[x][y]=t_est;
                    updated++;
                }/* 1D transmission towards future */
            }
            if((t_est=t[x+1][y0]+hs0*M_SQRT2)<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* diffraction */
            if(y_s0+future>=0) {
                hs2=hs[x][y_s0+future];
/* FUZZIFIED COMPARISON (11 July 2003) */
                if(hs2<hs0){
                    t_est=t[x+1][y]+hs2;
                    dt=t[x][y]-t_est;
                    if(dt>EPS_FUZZY*t[x][y]){
                        t[x][y]=t_est;
                        updated++;
                        if(!alert){
                            longhead++;
                            send_x_headwave(x,x_begin,
                                y,y_s0,y_s0+future);
                            alert=1;
                        }
                    }/* 1D transmission along the current side. */
                    /* if adopted, implies reverse propagation  */
/* added 090107: must try headwave in the reverse direction too ! */
                    else { 
                        alert=0;
                        t_est=t[x][y]+hs2;
                        dt=t[x+1][y]-t_est;
                        if(dt>EPS_FUZZY*t[x+1][y]){
                            t[x+1][y]=t_est;
                            updated++;
                            send_x_headwave(x+1,x_end,y,y_s0,y_s0+future);
                            longhead++;
                        }/* this happens only VERY rarely... but it DOES ! */
                    }/* 1D forward transmission along the current side. */
                }
            }
            x--;
        }

/* Proceed forwards until the next local time maximum (or x_end) */
        if((x=x_old)==x_end) break;
        x++;
        alert=0;
        while(x<=x_end && (dt=t[x][y0]-t[x-1][y0])>=0.0) {
            hs0=hs[x-1][y_s0];
            if(dt<hs0/M_SQRT2
                && (t_est=t[x][y0]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through horizontal interface */
            dt=t[x-1][y]-t[x-1][y0];
            if(dt>=0.0 && dt<hs0/M_SQRT2
                && (t_est=t[x-1][y]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through vertical interface */
            hs1=hs[x][y_s0];
            if((t_est=t[x][y0]+hs1)<t[x][y]){
                t[x][y]=t_est;
                updated++;
            }/* 1D transmission towards future */
            if((t_est=t[x-1][y0]+hs0*M_SQRT2)<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* diffraction */
            if(y_s0+future>=0) {
                hs2=hs[x-1][y_s0+future];
/* FUZZIFIED COMPARISON (11 July 2003) */
                if(hs2<hs0){
                    t_est=t[x-1][y]+hs2;
                    dt=t[x][y]-t_est;
                    if(dt>EPS_FUZZY*t[x][y]){
                        t[x][y]=t_est;
                        updated++;
                        if(!alert){
                            longhead++;
                            send_x_headwave(x,x_end,
                                y,y_s0,y_s0+future);
                            alert=1;
                        }
                    }/* 1D transmission along the current side. */
                    /* if adopted, implies reverse propagation  */
/* added 090107: must try headwave in the reverse direction too ! */
                    else { 
                        alert=0;
                        t_est=t[x][y]+hs2;
                        dt=t[x-1][y]-t_est;
                        if(dt>EPS_FUZZY*t[x-1][y]){
                            t[x-1][y]=t_est;
                            updated++;
                            send_x_headwave(x-1,x_begin,y,y_s0,y_s0+future);
                            longhead++;
                        }/* this happens only VERY rarely... but it DOES ! */
                    }/* 1D backward transmission along the current side. */
                }
            }
            x++;
        }
    }

/* times are now computed at every point on the current side.  */
/* Reverse propagation may be necessary if a headwave has been */
/* generated along this side...                                */

    if(longhead) {
        reverse_order++;
        if(VERBOSE) 
            printf("\nReverse-propagation(order %d) from row %d:",
                reverse_order,y);
        past= -future;
        for(y=y0;y!=current_side_limit;y+=past) {
            if(y<0 || y>nmesh_y) break; /* The Horrific Case ! (see Note) */
            if(VERBOSE) printf("\nrow #%d: ",y);
            if(x_side(y,past,x_begin,x_end)==0) break;
            if(VERBOSE) printf(" (%d):updated.",y);
        }
        reverse_order--;

    }

    return updated;

}


/*---------------------------------------------------send_x_headwave()------*/

static void send_x_headwave(int from,int to,int y,int y_s_now,int y_s_future)

/* enforces propagation of a headwave to the end of the current row.   */
/* (this is only needed in very severe models where such headwave may  */
/* provide first arrival at points located FAR from the region where   */
/* critical conditions were encountered...e.g., the slow disk example).*/
/* These headwaves may "jump" over local maxima of the past side.      */
/* In such cases, this information would be lost.                      */

{
    int    x;
    double  hsnow,hsmin,t_est;

    if(from<to) for(x=from;x<to;x++){
        hsnow=hs[x][y_s_now];
        if(y_s_future<0) hsmin=hsnow;
        else hsmin=min(hsnow,hs[x][y_s_future]);
        if((t_est=t[x][y]+hsmin)<t[x+1][y]) t[x+1][y]=t_est;
    }
    else for(x=from;x>to;x--){
        hsnow=hs[x-1][y_s_now];
        if(y_s_future<0) hsmin=hsnow;
        else hsmin=min(hsnow,hs[x-1][y_s_future]);
        if((t_est=t[x][y]+hsmin)<t[x-1][y]) t[x-1][y]=t_est;
    }
}
/*--------------------------------------------------END: TIME_2D------------*/

/* Scaling - Toby style */
void scaling(double* hsbuf, double* realxs, double* realys, double* S, double h, int size_S){
    int i;
	for ( i = 0; i < size_S-1; i++ ) {    /* model data required by time_2d is hs: grid_spacing*slowness*/
	    hsbuf[i] = S[i] * h;
	}
	*realxs = *realxs / h;           /* xs, ys must be provided in h unit*/
	*realys = *realys / h;
}   
