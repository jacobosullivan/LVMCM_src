/*
    Copyright (C) 2020  Jacob D. O'Sullivan, Axel G. Rossberg

    This file is part of pLVMCM

    pLVMCM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

//////////////////////////////////// PRE-RELEASE VERSION, PLEASE DO NOT DISTRIBUTE /////////////////////////////////////

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <iomanip>
#ifdef ON_SX5FSV
// for bzero on sx5
//#include <sys/types.h>
//#include <sys/ddi.h>
#include <strings.h> 
#endif
#include <stdlib.h>

/* EcoDyn headers: */
#include "ODE.h"
#include "error.h"

/* CVODE header files with a description of contents used in cvdx.c */
#include "sundials/sundials_types.h" 
                           /* definitions of types realtype and             */
                           /* integertype, and the constant FALSE           */
#include "cvode/cvode.h"   /* prototypes for CVodeMalloc, CVode, and        */
                           /* CVodeFree, constants OPT_SIZE, BDF, NEWTON,   */
                           /* SV, SUCCESS, NST,NFE,NSETUPS, NNI, NCFN, NETF */
#include "nvector/nvector_serial.h"
                           /* definitions of type N_Vector and macro        */
                           /* NV_Ith_S, prototypes for N_VNew, N_VFree      */
#include "cvode/cvode_dense.h"
                           /* prototype for CVDense, constant DENSE_NJE     */
#include "cvode/cvode_diag.h"
                           // diagonal solver (not so good??)
#include "cvode/cvode_spgmr.h"
#include "cvode/cvode_bandpre.h"


int newton_failure=0;

ODE_vector::ODE_vector(int length):
  the_length(length),
  the_elements(new realtype[length]),
  the_elements_are_mine(true)
{
    if(!the_elements) FATAL_ERROR("trouble allocating memory");
    for(int i=0;i<length;i++) the_elements[i]=0;
}
ODE_vector::ODE_vector(N_Vector vec):
  the_length(NV_LENGTH_S(vec)),
  the_elements(NV_DATA_S(vec)),
  the_elements_are_mine(false)
{
}
ODE_vector::ODE_vector(realtype * elements,int length):
  the_length(length),
  the_elements(elements),
  the_elements_are_mine(false)
{
}
ODE_vector::ODE_vector(const ODE_vector & other):
  the_length(other.the_length),
  the_elements(new realtype[the_length]),
  the_elements_are_mine(true){
  ASSERT(the_elements);
  memcpy(the_elements,other.the_elements,the_length*sizeof(*the_elements));
}
ODE_vector::~ODE_vector(){
    if(the_elements_are_mine)
      delete[] the_elements;
}

ODE_vector exp(ODE_vector & v){
  ODE_vector w(v.size());
  for(int i=v.size();i-->0;)
    w[i]=exp(v[i]);
  return w;
}

ODE_matrix::ODE_matrix(int length):the_length(length){
  typedef realtype * pointer_to_realtype;
  the_elements=new pointer_to_realtype[length];
  if(!the_elements) FATAL_ERROR("trouble allocating memory");
  for(int j=0;j<the_length;j++){
    typedef realtype tt;
    the_elements[j]=new tt[length];
    if(!the_elements[j]) FATAL_ERROR("trouble allocating memory");
  }
  the_elements_are_mine=true;
}
#if SUNDIALS_PRE_2_4_0
ODE_matrix::ODE_matrix(DenseMat mat){ // probably dead code
  ALWAYS_ASSERT( mat->M == mat->N );
  the_elements=(mat->data);
  the_length=(mat->M);
  the_elements_are_mine=false;
}
#else
ODE_matrix::ODE_matrix(DlsMat mat){ // probably dead code
  ALWAYS_ASSERT( mat->M == mat->N );
  the_elements=(mat->cols);
  the_length=(mat->M);
  the_elements_are_mine=false;
}
#endif
ODE_matrix::ODE_matrix(const ODE_matrix & other):
  the_elements_are_mine(true),
  the_length(other.the_length)
{
  int length=the_length;
  typedef realtype * pointer_to_realtype;
  the_elements=new pointer_to_realtype[length];
  if(!the_elements) FATAL_ERROR("trouble allocating memory");
  for(int j=0;j<the_length;j++){
    typedef realtype tt;
    the_elements[j]=new tt[length];
    if(!the_elements[j]) FATAL_ERROR("trouble allocating memory");
    memcpy(the_elements[j],other.the_elements[j],
	   length*sizeof(**the_elements));
  }
}

ODE_matrix::~ODE_matrix(){
  if(the_elements_are_mine){
    for(int j=0;j<the_length;j++)
      delete[] the_elements[j];
    delete[] the_elements;
  }
}


const  ODE_vector & ODE_vector::operator=(ODE_vector const &other){
  if(other.the_length!=the_length){
    if(the_elements_are_mine){
      delete[] the_elements;
    }
    the_length=other.the_length;
    the_elements=new realtype[the_length];
  }
  memcpy(the_elements, other.the_elements, the_length*sizeof(*the_elements));
  return *this;
}
ODE_vector& ODE_vector::operator+=(ODE_vector const &other){
  if(other.the_length!=the_length)
    FATAL_ERROR("ODE_vectors must be of the same size in Assignments");
  for(int i=0;i<the_length;i++)
    (*this)[i]+=other[i];
  return *this;
}
  
ODE_vector& ODE_vector::operator-=(ODE_vector const &other){
  if(other.the_length!=the_length)
    FATAL_ERROR("ODE_vectors must be of the same size in Assignments");
  for(int i=0;i<the_length;i++)
    (*this)[i]-=other[i];
  return *this;
}
  
ODE_vector&  ODE_vector::operator*=(double x){
  for(int i=0;i<the_length;i++)
    (*this)[i]*=x;
  return *this;
}

//  #include <string.h>


void ODE_vector::clear(){
  bzero(the_elements, the_length*sizeof(*the_elements));
}

void ODE_matrix::clear(){
  for(int i=0;i<the_length;i++){
    bzero(the_elements[i], the_length*sizeof(**the_elements));
  }
}

#if DEBUGGING
static int precision_of_vector_print=20;
#else
static int precision_of_vector_print=3;
#endif
extern double TolR, TolA;
static double default_relative_tolerance=TolR;
static double default_absolute_tolerance=TolA;
static int max_integrator_steps = 1<<15;
static double max_step_size = 1e15;
static double min_step_size = 0;
static int precondition = 1;
static int use_default_preconditioner = 0;
static double ode_initial_step_size = 1e-10;


using namespace std;

int ODE_dynamical_object::Jacobian(ODE_vector const & state,
				   ODE_vector const & dynamics,
				   ODE_matrix & jac){
  FATAL_ERROR("The Jacobian has not been implemented for this class.");
}

void ODE_dynamical_object::JTimes(ODE_vector const & state,
				  ODE_vector const & in,
				  ODE_vector & out){
  FATAL_ERROR(
"Jacobian multiplication has not been implemented for this class."
);
}

void ODE_dynamical_object::precondition(ODE_vector const & state,
					ODE_vector const & in,
					ODE_vector & out,
					realtype gamma,
					bool left_rather_than_right){
  FATAL_ERROR("No preconditioner defined!" << endl << \
	      "The preconditioner computes" << endl << 
	      "   out = P_i^(-1) in," << endl <<
	      "where i=1,2 and P=P_1*P_2 is an approximation of the Newton Matrix N " << endl <<
	      "   N = I - gamma J " << endl <<
	      "where J is the Jacobian of the dynamics f.");
}

realtype ODE_state::dummy_real=0;

void ODE_dynamical_object::test_Jacobian(){
  ODE_state state(this);
  int N=state.size();
  for(int i=0;i<N;i++){
    ODE_vector in(N),in2(N),inJ(N),out(N),out2(N),outJ(N);
    //derivative with respect to i
    in=state;
    in2=in;
    in2[i]+=0.001*state[i];
    inJ[i]=1;
    dynamics(in,out);
    dynamics(in2,out2);
    out*=-1;
    out2+=out;
    out2*=1/(in2[i]-in[i]);
    JTimes(state,inJ,outJ);
    for(int j=0;j<N;j++)
      cout << i << " " << j << " " << out2[j] << " " << outJ[j] << endl;
  }
}
    
/* Functions Called by the CVODE Solver */

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
#if SUNDIALS_PRE_2_4_0
static int Jac(long int N, DenseMat J, realtype t,
		N_Vector y, N_Vector fy, void *jac_data,
		N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
#else // sundials version > 2.4.0
#ifdef SUNDIALS_VERSION_2_4
static int Jac(int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, 
		void *jac_data,
		N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
#else // sundials version >= 2.5
static int Jac(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, 
		void *jac_data,
		N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
#endif
#endif
static int PrecondSetupFn(realtype t, N_Vector y, N_Vector fy, 
			  booleantype jok, booleantype * jcurPtr,
			  realtype gamma, void *p_data, 
			  N_Vector tmp1,N_Vector tmp2,N_Vector tmp3){
  *jcurPtr=FALSE; //we did not recompute Jacobian data here
  return 0; //0 means we did the setup successfully (there was nothing
            //to be done)
}
static int PrecondSolveFn(realtype t, N_Vector y, 
			  N_Vector fy, N_Vector r, N_Vector z,
			  realtype gamma, realtype delta,
			  int lr, void *p_data, N_Vector tmp);


void ODE_state::prepare_integrator(){
  const int NEQ=this->size();
  int flag;

  y = N_VMake_Serial(size(),&((*this)[0]));
  
#if SUNDIALS_PRE_2_4_0
  abstol=N_VNew_Serial(size());
  for(int i=0;i< NEQ; i++){
    NV_Ith_S(abstol,i)=default_absolute_tolerance;
  }
#else
  abstol=default_absolute_tolerance;
#endif

  // The recommended choices here are (CV_ADAMS, CV_FUNCTIONAL) for
  // nonstiff problems and (CV_BDF, CV_NEWTON) for stiff problems:
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (cvode_mem == NULL) { FATAL_ERROR("CVodeCreate failed."); }
#if SUNDIALS_PRE_2_4_0
  flag = CVodeMalloc(cvode_mem, f, 0, y, CV_SV, reltol, abstol);
  if (flag != CV_SUCCESS) { FATAL_ERROR("CVodeMalloc failed with code " <<
					flag ); }
#else
  flag = CVodeInit(cvode_mem, f, 0, y);
  if (flag != CV_SUCCESS) { FATAL_ERROR("CVodeInit failed with code " <<
					flag ); }
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (flag != CV_SUCCESS) { FATAL_ERROR("CVodeSStolerances failed with code " <<
					flag ); }
#endif
  
  if(precondition && use_default_preconditioner){
#if SUNDIALS_PRE_2_4_0
    bp_data = CVBandPrecAlloc(cvode_mem, size(), size()-1, size()-1);
#else
    flag = CVBand(cvode_mem, size(), size()-1, size()-1);
    if (flag != CV_SUCCESS) { 
      FATAL_ERROR("CVBand failed with code " <<
		  flag ); 
    }
#endif
  }


  if(
#if SUNDIALS_PRE_2_4_0
     CVodeSetFdata(cvode_mem, the_dynamics)!=CV_SUCCESS ||
#else // sundials >= 2.4.0
     CVodeSetUserData(cvode_mem, the_dynamics)!=CV_SUCCESS ||
#endif
     CVodeSetErrFile(cvode_mem, stdout)!=CV_SUCCESS ||
     CVodeSetMaxNumSteps(cvode_mem,max_integrator_steps) ||
     CVodeSetMaxStep(cvode_mem,max_step_size) ||
     CVodeSetMinStep(cvode_mem,min_step_size) ||
     false ){
    FATAL_ERROR("setting optinal CVODE parameters");
  }

  switch(linear_solver){
  case DIAG:
    flag = CVDiag(cvode_mem);
    if (flag != CVDIAG_SUCCESS) { printf("CVDiag failed.\n"); exit(1); }
    break;
  case DENSE:
    flag = CVDense(cvode_mem, size());
#if SUNDIALS_PRE_2_4_0
    if (flag != CVDENSE_SUCCESS) { printf("CVDense failed.\n"); exit(1); }
#else
    if (flag != CVDLS_SUCCESS) { printf("CVDense failed.\n"); exit(1); }
#endif
    if(the_dynamics->can_calculate_Jacobian()){
#if SUNDIALS_PRE_2_4_0
      flag = CVDenseSetJacFn(cvode_mem, Jac, NULL);
#else
      // If you encounter a problem with the type of the "Jac"
      // function argument below, please try passing a "-D
      // SUNDIALS_VERSION_2_4" compiler option.
      flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
#endif
      if (flag != 1/*SUCCESS*/){ 
	printf("Setting Jacobian failed.\n"); exit(1); 
      }
    }
    break;
  case CVSPGMR:
    //  pretype (int) specifies the preconditioning type and must be one
    //  of: PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH.
    {
      int pretype = 
	(precondition && (the_dynamics->has_preconditioner() || bp_data) ?
	 PREC_LEFT : PREC_NONE );
#if SUNDIALS_PRE_2_4_0
      if(bp_data){
	// using build-in preconditioner:
	flag = CVBPSpgmr(cvode_mem, pretype, 0, bp_data);
	if (flag != CVSPILS_SUCCESS) { 
	  printf("CVBPSpgmr failed.\n"); 
	  exit(1); 
	}
      }else{
	// using custom preconditioner (does not work well):
	flag = CVSpgmr(cvode_mem, pretype, 0);
	if (flag != CVSPILS_SUCCESS) { printf("CVSpgmr failed.\n"); exit(1); }
	flag = CVSpilsSetPreconditioner(cvode_mem, PrecondSetupFn, PrecondSolveFn, 
					the_dynamics);
	if (flag != CVSPILS_SUCCESS) 
	  { printf("setting preconditioner failed.\n"); exit(1); };
      }
#else
      flag = CVSpgmr(cvode_mem, pretype, 0);
      if (flag != CVSPILS_SUCCESS) { printf("CVSpgmr failed.\n"); exit(1); }
      flag = CVSpilsSetPreconditioner(cvode_mem, PrecondSetupFn, PrecondSolveFn);
      if (flag != CVSPILS_SUCCESS) 
	{ printf("setting preconditioner failed.\n"); exit(1); };
#endif
    }

    break;
  default:
    FATAL_ERROR("unknown linear_solver");
  } // switch(linear_solver)
}

void ODE_state::release_integrator(){
  N_VDestroy_Serial(y);
#if SUNDIALS_PRE_2_4_0
  N_VDestroy_Serial(abstol);
  if(bp_data)
    CVBandPrecFree(&bp_data);
#endif
  CVodeFree(&cvode_mem);        /* Free the CVODE problem memory */
}

ODE_state::ODE_state(ODE_dynamical_object * object) : 
    ODE_vector(object->number_of_variables()),
    the_time_since_start(0),
    the_start_time(object->the_start_time),
    reltol(default_absolute_tolerance),
    //  machEnv(M_EnvInit_Serial(this->size())),
    linear_solver(CVSPGMR), // fastest
    //linear_solver(DIAG),  // same speed as DENSE when I last tested
    //linear_solver(DENSE),
    the_dynamics(object),
    bp_data(0)
{
  object->prepare_for_integration();
  the_start_time=the_dynamics->current_time;
  object->write_state_to(*this);

  reltol = default_relative_tolerance;

  prepare_integrator();
  CVodeSetInitStep(cvode_mem,ode_initial_step_size);
}


ODE_state::~ODE_state()
{
  the_dynamics->read_state_from(*this);
  the_dynamics->current_time=the_start_time+the_time_since_start;
  the_dynamics->cleanup_after_integration();
  release_integrator();
}

realtype ODE_state::redo_with_shorter_step_size(){
  // (returns t from where we restart)
  realtype last_step_size;
  if(CVodeGetLastStep(cvode_mem,&last_step_size)!=CV_SUCCESS){
    FATAL_ERROR("Can't get last step size");
  }
  release_integrator();
  the_dynamics->write_state_to(*this);
  prepare_integrator();
  WARNING("reducing stepsize from " << 
	  last_step_size << " to " << last_step_size/2);
  CVodeSetInitStep(cvode_mem,last_step_size/2);
  //CVodeSetMaxStep(...);
  return the_start_time;
}

void ODE_state::restart(double t_start){
  release_integrator();
  the_start_time=t_start;
  the_dynamics->current_time=the_start_time;
  the_time_since_start=0;
  the_dynamics->read_state_from(*this);
  prepare_integrator();
  CVodeSetInitStep(cvode_mem,ode_initial_step_size);
}

int ODE_state::
integrate_until(realtype target_time){
  int flag;
  flag = CVode(cvode_mem, target_time-the_start_time, y, 
	       &the_time_since_start, CV_NORMAL);
  if (flag != CV_SUCCESS) {
    fprintf(stderr,"CVode failed, flag=%d.\n", flag); // <<<<<
    WARNING("integrator failed!");
    return 1;
    //    FATAL_ERROR("Exiting");
  }
  the_dynamics->current_time=the_start_time+the_time_since_start;
  return 0;
}

int ODE_state::
integrate_one_step(realtype target_time){
  int flag;
  flag = CVode(cvode_mem, target_time-the_start_time, y, 
	       &the_time_since_start, CV_ONE_STEP);
  the_dynamics->current_time=the_start_time+the_time_since_start;
  if (flag != CV_SUCCESS) {
    fprintf(stderr,"CVode failed, flag=%d.\n", flag);
    WARNING("integrator failed!");
    return 1;
    //    FATAL_ERROR("Exiting");
  }
  return 0;
}

int ODE_state::
current_time_derivative(ODE_vector & ddt){
  int flag;
  N_Vector ddt_N_Vector = N_VMake_Serial(size(),&(ddt[0]));
  flag = CVodeGetDky(cvode_mem, 
		     the_dynamics->current_time - 
		     the_time_since_start, 
		     1, 
		     ddt_N_Vector);
  if (flag != CV_SUCCESS) {
    fprintf(stderr,"CVodeGetDky failed, flag=%d.\n", flag);
    WARNING("Could not compute time derivative!");
    return 1;
    //    FATAL_ERROR("Exiting");
  }
  return 0;
}

    

  
int ODE_state::Yoshida_step(double dt){
  ODE_vector ODE_ydot(size());
  ODE_vector ODE_y(y);
  
  the_dynamics->dynamics(ODE_y,ODE_ydot);
  ODE_ydot*=dt;
  ODE_y+=ODE_ydot;
  the_time_since_start+=dt;
  the_dynamics->current_time=the_start_time+the_time_since_start;
  return 0;
}

void ODE_state::short_diagnosis(ostream &co){
  long int i;
  int j;
  realtype x;
  CVodeGetNumSteps(cvode_mem,&i);
  co << i << " steps" << endl;
}

void ODE_state::diagnosis(ostream &co){
  long int i;
  int j;
  realtype x;
  co << "Diagnosis" << endl;
  co << "current_time " << the_dynamics->current_time << endl;
  CVodeGetNumSteps(cvode_mem,&i);
  co << i << " internal ODE steps" << endl;
  CVodeGetNumRhsEvals(cvode_mem, &i);
  co << i << " function evaluations" << endl;
  if(bp_data){
    CVBandPrecGetNumRhsEvals(bp_data, &i);
    co << i << " function evaluations for preconditioner" << endl;
  }
  CVodeGetNumNonlinSolvIters(cvode_mem, &i);
  co << i << " nonlinear iterations" << endl;
  CVodeGetNumLinSolvSetups(cvode_mem, &i);
  co << i << " linear solver setups (since that's expensive)" << endl;
  if(linear_solver==CVSPGMR){
    CVSpilsGetNumPrecEvals(cvode_mem, &i);
    co << i << " preconditioner setups" << endl;
    CVSpilsGetNumPrecSolves(cvode_mem, &i);
    co << i << " preconditioner calls" << endl;
  }
  CVodeGetNumErrTestFails(cvode_mem, &i);
  co << i << " local error (accuracy) test failures" << endl;
  CVodeGetNumNonlinSolvConvFails(cvode_mem, &i);
  co << i << " nonlinear solver failures" << endl;
  CVodeGetLastOrder(cvode_mem, &j);
  co << j << "-" << (j==1?"st":j==2?"nd":j==3?"rd":"th") << 
    " order in final step" << endl;
  CVodeGetActualInitStep(cvode_mem, &x);
  co << x << " was the size of the first ODE step" << endl;
  CVodeGetLastStep(cvode_mem, &x);
  co << x << " was the size of the last ODE step" << endl;
}


//*********************************************************************
//*********************************************************************
//*********************************************************************

#include "kinsol/kinsol.h"
#include "kinsol/kinsol_dense.h"
#include "kinsol/kinsol_dense.h"
#include "kinsol/kinsol_spgmr.h"

#undef USE_HOMPACK90
static ODE_dynamical_object * HOMPACK90_object;
extern "C" {
  void __hompack90_MOD_fixpdf(int *N,double *Y,int *IFLAG,double *ARCTOL,
	       double *EPS,int *TRACE,double *A,int *NDIMA,
	       int *NFE,double *ARCLEN);
}


static int fixed_point_analyzer_f(N_Vector y, N_Vector ydot,
				  void *f_data);

const int fixed_point_analyzer::F_IS_NAN=-100;

fixed_point_analyzer::fixed_point_analyzer(ODE_dynamical_object * object,
					   double lower_bound_shift_request):
  ODE_vector(object->number_of_variables()),
  //linear_solver(DENSE), 
  linear_solver(SPGMR), 
  scaling_vector_D_u(N_VNew_Serial(size())),
  scaling_vector_D_F(N_VNew_Serial(size())),
  tmpl(N_VNew_Serial(size())),
  u(N_VNew_Serial(size())),
  constraints(N_VNew_Serial(size())),
  lower_bound(size()),
  lower_bound_shift(lower_bound_shift_request),
  the_dynamics(object),
  user_data(object,&lower_bound)
{
  prepare_solver();
  ODE_vector state(size());
  object->write_state_to(state);
  for(int i=size();i-->0;){
    NV_DATA_S(scaling_vector_D_u)[i]=1;
    NV_DATA_S(scaling_vector_D_F)[i]=1;
    lower_bound[i]=state[i]+lower_bound_shift;
    NV_DATA_S(constraints)[i]=1;
  }
  KINSetConstraints(kin_mem, constraints);
}

void fixed_point_analyzer::
set_time_scale_to(double dt){
  for(int i=size();i-->0;){
    NV_DATA_S(scaling_vector_D_F)[i]=dt/NV_DATA_S(scaling_vector_D_u)[i];
  }
}  


fixed_point_analyzer::~fixed_point_analyzer(){
  release_solver();
  N_VDestroy_Serial(scaling_vector_D_u);
  N_VDestroy_Serial(scaling_vector_D_F);
  N_VDestroy_Serial(tmpl);
  N_VDestroy_Serial(u);
  N_VDestroy_Serial(constraints);
}
    

void fixed_point_analyzer::
prepare_solver(){
  kin_mem = KINCreate();
  if (kin_mem == NULL) { FATAL_ERROR("KINCreate failed."); }
  KINSetPrintLevel(kin_mem, 1); //values 0,...,3 allowed
#if SUNDIALS_PRE_2_4_0  
  if(KIN_SUCCESS!=KINMalloc(kin_mem, fixed_point_analyzer_f, tmpl)){
    FATAL_ERROR("KINMalloc failed.");
  }
#else
  if(KIN_SUCCESS!=KINInit(kin_mem, fixed_point_analyzer_f, tmpl)){
    FATAL_ERROR("KINMalloc failed.");
  }
#endif
#if SUNDIALS_PRE_2_4_0
  KINSetFdata(kin_mem, &user_data);
#else
  KINSetUserData(kin_mem, &user_data);
#endif
  KINSetNumMaxIters(kin_mem,1<<15);
  KINSetFuncNormTol(kin_mem,1e-9); //
  //  KINSetScaledStepTol(kin_mem, 1e-100);
  switch(linear_solver){
  case DENSE:
#if SUNDIALS_PRE_2_4_0
    if(KINDENSE_SUCCESS!=KINDense(kin_mem, size())){
      FATAL_ERROR("Linear solver initialization failed");
    }
#else
    if(KINDLS_SUCCESS!=KINDense(kin_mem, size())){
      FATAL_ERROR("Linear solver initialization failed");
    }
#endif
    break;
  case SPGMR:
    if(KINSPILS_SUCCESS!=KINSpgmr(kin_mem, 0)){
      FATAL_ERROR("Linear solver initialization failed");
    }
    break;
  default: 
    FATAL_ERROR("Linear solver type not implemented, yet");
  }
}

void fixed_point_analyzer::
release_solver(){
  KINFree(&kin_mem);
}


int fixed_point_analyzer::
snap_to_fixed_point(){
#ifndef USE_HOMPACK90
  ODE_vector state(u);// state shares data with u
  the_dynamics->write_state_to(state);

  for(int i=state.size();i-->0;){
    state[i]-=lower_bound[i];
    
  }

  int retval;
  // Either the C++ catch/throw or the the C setjump/longjmp should
  // get us out of KINSol in case of trouble.  I do not know yet which
  // works better.
  try{
    if(setjmp(the_dynamics->the_sigjmp_buf)){
      retval=F_IS_NAN;
    }else{
      retval=
	KINSol(kin_mem, u, (1 ? KIN_LINESEARCH : KIN_NONE /* no line search */ ), 
	       scaling_vector_D_u,
	       scaling_vector_D_F);
    }
  }catch(fixed_point_analyzer::failure){
    retval=F_IS_NAN;
  }
    
  if(retval<0){
    WARNING("Could not find fixed point.");
    REPORT(retval);
  }
  
  for(int i=state.size();i-->0;){
    state[i]+=lower_bound[i];
  }

  the_dynamics->read_state_from(state);

  return retval;
#else // USE_HOMPACK90 is defined
  int N=the_dynamics->number_of_variables();
  ODE_vector Y(N+1); // =(lambda, state)
  ODE_vector state(&Y[1],N); 
  int IFLAG=-1; // means seach ODE fixed point
  double EPS=1e-15; // "local error allowed the ODE solver when very near
  double ARCTOL=1e-8; // tolerance of ODE solver
                //  the fixed point"
  int TRACE=0;  // nothing written to file
  ODE_vector A(N); // "parameter vector"
  int NDIMA=N; // size of A
  int NFE; // number of function evaluations
  double ARCLEN; // "length of path followed"
  
  HOMPACK90_object=the_dynamics;//pass this to F(x), etc.

  the_dynamics->prepare_for_integration();
  the_dynamics->write_state_to(state);Y[0]=0;
  do {
    __hompack90_MOD_fixpdf(&N,&Y[0],&IFLAG,&ARCTOL,&EPS,&TRACE,&A[0],
			   &NDIMA,&NFE,&ARCLEN);
    REPORT(ARCLEN);
    if(IFLAG == 5 && EPS > 1e-15) {
      EPS*=0.1;
      ARCTOL*=0.1;
      REPORT(EPS);
      IFLAG=-1;
    }
  } while(IFLAG == 2 || IFLAG == 3 || IFLAG == -1);
  the_dynamics->read_state_from(state);
  the_dynamics->cleanup_after_integration();  
  
  REPORT(IFLAG);
  REPORT(Y[0]);
  REPORT(NFE);
  REPORT(ARCLEN);

  switch(IFLAG){
  case 4:   
    WARNING("Jacobian matrix does not have full rank.  The algorithm\n"
"has failed (the zero curve of the homotopy map cannot be \n"
"followed any further).");
    break;
  case 6:   
    WARNING("I - DF(X)  is nearly singular at the fixed point \n"
"(DF(X) is nearly singular at the zero, or  D RHO(A,LAMBDA,X)/DX is\n"
"nearly singular at  LAMBDA = 1 ).  Answer may not be accurate.");
    break;
  case 7: FATAL_ERROR("Illegal input parameters, a fatal error.");
    break;
  case 8: WARNING("Memory allocation error.");
    break;
  default: 
    break;
  }
#endif
}



//*********************************************************************
// general use functions **********************************************
//*********************************************************************


/***************** Functions Called by the CVODE Solver ******************/

/* f routine. Compute f(t,y). */

static int f(realtype t, N_Vector y, N_Vector ydot,
              void *f_data)
{
  ODE_vector ODE_y(y),ODE_ydot(ydot);

  ODE_dynamical_object * obj=(ODE_dynamical_object*)f_data;
  obj->current_time=t+obj->the_start_time; 
  return obj->dynamics(ODE_y,ODE_ydot);
}

/* f routine. Compute f(y). */

static int fixed_point_analyzer_f(N_Vector y, N_Vector ydot,
				  void *f_data)
{
  fixed_point_analyzer::user_data_t * user_data=
    (fixed_point_analyzer::user_data_t *)f_data;
  ODE_dynamical_object * obj=user_data->first;
  int n=obj->number_of_variables();

  const ODE_vector & lower_bound=*(user_data->second);
  ODE_vector ODE_y(n),ODE_ydot(ydot);

  for(int i=n;i-->0;){
    ODE_y[i]=NV_DATA_S(y)[i]+lower_bound[i];
  }
  int retval = obj->dynamics(ODE_y,ODE_ydot);
  if(retval == 0){
    for(int i=n;i-->0;){
      if(!isfinite(ODE_ydot[i]) or isnan(ODE_ydot[i])){
	REPORT(i);
	REPORT(ODE_ydot[i]);
	WARNING("throw fixed_point_analyzer::failure();");
	//return 1;
	longjmp(obj->the_sigjmp_buf, 1);
      }
    }
  }
  
  if(obj->has_inherent_rates()){
    ODE_vector rates(n);
    obj->get_inherent_rates(rates);
    for(int i=n;i-->0;){
      ODE_ydot[i]/=rates[i];
    }
  }

  return retval;
}

/* Jacobian routine. Compute J(t,y). */

#if SUNDIALS_PRE_2_4_0
static int Jac(long int N, DenseMat J, realtype t,
		N_Vector y, N_Vector fy, void *jac_data,
		N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
#else // sundials version > 2.4.0
#ifdef SUNDIALS_VERSION_2_4
static int Jac(int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, 
		void *jac_data,
		N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
#else // sundials version >= 2.5
static int Jac(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, 
		void *jac_data,
		N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
#endif
#endif
 {
   ODE_vector ODE_y(y), ODE_ydot(fy);
   ODE_matrix ODE_J(J);

   ((ODE_dynamical_object*)jac_data)->Jacobian(ODE_y,ODE_ydot,ODE_J);
 }

static int PrecondSolveFn(realtype t, N_Vector y, 
			  N_Vector fy, N_Vector r, N_Vector z,
			  realtype gamma, realtype delta,
			  int lr, void *p_data, N_Vector tmp){
  // we assume some very simple preconditioning here, so we pass only
  // a minimum of all these parameters over to the dynamical object's
  // preconditioner.

  ODE_vector ODE_y(y), ODE_fy(fy);
  ODE_vector in(r), out(z);
  ODE_dynamical_object * obj=(ODE_dynamical_object*)p_data;
  obj->current_time=t+obj->the_start_time; 
  bool left_rather_than_right = (lr==1);
  obj->precondition(ODE_fy,in,out,gamma,left_rather_than_right);
  return 0;
}  

extern "C" 
{
  void hompack90_f_(double * x, double * v);
  void hompack90_fjac_(double * x, double * v, int * k);
}
void hompack90_f_(double * x, double * v)
{
  ODE_dynamical_object & obj = *HOMPACK90_object; //abbreviation

  // Wrap the data into ODE_vectors
  ODE_vector state(x,obj.number_of_variables());
  ODE_vector time_derivative(v,obj.number_of_variables());
  
  if(HOMPACK90_object->dynamics(state,time_derivative)){
    FATAL_ERROR("Problem computing fixed point");
  }
  for(int i=time_derivative.size();i-->0;){
    time_derivative[i]*=exp(state[i]);
  }
  time_derivative*=-1;

  static int count=0;
  if((count++)%100==0){
    cout << "*"; cout.flush();
  }
  // realtype sum=0;
  // for(int i=0;i<time_derivative.size();i++){
  //   sum+=time_derivative[i]*time_derivative[i];
  // }
  // REPORT(sum);
}
void hompack90_fjac_(double * x, double * v, int * k)
{
  ODE_dynamical_object & obj = *HOMPACK90_object; //abbreviation

  // Wrap the data into ODE_vectors
  ODE_vector state(x,obj.number_of_variables());
  ODE_vector time_derivative1(v,obj.number_of_variables());
  ODE_vector time_derivative2(obj.number_of_variables());
  
  const realtype eps=1e-5;
  const realtype x_k=x[*k];

  x[*k]=x_k+eps;
  hompack90_f_(&state[0],&time_derivative1[0]);

  x[*k]=x_k-eps;
  hompack90_f_(&state[0],&time_derivative2[0]);

  time_derivative1-=time_derivative2;
  time_derivative1*=1/(2*eps);

  x[*k]=x_k;
}


std::ostream & operator<<(std::ostream &stream, const ODE_vector  &vector){
  ios::fmtflags initial_flags=stream.flags();
  streamsize initial_prec=stream.precision();
  
  left(stream);scientific(stream);
  for(int i=0;i<vector.size();i++){
    stream << setprecision(precision_of_vector_print) 
	   << std::setw(precision_of_vector_print+7) << vector[i] << " ";
  }

  stream.precision(initial_prec);
  stream.flags(initial_flags);
  return stream;
}

std::ostream & operator<<(std::ostream &stream, 
			  const ODE_state & state){
//   left(stream);scientific(stream);
//   stream << setprecision(precision_of_vector_print);
  state.the_dynamics->line_print(state,stream);
  return stream;
}

void ODE_dynamical_object::line_print(ODE_vector const & state, ostream &co){
  for(int i=0;i<state.size();i++){
    co << state[i] << " ";
  }
}

// Local Variables:
// c-file-style: "stroustrup"
// End:
