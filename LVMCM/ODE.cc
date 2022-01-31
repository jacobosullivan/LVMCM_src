// $Id: ODE.cc 2519 2020-05-15 08:27:03Z axel $
// DENSE removed in cvodes v2.4.0
#include <math.h>
#include <stdio.h>
#include <iomanip>
#include <stdlib.h>

/* EcoDyn headers: */
//#define DEBUGGING 1
#include "ODE.h"
#include "error.h"

/* CVODE header files with a description of contents used in cvdx.c */
#include <sundials/sundials_types.h>
/* definitions of types realtype and             */
/* integertype, and the constant FALSE           */
#include <nvector/nvector_serial.h>
/* definitions of type N_Vector and macro        */
/* NV_Ith_S, prototypes for N_VNew, N_VFree      */
#include <cvode/cvode.h>   /* prototypes for CVodeMalloc, CVode, and        */
/* CVodeFree, constants OPT_SIZE, BDF, NEWTON,   */
/* SV, SUCCESS, NST,NFE,NSETUPS, NNI, NCFN, NETF */
#include <cvode/cvode_diag.h>
// diagonal solver (not so good??)
#include <cvode/cvode_spils.h>
#include <ida/ida.h>       /* prototypes for IDA fcts., consts.             */
#include <ida/ida_spils.h>             /* access to IDASpils interface         */
#include <sundials/sundials_types.h>   /* definition of type realtype          */


#ifndef CVLS_SUCCESS
#define SUNDIALS_3
#endif

#ifdef SUNDIALS_3
#define CVLS_SUCCESS CVSPILS_SUCCESS
#define IDALS_SUCCESS IDASPILS_SUCCESS
#define SUNLinSol_SPGMR SUNSPGMR
#endif

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
ODE_matrix::ODE_matrix(DlsMat mat){ // probably dead code
    ALWAYS_ASSERT( mat->M == mat->N );
    the_elements=(mat->cols);
    the_length=(mat->M);
    the_elements_are_mine=false;
}
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
static double default_relative_tolerance=0.0000;
static double default_absolute_tolerance=0.0001;
static int max_integrator_steps = 1<<15;
static double max_step_size = 1e15;
static double min_step_size = 0;
static int precondition = 1;
static int use_default_preconditioner = 0;
static double ode_initial_step_size = 1e-5;


using namespace std;

int ODE_dynamical_object::number_of_root_functions() const{
    return 0;
}


int ODE_dynamical_object::set_root_directions(int * directions) const{
    return -1; // non-zero return indicates that directions are not set
}

int ODE_dynamical_object::number_of_derived_quantities(){
    return 0;
}

int ODE_dynamical_object::root_functions(ODE_vector const & state,
                                         ODE_vector & gout){
    int d = number_of_derived_quantities();
    ODE_vector derived(&(state[number_of_variables()]),d);
    return root_functions(state,derived,gout);
}

int ODE_dynamical_object::root_functions(ODE_vector const & state,
                                         ODE_vector const & derived_quantities,
                                         ODE_vector & gout){
    FATAL_ERROR("root_functions has not been implemented for this class.");
    return 0;
}

int ODE_dynamical_object::dynamics(ODE_vector const & state,
                                   ODE_vector & time_derivative,
                                   ODE_vector & derived_quantities){
    if( number_of_derived_quantities() != 0){
        FATAL_ERROR("number_of_derived_quantities() is nonzero ..." << std::endl <<
                                                                    "but dynamics with derived quantities has not been defined.");
    }else{
        FATAL_ERROR("Did not really expect code to pass through here. Does it make sense?");
        return dynamics(state,time_derivative);
    }

    return 0;
}

int ODE_dynamical_object::Jacobian(ODE_vector const & state,
                                   ODE_vector const & dynamics,
                                   ODE_matrix & jac){
    FATAL_ERROR("The Jacobian has not been implemented for this class.");
    return 0;
}

void ODE_dynamical_object::JTimes(ODE_vector const & state,
                                  ODE_vector const & in,
                                  ODE_vector & out){
    FATAL_ERROR(
            "Jacobian multiplication has not been implemented for this class."
    );
}

int ODE_dynamical_object::precondition(ODE_vector const & state,
				       ODE_vector const & in,
				       ODE_vector & out,
				       realtype gamma,
				       bool left_rather_than_right) const{
  FATAL_ERROR("No preconditioner defined!" << endl << 
	      "The preconditioner computes" << endl <<
	      "   out = P_i^(-1) in," << endl <<
	      "where i=1,2 and P=P_1*P_2 is an approximation of the Newton Matrix N " << endl <<
	      "   N = I - gamma J " << endl <<
	      "where J is the Jacobian of the dynamics f.");
  return -1;
}

int ODE_dynamical_object::precondition_IDA(ODE_vector const & state,
					   ODE_vector const & derived,
					   ODE_vector const & in,
					   ODE_vector & out,
					   realtype alpha) const{
  FATAL_ERROR("No IDA preconditioner defined!" << endl <<
	      "The preconditioner computes " << endl <<
	      "  out = P^{-1} in," << endl <<
	      "where P is a coarse approximation of the IDA Jacobian" << endl <<
	      "  alpha I - J " << endl <<
	      "and J represents the r.h.s. of the original dynamics."
	      );
  return -1;
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
static int root_f(realtype t, N_Vector y, realtype * gout, void *f_data);
int f_IDA(realtype t, N_Vector y, N_Vector ydot,
          N_Vector residuals, void *f_data);
static int root_f_IDA(realtype t, N_Vector y, N_Vector ydot, realtype * gout,
                      void *f_data);
static int Jac(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J,
               void *jac_data,
               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
static int PrecondSetupFn(realtype t, N_Vector y, N_Vector fy,
                          booleantype jok, booleantype * jcurPtr,
                          realtype gamma, void *user_data){
    *jcurPtr=0; //we did not recompute Jacobian data here
    return 0; //0 means we did the setup successfully (there was nothing
    //to be done)
}
static int PrecondSolveFn(realtype t, N_Vector y,
                          N_Vector fy, N_Vector r, N_Vector z,
                          realtype gamma, realtype delta,
                          int lr, void *p_data);

static int PrecondSolveFn_IDA(realtype t, N_Vector y, N_Vector yp,
			      N_Vector res, N_Vector r, N_Vector z,
			      realtype alpha, realtype delta,
			      void *p_data);


ODE_state::ODE_state():
        the_start_time(dummy_real),the_number_of_derived_quantities(0){
    FATAL_ERROR("cannot use default constructor");
};

ODE_state::ODE_state(const ODE_state & other):
        the_start_time(dummy_real),the_number_of_derived_quantities(0){
    FATAL_ERROR("cannot use copy constructor");
};

void ODE_state::prepare_integrator_CVode(){
    const int NEQ=this->size();
    int flag;

    y = N_VMake_Serial(size(),&((*this)[0]));

    abstol=default_absolute_tolerance;

    // The recommended choices here are (CV_ADAMS, CV_FUNCTIONAL) for
    // nonstiff problems and (CV_BDF, CV_NEWTON) for stiff problems:
#define MULTISTEP_METHOD CV_ADAMS // try CV_ADAMS or CV_BDF
#ifdef SUNDIALS_3
    sundials_mem = CVodeCreate(MULTISTEP_METHOD, CV_NEWTON);
#else
    sundials_mem = CVodeCreate(MULTISTEP_METHOD);
#endif    
    if (sundials_mem == NULL) { FATAL_ERROR("CVodeCreate failed."); }
    flag = CVodeInit(sundials_mem, f, 0, y);
    if (flag != CV_SUCCESS) { FATAL_ERROR("CVodeInit failed with code " <<
                                                                        flag ); }
    flag = CVodeSStolerances(sundials_mem, reltol, abstol);
    if (flag != CV_SUCCESS) { FATAL_ERROR("CVodeSStolerances failed with code " <<
                                                                                flag ); }

    if(precondition && use_default_preconditioner){
        FATAL_ERROR("Default preconditioner not implemented.");
        //flag = CVBand(sundials_mem, size(), size()-1, size()-1);
        if (flag != CV_SUCCESS) {
            FATAL_ERROR("CVBand failed with code " <<
                                                   flag );
        }
    }


    if(
            CVodeSetUserData(sundials_mem, the_dynamics)!=CV_SUCCESS ||
            CVodeSetErrFile(sundials_mem, stdout)!=CV_SUCCESS ||
            CVodeSetMaxNumSteps(sundials_mem,max_integrator_steps) ||
            CVodeSetMaxStep(sundials_mem,max_step_size) ||
            CVodeSetMinStep(sundials_mem,min_step_size) ||
            false ){
        FATAL_ERROR("setting optional CVODE parameters");
    }

    N_Vector constraints = N_VNew_Serial(NEQ);
    {
      ODE_vector ODE_constraints(constraints);
      int ret=the_dynamics->set_inequality_constraints(ODE_constraints);
      if(ret==0){
	bool some_non_zero=false;
	for(int i=0;i<NEQ;i++){
	  if(ODE_constraints[i]!=0){
	    some_non_zero=true;
	    break;
	  }
	}
	if(some_non_zero){
	  FATAL_ERROR("Function CVodeSetConstraints is not available, yet");
	  //	  flag = CVodeSetConstraints(sundials_mem, constraints);
	  if (flag != CV_SUCCESS)
	    FATAL_ERROR("Setting constraints failed.");
	}
      }
    }
    N_VDestroy_Serial(constraints);

    switch(linear_solver){
        case DIAG:
            flag = CVDiag(sundials_mem);
            if (flag != CVDIAG_SUCCESS) { printf("CVDiag failed.\n"); exit(1); }
            break;
        case SPGMR:
            //  pretype (int) specifies the preconditioning type and must be one
            //  of: PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH.
        {
            int pretype =
                    (precondition && the_dynamics->has_preconditioner() ?
                     PREC_LEFT : PREC_NONE );
            LS = SUNLinSol_SPGMR(y, PREC_NONE, 0);
            flag = CVSpilsSetLinearSolver(sundials_mem, LS);
            if (flag != CVLS_SUCCESS) { printf("CVSpgmr failed.\n"); exit(1); }
            flag = CVSpilsSetPreconditioner(sundials_mem, PrecondSetupFn, PrecondSolveFn);
            if (flag != CVLS_SUCCESS)
            { printf("setting preconditioner failed.\n"); exit(1); };
        }

            break;
        default:
            FATAL_ERROR("unknown linear_solver");
    } // switch(linear_solver)

    const int nrtfn = the_dynamics->number_of_root_functions();

    if(nrtfn > 0){
        flag = CVodeRootInit(sundials_mem, nrtfn, root_f);

        {
            int root_directions[nrtfn];
            if(0==the_dynamics->set_root_directions(root_directions)){
                flag = CVodeSetRootDirection(sundials_mem, root_directions);
                if(flag != CV_SUCCESS)
                    FATAL_ERROR("setting root directions");
            }
        }
    }

    CVodeSetInitStep(sundials_mem,ode_initial_step_size);
}

void ODE_state::prepare_integrator_IDA(){
    // Number of equations form perspective of IDA:
    const int NEQ = size();
    const int n = the_dynamics->number_of_variables();
    int flag;

    y = N_VMake_Serial(NEQ,&((*this)[0]));
    ydot = N_VNew_Serial(NEQ);
    
    N_Vector constraints = N_VNew_Serial(NEQ);
    bool some_constraints_non_zero=false;
    {
      ODE_vector ODE_constraints(&(NV_DATA_S(constraints)[0]),n);
      int ret=the_dynamics->set_inequality_constraints(ODE_constraints);
      if(ret==0){
	for(int i=0;i<n;i++){
	  if(ODE_constraints[i]!=0){
	    some_constraints_non_zero = true;
	    if(int(ODE_constraints[i]) != ODE_constraints[i])
	      FATAL_ERROR("Unknown constraint " << ODE_constraints[i] << ".");
	    switch(int(ODE_constraints[i])){
	    case -2:
	      if(! (NV_DATA_S(y)[i] < 0.0) ){
		WARNING("Elements " << i << " is " << NV_DATA_S(y)[i] << " and therefore violates constraint " << ODE_constraints[i] << ".");
	      }
	      break;
	    case -1:
	      if(! (NV_DATA_S(y)[i] <= 0.0) ){
		WARNING("Elements " << i << " is " << NV_DATA_S(y)[i] << " and therefore violates constraint " << ODE_constraints[i] << ".");
		WARNING("Setting to zero.");
		NV_DATA_S(y)[i] = 0;
	      }
	      break;
	    case 1:
	      if(! (NV_DATA_S(y)[i] >= 0.0) ){
		WARNING("Elements " << i << " is " << NV_DATA_S(y)[i] << " and therefore violates constraint " << ODE_constraints[i] << ".");
		WARNING("Setting to zero.");
		NV_DATA_S(y)[i] = 0;
	      }
	      break;
	    case 2:
	      if(! (NV_DATA_S(y)[i] > 0.0) ){
		WARNING("Elements " << i << " is " << NV_DATA_S(y)[i] << " and therefore violates constraint " << ODE_constraints[i] << ".");
	      }
	      break;
	    default:
	      FATAL_ERROR("Unknown constraint " << ODE_constraints[i] << ".");
	    }
	  }
	}
      }
    }

    IDA_target_ydot = N_VNew_Serial(n);
    IDA_target_derived = N_VNew_Serial(the_number_of_derived_quantities);

    if(!y || !IDA_target_ydot || !IDA_target_derived){
        FATAL_ERROR("Problem allocating N_Vectors");
    }

    abstol=default_absolute_tolerance;

    sundials_mem = IDACreate();
    if (sundials_mem == NULL) { FATAL_ERROR("IDACreate failed."); }

    {
        ODE_vector ODE_ydot(ydot); // adaptor
        ODE_vector derived(&(*this)[n], the_number_of_derived_quantities);
        if(the_dynamics->dynamics(*this,ODE_ydot,derived))
            FATAL_ERROR("Call to dynamics during initialization failed");

	for(int i=n; i<NEQ; i++){
	  ODE_ydot[i]=0;
	}

        flag = IDAInit(sundials_mem, f_IDA, 0, y, ydot);
        if (flag != IDA_SUCCESS)
            FATAL_ERROR("IDAInit failed with code " <<  flag );
    }

    flag = IDASStolerances(sundials_mem, reltol, abstol);
    if (flag != IDA_SUCCESS)
        FATAL_ERROR("IDASStolerances failed with code " << flag );

    { // declare which variables are dynamics, which derived:
        ODE_vector id(NEQ);
        for(int i=NEQ;i-->0;){
            id[i] = (i < n);
        }
        N_Vector id_NV = N_VMake_Serial(NEQ,&(id[0]));
        flag = IDASetId(sundials_mem, id_NV);
        if (flag != IDA_SUCCESS)
            FATAL_ERROR("IDASetId failed with code " << flag);
        N_VDestroy_Serial(id_NV);
    }

    if(precondition && use_default_preconditioner){
        FATAL_ERROR("Default preconditioner not implemented.");
        //flag = IDABand(sundials_mem, size(), size()-1, size()-1);
        // if (flag != IDA_SUCCESS) {
        //     FATAL_ERROR("IDABand failed with code " <<
        //                                             flag );
        // }
    }


    if(
            IDASetUserData(sundials_mem, this)!=IDA_SUCCESS ||
            IDASetErrFile(sundials_mem, stdout)!=IDA_SUCCESS ||
            IDASetMaxNumSteps(sundials_mem,max_integrator_steps) ||
            IDASetMaxStep(sundials_mem,max_step_size) ||
            false ){
        FATAL_ERROR("setting optional IDA parameters");
    }
    
    if(some_constraints_non_zero){
      for(int i=n;i<NEQ;i++){
	NV_DATA_S(constraints)[i]=0;
      }
      flag = IDASetConstraints(sundials_mem, constraints);
      if (flag != IDA_SUCCESS)
	FATAL_ERROR("Setting constraints failed.");
    }
    N_VDestroy_Serial(constraints);

    switch(linear_solver){
    case SPGMR:
      //  pretype (int) specifies the preconditioning type and must be one
      //  of: PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH.
      {
	int pretype =
	  (precondition && the_dynamics->has_preconditioner() ?
	   PREC_LEFT : PREC_NONE );
	LS = SUNLinSol_SPGMR(y, pretype, 0);
	flag = IDASpilsSetLinearSolver(sundials_mem, LS);
	if (flag != IDALS_SUCCESS)
	  FATAL_ERROR("IDASpilsSetLinearSolver failed.");
	if(pretype != PREC_NONE){
	  flag = IDASpilsSetPreconditioner(sundials_mem, NULL, PrecondSolveFn_IDA);
	  if (flag != IDALS_SUCCESS)
	    FATAL_ERROR("Setting preconditioner failed.\n");
	}
      }
      
      break;
    default:
      FATAL_ERROR("unknown linear_solver");
    } // switch(linear_solver)
    
    // This might be useful, but I'm not sure what it does.
    // flag = IDASetNonlinConvCoefIC(sundials_mem, 1e-10);
    // if (flag != IDA_SUCCESS)
    //     FATAL_ERROR("IDASetNonlinConvCoefIC failed with code " << flag);
    
    flag = IDACalcIC(sundials_mem, IDA_YA_YDP_INIT, /*tout1*/ double(1));
    if (flag != IDA_SUCCESS)
        FATAL_ERROR("IDACalcIC failed with code " << flag);
    // Note that IDACalcIC will correct the values of y(t0 ) and ẏ(t0)
    // which were specified in the previous call to IDAInit.

    const int nrtfn = the_dynamics->number_of_root_functions();
    if(nrtfn > 0){
        flag = IDARootInit(sundials_mem, nrtfn, root_f_IDA);
    }
    {
      int root_directions[nrtfn];
      if(0==the_dynamics->set_root_directions(root_directions)){
	flag = IDASetRootDirection(sundials_mem, root_directions);
	if(flag != CV_SUCCESS)
	  FATAL_ERROR("setting root directions");
      }
    }


    IDASetInitStep(sundials_mem,ode_initial_step_size);
}

void ODE_state::release_integrator_CVode(){
    CVodeFree(&sundials_mem);        /* Free the CVODE problem memory */
    SUNLinSolFree(LS);            /* Free linear solver */
    N_VDestroy_Serial(y);
}

void ODE_state::release_integrator_IDA(){
    IDAFree(&sundials_mem);        /* Free the IDA problem memory */
    SUNLinSolFree(LS);            /* Free linear solver */
    N_VDestroy_Serial(y);
    N_VDestroy_Serial(ydot);
    N_VDestroy_Serial(IDA_target_ydot);
    N_VDestroy_Serial(IDA_target_derived);
}

void ODE_state::prepare_integrator(){
    if(the_number_of_derived_quantities == 0){
        prepare_integrator_CVode();
    }else{
        prepare_integrator_IDA();
    }
}

void ODE_state::release_integrator(){
    if(the_number_of_derived_quantities == 0){
        release_integrator_CVode();
    }else{
        release_integrator_IDA();
    }
}


ODE_state::ODE_state(ODE_dynamical_object * object) :
        ODE_vector(object->number_of_variables()+object->number_of_derived_quantities()),
        the_time_since_start(0),
        the_start_time(object->the_start_time),
        reltol(default_absolute_tolerance),
        //  machEnv(M_EnvInit_Serial(this->size())),
        linear_solver(SPGMR), // fastest
        //linear_solver(DIAG),  // same speed as DENSE when I last tested
        the_dynamics(object),
        the_number_of_derived_quantities(object->number_of_derived_quantities())
{
    if(size()==0) FATAL_ERROR("ODE_dynamical_object without variables");
    object->prepare_for_integration();
    the_start_time=the_dynamics->current_time;
    object->write_state_to(*this);

    reltol = default_relative_tolerance;
    prepare_integrator();
}


ODE_state::~ODE_state()
{
    the_dynamics->read_state_from(*this);
    the_dynamics->current_time=the_start_time+the_time_since_start;
    //    diagnosis();
    release_integrator();
    the_dynamics->cleanup_after_integration();
    //  M_EnvFree_Serial(machEnv);   /* Free the machine environment memory */
}

realtype ODE_state::redo_with_shorter_step_size(){
    // (returns t from where we restart)
    if(the_number_of_derived_quantities == 0){
        realtype last_step_size;
        if(CVodeGetLastStep(sundials_mem,&last_step_size)!=CV_SUCCESS){
            FATAL_ERROR("Can't get last step size");
        }
        release_integrator();
        the_dynamics->write_state_to(*this);
        prepare_integrator();
        WARNING("reducing stepsize from " <<
                                          last_step_size << " to " << last_step_size/2);
        CVodeSetInitStep(sundials_mem,last_step_size/2);
        //CVodeSetMaxStep(...);
    }else{ // IDA
        FATAL_ERROR("Please implement redo_with_shorter_step_size for IDA.");
    }
    return the_start_time;
}

void ODE_state::restart(double t_start){
    release_integrator();
    the_start_time=t_start;
    the_dynamics->current_time=the_start_time;
    the_time_since_start=0;
    the_dynamics->read_state_from(*this);
    prepare_integrator();
}

int ODE_state::
integrate_until(realtype target_time){
    int flag;
    if(the_number_of_derived_quantities==0){
        flag = CVode(sundials_mem, target_time-the_start_time, y,
                     &the_time_since_start, CV_NORMAL);
        if (flag != CV_SUCCESS) {
            fprintf(stderr,"CVode failed, flag=%d.\n", flag);
            WARNING("integrator failed!");
            return 1;
        }
    }else{
        flag = IDASolve(sundials_mem, target_time-the_start_time,
                        &the_time_since_start, y, ydot, IDA_NORMAL);
        if (flag != IDA_SUCCESS) {
            fprintf(stderr,"IDASolve failed, flag=%d.\n", flag);
            WARNING("integrator failed!");
            return 1;
        }
    }
    the_dynamics->current_time=the_start_time+the_time_since_start;
    return 0;
}

int ODE_state::
integrate_to_root(realtype target_time, root_indices_t & root_indices){
    bool RETURN_DIRECTION=true;
    int flag;
    if(the_number_of_derived_quantities==0){
        flag = CVode(sundials_mem, target_time-the_start_time, y,
                     &the_time_since_start, CV_NORMAL);
        the_dynamics->current_time=the_start_time+the_time_since_start;
        if (flag == CV_ROOT_RETURN){ // we found a root
            const int nrtfn = the_dynamics -> number_of_root_functions();
            int rootsfound[nrtfn];
            CVodeGetRootInfo(sundials_mem, rootsfound);
            for(int i=0; i < nrtfn; ++i){
                if(rootsfound[i] && !RETURN_DIRECTION) root_indices.push_back(i);
                // return index with sign indication direction of zero crossing
                // +1 needed due to zero indexing
                if(rootsfound[i] && RETURN_DIRECTION) {
                    root_indices.push_back(rootsfound[i]*(i+1));
                }
            }
            return 0;
        }
        if (flag != CV_SUCCESS) {
            fprintf(stderr,"CVode failed, flag=%d.\n", flag);
            WARNING("integrator failed!");
            return 1;
        }
    }else{
        bool print_y=false;
        if (print_y) {
            cout << ODE_vector(y) << endl;
            cout << ODE_vector(ydot) << endl;
        }
        flag = IDASolve(sundials_mem, target_time-the_start_time,
                        &the_time_since_start, y, ydot, IDA_NORMAL);
        the_dynamics->current_time=the_start_time+the_time_since_start;
        if (flag == IDA_ROOT_RETURN){ // we found a root
            const int nrtfn = the_dynamics -> number_of_root_functions();
            int rootsfound[nrtfn];
            IDAGetRootInfo(sundials_mem, rootsfound);
            for(int i=0; i < nrtfn; ++i){
                if(rootsfound[i] && !RETURN_DIRECTION) root_indices.push_back(i);
                // return index with sign indication direction of zero crossing
                // +1 needed due to zero indexing
                if(rootsfound[i] && RETURN_DIRECTION) {
                    root_indices.push_back(rootsfound[i]*(i+1));
                }
            }
            return 0;
        }
        if (flag != IDA_SUCCESS) {
            fprintf(stderr,"IDASolve failed, flag=%d.\n", flag);
            WARNING("integrator failed!");
            return 1;
        }
    }
    return 0;
}

int ODE_state::
integrate_one_step(realtype target_time){
    int flag;
    if(the_number_of_derived_quantities == 0){
        flag = CVode(sundials_mem, target_time-the_start_time, y,
                     &the_time_since_start, CV_ONE_STEP);
        the_dynamics->current_time=the_start_time+the_time_since_start;
        if (flag != CV_SUCCESS) {
            fprintf(stderr,"CVode failed, flag=%d.\n", flag);
            WARNING("integrator failed!");
            return 1;
            //    FATAL_ERROR("Exiting");
        }
    }else{ // IDA
        flag = IDASolve(sundials_mem, target_time-the_start_time,
                        &the_time_since_start, y, ydot, IDA_ONE_STEP);
        if (flag != IDA_SUCCESS) {
            fprintf(stderr,"IDASolve failed, flag=%d.\n", flag);
            WARNING("integrator failed!");
            return 1;
            //    FATAL_ERROR("Exiting");
        }
    }
    return 0;
}

int ODE_state::
current_time_derivative(ODE_vector & ddt){
    int flag;
    if(the_number_of_derived_quantities == 0){
        N_Vector ddt_N_Vector = N_VMake_Serial(size(),&(ddt[0]));
        flag = CVodeGetDky(sundials_mem,
                           the_dynamics->current_time -
                           the_start_time,
                           1,
                           ddt_N_Vector);
        N_VDestroy(ddt_N_Vector);
        if (flag != CV_SUCCESS) {
            fprintf(stderr,"CVodeGetDky failed, flag=%d.\n", flag);
            WARNING("Could not compute time derivative!");
            return 1;
            //    FATAL_ERROR("Exiting");
        }
    }else{
        N_Vector ddt_N_Vector = N_VMake_Serial(size(),&(ddt[0]));
        flag = IDAGetDky(sundials_mem,
                         the_dynamics->current_time -
                         the_start_time,
                         1,
                         ddt_N_Vector);
        N_VDestroy(ddt_N_Vector);
        if (flag != IDA_SUCCESS) {
            fprintf(stderr,"IDAGetDky failed, flag=%d.\n", flag);
            WARNING("Could not compute time derivative!");
            return 1;
            //    FATAL_ERROR("Exiting");
        }
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
    CVodeGetNumSteps(sundials_mem,&i);
    co << i << " steps" << endl;
}

void ODE_state::diagnosis(ostream &co){
    long int i;
    int j;
    realtype x;
    co << "Diagnosis" << endl;
    if(the_number_of_derived_quantities == 0){
        co << "current_time " << the_dynamics->current_time << endl;
        CVodeGetNumSteps(sundials_mem,&i);
        co << i << " internal ODE steps" << endl;
        CVodeGetNumRhsEvals(sundials_mem, &i);
        co << i << " function evaluations" << endl;
        CVodeGetNumNonlinSolvIters(sundials_mem, &i);
        co << i << " nonlinear iterations" << endl;
        CVodeGetNumLinSolvSetups(sundials_mem, &i);
        co << i << " linear solver setups (since that's expensive)" << endl;
        if(linear_solver==SPGMR){
            CVSpilsGetNumPrecEvals(sundials_mem, &i);
            co << i << " preconditioner setups" << endl;
            CVSpilsGetNumPrecSolves(sundials_mem, &i);
            co << i << " preconditioner calls" << endl;
        }
        CVodeGetNumErrTestFails(sundials_mem, &i);
        co << i << " local error (accuracy) test failures" << endl;
        CVodeGetNumNonlinSolvConvFails(sundials_mem, &i);
        co << i << " nonlinear solver failures" << endl;
        CVodeGetLastOrder(sundials_mem, &j);
        co << j << "-" << (j==1?"st":j==2?"nd":j==3?"rd":"th") <<
           " order in final step" << endl;
        CVodeGetActualInitStep(sundials_mem, &x);
        co << x << " was the size of the first ODE step" << endl;
        CVodeGetLastStep(sundials_mem, &x);
        co << x << " was the size of the last ODE step" << endl;
    }else{ // IDA:
        co << "current_time " << the_dynamics->current_time << endl;
        IDAGetNumSteps(sundials_mem,&i);
        co << i << " internal ODE steps" << endl;
        IDAGetNumResEvals(sundials_mem, &i);
        co << i << " function evaluations" << endl;
        IDAGetNumNonlinSolvIters(sundials_mem, &i);
        co << i << " nonlinear iterations" << endl;
        IDAGetNumLinSolvSetups(sundials_mem, &i);
        co << i << " linear solver setups (since that's expensive)" << endl;
        if(linear_solver==SPGMR){
            IDASpilsGetNumPrecEvals(sundials_mem, &i);
            co << i << " preconditioner setups" << endl;
            IDASpilsGetNumPrecSolves(sundials_mem, &i);
            co << i << " preconditioner calls" << endl;
        }
        IDAGetNumErrTestFails(sundials_mem, &i);
        co << i << " local error (accuracy) test failures" << endl;
        IDAGetNumNonlinSolvConvFails(sundials_mem, &i);
        co << i << " nonlinear solver failures" << endl;
        IDAGetLastOrder(sundials_mem, &j);
        co << j << "-" << (j==1?"st":j==2?"nd":j==3?"rd":"th") <<
           " order in final step" << endl;
        IDAGetActualInitStep(sundials_mem, &x);
        co << x << " was the size of the first ODE step" << endl;
        IDAGetLastStep(sundials_mem, &x);
        co << x << " was the size of the last ODE step" << endl;
    }
}


//*********************************************************************
//*********************************************************************
//*********************************************************************

#include "kinsol/kinsol.h"

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
    if(KIN_SUCCESS!=KINInit(kin_mem, fixed_point_analyzer_f, tmpl)){
        FATAL_ERROR("KINMalloc failed.");
    }
    KINSetUserData(kin_mem, &user_data);
    KINSetNumMaxIters(kin_mem,1<<15);
    KINSetFuncNormTol(kin_mem,1e-9); //
    //  KINSetScaledStepTol(kin_mem, 1e-100);
    switch(linear_solver){
        case SPGMR:
            FATAL_ERROR("KINSpgmr disappeared from SUNDIALS (since v2.4.0?), interface needs reimplementing");
            // if(KINSPILS_SUCCESS!=KINSpgmr(kin_mem, 0)){
            //   FATAL_ERROR("Linear solver initialization failed");
            // }
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


int f_IDA(realtype t, N_Vector y, N_Vector ydot,
          N_Vector residuals, void *f_data)
{
    ODE_state * state=(ODE_state *)f_data;
    ODE_dynamical_object * obj=state->the_dynamics;

    ODE_vector ODE_y(y),ODE_ydot(ydot),ODE_residuals(residuals);
    ODE_vector ODE_target_ydot(state->IDA_target_ydot);
    ODE_vector ODE_target_derived(state->IDA_target_derived);

    obj->current_time=t+obj->the_start_time;
    int retval = obj->dynamics(ODE_y,ODE_target_ydot,ODE_target_derived);
    const int d = state->the_number_of_derived_quantities;
    const int n = state->size()-d;
    int i=0;
    for(; i < n; ++i){
        ODE_residuals[i] = ODE_ydot[i] - ODE_target_ydot[i];
    }
    for(; i < d+n; ++i){
        ODE_residuals[i] = ODE_y[i] - ODE_target_derived[i-n];
    }
    return retval;
}

static int PrecondSolveFn_IDA(realtype t, N_Vector y, N_Vector yp,
			      N_Vector res, N_Vector r, N_Vector z,
			      realtype alpha, realtype delta,
			      void *f_data){
  // we assume some very simple preconditioning here, so we pass only
  // a minimum of all these parameters over to the dynamical object's
  // preconditioner.
  ODE_state * state=(ODE_state *)f_data;
  ODE_dynamical_object * obj=state->the_dynamics;
  const int n = obj->number_of_variables();
  const int d = obj->number_of_derived_quantities();
  const ODE_vector ODE_y(&(NV_DATA_S(y)[0]),n);
  const ODE_vector ODE_derived(&(NV_DATA_S(y)[n]),d);
  const ODE_vector in(r);
  ODE_vector out(z);
  ODE_vector out_vars(&(out[0]),n);
  obj->current_time=t+obj->the_start_time;
  int retval =
    obj->precondition_IDA(ODE_y,ODE_derived,ODE_vector(&in[0],n),out_vars,alpha);
  if(retval != 0)
    return retval;
  // Jacobian for derived quantities is approximately the identity matrix:
  for(int i=n+d;i-->n;){
    out[i]=in[i];
  }
  return 0;
}


/* root function, with zeros indicatin times to stop */
static int root_f(realtype t, N_Vector y, realtype * gout,
                  void *f_data)
{
    ODE_dynamical_object * obj=(ODE_dynamical_object*)f_data;
    ODE_vector ODE_y(y),ODE_gout(gout, obj->number_of_root_functions());

    obj->current_time=t+obj->the_start_time;
    return obj->root_functions(ODE_y, ODE_gout);
}

static int root_f_IDA(realtype t, N_Vector y, N_Vector ydot, realtype * gout,
                      void *f_data)
{
    ODE_state * state=(ODE_state *)f_data;
    return root_f(t, y, gout, state->the_dynamics);
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

static int Jac(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J,
               void *jac_data,
               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
    ODE_vector ODE_y(y), ODE_ydot(fy);
    ODE_matrix ODE_J(J);

    return ((ODE_dynamical_object*)jac_data)->Jacobian(ODE_y,ODE_ydot,ODE_J);
}

static int PrecondSolveFn(realtype t, N_Vector y,
                          N_Vector fy, N_Vector r, N_Vector z,
                          realtype gamma, realtype delta,
                          int lr, void *p_data){
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

//*********************************************************************
//*********************************************************************
//*********************************************************************

// Manage adjustable parameters:
//#include "cfgList.h"
//static cfgStruct cfg[] =
//{
//       /* parameter,          type,       address of variable */
//        {"PRECISION_OF_VECTOR_PRINT", CFG_INT, &precision_of_vector_print},
//        {"DEFAULT_ABSOLUTE_TOLERANCE", CFG_DOUBLE,&default_absolute_tolerance},
//        {"DEFAULT_RELATIVE_TOLERANCE", CFG_DOUBLE,&default_relative_tolerance},
//	{"MAX_INTEGRATOR_STEPS",CFG_INT ,&max_integrator_steps},
//	{"MAX_STEP_SIZE",CFG_DOUBLE ,&max_step_size},
//	{"MIN_STEP_SIZE",CFG_DOUBLE ,&min_step_size},
//	CFGINT(precondition),
//	CFGINT(use_default_preconditioner),
//	CFGDOUBLE(ode_initial_step_size),
//	{0, CFG_END, 0}   /* no more parameters */
//};
//static cfg_add dummy(cfg);

