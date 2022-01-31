// -*- mode: c++ -*-
// $Id: ODE.h 2519 2020-05-15 08:27:03Z axel $

/// \file: Interface to a stiff ODE solver (here sundials):

#ifndef __ODE_H__
#define __ODE_H__

#include <iostream>
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver      */
#include <sundials/sundials_dense.h>
#include <utility>
#include <vector>
#include <math.h>
#include <setjmp.h>
#include "error.h"

/* for default preconditioner: */

extern int newton_failure; // would ideally be local to each ODE_state

// Tests for sundials < 2.4.0
#define SUNDIALS_PRE_2_4_0 (!defined(SUNDIALS_EXPORT))
//#define SUNDIALS_VERSION_2_4

/// Encapsulates the vector type of the ODE solver
class  ODE_vector {
    int the_length;
    // this is usually double:
    realtype * the_elements;
    bool the_elements_are_mine;
public:
    ODE_vector(int length=0);
    ODE_vector(N_Vector vec); // N_Vector is a pointer!
    ODE_vector(realtype * elements,int length);
    ODE_vector(const ODE_vector & other);
    ~ODE_vector();
    inline double& operator[](int i){
        ASSERT(0 <= i && i < the_length);
        return the_elements[i];
    }
    inline double& operator[](int i) const{
        ASSERT(0 <= i && i < the_length);
        return the_elements[i];
    }
    int size() const{
        return the_length;
    }
    void clear();
    // this is only for exceptional use:
    const ODE_vector& operator=(ODE_vector const &other);
    ODE_vector& operator+=(ODE_vector const &other);
    ODE_vector& operator-=(ODE_vector const &other);
    ODE_vector& operator*=(double x);
};

ODE_vector exp(ODE_vector & v);

/// Encapsulates the matrix type of the ODE solver
class  ODE_matrix {
    // this is usually double:
    realtype ** the_elements;
    bool the_elements_are_mine;
    int the_length;
private:
    //this is only for exeptional use:
    const ODE_matrix& operator=(ODE_matrix const &other){
        FATAL_ERROR("cannot assign ODE_matrixes");
        return *this;
    }
public:
    ODE_matrix(int length=0);
    ODE_matrix(DlsMat vec);
    ODE_matrix(const ODE_matrix & other);
    ~ODE_matrix();
/*   double* & operator[](int i) const{ */
/*     return the_elements[i]; */
/*   } */
    double* operator[](int i) const{
        return the_elements[i];
    }
    int size(){
        return the_length;
    }
    void clear();
};

/// Abstract base class for describing the system to be simulated.
class ODE_dynamical_object
{
public:
    virtual int dynamics(ODE_vector const & state,
                         ODE_vector & time_derivative)=0;
    virtual int dynamics(ODE_vector const & state,
                         ODE_vector & time_derivative,
                         ODE_vector & derived_quantities);
    virtual int Jacobian(ODE_vector const & state,
                         ODE_vector const & dynamics,
                         ODE_matrix & jac);
    virtual int root_functions(ODE_vector const & state,
                               ODE_vector & gout);
    virtual int root_functions(ODE_vector const & state,
                               ODE_vector const & derived_quantities,
                               ODE_vector & gout);
    virtual void JTimes(ODE_vector const & state,
                        ODE_vector const & in,
                        ODE_vector & out);
    virtual int precondition(ODE_vector const & state,
			     ODE_vector const & in,
			     ODE_vector & out,
			     realtype gamma,
			     bool left_rather_than_right) const;
    virtual int precondition_IDA(ODE_vector const & state,
				 ODE_vector const & derived,
				 ODE_vector const & in,
				 ODE_vector & out,
				 realtype alpha) const;
  // constraint value of -2, -1, 0, 1, 2 mean y<0, y<=0, any y, y>=0, y>0
    virtual int set_inequality_constraints(ODE_vector & constraints) const {return -1;};
    virtual bool can_calculate_Jacobian(){return false;};
    virtual bool has_preconditioner(){return false;};
    virtual bool has_inherent_rates(){return false;};
    virtual void get_inherent_rates(ODE_vector & rates){};//hook
    virtual void write_state_to(ODE_vector & state) const=0;
    virtual void read_state_from(const ODE_vector & state)=0;
    virtual int number_of_variables() const =0;
    virtual int number_of_root_functions() const;
    virtual int set_root_directions(int * directions) const;
    virtual int number_of_derived_quantities();
    virtual void line_print(ODE_vector const & state,std::ostream &co);
    virtual void prepare_for_integration(){};  //hook
    virtual void cleanup_after_integration(){}; //hook
    template<typename MATRIX> int
    numerical_Jacobian(MATRIX & Jac,double dx=sqrt(DBL_EPSILON));
    void test_Jacobian();
    double current_time; //set dynamically
    double the_start_time; //set indirectly via ODE_state
    ODE_dynamical_object():current_time(0){};
    sigjmp_buf the_sigjmp_buf; // This is for "setjmp / longjump" trickery.
};

template<typename MATRIX>
int ODE_dynamical_object::
numerical_Jacobian(MATRIX & Jac,double dx){
    // use a central differences scheme
    const double norm=1/(2*dx);
    const int n=number_of_variables();
    ODE_vector state(n);
    write_state_to(state);

    prepare_for_integration();
    ODE_vector state1(n),time_derivative1(n);
    ODE_vector time_derivative2(n);
    state1=state;
    for(int i=n;i-->0;){
        state1[i]=state[i]+dx;
        dynamics(state1,time_derivative1);
        state1[i]=state[i]-dx;
        dynamics(state1,time_derivative2);
        state1[i]=state[i];
        for(int k=n;k-->0;){
            Jac[k][i]=(time_derivative1[k]-time_derivative2[k])*norm;
        }
    }

    cleanup_after_integration();
    return 1;
}



/// Encapsulates the ODE solver state and operations on it.
/** The constructor takes an ODE_dynamical_object as argument.  While
    the ODE_state exists, it "owns" this ODE_dynamical_object, keeping
    the values of the dependent variables in an ODE_vector.  At the
    moment where the ODE_state is destroyed, these dependent variables
    are written back to the ODE_dynamical_object, which can then be
    analyzed for the structure of the final state.*/
class ODE_state : public ODE_vector{
public:
    ODE_dynamical_object * the_dynamics;
private:
    ODE_state(); // outputs error message
    ODE_state(const ODE_state & other); // outputs error message
public:
    ODE_state(ODE_dynamical_object * object);
    ~ODE_state();
    int integrate_until(realtype target_time);
    int integrate_one_step(realtype target_time);
    typedef std::vector<int> root_indices_t;
    int integrate_to_root(realtype target_time, root_indices_t & root_indices);
    int current_time_derivative(ODE_vector & ddt);
    int Yoshida_step(double dt);
    void diagnosis(std::ostream &co=std::cout);
    void short_diagnosis(std::ostream &co=std::cout);
    realtype redo_with_shorter_step_size(); // returns t from where we restart;
    void restart(double t_start); // re-initialize integrator
    realtype time_since_start(){return the_time_since_start;}
    realtype start_time(){return the_start_time;}
private:
    static realtype dummy_real;
    realtype the_time_since_start;
    realtype & the_start_time;
    realtype reltol;
    void *sundials_mem;
    SUNLinearSolver LS;
    N_Vector y;  // will be identified with base ODE_vector
    N_Vector ydot;  // for IDA only
#if SUNDIALS_PRE_2_4_0
    N_Vector abstol;
#else
    realtype abstol;
#endif
    const int the_number_of_derived_quantities;
    N_Vector IDA_target_ydot; // memory allocated for passing over target ydot
    N_Vector IDA_target_derived; // memory allocated for passing over derived quantities
    typedef enum {DIAG,SPGMR} linear_solver_t;
    linear_solver_t linear_solver;
    void prepare_integrator_CVode();
    void release_integrator_CVode();
    void prepare_integrator_IDA();
    void release_integrator_IDA();
    void prepare_integrator();
    void release_integrator();
    friend int f_IDA(realtype t, N_Vector y, N_Vector ydot,
                     N_Vector residuals, void *f_data);
};

/// Encapsulates a nonlinear equation solver for finding fixed points of ODEs.
/** The constructor takes an ODE_dynamical_object as argument.  While
    the fixed_point_analyzer exists, it "owns" this
    ODE_dynamical_object, keeping the values of the dependent
    variables in an ODE_vector.  At the moment where the
    fixed_point_analyzer is destroyed, these dependent variables are
    written back to the ODE_dynamical_object, which can then be
    analyzed for the structure of the final state.  The class can be
    extended to do more analyses of the fixed point. */
class fixed_point_analyzer : public ODE_vector{
    static realtype dummy_real;
//   realtype ftol;
//   realtype steptol;
    void *kin_mem;
    N_Vector scaling_vector_D_u;
    N_Vector scaling_vector_D_F;
    N_Vector tmpl; // template vector, needed by KINMalloc
    N_Vector u;
    N_Vector constraints;
    typedef enum {SPGMR,SPBCG,SPTFQMR} linear_solver_t;
    linear_solver_t linear_solver;
    ODE_vector lower_bound;
    double lower_bound_shift;
public:
    ODE_dynamical_object * the_dynamics;
    typedef std::pair<ODE_dynamical_object *,ODE_vector *> user_data_t;
private:
    user_data_t user_data;
    fixed_point_analyzer(){FATAL_ERROR("cannot use default constructor");};
    fixed_point_analyzer(const fixed_point_analyzer & other){FATAL_ERROR("cannot use copy constructor");};
public:
    fixed_point_analyzer(ODE_dynamical_object * object,
                         double lower_bound_shift_request=-5.0*log(10.0));
    ~fixed_point_analyzer();
    void set_time_scale_to(double dt);
    int snap_to_fixed_point();
//   void diagnosis(std::ostream &co=std::cout);
//   void short_diagnosis(std::ostream &co=std::cout);
private:
    void prepare_solver();
    void release_solver();
public:
    class failure {};
    static const int F_IS_NAN;
};

#include <iosfwd>
std::ostream & operator<<(std::ostream &stream, const ODE_state & av);
std::ostream & operator<<(std::ostream &stream, const ODE_vector & vector);

#endif // __ODE_H__
