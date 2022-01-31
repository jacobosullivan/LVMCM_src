#include "ODE.h"
#include <armadillo>  // for integrating_dynamical_object
#include <iomanip>      // std::setw
#include <fenv.h>

using namespace arma;
using namespace std;

class testODE : public ODE_dynamical_object {
public: // naked private:
    double x;
    double v;
    double bounce_factor;
    const double gravity = 0.2;
public:
    virtual int dynamics(ODE_vector const & state,
                         ODE_vector & time_derivative);
    virtual int dynamics(ODE_vector const & state,
                         ODE_vector & time_derivative,
                         ODE_vector & derived_quantities);
    virtual void write_state_to(ODE_vector & state) const;
    virtual void read_state_from(const ODE_vector & state);
    virtual int number_of_variables() const {return 2;};
    virtual int number_of_derived_quantities() {return 1;}; //comment to switch
    testODE(double x,double v, double bounce_fric=double(1));
    // root finding stuff:
    virtual int number_of_root_functions(){return 1;}
    virtual int root_functions(ODE_vector const & state,
                               ODE_vector const & derived,
                               ODE_vector & gout){
        gout[0] = state[0]; // find root if x hits zero;
        return 0;
    }
    void react_to_roots(ODE_state::root_indices_t & root_indices){
        if( !root_indices.empty() && root_indices.back()==0){
            v = bounce_factor * abs(v); // bounce upwards;
            x = 0;
            root_indices.pop_back();
        }
    }
    double energy() const{
        return x*gravity + v*v/2;
    }
};

testODE::testODE(double xx, double vv, double ff):
        x(xx),
        v(vv),
        bounce_factor(1-ff)
{}

void testODE::read_state_from(const ODE_vector & state){
    x = state[0];
    v = state[1];
}

void testODE::write_state_to(ODE_vector & state) const {
    state[0] = x;
    state[1] = v;
}

int testODE::dynamics(ODE_vector const & state,
                      ODE_vector & time_derivative){
    ODE_vector derived(1);  // memory allocation, would be better outside the loop
    int retval = dynamics(state,time_derivative,derived);
    return retval;
}

int testODE::dynamics(ODE_vector const & state,
                      ODE_vector & time_derivative,
                      ODE_vector & derived_quantities){
    time_derivative[0] = state[1];
    time_derivative[1] = - gravity; // acceleration by gravity
    derived_quantities[0] = state[1]*state[1]/2; // kinetic energy
    return 0;
}

class integrating_dynamical_object : public ODE_dynamical_object {
    ODE_dynamical_object & main;
    vec integrals;
    integrating_dynamical_object();
public:
    int dynamics(ODE_vector const & state,
                 ODE_vector & time_derivative) {
        int retval = main.dynamics(state,time_derivative);
        int n = main.number_of_variables();
        for(int i = 0; i < n; ++i){
            time_derivative[n+i]=state[i];
        }
        return(retval);
    }
    void write_state_to(ODE_vector & state) const {
        main.write_state_to(state);
        int n = main.number_of_variables();
        for(int i = 0; i < n; ++i){
            state[n+i]=integrals[i];
        }
    }
    void read_state_from(const ODE_vector & state) {
        main.read_state_from(state);
        int n = main.number_of_variables();
        for(int i = 0; i < n; ++i){
            integrals[i]=state[n+i];
        }
    }
    int number_of_variables() const {
        return(2*main.number_of_variables());
    }
    vec get_integrals() const {return integrals;}
    integrating_dynamical_object(ODE_dynamical_object & other) :
            main(other) , integrals(other.number_of_variables(),fill::zeros) {
    }
    virtual void prepare_for_integration(){
        main.prepare_for_integration();
    };
    virtual void cleanup_after_integration(){
        main.cleanup_after_integration();
    };
    virtual int number_of_root_functions(){
        return main.number_of_root_functions();
    }
    virtual int root_functions(ODE_vector const & state,
                               ODE_vector & gout){
        return main.root_functions(state, gout);
    }
};

#define COMPUTE_INTEGRALS 0

int main(){

    if (0) {

        feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

        testODE myball_raw(/*height: */15,/*speed: */0,/*bounce friction: */0.0);
#if COMPUTE_INTEGRALS
        integrating_dynamical_object myball(myball_raw);
#else
        testODE &myball = myball_raw;
#endif

        const double tMax = 100;

        cout << "Simulate a bouncing ball" << endl;
        cout << setprecision(3) << setw(5) << "time" << " ";
        cout << setprecision(5) << setw(11) << "height" << " speed";
        cout << endl;

        while (myball.current_time < tMax && myball_raw.energy() > 1e-9) {
            ODE_state::root_indices_t root_indices;
            {
                ODE_state state(&myball);

                cout << setprecision(3) << setw(5) << myball.current_time << " ";
                cout << setprecision(5) << setw(11) << state;
                cout << endl;

                for (double tNext = ceil(myball.current_time + 0.1);
                     tNext <= tMax;
                     tNext = ceil(myball.current_time + 0.1)) {
                    state.integrate_to_root(tNext, root_indices);
                    cout << setprecision(3) << setw(5) << myball.current_time << " ";
                    cout << setprecision(5) << setw(11) << state;
                    // See if kinetic energy is reproduced:
                    cout << setprecision(5) << setw(11) << state[1] * state[1] / 2;
                    cout << endl;
                    ODE_vector ydot(myball.number_of_variables() +
                                    myball.number_of_derived_quantities());
                    state.current_time_derivative(ydot);
                    cout << ydot << endl;
                    if (!root_indices.empty()) break; // root was found
                }
//	    state.diagnosis(cout);  // use to get number of function evals
            }
            myball_raw.react_to_roots(root_indices);  // <- outside of ODE_state scope
        }
    } else { // testing ground
        uvec A = regspace<uvec>(4, 1);  // 4, 3, 2, 1
        uvec B = regspace<uvec>(3, 6);  // 3, 4, 5, 6
        A.print("A");
        B.print("B");

        uvec C = intersect(A,B);       // 3, 4
        C.print("C");

        uvec CC;
        uvec iA;
        uvec iB;

        intersect(CC, iA, iB, A, B);

        CC.print("CC");
        iA.print("ia");
        iB.print("iB");

    }

    exit(0);


}

// Local Variables:
// c-file-style: "stroustrup"
// End:
