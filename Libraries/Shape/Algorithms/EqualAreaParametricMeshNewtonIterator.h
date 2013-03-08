// HEADER FILE EqualAreaParametricMeshEqualAreaParametricMeshNewtonIterator.h
// last modified: August 05  by Martin Styner
//

#ifndef __iterator_hh_

struct    IteratorSurfaceVertex    { int x, y, z, count, neighb[14]; };

struct    IteratorSurfaceNet    { int nvert, nface, *face; IteratorSurfaceVertex *vert; };

class EqualAreaParametricMeshSparseMatrix
{
private:
  int max_col, max_row, max_nonzero;
public:
  int     n_row, n_col, *iaT, *jaT, *rowT;
  int *   ia, *ja;
  double *a;

  EqualAreaParametricMeshSparseMatrix(int, int, int);
  ~EqualAreaParametricMeshSparseMatrix();
  void from_net(const IteratorSurfaceNet &, int n_active,  int *active_scatter);

  void invTables();                  // set up inverted tables

  void mult(double *vec, double *result);       // multiply this.vec

  void multT(double *vec, double *result);      // multiply this^T.vec

  void solve(int change, double *rhs, double *x); // solve this.x == rhs

  void set_aTa(const EqualAreaParametricMeshSparseMatrix& aT);          // set this = aT . a

  void print(const char* name, const int append);      // writes matrix to cout

  void test_matrix();

};

typedef struct
  {
  int
    max_active,                     // max. #of active inequalities
    print_itn;                      // iter #for debug output; <0 : never
  //    reread_param;                        // bool: re-read every 20 itns
  double
    delta,                                    // step for finite differences
    constr_tol,                               // -tolerable constraint violation
    line_tol,                                 // tolerance in line search
    ineq_low,                                 // c > -ineq_tol is ok ...
    ineq_init,                                // initial value for ineq_high
    ineq_final,                               // gradually reduce ineq_high to this
    ineq_slack,                               // but not too tightly (slack > 1)
    newton_tol,                               // max badness accepted by Newton step
    rho_init, c0rho, c1rho, c2rho, rho_limit, // for adapting rho in line_search
    step_small, step_large;                   // factors adjust step in line_search

  } EqualAreaParametricMeshParameter;

class EqualAreaParametricMeshNewtonIterator
{
private:
  EqualAreaParametricMeshParameter & par;
  double                             ineq_high; // ... 2 values for hysteresis
  double                             min_ineq;  // value of worst inactive inequality
  // dangerous: file scope variable!

  double *                            m_x;
  double *                            m_x_try;
  double *                            m_newton_dir;
  double *                            m_dx;
  double *                            m_proj_dx;
  double *                            m_lambda;
  double                              m_rho;
  double                              m_alpha_step;
  double *                            gCG;
  double *                            hCG;
  double *                            constr_ineq, *grad, *gradY, *aTgrad;
  int                                 n_active;
  int *                               active;
  double *                            c_hat;
  double *                            c_hat_l;
  double *                            c_hat_r;
  int                                 count;
  int                                 last_complete;
  const IteratorSurfaceNet &          net;
  EqualAreaParametricMeshSparseMatrix jacobi_aT;
  EqualAreaParametricMeshSparseMatrix jacobi_aTa;
  char *                              activity;

  int     line_search(double & step, double & c_sqr_sum);

  int     step_and_check(double step, double *c_hat_try, double & badness);

  double aug_lagrangian(double step, int & activate, double & badness, double & c_sqr_sum, double *c_hat_try);

  int    check_constraints(int n_equal, int n_inequal, double *inequal, double* c_hat_old, double & badness);

  void   activate_initial(const int n_equal, const int n_inequal, const double *inequal);

  void     calc_gradient();

  void     jacobian(const EqualAreaParametricMeshSparseMatrix &);

  void     estimate_gradient();

  void     estimate_jacobian();

  int     activate(int act, const char *); // returns 0 if no_activation, else 1

  int     inactivate(int);              // returns 0 if delayed, else 1

  void   constraints(const IteratorSurfaceNet& net, const double *x, double *equal, double *inequal);

  double one_inequality(const IteratorSurfaceNet & net, const double *x, int which);

  void write_YSMP(const char* name, int n_row, int n_col, int* ia, int* ja, double* a, int append_flag);

  void write_vector(const char* name, int n, double* vector, int append_flag);

  // computes the initial parametrization using the heat equation stuff
  void start_values(const IteratorSurfaceNet &, double *);

  void set_longi_rhs(const IteratorSurfaceNet& net, double * rhs, double * lati);

  int * cycle2(int* * pos, IteratorSurfaceVertex* vert);

  void set_lati_rhs(const IteratorSurfaceNet& net,   double* rhs);

  void modify_matrix(struct CompRows *mat, const IteratorSurfaceNet& net);

  void generate_matrix(struct CompRows *mat, const IteratorSurfaceNet& net);

  // utility functions
  double    dotproduct3(const double *, const double *);

  void    normalize(const int nvectors, const int dimension, double *x);

  double spher_area4(const double *x, const int corner[4], double angle[4]);
  void    spher_step(double step, double *src, double *vec,  int, double *dest);

  void    copy_vector(double *dest, const double *src, const  int length);

  double det3(const double *, const double *, const double *);

  inline double sqr(double val)
  {
    return val * val;
  }

public:
  EqualAreaParametricMeshNewtonIterator(const IteratorSurfaceNet& netref, EqualAreaParametricMeshParameter& par);
  EqualAreaParametricMeshNewtonIterator(const IteratorSurfaceNet& netref, EqualAreaParametricMeshParameter& par,
                                        double *initPar);
  ~EqualAreaParametricMeshNewtonIterator();
  double iterate();

  double goal_func();

  void   get_solution(double *);

  void   set_parameters(EqualAreaParametricMeshParameter *para);          // call once before making iterators

};

#endif
