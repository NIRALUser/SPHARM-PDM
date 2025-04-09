#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include "EqualAreaParametricMeshNewtonIterator.h"

// sparse matrix library stuff
#include "cg.h"           // IML++ Conjugate gradient template
#include "icpre_double.h" // IC preconditioner

const int no_activation = -1;
const double quite_bad = 1000;

#ifndef M_PI
#define M_PI 3.1415926535897932
#endif

EqualAreaParametricMeshSparseMatrix::EqualAreaParametricMeshSparseMatrix() {}

EqualAreaParametricMeshSparseMatrix::~EqualAreaParametricMeshSparseMatrix() {}

/**
 * Define the structure of the matrix given the surface net and active constraints. Does not set any values.
 */
void EqualAreaParametricMeshSparseMatrix::from_net(const IteratorSurfaceNet &net, int n_active, int *active) {
  int n_row = net.nface - 1 + n_active;
  int n_col = 3 * net.nvert;

  // Reuse the same storage and avoid redundant allocations.
  static std::vector<Triplet> triplets;
  triplets.resize(0);

  int i;
  for (i = 0; i < net.nface - 1; i++) {
    int row = i;
    for (int corner = 0; corner < 4; corner++) {
      int vertNr = net.face[4 * i + corner];
      for (int coord = 0; coord < 3; coord++) {
        int col = 3 * vertNr + coord;
        triplets.emplace_back(row, col, 1);
      }
    }
  }

  for (i = 0; i < n_active; i++) {
    int row = i + net.nface - 1;
    int faceNr = active[i] / 4;
    for (int corner = 0; corner < 4; corner++) {
      if ((active[i] + 4 - corner) % 4 != 2) // opposite corner has no influence
      {
        int vertNr = net.face[4 * faceNr + corner];
        for (int coord = 0; coord < 3; coord++) {
          int col = 3 * vertNr + coord;
          triplets.emplace_back(row, col, 1);
        }
      }
    }
  }

  mat.resize(n_row, n_col);
  mat.setFromTriplets(triplets.begin(), triplets.end());
}

// void EqualAreaParametricMeshNewtonIterator::new_jacobian(EqualAreaParametricMeshSparseMatrix &wrapper) {
//   const double desired_area = 4 * M_PI / net.nface;
//
//   using Matrix = EqualAreaParametricMeshSparseMatrix::Matrix;
//   using Iterator = EqualAreaParametricMeshSparseMatrix::Iterator;
//
//   Matrix &A = wrapper.mat;
//   Eigen::Map<Eigen::Matrix3Xd> x(this->m_x, 3, net.nvert);
//
//   // Reuse the same storage for temporary vector data throughout.
//   static Eigen::Matrix3Xd tmpx;
//   tmpx = x;
//
//   assert(A.cols() == x.size());
//
//   for (int idx = 0; idx < A.cols(); idx++) {
//     // Remember that `idx` here corresponds to a _coordinate_ in `x`. Divide to find the current _vector_.
//
//     int xcol = (idx / 3);
//     tmpx(idx) += par.delta;
//     tmpx.col(xcol).normalize(); // project back to sphere
//
//     // finite differences in this coordinate
//     tmpx(idx) += par.delta;
//
//     Iterator area_it(A, idx);
//     Iterator ineq_it(A, idx);
//
//     // scan forward to find the inequalities part of this column
//     for (; ineq_it and ineq_it.row() < net.nface - 1; ++ineq_it) {}
//
//     // iterate over the areas part of this column
//     for (; area_it and area_it.row() < net.nface - 1; ++area_it) {
//       int row = area_it.row();
//       double sines[4];
//       double area_c = spher_area4(tmpx.data(), net.face + 4 * row, sines) - desired_area;
//       area_it.valueRef() = (area_c - c_hat[row]) / par.delta;
//
//       // find the inequalities that correspond to this face.
//       for (; ineq_it and (ineq_it.row() - (net.nface - 1)) / 4 == row; ++ineq_it) {
//         ineq_it.valueRef() =
//             (sines[(ineq_it.row() - (net.nface - 1)) % 4] - c_hat[ineq_it.row()]) / par.delta;
//       }
//     }
//
//     // Extra inequalities. Face `net.nface - 1` does not have any `spher_area4` constraint so it cannot be
//     // handled in the loop above, but it may have inequality constraints, so populate those here.
//     if (ineq_it) {
//       double sines[4];
//       spher_area4(tmpx.data(), net.face + 4 * (net.nface - 1), sines);
//       for (; ineq_it; ++ineq_it) {
//         ineq_it.valueRef() = (sines[(ineq_it.row() - (net.nface - 1)) % 4] - c_hat[ineq_it.row()]) / par.delta;
//       }
//     }
//
//     // undo finite difference step
//     tmpx.col(xcol) = x.col(xcol);
//   }
//
//   // verify that mat has correct data.
// #ifndef NDEBUG
//   Eigen::SparseMatrix<double, Eigen::ColMajor> leg(wrapper.n_row, wrapper.n_col);
//   leg.reserve(wrapper.ia[wrapper.n_row]);
//   for (int row = 0; row < wrapper.n_row; ++row) {
//     for (int i = wrapper.ia[row]; i < wrapper.ia[row + 1]; ++i) {
//       int col = wrapper.ja[i];
//       leg.insert(row, col) = wrapper.a[i];
//     }
//   }
//   leg.makeCompressed();
//
//   std::cout << wrapper.mat - leg << std::endl;
//   assert(wrapper.mat.isApprox(leg));
// #endif
// }

/**
 * Define the values of the matrix to be the jacobian by finite differences. Assume the structure has already
 * been defined with ::from_net.
 */
void EqualAreaParametricMeshNewtonIterator::jacobian(EqualAreaParametricMeshSparseMatrix &A) {
  const double desired_area = 4 * M_PI / net.nface;

  Eigen::Block<Eigen::SparseMatrix<double, Eigen::ColMajor>> top = A.mat.topRows(net.nface - 1);
  Eigen::Block<Eigen::SparseMatrix<double, Eigen::ColMajor>> bot = A.mat.bottomRows(n_active);

  for (int col = 0; col < A.mat.cols(); col++) // assume x == this->m_x_try
  {
    this->m_x_try[col] = this->m_x[col] + par.delta; // go a finite step
    int col_0 = 3 * (col / 3);                       // column rounded down: x-component
    normalize(1, 3, this->m_x_try + col_0);

    Eigen::InnerIterator itop(top, col);
    Eigen::InnerIterator ibot(bot, col);

    double sines[4];
    for (; itop; ++itop) {
      int trow = itop.row();

      double area_c = spher_area4(this->m_x_try, net.face + 4 * trow, sines) - desired_area;

      top.coeffRef(trow, col) = (area_c - c_hat[trow]) / par.delta;

      for (; ibot; ++ibot) {
        int brow = ibot.row();
        int cid = active[brow];
        if ((cid / 4) != trow)
          break;

        bot.coeffRef(brow, col) = (sines[cid % 4] - c_hat[brow + top.rows()]) / par.delta;
      }
    }

    for (; ibot; ++ibot) {
      int brow = ibot.row();
      int cid = active[brow];
      spher_area4(this->m_x_try, net.face + 4 * (net.nface - 1), sines) - desired_area;
      bot.coeffRef(brow, col) = (sines[cid % 4] - c_hat[brow + net.nface - 1]) / par.delta;
    }

    this->m_x_try[col_0] = this->m_x[col_0]; // back up the finite step
    this->m_x_try[col_0 + 1] = this->m_x[col_0 + 1];
    this->m_x_try[col_0 + 2] = this->m_x[col_0 + 2];
  }
}

void EqualAreaParametricMeshSparseMatrix::mult(double *vec, double *result) {
  using Map = Eigen::Map<Eigen::VectorXd>;
  Map(result, mat.rows()) = mat * Map(vec, mat.cols());
}

void EqualAreaParametricMeshSparseMatrix::multT(double *vec, double *result) {
  using Map = Eigen::Map<Eigen::VectorXd>;
  Map(result, mat.cols()) = mat.transpose() * Map(vec, mat.rows());
}

void EqualAreaParametricMeshSparseMatrix::solve(int /* structure_change */, double *b, double *x) {
  using namespace Eigen;
  using Map = Eigen::Map<Eigen::VectorXd>;

  ConjugateGradient<SparseMatrix<double>> cg;
  cg.setMaxIterations(5000);
  cg.setTolerance(1e-14);
  cg.compute(mat);

  Map(x, mat.rows()) = cg.solve(Map(b, mat.cols()));

  std::cout << "i" << cg.iterations();
}

void EqualAreaParametricMeshSparseMatrix::set_aTa(const EqualAreaParametricMeshSparseMatrix &aT) {
  mat = aT.mat * aT.mat.transpose();
}

void EqualAreaParametricMeshNewtonIterator::set_parameters(EqualAreaParametricMeshParameter *par_in) {

  par.max_active = par_in->max_active;
  par.print_itn = par_in->print_itn;
  // par.reread_param = par_in->reread_param;
  par.delta = par_in->delta;
  par.constr_tol = par_in->constr_tol;
  par.line_tol = par_in->line_tol;
  par.ineq_low = par_in->ineq_low;
  par.ineq_init = par_in->ineq_init;
  par.ineq_final = par_in->ineq_final;
  par.ineq_slack = par_in->ineq_slack;
  par.newton_tol = par_in->newton_tol;
  par.rho_init = par_in->rho_init;
  par.c0rho = par_in->c0rho;
  par.c1rho = par_in->c1rho;
  par.c2rho = par_in->c2rho;
  par.rho_limit = par_in->rho_limit;
  par.step_small = par_in->step_small;
  par.step_large = par_in->step_large;
}

int EqualAreaParametricMeshNewtonIterator::check_constraints(int n_equal, int n_inequal, /*int *active,*/
                                                             double *inequal,
                                                             /* double* c_hat,*/
                                                             double *c_hat_try, double &badness) {
  int act_i = 0;

  min_ineq = 1;     // init with max any ineq could reach
  int min_pos = -1; // where the minimum was found
  for (int i = 0; i < n_inequal; i++) {
    if (active[act_i] == i) {
      c_hat_try[n_equal + act_i++] = inequal[i];
    } else if (inequal[i] < min_ineq) {
      min_ineq = inequal[i];
      min_pos = i;
    }
  }
  if (min_ineq < -ineq_high) {
    return min_pos;
  }
  // assert(active[act_i] == -1, "there must be a sentinel.");
  // assert(act_i == n_active, "all actives must be visited.");

  badness = 0;
  // int worst = -1;
  for (act_i = n_equal + n_active; act_i-- > 0;) // also over positive inequalities(?)
  {
    double increase = sqr(c_hat_try[act_i]) - sqr(c_hat[act_i]);
    if (increase > badness) {
      badness = increase;
      // worst = act_i;
    }
  }
  // if (badness > 1)
  //  std::cout << " worst: " << worst << "(" << badness << ") ";
  return no_activation;
}

void EqualAreaParametricMeshNewtonIterator::activate_initial(const int n_equal, const int n_inequal,
                                                             const double *inequal) {
  std::cout << "initially active: {";
  n_active = 0;
  for (int i = 0; i < n_inequal; i++) {
    if (inequal[i] < -ineq_high) {
      std::cout << (n_active > 0 ? ", " : "") << i;
      active[n_active] = i;
      c_hat[n_equal + n_active++] = inequal[i];
      activity[i] = 6;
    } else {
      activity[i] = 0;
    }
  }
  active[-1] = active[n_active] = no_activation; // mark start and end with a sentinel
  std::cout << "}.\n";
}

void EqualAreaParametricMeshNewtonIterator::constraints(const IteratorSurfaceNet & /* net */,
                                                        const double *_x /* x */, double *equal,
                                                        double *inequal) {
  const double desired_area = 4 * M_PI / net.nface;
  for (int face = 0; face < net.nface; face++) {
    // determine inequal through side effect
    equal[face] = spher_area4(_x, net.face + 4 * face, inequal + 4 * face) - desired_area;
  }
}

double EqualAreaParametricMeshNewtonIterator::one_inequality(const IteratorSurfaceNet & /* net */,
                                                             const double *x_ /* x */, int which) {
  double sines[4];
  (void)spher_area4(x_, net.face + which / 4 * 4, sines);
  return sines[which % 4];
}

EqualAreaParametricMeshNewtonIterator::EqualAreaParametricMeshNewtonIterator(
    const IteratorSurfaceNet &netR, EqualAreaParametricMeshParameter &par_in)
    : par(par_in), net(netR) {

  count = 0;
  last_complete = -1; // any value < count
  this->m_rho = par.rho_init;
  this->m_alpha_step = 1.0;

  this->m_x = new double[3 * net.nvert];
  this->m_x_try = new double[3 * net.nvert];
  this->m_newton_dir = new double[3 * net.nvert];
  grad = new double[3 * net.nvert];
  gradY = new double[3 * net.nvert];
  gCG = new double[3 * net.nvert];
  hCG = new double[3 * net.nvert];
  constr_ineq = new double[4 * net.nface]; // values of the inequality constraints

  const int max_eq_act = net.nface - 1 + par.max_active;

  this->m_lambda = new double[max_eq_act]; // Lagrange multipliers
  this->m_proj_dx = new double[max_eq_act];
  aTgrad = new double[max_eq_act];
  c_hat = new double[max_eq_act]; // values of eql. & active constraints
  c_hat_l = new double[max_eq_act];
  c_hat_r = new double[max_eq_act];
  active = new int[max_eq_act + 2] + 1; // first and last for sentinels
  activity = new char[4 * net.nface];   // activity levels of inequalities

  start_values(net, this->m_x);
  copy_vector(gCG, this->m_x, 3 * net.nvert);
  copy_vector(this->m_x_try, this->m_x, 3 * net.nvert);
  for (int i = 0; i < 3 * net.nvert; hCG[i++] = 0.0) {
    ;
  }
  this->m_dx = hCG; // dx is set to this->m_newton_dir or hCG

  constraints(net, this->m_x, c_hat, constr_ineq);
  ineq_high = par.ineq_init;
  activate_initial(net.nface - 1, 4 * net.nface, constr_ineq);

  std::cout << "fromC[expr_] := expr /. (mant_Real)*e + (expon_Integer) :> mant*10^expon;\n";
  std::cout.flush();
} /* EqualAreaParametricMeshNewtonIterator (constructor) */

EqualAreaParametricMeshNewtonIterator ::EqualAreaParametricMeshNewtonIterator(
    const IteratorSurfaceNet &netR, EqualAreaParametricMeshParameter &par_in, double *initPar)
    : par(par_in), net(netR) {

  count = 0;
  last_complete = -1; // any value < count
  this->m_rho = par.rho_init;
  this->m_alpha_step = 1.0;

  this->m_x = new double[3 * net.nvert];
  this->m_x_try = new double[3 * net.nvert];
  this->m_newton_dir = new double[3 * net.nvert];
  grad = new double[3 * net.nvert];
  gradY = new double[3 * net.nvert];
  gCG = new double[3 * net.nvert];
  hCG = new double[3 * net.nvert];
  constr_ineq = new double[4 * net.nface]; // values of the inequality constraints

  const int max_eq_act = net.nface - 1 + par.max_active;

  this->m_lambda = new double[max_eq_act]; // Lagrange multipliers
  this->m_proj_dx = new double[max_eq_act];
  aTgrad = new double[max_eq_act];
  c_hat = new double[max_eq_act]; // values of eql. & active constraints
  c_hat_l = new double[max_eq_act];
  c_hat_r = new double[max_eq_act];
  active = new int[max_eq_act + 2] + 1; // first and last for sentinels
  activity = new char[4 * net.nface];   // activity levels of inequalities

  copy_vector(this->m_x, initPar, 3 * net.nvert);
  copy_vector(gCG, this->m_x, 3 * net.nvert);
  copy_vector(this->m_x_try, this->m_x, 3 * net.nvert);
  for (int i = 0; i < 3 * net.nvert; hCG[i++] = 0.0) {
    ;
  }
  this->m_dx = hCG; // dx is set to this->m_newton_dir or hCG

  constraints(net, this->m_x, c_hat, constr_ineq);
  ineq_high = par.ineq_init;
  activate_initial(net.nface - 1, 4 * net.nface, constr_ineq);

  std::cout << "fromC[expr_] := expr /. (mant_Real)*e + (expon_Integer) :> mant*10^expon;\n";
  std::cout.flush();
} /* EqualAreaParametricMeshNewtonIterator (constructor) */

EqualAreaParametricMeshNewtonIterator::~EqualAreaParametricMeshNewtonIterator() {
  delete[] m_x;
  delete[] this->m_x_try;
  delete[] this->m_dx;
  delete[] this->m_lambda;
  delete[] this->m_proj_dx;
  delete[] this->m_newton_dir;
  delete[] grad;
  delete[] gradY;
  delete[] gCG;
  delete[] hCG;
  delete[] constr_ineq;
  delete[] aTgrad;
  delete[] c_hat;
  delete[] c_hat_l;
  delete[] c_hat_r;
  delete[] active;
  delete[] activity;
}

double EqualAreaParametricMeshNewtonIterator::iterate() {
  char form[1 << 10];

  sprintf(form, "count=%3d ", count);
  std::cout << form; // print iteration number

  int struct_change = last_complete != count++;

  // ------------------------------------------------ Newton: satisfy constraints -------
  if (struct_change) {
    jacobi_aT.from_net(net, n_active, active);
  }
  jacobian(jacobi_aT);
  jacobi_aTa.set_aTa(jacobi_aT);
  jacobi_aTa.solve(struct_change, c_hat, this->m_proj_dx);
  jacobi_aT.multT(this->m_proj_dx, this->m_newton_dir);

  int act, act_keep = no_activation;
  double badness;

  this->m_dx = this->m_newton_dir;
  double newtonStep;
  for (newtonStep = -1;
       (act = step_and_check(newtonStep, c_hat_l, badness)) != no_activation || badness > par.newton_tol;
       newtonStep /= 2) {
    act_keep = act;
  }

  sprintf(form, " step=%6.0e act_keep=%6d", newtonStep, act_keep);
  std::cout << form; // debug
  copy_vector(this->m_x, this->m_x_try, 3 * net.nvert); // accept new x
  // // write_vector("x", 3*net.nvert, x, 0);      // debug
  double *help = c_hat;
  c_hat = c_hat_l;
  c_hat_l = help;

  // ------------------------------------------------ Lagrange multipliers --------------
  calc_gradient();
  jacobian(jacobi_aT); // re-calculate after Newton-step
  jacobi_aT.mult(grad, aTgrad);
  jacobi_aTa.set_aTa(jacobi_aT); // re-calculate as well
  jacobi_aTa.solve(0, aTgrad, this->m_lambda);

  int inactivated = 0;
  int i;
  // for( i = 0; i < n_active; i++ )
  for (i = n_active - 1; i >= 0; i--) {
    if (this->m_lambda[net.nface - 1 + i] < 0 && c_hat[net.nface - 1 + i] > -par.ineq_low) {
      inactivated += inactivate(i);
    }
  }
  if (inactivated > 0) {
    activate(act_keep, "in_out"); // activate Newton's stopper anyway
    if (act_keep != no_activation) {
      std::cout << "\n";
    }
    return quite_bad;
  }

  // ------------------------------------------------ Conjugate gradients ---------------
  jacobi_aT.multT(this->m_lambda, gradY);
  double numerator = 0.0, denominator = 0.0, norm2gradZ = 0.0;
  for (i = 0; i < 3 * net.nvert; i++) {
    double gradZ = gradY[i] - grad[i]; // only one component is kept at a time
    numerator += (gradZ - gCG[i]) * gradZ;
    denominator += sqr(gCG[i]);
    norm2gradZ += sqr(gradZ);
    gCG[i] = gradZ;
  }
  // // write_vector("gCG", 3*net.nvert, gCG, 0);  // debug
  double gamma = numerator / denominator;
  for (i = 0; i < 3 * net.nvert; i++) {
    hCG[i] = gCG[i] + gamma * hCG[i];
  }
  // // write_vector("hCG", 3*net.nvert, hCG, 0);  // debug

  // ------------------------------------------------ Linear minimization ---------------
  double c_sqr_sum;
  this->m_dx = hCG;
  // // estimate_jacobian();
  sprintf(form, " gamma=%6.3f |gradZ|=%8.1e", gamma, sqrt((double)norm2gradZ));
  std::cout << form;

  int activated = line_search(this->m_alpha_step, c_sqr_sum);
  activated += activate(act_keep, "newton");
  if (activated == 0) {
    last_complete = count;
  }

  if (-min_ineq * par.ineq_slack < ineq_high) {
    ineq_high = -min_ineq * par.ineq_slack;
  }

  if (par.ineq_final > ineq_high) {
    ineq_high = par.ineq_final;
  }

  double cost = sqrt((double)c_sqr_sum) + norm2gradZ;
  sprintf(form, " |c_hat|=%10.3e cost=%10.3e ineq_high=%10.3e\n", sqrt((double)c_sqr_sum), cost, ineq_high);
  std::cout << form;
  return cost;
}

void EqualAreaParametricMeshNewtonIterator::get_solution(double *copy) {
  copy_vector(copy, this->m_x, 3 * net.nvert);
}

double EqualAreaParametricMeshNewtonIterator::goal_func() {
  double goal = 0.0;

  for (int v = 0; v < net.nvert; v++) {
    IteratorSurfaceVertex *vertex = net.vert + v;
    for (int nb = 0; nb < vertex->count; nb += 2) {
      goal += 1.0 - dotproduct3(this->m_x_try + 3 * v, this->m_x_try + 3 * vertex->neighb[nb]);
    }
  }
  return goal / 2.0;
}

void EqualAreaParametricMeshNewtonIterator::calc_gradient() {
  for (int v = 0; v < net.nvert; v++) {
    double nbsum[3];
    nbsum[0] = nbsum[1] = nbsum[2] = 0.0;
    IteratorSurfaceVertex *vertex = net.vert + v;
    for (int nb = 0; nb < vertex->count; nb += 2) {
      double *neighbor = this->m_x + 3 * vertex->neighb[nb];
      nbsum[0] += neighbor[0];
      nbsum[1] += neighbor[1];
      nbsum[2] += neighbor[2];
    }
    double prod = dotproduct3(this->m_x + 3 * v, nbsum);
    grad[3 * v] = prod * this->m_x[3 * v] - nbsum[0];
    grad[3 * v + 1] = prod * this->m_x[3 * v + 1] - nbsum[1];
    grad[3 * v + 2] = prod * this->m_x[3 * v + 2] - nbsum[2];
  }
}

int EqualAreaParametricMeshNewtonIterator::activate(int act, const char *info) {
  if (act == no_activation) {
    return 0;
  }
  // assert(act >= 0 && act < 4*net.nvert, "activate: act out of range");
  std::cout << "\n]]] " << info << " : constraint " << act << " -> level " << activity[act] + 1
            << "                            ";
  if (++activity[act] != 3) {
    return 0;
  }
  activity[act] = 5;
  int i;
  for (i = n_active; active[i - 1] > act; i--) {
    active[i] = active[i - 1];
    c_hat[net.nface - 1 + i] = c_hat[net.nface - 1 + i - 1];
  }
  c_hat[net.nface - 1 + i] = one_inequality(net, this->m_x_try, act);
  active[i] = act;
  active[++n_active] = -1; // sentinel for check_constraints
  std::cout << ">>> " << info << " activates constraint " << act << " (pos " << i
            << "); now active: " << n_active
            << "\n                                                                   ";
  // assert(active[i-1] < act, "only an inactive constraint may be activated");
  // assert(c_hat[net.nface-1+i] >= -ineq_high,
  //     "the activated constraint must not be violated");
  // assert(n_active < par.max_active, "too many constraints are active at the same time");
  return 1;
}

int EqualAreaParametricMeshNewtonIterator::inactivate(int pos) {
  char form[1000];

  // assert(pos >= 0 && pos < n_active, "inactivate: pos out of range");
  sprintf(form, "%6d", active[pos]);
  std::cout << "\n[[[ constraint " << form << " -> level " << activity[active[pos]] - 1 << ";          ";
  if (--activity[active[pos]] != 3) {
    return 0;
  }
  activity[active[pos]] = 2;
  n_active--;
  std::cout << "<<< inactivate constraint " << active[pos] << " at postition " << pos
            << "; now active: " << n_active << "\n";
  for (int i = pos; i < n_active; i++) {
    active[i] = active[i + 1];
    c_hat[net.nface - 1 + i] = c_hat[net.nface - 1 + i + 1];
  }
  active[n_active] = -1; // sentinel for check_constraints
  return 1;
}

int EqualAreaParametricMeshNewtonIterator::line_search(double &fullStep, double &c_sqr_sum) {
  char form[1000];
  double m = 0.0;
  double f_m, f_l, f_r, c2s_l, c2s_r, bad_m, bad_l, bad_r;
  int viol_l, viol_r, act_l, act_m, act_r;
  int act_keep_l = no_activation, act_keep_r = no_activation;

  f_m = aug_lagrangian(m, act_m, bad_m, c_sqr_sum, c_hat);
  // assert (act_m == no_activation, "line_search: activation at m = 0");
  // assert (bad_m == 0, "line_search: bad_m should be 0.");
  for (double stepSize = par.step_large; stepSize > par.line_tol; stepSize /= 2.0) {
    double step = fullStep * stepSize;
    f_l = aug_lagrangian(m - step, act_l, bad_l, c2s_l, c_hat_l);
    f_r = aug_lagrangian(m + step, act_r, bad_r, c2s_r, c_hat_r);
    viol_l = act_l != no_activation || bad_l > par.constr_tol;
    viol_r = act_r != no_activation || bad_r > par.constr_tol;
    if (!viol_l && f_l < f_m && (viol_r || f_l < f_r)) {
      f_m = f_l;
      bad_m = bad_l;
      m -= step;
      act_keep_r = no_activation;
    } else if (!viol_r && f_r < f_m && (viol_l || f_r < f_l)) {
      f_m = f_r;
      bad_m = bad_r;
      m += step;
      act_keep_l = no_activation;
    } else {
      act_keep_r = act_r;
      act_keep_l = act_l;
    }
  }

  // not the most efficient solution:
  double dum_bad;
  f_m = aug_lagrangian(m, act_m, dum_bad, c_sqr_sum, c_hat); // update also this->m_x_try
  // assert (act_m == no_activation, "line_search: no activation must occur at the end");
  // assert (dum_bad == 0, "line_search: dummy badness should be 0.");

  copy_vector(this->m_x, this->m_x_try, 3 * net.nvert); // accept new x
  if (m != 0.0) {
    this->m_rho =
        this->m_rho * (par.constr_tol * par.c1rho / (par.constr_tol * par.c0rho - bad_m) + par.c2rho) +
        par.rho_limit;
  }

  sprintf(form, " %6.0e %10.6f %8.1e %8.1e", bad_m, f_m, this->m_rho, m);
  std::cout << form;

  if (activate(act_keep_l, "left") | activate(act_keep_r, "right")) {
    return 1;
  }
  // else adjust stepsize
  fullStep *= par.step_small;
  if (fabs(m) > fullStep) {
    fullStep = fabs(m);
  }

  return 0;
}

int EqualAreaParametricMeshNewtonIterator::step_and_check(double step, double *c_hat_try, double &badness) {
  spher_step(step, this->m_x, this->m_dx, net.nvert, this->m_x_try);
  constraints(net, this->m_x_try, c_hat_try, constr_ineq);
  return check_constraints(net.nface - 1, 4 * net.nface, /*active,*/ constr_ineq,
                           /*c_hat,*/ c_hat_try, badness);
}

double EqualAreaParametricMeshNewtonIterator::aug_lagrangian(double step, int &_activate, double &badness,
                                                             double &c_sqr_sum, double *c_hat_try) {
  _activate = step_and_check(step, c_hat_try, badness);
  if (_activate != no_activation) {
    return 0; // new active inequalities
  }
  double lagr = goal_func();
  c_sqr_sum = 0;
  for (int i = 0; i < net.nface - 1 + n_active; i++) {
    lagr -= this->m_lambda[i] * c_hat_try[i];
    if (i < net.nface - 1 || c_hat_try[i] < 0) {
      c_sqr_sum += sqr(c_hat_try[i]);
    }
  }
  return lagr + this->m_rho * c_sqr_sum;
}

// void EqualAreaParametricMeshNewtonIterator::estimate_jacobian() {
//   char form[1000];
//   const double stepSize = 0.001;
//   int i;
//   double *old_dx = this->m_dx, dum_bad, *a_row = new double[net.nface - 1 + n_active];
//
//   this->m_dx = new double[3 * net.nvert];
//   std::cout << "a = Table[,{" << 3 * net.nvert << "}];\n";
//   copy_vector(this->m_x_try, this->m_x, 3 * net.nvert);
//   normalize(net.nvert, 3, this->m_x_try); // almost no effect; for exact equality
//   for (i = 0; i < 3 * net.nvert; this->m_dx[i++] = 0.0) {}
//
//   for (i = 0; i < 3 * net.nvert; i++) {
//     this->m_dx[i] = 1.0;
//     int act = step_and_check(stepSize, c_hat_l, dum_bad);
//     if (act == no_activation) {
//       for (int j = 0; j < net.nface - 1 + n_active; j++) {
//         a_row[j] = (c_hat_l[j] - c_hat[j]) / stepSize;
//       }
//       sprintf(form, "a[[%d]]", i + 1);
//       write_vector(form, net.nface - 1 + n_active, a_row, 0);
//     } else {
//       std::cout << "a[[" << i + 1 << "]] = activate[" << act << "];\n";
//     }
//     this->m_dx[i] = 0.0;
//   }
//   copy_vector(this->m_x_try, this->m_x, 3 * net.nvert); // return to initial position
//   delete this->m_dx;
//   delete[] a_row;
//   this->m_dx = old_dx;
// }
//
// void EqualAreaParametricMeshNewtonIterator::estimate_gradient() {
//   const double stepSize = 0.00001;
//   double *old_dx = this->m_dx, *gradEstim = new double[3 * net.nvert];
//
//   this->m_dx = new double[3 * net.nvert];
//   int i;
//   for (i = 0; i < 3 * net.nvert; this->m_dx[i++] = 0.0) {
//     ;
//   }
//   spher_step(0, this->m_x, this->m_dx, net.nvert, this->m_x_try);
//   double goal = goal_func();
//   for (i = 0; i < 3 * net.nvert; i++) {
//     this->m_dx[i] = 1.0;
//     spher_step(stepSize, this->m_x, this->m_dx, net.nvert, this->m_x_try);
//     this->m_dx[i] = 0.0;
//     gradEstim[i] = (goal_func() - goal) / stepSize;
//   }
//   copy_vector(this->m_x_try, this->m_x, 3 * net.nvert); // return to initial position
//   write_vector("gradEstim", 3 * net.nvert, gradEstim, 0);
//   delete this->m_dx;
//   this->m_dx = old_dx;
// }

struct CompRows {
  int n;
  int *ia;
  int *ja;
  double *a;
};

void EqualAreaParametricMeshNewtonIterator::generate_matrix(struct CompRows *mat,
                                                            const IteratorSurfaceNet &_net) {
  const int north_pole = 0, south_pole = _net.nvert - 1;
  int n_nonzero = 0;
  int v;

  for (v = 1; v < south_pole; v++) {
    n_nonzero += _net.vert[v].count / 2 + 1; // + 1 diagonal element
  }
  mat->n = _net.nvert - 2;
  mat->ja = new int[n_nonzero];
  mat->ia = new int[mat->n + 1];
  mat->a = new double[n_nonzero];
  int n_used = 0;
  for (v = 1; v < south_pole; v++) {
    IteratorSurfaceVertex *vertex = _net.vert + v;
    mat->ia[v - 1] = n_used; // _net_index = matrix_index+1
    mat->ja[n_used] = v - 1;
    mat->a[n_used++] = vertex->count / 2; // put diagonal
    for (int nb = 0; nb < vertex->count; nb += 2) {
      int neighbour = vertex->neighb[nb];
      if (neighbour != north_pole && neighbour != south_pole) {
        mat->ja[n_used] = neighbour - 1;
        mat->a[n_used++] = -1.0;
      }
    }
  }
  mat->ia[mat->n] = n_used;
}

void EqualAreaParametricMeshNewtonIterator::modify_matrix(struct CompRows *mat,
                                                          const IteratorSurfaceNet &_net)
// modify the matrix for longitude computation by removing the
// connections to the poles and increasing the diagonal element of an
// arbitrary row (we select the first) by an arbitrary value (we select
// 2.0) that makes the matrix nonsingular. This allows us to specify a
// desired value for the longitude of vertex 1. Changing only the diagonal
// element preserves the symmetry of the matrix.
{
  for (int p = 0; p < 2; p++) {
    IteratorSurfaceVertex *pole = _net.vert + p * (_net.nvert - 1); // north (p=0) or south (p=1)
    for (int nb = 0; nb < pole->count; nb += 2) {
      mat->a[mat->ia[pole->neighb[nb] - 1]]--;
    }
  }
  mat->a[0] += 2.0; // make matrix regular: longi[1] = 0
}

void EqualAreaParametricMeshNewtonIterator::set_lati_rhs(const IteratorSurfaceNet &_net, double *rhs) {
  for (int i = 0; i < _net.nvert; rhs[i++] = 0.0) {
    ;
  }
  IteratorSurfaceVertex *south = _net.vert + _net.nvert - 1;
  for (int nb = 0; nb < south->count; nb += 2) {
    rhs[south->neighb[nb]] = M_PI;
  }
}

int *EqualAreaParametricMeshNewtonIterator::cycle2(int **pos, IteratorSurfaceVertex *vert)
// Advance a pointer to any of vert's neighbour cursors
// cyclically by 2 positions.
{
  *pos += 2;
  if (*pos >= vert->neighb + vert->count) {
    *pos -= vert->count;
  }
  return *pos;
}

void EqualAreaParametricMeshNewtonIterator::set_longi_rhs(const IteratorSurfaceNet &_net, double *rhs,
                                                          double *lati) {
  for (int i = 0; i < _net.nvert; i++) {
    rhs[i] = 0;
  }
  int prev = 0; // northpole
  int here = 1; // any neighbor of northpole
  double maximum = 0.0;
  while (here != _net.nvert - 1) {
    IteratorSurfaceVertex *hereV = _net.vert + here;
    int *maxpos = 0, *prevpos = 0;
    int *nb;
    for (nb = hereV->neighb; nb < hereV->neighb + hereV->count; nb += 2) {
      if (lati[*nb] > maximum) {
        maximum = lati[*nb];
        maxpos = nb;
      } else if (*nb == prev) {
        prevpos = nb;
      }
    }
    for (nb = prevpos; cycle2(&nb, hereV) != maxpos;) {
      rhs[*nb] += 2.0 * M_PI;
      rhs[here] -= 2.0 * M_PI; // totally: 2 pi * #east_neighbours
    }
    prev = here;
    here = *maxpos;
  }
}

void EqualAreaParametricMeshNewtonIterator::start_values(const IteratorSurfaceNet &_net, double *cartesian) {
  struct CompRows mat;

  double *rhs = new double[_net.nvert];
  double *lati = new double[_net.nvert];
  double *longi = new double[_net.nvert];
  int iter;
  double tol = 1e-14;
  int n_nonzero;

  generate_matrix(&mat, _net);
  set_lati_rhs(_net, rhs);
  n_nonzero = mat.ia[mat.n];

  CompRow_Mat_double A(mat.n, mat.n, n_nonzero, mat.a, mat.ia, mat.ja);
  // call of compressed row matrix = nxn with n_nonzero items that are non-zero, mat.a = value vector,
  // mat.ia = row pointer, mat.ja = column index
  // Build up the structure from a grid in 0-relative compressed row format

  ICPreconditioner_double D(A);
  // CompRow_ILUPreconditioner_double D(A);
  // DiagPreconditioner_double D(A);

  iter = 5000; // max iter, actual iterations returned in this variable

  VECTOR_double b(mat.n, 0.0);
  VECTOR_double x(mat.n, 0.0);
  int i;
  for (i = 0; i < mat.n; i++) {
    b[i] = rhs[i + 1]; // copy rhs
  }
  int result = CG(A, x, b, D, iter, tol); // A * x = b, max it, tol  6
  // int result = CGS(A, x, b, D ,iter, tol); // A * x = b, max it, tol BAD
  // int result = IR(A, x, b, D ,iter, tol); // A * x = b, max it, tol
  // int result = BiCGSTAB(A, x, b, D ,iter, tol); // A * x = b, max it, tol 9
  // int result = BiCG(A, x, b, D ,iter, tol); // A * x = b, max it, tol 12

  std::cout << "flag = " << result << " iterations performed: " << iter << std::endl;
  for (i = 0; i < mat.n; i++) {
    lati[i + 1] = x[i]; // copy result
  }
  lati[0] = 0.0;
  lati[_net.nvert - 1] = M_PI;

  modify_matrix(&mat, _net);
  set_longi_rhs(_net, rhs, lati);
  n_nonzero = mat.ia[mat.n];
  CompRow_Mat_double A2(mat.n, mat.n, n_nonzero, mat.a, mat.ia, mat.ja);

  ICPreconditioner_double D2(A2);
  // CompRow_ILUPreconditioner_double D2(A2);
  // DiagPreconditioner_double D2(A2);

  tol = 1e-14;
  iter = 5000; // max iter, actual iterations returned in this variable
  for (i = 0; i < mat.n; i++) {
    b[i] = rhs[i + 1]; // copy rhs
  }
  int result2 = CG(A2, x, b, D2, iter, tol); // A2 * x = b, max it, tol
  // int result2 = CGS(A2, x, b, D2 ,iter, tol); // A2 * x = b, max it, tol
  // int result2 = IR(A2, x, b, D2 ,iter, tol); // A2 * x = b, max it, tol
  // int result2 = BiCGSTAB(A2, x, b, D2 ,iter, tol); // A2 * x = b, max it, tol
  // int result2 = BiCG(A2, x, b, D2 ,iter, tol); // A2 * x = b, max it, tol

  std::cout << "flag = " << result2 << " iterations performed: " << iter << std::endl;
  for (i = 0; i < mat.n; i++) {
    longi[i + 1] = x[i]; // copy result
  }
  longi[0] = longi[_net.nvert - 1] = 0.0;
  for (i = 0; i < _net.nvert; i++) {
    cartesian[i * 3 + 0] = sin(lati[i]) * cos(longi[i]);
    cartesian[i * 3 + 1] = sin(lati[i]) * sin(longi[i]);
    cartesian[i * 3 + 2] = cos(lati[i]);
  }
}

double EqualAreaParametricMeshNewtonIterator::det3(const double *a, const double *b, const double *c) {
  return a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) +
         a[2] * (b[0] * c[1] - b[1] * c[0]);
}

double EqualAreaParametricMeshNewtonIterator::dotproduct3(const double *a, const double *b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void EqualAreaParametricMeshNewtonIterator::normalize(const int nvectors, const int dim, double *x) {
  for (int v = 0; v < nvectors; v++) {
    double length2 = sqr(x[v * dim]);
    int c;
    for (c = 1; c < dim; c++) {
      length2 += sqr(x[v * dim + c]);
    }
    double length = sqrt((double)length2);
    for (c = 0; c < dim; c++) {
      x[v * dim + c] /= length;
    }
  }
}

double EqualAreaParametricMeshNewtonIterator::spher_area4(const double *x, const int corner[4],
                                                          double spat[4]) {
  const double *a = x + 3 * corner[0];
  const double *b = x + 3 * corner[1];
  const double *c = x + 3 * corner[2];
  const double *d = x + 3 * corner[3];

  const double ab = dotproduct3(a, b);
  const double ac = dotproduct3(a, c);
  const double ad = dotproduct3(a, d);
  const double bc = dotproduct3(b, c);
  const double bd = dotproduct3(b, d);
  const double cd = dotproduct3(c, d);

  const double Ca = ad * ab - bd;
  const double Cb = ab * bc - ac;
  const double Cc = bc * cd - bd;
  const double Cd = cd * ad - ac;

  spat[0] = det3(d, a, b);
  spat[1] = det3(a, b, c);
  spat[2] = det3(b, c, d);
  spat[3] = det3(c, d, a);

  double area = std::atan2(Ca, spat[0]) + std::atan2(Cb, spat[1]) + std::atan2(Cc, spat[2]) + std::atan2(Cd, spat[3]);
  return std::remainder(area, M_PI);
}

void EqualAreaParametricMeshNewtonIterator::spher_step(double step, double *src, double *vec, int nvect,
                                                       double *dest) {
  for (int i = 3 * nvect; i--;) {
    dest[i] = src[i] + step * vec[i];
  }
  normalize(nvect, (long int)3, dest);
}

void EqualAreaParametricMeshNewtonIterator::copy_vector(double *dest, const double *src, const int length) {
  for (int i = 0; i < length; i++) {
    dest[i] = src[i];
  }
}
