// In-process IPOPT solver exposed to Python through nanobind.
//
// One function is exported: solve_nlp(...). It receives the problem data (bounds,
// starting point, Jacobian/Hessian sparsity) as numpy arrays and the evaluation
// callbacks as Python callables, drives Ipopt::IpoptApplication over a TNLP that
// trampolines every evaluation back into Python, and returns the solution as a dict.
//
// Threading: everything runs on the calling thread and the GIL is held throughout,
// so the Python callbacks need no locking. A Python exception raised inside a
// callback is captured, the evaluation reports failure to IPOPT, and the exception
// is re-raised to the caller once IPOPT returns.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <IpIpoptApplication.hpp>
#include <IpIpoptData.hpp>
#include <IpSolveStatistics.hpp>
#include <IpTNLP.hpp>

#include <cstring>
#include <exception>
#include <string>
#include <vector>

namespace nb = nanobind;

using DArray = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using IArray = nb::ndarray<int32_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

namespace {

std::vector<double> to_vector(const DArray &a) {
    return std::vector<double>(a.data(), a.data() + a.size());
}

std::vector<int32_t> to_vector_i(const IArray &a) {
    return std::vector<int32_t>(a.data(), a.data() + a.size());
}

// A fresh numpy array holding a copy of `x`, handed to the Python callbacks.
nb::object make_numpy(const double *x, size_t k) {
    double *buf = new double[k];
    std::memcpy(buf, x, k * sizeof(double));
    nb::capsule owner(buf, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    return nb::cast(nb::ndarray<nb::numpy, double, nb::ndim<1>>(buf, {k}, owner));
}

// Copy a callback's numpy result into an IPOPT buffer, checking the length.
void copy_result(const nb::object &res, double *out, size_t expected, const char *what) {
    DArray arr = nb::cast<DArray>(res);
    if (arr.size() != expected)
        throw std::runtime_error(std::string(what) + ": expected " + std::to_string(expected) +
                                 " values, got " + std::to_string(arr.size()));
    std::memcpy(out, arr.data(), expected * sizeof(double));
}

const char *status_name(Ipopt::ApplicationReturnStatus s) {
    using S = Ipopt::ApplicationReturnStatus;
    switch (s) {
        case S::Solve_Succeeded: return "Solve_Succeeded";
        case S::Solved_To_Acceptable_Level: return "Solved_To_Acceptable_Level";
        case S::Infeasible_Problem_Detected: return "Infeasible_Problem_Detected";
        case S::Search_Direction_Becomes_Too_Small: return "Search_Direction_Becomes_Too_Small";
        case S::Diverging_Iterates: return "Diverging_Iterates";
        case S::User_Requested_Stop: return "User_Requested_Stop";
        case S::Feasible_Point_Found: return "Feasible_Point_Found";
        case S::Maximum_Iterations_Exceeded: return "Maximum_Iterations_Exceeded";
        case S::Restoration_Failed: return "Restoration_Failed";
        case S::Error_In_Step_Computation: return "Error_In_Step_Computation";
        case S::Maximum_CpuTime_Exceeded: return "Maximum_CpuTime_Exceeded";
        case S::Maximum_WallTime_Exceeded: return "Maximum_WallTime_Exceeded";
        case S::Not_Enough_Degrees_Of_Freedom: return "Not_Enough_Degrees_Of_Freedom";
        case S::Invalid_Problem_Definition: return "Invalid_Problem_Definition";
        case S::Invalid_Option: return "Invalid_Option";
        case S::Invalid_Number_Detected: return "Invalid_Number_Detected";
        case S::Unrecoverable_Exception: return "Unrecoverable_Exception";
        case S::NonIpopt_Exception_Thrown: return "NonIpopt_Exception_Thrown";
        case S::Insufficient_Memory: return "Insufficient_Memory";
        case S::Internal_Error: return "Internal_Error";
        default: return "Unknown";
    }
}

class CallbackNLP : public Ipopt::TNLP {
public:
    // problem data
    Ipopt::Index n_ = 0, m_ = 0;
    std::vector<double> x_l_, x_u_, g_l_, g_u_, x0_;
    std::vector<int32_t> jac_rows_, jac_cols_, hess_rows_, hess_cols_;
    nb::object eval_f_, eval_grad_, eval_g_, eval_jac_, eval_h_;
    bool have_hess_ = false;

    // captured Python error from a callback
    std::exception_ptr error_;

    // solution
    bool solved_ = false;
    double obj_ = 0.0;
    std::vector<double> sol_x_, sol_g_, sol_mult_g_, sol_z_l_, sol_z_u_;

    bool get_nlp_info(Ipopt::Index &n, Ipopt::Index &m, Ipopt::Index &nnz_jac_g,
                      Ipopt::Index &nnz_h_lag, IndexStyleEnum &index_style) override {
        n = n_;
        m = m_;
        nnz_jac_g = static_cast<Ipopt::Index>(jac_rows_.size());
        nnz_h_lag = static_cast<Ipopt::Index>(hess_rows_.size());
        index_style = C_STYLE;
        return true;
    }

    bool get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l, Ipopt::Number *x_u,
                         Ipopt::Index m, Ipopt::Number *g_l, Ipopt::Number *g_u) override {
        std::memcpy(x_l, x_l_.data(), n * sizeof(double));
        std::memcpy(x_u, x_u_.data(), n * sizeof(double));
        std::memcpy(g_l, g_l_.data(), m * sizeof(double));
        std::memcpy(g_u, g_u_.data(), m * sizeof(double));
        return true;
    }

    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number *x, bool init_z,
                            Ipopt::Number * /*z_L*/, Ipopt::Number * /*z_U*/, Ipopt::Index /*m*/,
                            bool init_lambda, Ipopt::Number * /*lambda*/) override {
        if (init_z || init_lambda)
            return false;
        if (init_x)
            std::memcpy(x, x0_.data(), n * sizeof(double));
        return true;
    }

    bool eval_f(Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Number &obj_value) override {
        try {
            obj_value = nb::cast<double>(eval_f_(make_numpy(x, n)));
            return true;
        } catch (...) {
            error_ = std::current_exception();
            return false;
        }
    }

    bool eval_grad_f(Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Number *grad_f) override {
        try {
            copy_result(eval_grad_(make_numpy(x, n)), grad_f, n, "gradient");
            return true;
        } catch (...) {
            error_ = std::current_exception();
            return false;
        }
    }

    bool eval_g(Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number *g) override {
        try {
            copy_result(eval_g_(make_numpy(x, n)), g, m, "constraints");
            return true;
        } catch (...) {
            error_ = std::current_exception();
            return false;
        }
    }

    bool eval_jac_g(Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Index /*m*/,
                    Ipopt::Index nele_jac, Ipopt::Index *iRow, Ipopt::Index *jCol,
                    Ipopt::Number *values) override {
        if (values == nullptr) {
            for (Ipopt::Index k = 0; k < nele_jac; ++k) {
                iRow[k] = jac_rows_[k];
                jCol[k] = jac_cols_[k];
            }
            return true;
        }
        try {
            copy_result(eval_jac_(make_numpy(x, n)), values, nele_jac, "jacobian");
            return true;
        } catch (...) {
            error_ = std::current_exception();
            return false;
        }
    }

    bool eval_h(Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Number obj_factor,
                Ipopt::Index m, const Ipopt::Number *lambda, bool /*new_lambda*/, Ipopt::Index nele_hess,
                Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values) override {
        if (!have_hess_)
            return false;
        if (values == nullptr) {
            for (Ipopt::Index k = 0; k < nele_hess; ++k) {
                iRow[k] = hess_rows_[k];
                jCol[k] = hess_cols_[k];
            }
            return true;
        }
        try {
            copy_result(eval_h_(make_numpy(x, n), obj_factor, make_numpy(lambda, m)), values, nele_hess,
                        "hessian");
            return true;
        } catch (...) {
            error_ = std::current_exception();
            return false;
        }
    }

    void finalize_solution(Ipopt::SolverReturn /*status*/, Ipopt::Index n, const Ipopt::Number *x,
                           const Ipopt::Number *z_L, const Ipopt::Number *z_U, Ipopt::Index m,
                           const Ipopt::Number *g, const Ipopt::Number *lambda, Ipopt::Number obj_value,
                           const Ipopt::IpoptData * /*ip_data*/,
                           Ipopt::IpoptCalculatedQuantities * /*ip_cq*/) override {
        solved_ = true;
        obj_ = obj_value;
        sol_x_.assign(x, x + n);
        sol_g_.assign(g, g + m);
        sol_mult_g_.assign(lambda, lambda + m);
        sol_z_l_.assign(z_L, z_L + n);
        sol_z_u_.assign(z_U, z_U + n);
    }
};

nb::object vec_to_numpy(const std::vector<double> &v) {
    return make_numpy(v.data(), v.size());
}

nb::dict solve_nlp(int n, int m, DArray x_l, DArray x_u, DArray g_l, DArray g_u, DArray x0,
                   IArray jac_rows, IArray jac_cols, nb::object eval_f, nb::object eval_grad,
                   nb::object eval_g, nb::object eval_jac, nb::object hess_rows, nb::object hess_cols,
                   nb::object eval_h, nb::dict options) {
    if ((int)x_l.size() != n || (int)x_u.size() != n || (int)x0.size() != n)
        throw std::invalid_argument("variable arrays must have length n");
    if ((int)g_l.size() != m || (int)g_u.size() != m)
        throw std::invalid_argument("constraint bound arrays must have length m");
    if (jac_rows.size() != jac_cols.size())
        throw std::invalid_argument("jac_rows and jac_cols must have the same length");

    // IPOPT indexes its matrices with these values without any bounds checking of its
    // own, so an out-of-range entry coming from Python would corrupt memory
    auto check_indices = [](const IArray &idx, int limit, const char *what) {
        for (size_t k = 0; k < idx.size(); ++k)
            if (idx.data()[k] < 0 || idx.data()[k] >= limit)
                throw std::invalid_argument(std::string(what) + "[" + std::to_string(k) + "] = " +
                                            std::to_string(idx.data()[k]) + " is out of range [0, " +
                                            std::to_string(limit) + ")");
    };
    check_indices(jac_rows, m, "jac_rows");
    check_indices(jac_cols, n, "jac_cols");

    Ipopt::SmartPtr<CallbackNLP> nlp = new CallbackNLP();
    nlp->n_ = n;
    nlp->m_ = m;
    nlp->x_l_ = to_vector(x_l);
    nlp->x_u_ = to_vector(x_u);
    nlp->g_l_ = to_vector(g_l);
    nlp->g_u_ = to_vector(g_u);
    nlp->x0_ = to_vector(x0);
    nlp->jac_rows_ = to_vector_i(jac_rows);
    nlp->jac_cols_ = to_vector_i(jac_cols);
    nlp->eval_f_ = eval_f;
    nlp->eval_grad_ = eval_grad;
    nlp->eval_g_ = eval_g;
    nlp->eval_jac_ = eval_jac;

    if (!eval_h.is_none()) {
        IArray hr = nb::cast<IArray>(hess_rows);
        IArray hc = nb::cast<IArray>(hess_cols);
        if (hr.size() != hc.size())
            throw std::invalid_argument("hess_rows and hess_cols must have the same length");
        check_indices(hr, n, "hess_rows");
        check_indices(hc, n, "hess_cols");
        nlp->hess_rows_ = to_vector_i(hr);
        nlp->hess_cols_ = to_vector_i(hc);
        nlp->eval_h_ = eval_h;
        nlp->have_hess_ = true;
    }

    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    app->RethrowNonIpoptException(false);

    for (auto item : options) {
        std::string key = nb::cast<std::string>(item.first);
        nb::handle value = item.second;
        if (nb::isinstance<nb::str>(value)) {
            app->Options()->SetStringValue(key, nb::cast<std::string>(value));
        } else if (nb::isinstance<nb::bool_>(value)) {
            app->Options()->SetStringValue(key, nb::cast<bool>(value) ? "yes" : "no");
        } else if (nb::isinstance<nb::int_>(value)) {
            app->Options()->SetIntegerValue(key, nb::cast<int>(value));
        } else if (nb::isinstance<nb::float_>(value)) {
            app->Options()->SetNumericValue(key, nb::cast<double>(value));
        } else {
            throw std::invalid_argument("option '" + key + "' must be str, bool, int or float");
        }
    }

    // without a Hessian callback IPOPT must run quasi-Newton, or it aborts at startup
    if (!nlp->have_hess_)
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    Ipopt::ApplicationReturnStatus status = app->Initialize();
    if (status == Ipopt::Solve_Succeeded)
        status = app->OptimizeTNLP(Ipopt::SmartPtr<Ipopt::TNLP>(Ipopt::GetRawPtr(nlp)));

    // A failed evaluation is recoverable for IPOPT (it cuts the step and retries), so a
    // captured Python exception only surfaces when the solve actually failed — a solve
    // that converged despite a transient callback error keeps its solution.
    if (nlp->error_ && status != Ipopt::Solve_Succeeded && status != Ipopt::Solved_To_Acceptable_Level)
        std::rethrow_exception(nlp->error_);

    nb::dict out;
    out["status"] = static_cast<int>(status);
    out["status_name"] = status_name(status);
    out["obj"] = nlp->obj_;
    out["x"] = vec_to_numpy(nlp->solved_ ? nlp->sol_x_ : nlp->x0_);
    out["g"] = vec_to_numpy(nlp->sol_g_);
    out["mult_g"] = vec_to_numpy(nlp->sol_mult_g_);
    out["mult_x_l"] = vec_to_numpy(nlp->sol_z_l_);
    out["mult_x_u"] = vec_to_numpy(nlp->sol_z_u_);
    if (Ipopt::IsValid(app->Statistics()))
        out["iterations"] = app->Statistics()->IterationCount();
    return out;
}

}  // namespace

NB_MODULE(_core, mod) {
    mod.doc() = "In-process IPOPT solver for compas_cra";
    mod.def("solve_nlp", &solve_nlp, nb::arg("n"), nb::arg("m"), nb::arg("x_l"), nb::arg("x_u"),
            nb::arg("g_l"), nb::arg("g_u"), nb::arg("x0"), nb::arg("jac_rows"), nb::arg("jac_cols"),
            nb::arg("eval_f"), nb::arg("eval_grad"), nb::arg("eval_g"), nb::arg("eval_jac"),
            nb::arg("hess_rows").none(), nb::arg("hess_cols").none(), nb::arg("eval_h").none(),
            nb::arg("options"),
            "Solve a sparse NLP with IPOPT; see compas_cra.nlp for the friendly wrapper.");
    mod.attr("IPOPT_VERSION") = IPOPT_VERSION;
}
