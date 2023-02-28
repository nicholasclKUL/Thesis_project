#pragma once

#include <alpaqa/config/config.hpp>
#include <Thesis/para-box.hpp>
#include <alpaqa/problem/problem-counters.hpp>
#include <alpaqa/util/alloc-check.hpp>
#include <alpaqa/util/check-dim.hpp>
#include <alpaqa/util/not-implemented.hpp>
#include <alpaqa/util/required-method.hpp>
#include <alpaqa/util/type-erasure.hpp>
#include <chrono>
#include <stdexcept>
#include <type_traits>
#include <algorithm>

namespace alpaqa {

template <Config Conf>
struct ParaProblemVTable : util::BasicVTable {
    USING_ALPAQA_CONFIG(Conf);
    using Box = alpaqa::Box<config_t>;

    template <class F>
    using optional_function_t = util::BasicVTable::optional_function_t<F, ParaProblemVTable>;
    template <class F>
    using optional_const_function_t =
        util::BasicVTable::optional_const_function_t<F, ParaProblemVTable>;

    // clang-format off

    // Required
    required_const_function_t<void(crvec z, rvec p, unsigned int i)>
        eval_proj_diff_g;
    required_const_function_t<void(rvec y, real_t M, index_t penalty_alm_split, unsigned int i)>
        eval_proj_multipliers;
    required_const_function_t<real_t(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p, unsigned int i)>
        eval_prox_grad_step;
    required_const_function_t<real_t(crvec x)>
        eval_f;
    required_const_function_t<void(crvec x, rvec grad_fx, unsigned int i)>
        eval_grad_f;
    required_const_function_t<void(crvec x, rvec gxn, unsigned int i)>
        eval_g;
    required_const_function_t<void(crvec x, crvec y, rvec grad_gxy, unsigned int i)>
        eval_grad_g_prod;

    // Second order
    optional_const_function_t<void(crvec x, index_t i, rvec grad_gi)>
        eval_grad_gi = &default_eval_grad_gi;
    optional_const_function_t<void(crvec x, crvec y, crvec v, rvec Hv)>
        eval_hess_L_prod = &default_eval_hess_L_prod;
    optional_const_function_t<void(crvec x, crvec y, rmat H)>
        eval_hess_L = &default_eval_hess_L;
    optional_const_function_t<void(crvec x, crvec y, crvec Σ, crvec v, rvec Hv)>
        eval_hess_ψ_prod = &default_eval_hess_ψ_prod;
    optional_const_function_t<void(crvec x, crvec y, crvec Σ, rmat H)>
        eval_hess_ψ = &default_eval_hess_ψ;

    // Combined evaluations
    optional_const_function_t<real_t(crvec x, rvec grad_fx)>
        eval_f_grad_f = &default_eval_f_grad_f;
    optional_const_function_t<real_t(crvec x, rvec g)>
        eval_f_g = &default_eval_f_g;
    optional_const_function_t<real_t(crvec x, rvec grad_fx, rvec g)>
        eval_f_grad_f_g = &default_eval_f_grad_f_g;
    optional_const_function_t<void(crvec x, crvec y, rvec grad_f, rvec grad_gxy)>
        eval_grad_f_grad_g_prod = &default_eval_grad_f_grad_g_prod;

    // Lagrangian and augmented lagrangian evaluations
    optional_const_function_t<void(crvec x, crvec y, rvec grad_L, rvec work_n)>
        eval_grad_L = &default_eval_grad_L;
    optional_const_function_t<real_t(crvec x, crvec y, crvec Σ, rvec ŷ)>
        eval_ψ = &default_eval_ψ;
    optional_const_function_t<void(crvec x, crvec ŷ, rvec grad_ψ, rvec work_n)>
        eval_grad_ψ_from_ŷ = &default_eval_grad_ψ_from_ŷ;
    optional_const_function_t<void(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m)>
        eval_grad_ψ = &default_eval_grad_ψ;
    optional_const_function_t<real_t(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m)>
        eval_ψ_grad_ψ = &default_eval_ψ_grad_ψ;

    // Constraint sets
    optional_const_function_t<const Box &()>
        get_box_C = &default_get_box_C;
    optional_const_function_t<const Box &()>
        get_box_D = &default_get_box_D;

    // Check
    required_const_function_t<void()>
        check;

    // clang-format on

    length_t n, m;

    static real_t calc_ŷ_dᵀŷ(const void *self, rvec g_ŷ, crvec y, crvec Σ,
                             const ParaProblemVTable &vtable, unsigned int i) {
        if (Σ.size() == 1) {
            // ζ = g(x) + Σ⁻¹y
            g_ŷ[i] = (1 / Σ(0)) * y[i];
            // d = ζ - Π(ζ, D)
            vtable.eval_proj_diff_g(self, g_ŷ[i], g_ŷ[i]);
            // dᵀŷ, ŷ = Σ d
            real_t dᵀŷ = Σ(0) * g_ŷ[i].dot(g_ŷ[i]);
            g_ŷ[i] *= Σ(0);
            return dᵀŷ;
        } else {
            // ζ = g(x) + Σ⁻¹y
            g_ŷ[i] = Σ[i] * y[i];
            // d = ζ - Π(ζ, D)
            vtable.eval_proj_diff_g(self, g_ŷ[i], g_ŷ[i]);
            // dᵀŷ, ŷ = Σ d
            real_t dᵀŷ = 0;
            dᵀŷ = g_ŷ[i] * Σ[i] * g_ŷ[i];
            g_ŷ[i] = Σ[i] * g_ŷ[i];
            return dᵀŷ;
        }
    }

    static void default_eval_grad_gi(const void *, crvec, index_t, rvec, const ParaProblemVTable &) {
        throw not_implemented_error("eval_grad_gi");
    }
    static void default_eval_hess_L_prod(const void *, crvec, crvec, crvec, rvec,
                                         const ParaProblemVTable &) {
        throw not_implemented_error("eval_hess_L_prod");
    }
    static void default_eval_hess_L(const void *, crvec, crvec, rmat, const ParaProblemVTable &) {
        throw not_implemented_error("eval_hess_L");
    }
    static void default_eval_hess_ψ_prod(const void *self, crvec x, crvec y, crvec, crvec v,
                                         rvec Hv, const ParaProblemVTable &vtable) {
        if (y.size() == 0 && vtable.eval_hess_L_prod != default_eval_hess_L_prod)
            return vtable.eval_hess_L_prod(self, x, y, v, Hv, vtable);
        throw not_implemented_error("eval_hess_ψ_prod");
    }
    static void default_eval_hess_ψ(const void *self, crvec x, crvec y, crvec, rmat H,
                                    const ParaProblemVTable &vtable) {
        if (y.size() == 0 && vtable.eval_hess_L != default_eval_hess_L)
            return vtable.eval_hess_L(self, x, y, H, vtable);
        throw not_implemented_error("eval_hess_ψ");
    }
    //Need review
    static real_t default_eval_f_grad_f(const void *self, crvec x, rvec grad_fx,
                                        const ParaProblemVTable &vtable, unsigned int i) {
        vtable.eval_grad_f(self, x, grad_fx[i]);
        return vtable.eval_f(self, x);
    }
    //Need review
    static real_t default_eval_f_g(const void *self, crvec x, rvec g, const ParaProblemVTable &vtable,
                                   unsigned int i) {
        vtable.eval_g(self, x, g[i]);
        return vtable.eval_f(self, x);
    }
    //Need review
    static real_t default_eval_f_grad_f_g(const void *self, crvec x, rvec grad_fx, rvec g,
                                          const ParaProblemVTable &vtable, unsigned int i) {
        vtable.eval_g(self, x, g[i]);
        return vtable.eval_f_grad_f(self, x, grad_fx[i], vtable, i);
    }
    //Need review
    static void default_eval_grad_f_grad_g_prod(const void *self, crvec x, crvec y, rvec grad_f,
                                                rvec grad_gxy, const ParaProblemVTable &vtable, 
                                                unsigned int i) {
        vtable.eval_grad_f(self, x, grad_f[i]);
        vtable.eval_grad_g_prod(self, x, y, grad_gxy[i]);
    }
    static void default_eval_grad_L(const void *self, crvec x, crvec y, rvec grad_L, rvec work_n,
                                    const ParaProblemVTable &vtable, unsigned int i) {
        vtable.eval_grad_f_grad_g_prod(self, x, y, grad_L, work_n, vtable, i);
        grad_L[i] = work_n;
    }
    static real_t default_eval_ψ(const void *self, crvec x, crvec y, crvec Σ, rvec ŷ,
                                 const ParaProblemVTable &vtable, unsigned int i) {
        if (y.size() == 0) /* [[unlikely]] */
            return vtable.eval_f(self, x);

        real_t f   = vtable.eval_f_g(self, x, ŷ, vtable, i);
        real_t dᵀŷ = calc_ŷ_dᵀŷ(self, ŷ, y, Σ, vtable, i);
        // ψ(x) = f(x) + ½ dᵀŷ
        real_t ψ = f + real_t(0.5) * dᵀŷ;
        return ψ;
    }
    static void default_eval_grad_ψ_from_ŷ(const void *self, crvec x, crvec ŷ, rvec grad_ψ,
                                           rvec work_n, const ParaProblemVTable &vtable, unsigned int i) {
        if (ŷ.size() == 0) /* [[unlikely]] */
            vtable.eval_grad_f(self, x, grad_ψ[i]);
        else
            vtable.eval_grad_L(self, x, ŷ, grad_ψ, work_n, vtable, i);
    }
    static void default_eval_grad_ψ(const void *self, crvec x, crvec y, crvec Σ, rvec grad_ψ,
                                    rvec work_n, rvec work_m, const ParaProblemVTable &vtable, unsigned int i) {
        if (y.size() == 0) /* [[unlikely]] */ {
            vtable.eval_grad_f(self, x, grad_ψ[i]);
        } else {
            vtable.eval_g(self, x, work_m[i]);
            (void)calc_ŷ_dᵀŷ(self, work_m, y, Σ, vtable, i);
            vtable.eval_grad_ψ_from_ŷ(self, x, work_m, grad_ψ, work_n, vtable, i);
        }
    }
    static real_t default_eval_ψ_grad_ψ(const void *self, crvec x, crvec y, crvec Σ, rvec grad_ψ,
                                        rvec work_n, rvec work_m, const ParaProblemVTable &vtable, unsigned int i) {
        if (y.size() == 0) /* [[unlikely]] */
            return vtable.eval_f_grad_f(self, x, grad_ψ, vtable, i);

        auto &ŷ = work_m;
        // ψ(x) = f(x) + ½ dᵀŷ
        real_t f   = vtable.eval_f_g(self, x, ŷ, vtable, i);
        real_t dᵀŷ = calc_ŷ_dᵀŷ(self, ŷ, y, Σ, vtable, i);
        real_t ψ   = f + real_t(0.5) * dᵀŷ;
        // ∇ψ(x) = ∇f(x) + ∇g(x) ŷ
        vtable.eval_grad_L(self, x, ŷ, grad_ψ, work_n, vtable, i);
        return ψ;
    }
    static const Box &default_get_box_C(const void *, const ParaProblemVTable &) {
        throw not_implemented_error("get_box_C");
    }
    static const Box &default_get_box_D(const void *, const ParaProblemVTable &) {
        throw not_implemented_error("get_box_D");
    }

    template <class P>
    ParaProblemVTable(util::VTableTypeTag<P> t) : util::BasicVTable{t} {
        auto &vtable = *this;
        assert(t.t);

        // Initialize all methods

        // Required
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_proj_diff_g);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_proj_multipliers);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_prox_grad_step);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_f);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_grad_f);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_g);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_grad_g_prod);
        // Second order
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_gi, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_hess_L_prod, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_hess_L, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_hess_ψ_prod, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_hess_ψ, t.t);
        // Combined evaluations
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_f_grad_f, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_f_g, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_f_grad_f_g, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_f_grad_g_prod, t.t);
        // Lagrangian and augmented lagrangian evaluations
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_L, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_ψ, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_ψ_from_ŷ, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_ψ, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_ψ_grad_ψ, t.t);
        // Constraint set
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, get_box_C, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, get_box_D, t.t);
        // Check
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, check);

        // Dimensions
        vtable.n = t.t->get_n();
        vtable.m = t.t->get_m();
    }
    ParaProblemVTable() = default;
};

/// @addtogroup grp_Problems
/// @{

template <Config Conf = DefaultConfig, class Allocator = std::allocator<std::byte>>
class ParaTypeErasedProblem : public util::TypeErased<ParaProblemVTable<Conf>, Allocator> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using Box            = alpaqa::Box<config_t>;
    using VTable         = ParaProblemVTable<config_t>;
    using allocator_type = Allocator;
    using TypeErased     = util::TypeErased<VTable, allocator_type>;
    using TypeErased::TypeErased;

  protected:
    using TypeErased::call;
    using TypeErased::self;
    using TypeErased::vtable;

  public:
    template <class T, class... Args>
    static ParaTypeErasedProblem make(Args &&...args) {
        return TypeErased::template make<ParaTypeErasedProblem, T>(std::forward<Args>(args)...);
    }

    /// @name Problem dimensions
    /// @{

    /// **[Required]**
    /// Number of decision variables.
    length_t get_n() const;
    /// **[Required]**
    /// Number of constraints.
    length_t get_m() const;

    /// @}

    /// @name Required cost and constraint functions
    /// @{

    /// **[Required]**
    /// Function that evaluates the cost, @f$ f(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    real_t eval_f(crvec x, unsigned int i) const;
    /// **[Required]**
    /// Function that evaluates the gradient of the cost, @f$ \nabla f(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [out] grad_fx
    ///         Gradient of cost function @f$ \nabla f(x) \in \R^n @f$
    void eval_grad_f(crvec x, rvec grad_fx, unsigned int i) const;
    /// **[Required]**
    /// Function that evaluates the constraints, @f$ g(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [out] gx
    ///         Value of the constraints @f$ g(x) \in \R^m @f$
    void eval_g(crvec x, rvec gx, unsigned int i) const;
    /// **[Required]**
    /// Function that evaluates the gradient of the constraints times a vector,
    /// @f$ \nabla g(x)\,y = \tp{\jac_g(x)}y @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] y
    ///         Vector @f$ y \in \R^m @f$ to multiply the gradient by
    /// @param  [out] grad_gxy
    ///         Gradient of the constraints
    ///         @f$ \nabla g(x)\,y \in \R^n @f$
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy, unsigned int i) const;

    /// @}

    /// @name Projections onto constraint sets and proximal mappings
    /// @{

    /// **[Required]**
    /// Function that evaluates the difference between the given point @f$ z @f$
    /// and its projection onto the constraint set @f$ D @f$.
    /// @param  [in] z
    ///         Slack variable, @f$ z \in \R^m @f$
    /// @param  [out] p
    ///         The difference relative to its projection,
    ///         @f$ p = z - \Pi_D(z) \in \R^m @f$
    /// @note   @p z and @p p can refer to the same vector.
    void eval_proj_diff_g(crvec z, rvec p, unsigned int i) const;
    /// **[Required]**
    /// Function that projects the Lagrange multipliers for ALM.
    /// @param  [inout] y
    ///         Multipliers, @f$ y \leftarrow \Pi_Y(y) \in \R^m @f$
    /// @param  [in] M
    ///         The radius/size of the set @f$ Y @f$. See @ref ALMParams::M.
    /// @param  [in] penalty_alm_split
    ///         See @ref ALMParams::penalty_alm_split.
    void eval_proj_multipliers(rvec y, real_t M, index_t penalty_alm_split) const;
    /// **[Required]**
    /// Function that evaluates a proximal gradient step.
    /// @param  [in] γ
    ///         Step size, @f$ \gamma \in \R_{>0} @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] grad_ψ
    ///         Gradient of the subproblem cost, @f$ \nabla\psi(x) \in \R^n @f$
    /// @param  [out] x̂
    ///         Next proximal gradient iterate, @f$ \hat x = T_\gamma(x) =
    ///         \prox_{\gamma h}(x - \gamma\nabla\psi(x)) \in \R^n @f$
    /// @param  [out] p
    ///         The proximal gradient step,
    ///         @f$ p = \hat x - x \in \R^n @f$
    /// @return The nonsmooth function evaluated at x̂,
    ///         @f$ h(\hat x) @f$.
    /// @note   The vector @f$ p @f$ is often used in stopping criteria, so its
    ///         numerical accuracy is more important than that of @f$ \hat x @f$.
    real_t eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p, unsigned int i) const;

    /// @}

    /// @name Constraint sets
    /// @{

    /// **[Optional]**
    /// Get the rectangular constraint set of the decision variables,
    /// @f$ x \in C @f$.
    const Box &get_box_C() const;
    /// **[Optional]**
    /// Get the rectangular constraint set of the general constraint function,
    /// @f$ g(x) \in D @f$.
    const Box &get_box_D() const;

    /// @}

    /// @name Functions for second-order solvers
    /// @{

    /// **[Optional]**
    /// Function that evaluates the gradient of one specific constraint,
    /// @f$ \nabla g_i(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] i
    ///         Which constraint @f$ 0 \le i \lt m @f$
    /// @param  [out] grad_gi
    ///         Gradient of the constraint
    ///         @f$ \nabla g_i(x) \in \R^n @f$
    ///
    /// Required for second-order solvers only.
    void eval_grad_gi(crvec x, index_t i, rvec grad_gi) const;
    /// **[Optional]**
    /// Function that evaluates the Hessian of the Lagrangian multiplied by a
    /// vector,
    /// @f$ \nabla_{xx}^2L(x, y)\,v @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] y
    ///         Lagrange multipliers @f$ y \in \R^m @f$
    /// @param  [in] v
    ///         Vector to multiply by @f$ v \in \R^n @f$
    /// @param  [out] Hv
    ///         Hessian-vector product
    ///         @f$ \nabla_{xx}^2 L(x, y)\,v \in \R^{n} @f$
    ///
    /// Required for second-order solvers only.
    void eval_hess_L_prod(crvec x, crvec y, crvec v, rvec Hv) const;
    /// **[Optional]**
    /// Function that evaluates the Hessian of the Lagrangian,
    /// @f$ \nabla_{xx}^2L(x, y) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] y
    ///         Lagrange multipliers @f$ y \in \R^m @f$
    /// @param  [out] H
    ///         Hessian @f$ \nabla_{xx}^2 L(x, y) \in \R^{n\times n} @f$
    ///
    /// Required for second-order solvers only.
    void eval_hess_L(crvec x, crvec y, rmat H) const;
    /// **[Optional]**
    /// Function that evaluates the Hessian of the augmented Lagrangian
    /// multiplied by a vector,
    /// @f$ \nabla_{xx}^2L_\Sigma(x, y)\,v @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] y
    ///         Lagrange multipliers @f$ y \in \R^m @f$
    /// @param  [in] Σ
    ///         Penalty weights @f$ \Sigma @f$
    /// @param  [in] v
    ///         Vector to multiply by @f$ v \in \R^n @f$
    /// @param  [out] Hv
    ///         Hessian-vector product
    ///         @f$ \nabla_{xx}^2 L_\Sigma(x, y)\,v \in \R^{n} @f$
    ///
    /// Required for second-order solvers only.
    void eval_hess_ψ_prod(crvec x, crvec y, crvec Σ, crvec v, rvec Hv) const;
    /// **[Optional]**
    /// Function that evaluates the Hessian of the augmented Lagrangian,
    /// @f$ \nabla_{xx}^2L_\Sigma(x, y) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] y
    ///         Lagrange multipliers @f$ y \in \R^m @f$
    /// @param  [in] Σ
    ///         Penalty weights @f$ \Sigma @f$
    /// @param  [out] H
    ///         Hessian @f$ \nabla_{xx}^2 L_\Sigma(x, y) \in \R^{n\times n} @f$
    ///
    /// Required for second-order solvers only.
    void eval_hess_ψ(crvec x, crvec y, crvec Σ, rmat H) const;

    /// @}

    /// @name Combined evaluations
    /// @{

    /// **[Optional]**
    /// Evaluate both @f$ f(x) @f$ and its gradient, @f$ \nabla f(x) @f$.
    real_t eval_f_grad_f(crvec x, rvec grad_fx, unsigned int i) const;
    /// **[Optional]**
    /// Evaluate both @f$ f(x) @f$ and @f$ g(x) @f$.
    real_t eval_f_g(crvec x, rvec g, unsigned int i) const;
    /// **[Optional]**
    /// Evaluate @f$ f(x) @f$, its gradient @f$ \nabla f(x) @f$ and @f$ g(x) @f$.
    real_t eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g, unsigned int i) const;
    /// **[Optional]**
    /// Evaluate both @f$ \nabla f(x) @f$ and @f$ \nabla g(x)\,y @f$.
    void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f, rvec grad_gxy, unsigned int i) const;
    /// **[Optional]**
    /// Evaluate the gradient of the Lagrangian
    /// @f$ \nabla_x L(x, y) = \nabla f(x) + \nabla g(x)\,y @f$
    void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n, unsigned int i) const;

    /// @}

    /// @name Augmented Lagrangian
    /// @{

    /// **[Optional]**
    /// Calculate both ψ(x) and the vector ŷ that can later be used to compute
    /// ∇ψ.
    /// @f[ \psi(x) = f(x) + \tfrac{1}{2}
    ///   \text{dist}_\Sigma^2\left(g(x) + \Sigma^{-1}y,\;D\right) @f]
    /// @f[ \hat y = \Sigma\, \left(g(x) + \Sigma^{-1}y - \Pi_D\left(g(x)
    ///   + \Sigma^{-1}y\right)\right) @f]
    real_t eval_ψ(crvec x, ///< [in]  Decision variable @f$ x @f$
                  crvec y, ///< [in]  Lagrange multipliers @f$ y @f$
                  crvec Σ, ///< [in]  Penalty weights @f$ \Sigma @f$
                  rvec ŷ,  ///< [out] @f$ \hat y @f$
                  unsigned int i
    ) const;
    /// **[Optional]**
    /// Calculate ∇ψ(x) using ŷ.
    void eval_grad_ψ_from_ŷ(crvec x,     ///< [in]  Decision variable @f$ x @f$
                            crvec ŷ,     ///< [in]  @f$ \hat y @f$
                            rvec grad_ψ, ///< [out] @f$ \nabla \psi(x) @f$
                            rvec work_n, ///<       Dimension @f$ n @f$
                            unsigned int i
    ) const;
    /// **[Optional]**
    /// Calculate the gradient ∇ψ(x).
    /// @f[ \nabla \psi(x) = \nabla f(x) + \nabla g(x)\,\hat y(x) @f]
    void eval_grad_ψ(crvec x,     ///< [in]  Decision variable @f$ x @f$
                     crvec y,     ///< [in]  Lagrange multipliers @f$ y @f$
                     crvec Σ,     ///< [in]  Penalty weights @f$ \Sigma @f$
                     rvec grad_ψ, ///< [out] @f$ \nabla \psi(x) @f$
                     rvec work_n, ///<       Dimension @f$ n @f$
                     rvec work_m, ///<       Dimension @f$ m @f$
                     unsigned int i
    ) const;
    /// **[Optional]**
    /// Calculate both ψ(x) and its gradient ∇ψ(x).
    /// @f[ \psi(x) = f(x) + \tfrac{1}{2}
    /// \text{dist}_\Sigma^2\left(g(x) + \Sigma^{-1}y,\;D\right) @f]
    /// @f[ \nabla \psi(x) = \nabla f(x) + \nabla g(x)\,\hat y(x) @f]
    real_t eval_ψ_grad_ψ(crvec x,     ///< [in]  Decision variable @f$ x @f$
                         crvec y,     ///< [in]  Lagrange multipliers @f$ y @f$
                         crvec Σ,     ///< [in]  Penalty weights @f$ \Sigma @f$
                         rvec grad_ψ, ///< [out] @f$ \nabla \psi(x) @f$
                         rvec work_n, ///<       Dimension @f$ n @f$
                         rvec work_m, ///<       Dimension @f$ m @f$
                         unsigned int i
    ) const;

    /// @}

    /// @name Checks
    /// @{

    /// **[Required]**
    /// Check that the problem formulation is well-defined, the dimensions match,
    /// etc. Throws an exception if this is not the case.
    void check() const;

    /// @}

    /// @name Querying specialized implementations
    /// @{

    /// Returns true if the problem provides an implementation of
    /// @ref eval_grad_gi.
    bool provides_eval_grad_gi() const {
        return vtable.eval_grad_gi != vtable.default_eval_grad_gi;
    }
    /// Returns true if the problem provides an implementation of
    /// @ref eval_hess_L_prod.
    bool provides_eval_hess_L_prod() const {
        return vtable.eval_hess_L_prod != vtable.default_eval_hess_L_prod;
    }
    /// Returns true if the problem provides an implementation of
    /// @ref eval_hess_L.
    bool provides_eval_hess_L() const { return vtable.eval_hess_L != vtable.default_eval_hess_L; }
    /// Returns true if the problem provides an implementation of
    /// @ref eval_hess_ψ_prod.
    bool provides_eval_hess_ψ_prod() const {
        return vtable.eval_hess_ψ_prod != vtable.default_eval_hess_ψ_prod;
    }
    /// Returns true if the problem provides an implementation of
    /// @ref eval_hess_ψ.
    bool provides_eval_hess_ψ() const { return vtable.eval_hess_ψ != vtable.default_eval_hess_ψ; }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_f_grad_f, false if it uses the default implementation.
    bool provides_eval_f_grad_f() const {
        return vtable.eval_f_grad_f != vtable.default_eval_f_grad_f;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_f_g, false if it uses the default implementation.
    bool provides_eval_f_g() const { return vtable.eval_f_g != vtable.default_eval_f_g; }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_f_grad_f_g, false if it uses the default implementation.
    bool provides_eval_f_grad_f_g() const {
        return vtable.eval_f_grad_f_g != vtable.default_eval_f_grad_f_g;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_grad_f_grad_g_prod, false if it uses the default implementation.
    bool provides_eval_grad_f_grad_g_prod() const {
        return vtable.eval_grad_f_grad_g_prod != vtable.default_eval_grad_f_grad_g_prod;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_grad_L, false if it uses the default implementation.
    bool provides_eval_grad_L() const { return vtable.eval_grad_L != vtable.default_eval_grad_L; }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_ψ, false if it uses the default implementation.
    bool provides_eval_ψ() const { return vtable.eval_ψ != vtable.default_eval_ψ; }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_grad_ψ_from_ŷ, false if it uses the default implementation.
    bool provides_eval_grad_ψ_from_ŷ() const {
        return vtable.eval_grad_ψ_from_ŷ != vtable.default_eval_grad_ψ_from_ŷ;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_grad_ψ, false if it uses the default implementation.
    bool provides_eval_grad_ψ() const { return vtable.eval_grad_ψ != vtable.default_eval_grad_ψ; }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_ψ_grad_ψ, false if it uses the default implementation.
    bool provides_eval_ψ_grad_ψ() const {
        return vtable.eval_ψ_grad_ψ != vtable.default_eval_ψ_grad_ψ;
    }
    /// Returns true if the problem provides an implementation of
    /// @ref get_box_C.
    bool provides_get_box_C() const { return vtable.get_box_C != vtable.default_get_box_C; }
    /// Returns true if the problem provides an implementation of
    /// @ref get_box_D.
    bool provides_get_box_D() const { return vtable.get_box_D != vtable.default_get_box_D; }

    /// @}

    /// @name Helpers
    /// @{

    /// Given g(x), compute the intermediate results ŷ and dᵀŷ that can later be
    /// used to compute ψ(x) and ∇ψ(x).
    ///
    /// Computes the result using the following algorithm:
    /// @f[ \begin{aligned}
    ///     \zeta &= g(x) + \Sigma^{-1} y \\[]
    ///     d &= \zeta - \Pi_D(\zeta)
    ///        = \operatorname{eval\_proj\_diff\_g}(\zeta, \zeta) \\[]
    ///     \hat y &= \Sigma d \\[]
    /// \end{aligned} @f]
    /// @see @ref page_math
    ///
    /// @param[inout]   g_ŷ
    ///                 Input @f$ g(x) @f$, outputs @f$ \hat y @f$
    /// @param[in]      y
    ///                 Lagrange multipliers @f$ y @f$
    /// @param[in]      Σ
    ///                 Penalty weights @f$ \Sigma @f$
    /// @return The inner product @f$ d^\top \hat y @f$
    real_t calc_ŷ_dᵀŷ(rvec g_ŷ, crvec y, crvec Σ, unsigned int i) const;

    /// @}
};

/// @}

#ifndef DOXYGEN
template <class Tref>
explicit ParaTypeErasedProblem(Tref &&d)
    -> ParaTypeErasedProblem<typename std::remove_cvref_t<Tref>::config_t>;

template <class Tref, class Allocator>
explicit ParaTypeErasedProblem(Tref &&d, Allocator alloc)
    -> ParaTypeErasedProblem<typename std::remove_cvref_t<Tref>::config_t, Allocator>;
#endif

template <Config Conf, class Allocator>
auto ParaTypeErasedProblem<Conf, Allocator>::get_n() const -> length_t {
    return vtable.n;
}
template <Config Conf, class Allocator>
auto ParaTypeErasedProblem<Conf, Allocator>::get_m() const -> length_t {
    return vtable.m;
}

template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::eval_proj_diff_g(crvec z, rvec p, unsigned int i) const {
    return call(vtable.eval_proj_diff_g, z, p, i);
}
template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::eval_proj_multipliers(rvec y, real_t M,
                                                               index_t penalty_alm_split) const {
    return call(vtable.eval_proj_multipliers, y, M, penalty_alm_split);
}
template <Config Conf, class Allocator>
auto ParaTypeErasedProblem<Conf, Allocator>::eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ,
                                                             rvec x̂, rvec p, unsigned int i) const -> real_t {
    return call(vtable.eval_prox_grad_step, γ, x, grad_ψ, x̂, p, i);
}
template <Config Conf, class Allocator>
auto ParaTypeErasedProblem<Conf, Allocator>::eval_f(crvec x, unsigned int i) const -> real_t {
    return call(vtable.eval_f, x, i);
}
template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::eval_grad_f(crvec x, rvec grad_fx, unsigned int i) const {
    return call(vtable.eval_grad_f, x, grad_fx, i);
}
template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::eval_g(crvec x, rvec gx, unsigned int i) const {
    return call(vtable.eval_g, x, gx, i);
}
template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy, unsigned int i) const {
    return call(vtable.eval_grad_g_prod, x, y, grad_gxy, i);
}
template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::eval_grad_gi(crvec x, index_t i, rvec grad_gi) const {
    return call(vtable.eval_grad_gi, x, i, grad_gi);
}
template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::eval_hess_L_prod(crvec x, crvec y, crvec v,
                                                          rvec Hv) const {
    return call(vtable.eval_hess_L_prod, x, y, v, Hv);
}
template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::eval_hess_L(crvec x, crvec y, rmat H) const {
    return call(vtable.eval_hess_L, x, y, H);
}
template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::eval_hess_ψ_prod(crvec x, crvec y, crvec Σ, crvec v,
                                                          rvec Hv) const {
    return call(vtable.eval_hess_ψ_prod, x, y, Σ, v, Hv);
}
template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::eval_hess_ψ(crvec x, crvec y, crvec Σ, rmat H) const {
    return call(vtable.eval_hess_ψ, x, y, Σ, H);
}
template <Config Conf, class Allocator>
auto ParaTypeErasedProblem<Conf, Allocator>::eval_f_grad_f(crvec x, rvec grad_fx, unsigned int i) const -> real_t {
    return call(vtable.eval_f_grad_f, x, grad_fx, i);
}
template <Config Conf, class Allocator>
auto ParaTypeErasedProblem<Conf, Allocator>::eval_f_g(crvec x, rvec g, unsigned int i) const -> real_t {
    return call(vtable.eval_f_g, x, g, i);
}
template <Config Conf, class Allocator>
auto ParaTypeErasedProblem<Conf, Allocator>::eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g, unsigned int i) const
    -> real_t {
    return call(vtable.eval_f_grad_f_g, x, grad_fx, g, i);
}
template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f,
                                                                 rvec grad_gxy, unsigned int i) const {
    return call(vtable.eval_grad_f_grad_g_prod, x, y, grad_f, grad_gxy, i);
}
template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::eval_grad_L(crvec x, crvec y, rvec grad_L,
                                                     rvec work_n, unsigned int i) const {
    return call(vtable.eval_grad_L, x, y, grad_L, work_n, i);
}
template <Config Conf, class Allocator>
auto ParaTypeErasedProblem<Conf, Allocator>::eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ,
                                                unsigned int i) const -> real_t {
    return call(vtable.eval_ψ, x, y, Σ, ŷ, i);
}
template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::eval_grad_ψ_from_ŷ(crvec x, crvec ŷ, rvec grad_ψ,
                                                            rvec work_n, unsigned int i) const {
    return call(vtable.eval_grad_ψ_from_ŷ, x, ŷ, grad_ψ, work_n, i);
}
template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ,
                                                     rvec work_n, rvec work_m, unsigned int i) const {
    return call(vtable.eval_grad_ψ, x, y, Σ, grad_ψ, work_n, work_m, i);
}
template <Config Conf, class Allocator>
auto ParaTypeErasedProblem<Conf, Allocator>::eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ,
                                                       rvec work_n, rvec work_m, unsigned int i) const -> real_t {
    return call(vtable.eval_ψ_grad_ψ, x, y, Σ, grad_ψ, work_n, work_m, i);
}
template <Config Conf, class Allocator>
auto ParaTypeErasedProblem<Conf, Allocator>::calc_ŷ_dᵀŷ(rvec g_ŷ, crvec y, crvec Σ, unsigned int i) const -> real_t {
    return call(vtable.calc_ŷ_dᵀŷ, g_ŷ, y, Σ, i);
}
template <Config Conf, class Allocator>
auto ParaTypeErasedProblem<Conf, Allocator>::get_box_C() const -> const Box & {
    return call(vtable.get_box_C);
}
template <Config Conf, class Allocator>
auto ParaTypeErasedProblem<Conf, Allocator>::get_box_D() const -> const Box & {
    return call(vtable.get_box_D);
}
template <Config Conf, class Allocator>
void ParaTypeErasedProblem<Conf, Allocator>::check() const {
    return call(vtable.check);
}

/// @addtogroup grp_Problems
/// @{

/// Problem wrapper that keeps track of the number of evaluations and the run
/// time of each function.
/// You probably want to use @ref problem_with_counters or
/// @ref problem_with_counters_ref instead of instantiating this class directly.
/// @note   The evaluation counters are stored using a `std::shared_pointers`,
///         which means that different copies of a @ref ProblemWithCounters
///         instance all share the same counters. To opt out of this behavior,
///         you can use the @ref decouple_evaluations function.
template <class Problem>
struct ParaProblemWithCounters {
    USING_ALPAQA_CONFIG_TEMPLATE(std::remove_cvref_t<Problem>::config_t);
    using Box = typename ParaTypeErasedProblem<config_t>::Box;

    // clang-format off
    void eval_proj_diff_g(crvec z, rvec p, unsigned int i) const { ++evaluations->proj_diff_g; return timed(evaluations->time.proj_diff_g, std::bind(&std::remove_cvref_t<Problem>::eval_proj_diff_g, &problem, z, p)); }
    void eval_proj_multipliers(rvec y, real_t M, index_t penalty_alm_split) const { ++evaluations->proj_multipliers; return timed(evaluations->time.proj_multipliers, std::bind(&std::remove_cvref_t<Problem>::eval_proj_multipliers, &problem, y, M, penalty_alm_split)); }
    real_t eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p, unsigned int i) const { ++evaluations->prox_grad_step; return timed(evaluations->time.prox_grad_step, std::bind(&std::remove_cvref_t<Problem>::eval_prox_grad_step, &problem, γ, x, grad_ψ, x̂, p)); }
    real_t eval_f(crvec x, unsigned int i) const { ++evaluations->f; return timed(evaluations->time.f, std::bind(&std::remove_cvref_t<Problem>::eval_f, &problem, x)); }
    void eval_grad_f(crvec x, rvec grad_fx, unsigned int i) const { ++evaluations->grad_f; return timed(evaluations->time.grad_f, std::bind(&std::remove_cvref_t<Problem>::eval_grad_f, &problem, x, grad_fx)); }
    void eval_g(crvec x, rvec gx, unsigned int i) const { ++evaluations->g; return timed(evaluations->time.g, std::bind(&std::remove_cvref_t<Problem>::eval_g, &problem, x, gx)); }
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy, unsigned int i) const { ++evaluations->grad_g_prod; return timed(evaluations->time.grad_g_prod, std::bind(&std::remove_cvref_t<Problem>::eval_grad_g_prod, &problem, x, y, grad_gxy)); }
    void eval_grad_gi(crvec x, index_t i, rvec grad_gi, unsigned int i) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_gi; } { ++evaluations->grad_gi; return timed(evaluations->time.grad_gi, std::bind(&std::remove_cvref_t<Problem>::eval_grad_gi, &problem, x, i, grad_gi)); }
    void eval_hess_L_prod(crvec x, crvec y, crvec v, rvec Hv) const requires requires { &std::remove_cvref_t<Problem>::eval_hess_L_prod; } { ++evaluations->hess_L_prod; return timed(evaluations->time.hess_L_prod, std::bind(&std::remove_cvref_t<Problem>::eval_hess_L_prod, &problem, x, y, v, Hv)); }
    void eval_hess_L(crvec x, crvec y, rmat H) const requires requires { &std::remove_cvref_t<Problem>::eval_hess_L; } { ++evaluations->hess_L; return timed(evaluations->time.hess_L, std::bind(&std::remove_cvref_t<Problem>::eval_hess_L, &problem, x, y, H)); }
    void eval_hess_ψ_prod(crvec x, crvec y, crvec Σ, crvec v, rvec Hv) const requires requires { &std::remove_cvref_t<Problem>::eval_hess_ψ_prod; } { ++evaluations->hess_ψ_prod; return timed(evaluations->time.hess_ψ_prod, std::bind(&std::remove_cvref_t<Problem>::eval_hess_ψ_prod, &problem, x, y, Σ, v, Hv)); }
    void eval_hess_ψ(crvec x, crvec y, crvec Σ, rmat H) const requires requires { &std::remove_cvref_t<Problem>::eval_hess_ψ; } { ++evaluations->hess_ψ; return timed(evaluations->time.hess_ψ, std::bind(&std::remove_cvref_t<Problem>::eval_hess_ψ, &problem, x, y, Σ, H)); }
    real_t eval_f_grad_f(crvec x, rvec grad_fx, unsigned int i) const requires requires { &std::remove_cvref_t<Problem>::eval_f_grad_f; } { ++evaluations->f_grad_f; return timed(evaluations->time.f_grad_f, std::bind(&std::remove_cvref_t<Problem>::eval_f_grad_f, &problem, x, grad_fx)); }
    real_t eval_f_g(crvec x, rvec g, unsigned int i) const requires requires { &std::remove_cvref_t<Problem>::eval_f_g; } { ++evaluations->f_g; return timed(evaluations->time.f_g, std::bind(&std::remove_cvref_t<Problem>::eval_f_g, &problem, x, g)); }
    real_t eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g, unsigned int i) const requires requires { &std::remove_cvref_t<Problem>::eval_f_grad_f_g; } { ++evaluations->f_grad_f_g; return timed(evaluations->time.f_grad_f_g, std::bind(&std::remove_cvref_t<Problem>::eval_f_grad_f_g, &problem, x, grad_fx, g)); }
    void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f, rvec grad_gxy, unsigned int i) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_f_grad_g_prod; } { ++evaluations->grad_f_grad_g_prod; return timed(evaluations->time.grad_f_grad_g_prod, std::bind(&std::remove_cvref_t<Problem>::eval_grad_f_grad_g_prod, &problem, x, y, grad_f, grad_gxy)); }
    void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n, unsigned int i) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_L; } { ++evaluations->grad_L; return timed(evaluations->time.grad_L, std::bind(&std::remove_cvref_t<Problem>::eval_grad_L, &problem, x, y, grad_L, work_n)); }
    real_t eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ, unsigned int i) const requires requires { &std::remove_cvref_t<Problem>::eval_ψ; } { ++evaluations->ψ; return timed(evaluations->time.ψ, std::bind(&std::remove_cvref_t<Problem>::eval_ψ, &problem, x, y, Σ, ŷ)); }
    void eval_grad_ψ_from_ŷ(crvec x, crvec ŷ, rvec grad_ψ, rvec work_n, unsigned int i) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_ψ_from_ŷ; } { ++evaluations->grad_ψ_from_ŷ; return timed(evaluations->time.grad_ψ_from_ŷ, std::bind(&std::remove_cvref_t<Problem>::eval_grad_ψ_from_ŷ, &problem, x, ŷ, grad_ψ, work_n)); }
    void eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m, unsigned int i) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_ψ; } { ++evaluations->grad_ψ; return timed(evaluations->time.grad_ψ, std::bind(&std::remove_cvref_t<Problem>::eval_grad_ψ, &problem, x, y, Σ, grad_ψ, work_n, work_m)); }
    real_t eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m, unsigned int i) const requires requires { &std::remove_cvref_t<Problem>::eval_ψ_grad_ψ; } { ++evaluations->ψ_grad_ψ; return timed(evaluations->time.ψ_grad_ψ, std::bind(&std::remove_cvref_t<Problem>::eval_ψ_grad_ψ, &problem, x, y, Σ, grad_ψ, work_n, work_m)); }
    const Box &get_box_C() const requires requires { &std::remove_cvref_t<Problem>::get_box_C; } { return problem.get_box_C(); }
    const Box &get_box_D() const requires requires { &std::remove_cvref_t<Problem>::get_box_D; } { return problem.get_box_D(); }
    void check() const { problem.check(); }

    [[nodiscard]] bool provides_eval_grad_gi() const requires requires (Problem p) { { p.provides_eval_grad_gi() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_gi(); }
    [[nodiscard]] bool provides_eval_hess_L_prod() const requires requires (Problem p) { { p.provides_eval_hess_L_prod() } -> std::convertible_to<bool>; } { return problem.provides_eval_hess_L_prod(); }
    [[nodiscard]] bool provides_eval_hess_L() const requires requires (Problem p) { { p.provides_eval_hess_L() } -> std::convertible_to<bool>; } { return problem.provides_eval_hess_L(); }
    [[nodiscard]] bool provides_eval_hess_ψ() const requires requires (Problem p) { { p.provides_eval_hess_ψ() } -> std::convertible_to<bool>; } { return problem.provides_eval_hess_ψ(); }
    [[nodiscard]] bool provides_eval_f_grad_f() const requires requires (Problem p) { { p.provides_eval_f_grad_f() } -> std::convertible_to<bool>; } { return problem.provides_eval_f_grad_f(); }
    [[nodiscard]] bool provides_eval_f_g() const requires requires (Problem p) { { p.provides_eval_f_g() } -> std::convertible_to<bool>; } { return problem.provides_eval_f_g(); }
    [[nodiscard]] bool provides_eval_f_grad_f_g() const requires requires (Problem p) { { p.provides_eval_f_grad_f_g() } -> std::convertible_to<bool>; } { return problem.provides_eval_f_grad_f_g(); }
    [[nodiscard]] bool provides_eval_grad_f_grad_g_prod() const requires requires (Problem p) { { p.provides_eval_grad_f_grad_g_prod() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_f_grad_g_prod(); }
    [[nodiscard]] bool provides_eval_grad_L() const requires requires (Problem p) { { p.provides_eval_grad_L() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_L(); }
    [[nodiscard]] bool provides_eval_ψ() const requires requires (Problem p) { { p.provides_eval_ψ() } -> std::convertible_to<bool>; } { return problem.provides_eval_ψ(); }
    [[nodiscard]] bool provides_eval_grad_ψ_from_ŷ() const requires requires (Problem p) { { p.provides_eval_grad_ψ_from_ŷ() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_ψ_from_ŷ(); }
    [[nodiscard]] bool provides_eval_grad_ψ() const requires requires (Problem p) { { p.provides_eval_grad_ψ() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_ψ(); }
    [[nodiscard]] bool provides_eval_ψ_grad_ψ() const requires requires (Problem p) { { p.provides_eval_ψ_grad_ψ() } -> std::convertible_to<bool>; } { return problem.provides_eval_ψ_grad_ψ(); }
    [[nodiscard]] bool provides_get_box_C() const requires requires (Problem p) { { p.provides_get_box_C() } -> std::convertible_to<bool>; } { return problem.provides_get_box_C(); }
    [[nodiscard]] bool provides_get_box_D() const requires requires (Problem p) { { p.provides_get_box_D() } -> std::convertible_to<bool>; } { return problem.provides_get_box_D(); }
    // clang-format on

    [[nodiscard]] length_t get_n() const { return problem.get_n(); }
    [[nodiscard]] length_t get_m() const { return problem.get_m(); }

    std::shared_ptr<EvalCounter> evaluations = std::make_shared<EvalCounter>();
    Problem problem;

    ParaProblemWithCounters(const Problem &problem) : problem(problem) {}
    ParaProblemWithCounters(Problem &&problem)
        requires(!std::is_lvalue_reference_v<Problem>)
        : problem(std::forward<Problem>(problem)) {}

    /// Reset all evaluation counters and timers to zero. Affects all instances
    /// that share the same evaluations. If you only want to reset the counters
    /// of this instance, use @ref decouple_evaluations first.
    void reset_evaluations() { evaluations.reset(); }
    /// Give this instance its own evaluation counters and timers, decoupling
    /// it from any other instances they might have previously been shared with.
    /// The evaluation counters and timers are preserved (a copy is made).
    void decouple_evaluations() { evaluations = std::make_shared<EvalCounter>(*evaluations); }

  private:
    template <class TimeT, class FunT>
    static decltype(auto) timed(TimeT &time, FunT &&f) {
        detail::Timed timed{time};
        return std::forward<FunT>(f)();
    }
};

/// Wraps the given problem into a @ref ProblemWithCounters and keeps track of
/// how many times each function is called, and how long these calls took.
/// The wrapper has its own copy of the given problem. Making copies of the
/// wrapper also copies the underlying problem, but does not copy the evaluation
/// counters, all copies share the same counters.
template <class Problem>
[[nodiscard]] auto para_problem_with_counters(Problem &&p) {
    using Prob        = std::remove_cvref_t<Problem>;
    using ProbWithCnt = ParaProblemWithCounters<Prob>;
    return ProbWithCnt{std::forward<Problem>(p)};
}

/// Wraps the given problem into a @ref ProblemWithCounters and keeps track of
/// how many times each function is called, and how long these calls took.
/// The wrapper keeps only a reference to the given problem, it is the
/// responsibility of the caller to make sure that the wrapper does not outlive
/// the original problem. Making copies of the wrapper does not copy the
/// evaluation counters, all copies share the same counters.
template <class Problem>
[[nodiscard]] auto para_problem_with_counters_ref(Problem &p) {
    using Prob        = std::remove_cvref_t<Problem>;
    using ProbWithCnt = ParaProblemWithCounters<const Prob &>;
    return ProbWithCnt{p};
}

template <Config Conf>
class ParaBoxConstrProblem {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using Box = alpaqa::Box<config_t>;

    /// Number of decision variables, dimension of x
    length_t n;
    /// Number of constraints, dimension of g(x) and z
    length_t m;

    ParaBoxConstrProblem(length_t n, ///< Number of decision variables
                     length_t m) ///< Number of constraints
        : n{n}, m{m} {}

    ParaBoxConstrProblem(Box C, Box D)
        : n{C.lowerbound.size()}, m{D.lowerbound.size()}, C{std::move(C)}, D{std::move(D)} {}

    ParaBoxConstrProblem(const ParaBoxConstrProblem &)                = default;
    ParaBoxConstrProblem &operator=(const ParaBoxConstrProblem &)     = default;
    ParaBoxConstrProblem(ParaBoxConstrProblem &&) noexcept            = default;
    ParaBoxConstrProblem &operator=(ParaBoxConstrProblem &&) noexcept = default;

    /// Constraints of the decision variables, @f$ x \in C @f$
    Box C{this->n};
    /// Other constraints, @f$ g(x) \in D @f$
    Box D{this->m};

    /// Number of decision variables, @ref n
    length_t get_n() const { return n; }
    /// Number of constraints, @ref m
    length_t get_m() const { return m; }

    /// @f$ \hat x = \Pi_C(x - \gamma\nabla\psi(x)) @f$
    /// @f$ p = \hat x - x @f$
    static real_t eval_proj_grad_step_box(const Box &C, real_t γ, crvec x, crvec grad_ψ, rvec x̂,
                                          rvec p, unsigned int i) {
        p[i] = std::min(std::max(-γ * grad_ψ[i], C.upperbound[i]-x[i]), C.lowerbound[i]-x[i]);
        x̂[i] = x[i] + p[i];
        return real_t(0);
    }

    /// @see @ref TypeErasedProblem::eval_prox_grad_step
    real_t eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p, unsigned int i) const {
        return eval_proj_grad_step_box(C, γ, x, grad_ψ, x̂, p, i);
    }

    /// @see @ref TypeErasedProblem::eval_proj_diff_g
    void eval_proj_diff_g(crvec z, rvec p, unsigned int i) const { p[i] = alpaqa::para_projecting_difference(z, D, i); }

    static void eval_proj_multipliers_box(const Box &D, rvec y, real_t M,
                                          index_t penalty_alm_split) {
        auto max_lb = [M](real_t y, real_t z_lb) {
            real_t y_lb = z_lb == -alpaqa::inf<config_t> ? 0 : -M;
            return std::max(y, y_lb);
        };
        auto min_ub = [M](real_t y, real_t z_ub) {
            real_t y_ub = z_ub == alpaqa::inf<config_t> ? 0 : M;
            return std::min(y, y_ub);
        };
        auto num_alm    = y.size() - penalty_alm_split;
        auto &&y_alm    = y.bottomRows(num_alm);
        auto &&z_alm_lb = D.lowerbound.bottomRows(num_alm);
        auto &&z_alm_ub = D.upperbound.bottomRows(num_alm);
        y_alm           = y_alm.binaryExpr(z_alm_lb, max_lb).binaryExpr(z_alm_ub, min_ub);
    }

    /// @see @ref TypeErasedProblem::eval_proj_multipliers
    void eval_proj_multipliers(rvec y, real_t M, index_t penalty_alm_split) const {
        eval_proj_multipliers_box(D, y, M, penalty_alm_split);
    }

    /// @see @ref TypeErasedProblem::get_box_C
    const Box &get_box_C() const { return C; }
    /// @see @ref TypeErasedProblem::get_box_D
    const Box &get_box_D() const { return D; }

    /// @see @ref TypeErasedProblem::check
    void check() const {
        util::check_dim_msg<config_t>(
            C.lowerbound, n,
            "Length of problem.C.lowerbound does not match problem size problem.n");
        util::check_dim_msg<config_t>(
            C.upperbound, n,
            "Length of problem.C.upperbound does not match problem size problem.n");
        util::check_dim_msg<config_t>(
            D.lowerbound, m,
            "Length of problem.D.lowerbound does not match problem size problem.m");
        util::check_dim_msg<config_t>(
            D.upperbound, m,
            "Length of problem.D.upperbound does not match problem size problem.m");
    }
};

/// Problem class that allows specifying the basic functions as C++
/// `std::function`s.
/// @ingroup grp_Problems
template <Config Conf = DefaultConfig>
class ParaFunctionalProblem : public ParaBoxConstrProblem<Conf> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using ParaBoxConstrProblem<Conf>::ParaBoxConstrProblem;

    std::function<real_t(crvec)> f;
    std::function<void(crvec, rvec)> grad_f;
    std::function<void(crvec, rvec)> g;
    std::function<void(crvec, crvec, rvec)> grad_g_prod;
    std::function<void(crvec, index_t, rvec)> grad_gi;
    std::function<void(crvec, crvec, crvec, rvec)> hess_L_prod;
    std::function<void(crvec, crvec, rmat)> hess_L;
    std::function<void(crvec, crvec, crvec, crvec, rvec)> hess_ψ_prod;
    std::function<void(crvec, crvec, crvec, rmat)> hess_ψ;

    // clang-format off
    real_t eval_f(crvec x) const { ScopedMallocAllower ma; return f(x); }
    void eval_grad_f(crvec x, rvec grad_fx) const { ScopedMallocAllower ma; grad_f(x, grad_fx); }
    void eval_g(crvec x, rvec gx) const { ScopedMallocAllower ma; g(x, gx); }
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const { ScopedMallocAllower ma; grad_g_prod(x, y, grad_gxy); }
    void eval_grad_gi(crvec x, index_t i, rvec grad_i) const { ScopedMallocAllower ma; grad_gi(x, i, grad_i); }
    void eval_hess_L_prod(crvec x, crvec y, crvec v, rvec Hv) const { ScopedMallocAllower ma; hess_L_prod(x, y, v, Hv); }
    void eval_hess_L(crvec x, crvec y, rmat H) const { ScopedMallocAllower ma; hess_L(x, y, H); }
    void eval_hess_ψ_prod(crvec x, crvec y, crvec Σ, crvec v, rvec Hv) const { ScopedMallocAllower ma; hess_ψ_prod(x, y, Σ, v, Hv); }
    void eval_hess_ψ(crvec x, crvec y, crvec Σ, rmat H) const { ScopedMallocAllower ma; hess_ψ(x, y, Σ, H); }
    // clang-format on

    /// @see @ref TypeErasedProblem::provides_eval_grad_gi
    bool provides_eval_grad_gi() const { return bool{grad_gi}; }
    /// @see @ref TypeErasedProblem::provides_eval_hess_L_prod
    bool provides_eval_hess_L_prod() const { return bool{hess_L_prod}; }
    /// @see @ref TypeErasedProblem::provides_eval_hess_L
    bool provides_eval_hess_L() const { return bool{hess_L}; }
    /// @see @ref TypeErasedProblem::provides_eval_hess_ψ_prod
    bool provides_eval_hess_ψ_prod() const { return bool{hess_ψ_prod}; }
    /// @see @ref TypeErasedProblem::provides_eval_hess_ψ
    bool provides_eval_hess_ψ() const { return bool{hess_ψ}; }

    ParaFunctionalProblem(const ParaFunctionalProblem &)                = default;
    ParaFunctionalProblem &operator=(const ParaFunctionalProblem &)     = default;
    ParaFunctionalProblem(ParaFunctionalProblem &&) noexcept            = default;
    ParaFunctionalProblem &operator=(ParaFunctionalProblem &&) noexcept = default;
};

template <Config Conf>
void para_print_provided_functions(std::ostream &os, const ParaTypeErasedProblem<Conf> &problem) {
    os << "           eval_grad_gi: " << problem.provides_eval_grad_gi() << '\n'
       << "       eval_hess_L_prod: " << problem.provides_eval_hess_L_prod() << '\n'
       << "            eval_hess_L: " << problem.provides_eval_hess_L() << '\n'
       << "       eval_hess_ψ_prod: " << problem.provides_eval_hess_ψ_prod() << '\n'
       << "            eval_hess_ψ: " << problem.provides_eval_hess_ψ() << '\n'
       << "          eval_f_grad_f: " << problem.provides_eval_f_grad_f() << '\n'
       << "               eval_f_g: " << problem.provides_eval_f_g() << '\n'
       << "        eval_f_grad_f_g: " << problem.provides_eval_f_grad_f_g() << '\n'
       << "eval_grad_f_grad_g_prod: " << problem.provides_eval_grad_f_grad_g_prod() << '\n'
       << "            eval_grad_L: " << problem.provides_eval_grad_L() << '\n'
       << "                 eval_ψ: " << problem.provides_eval_ψ() << '\n'
       << "     eval_grad_ψ_from_ŷ: " << problem.provides_eval_grad_ψ_from_ŷ() << '\n'
       << "            eval_grad_ψ: " << problem.provides_eval_grad_ψ() << '\n'
       << "          eval_ψ_grad_ψ: " << problem.provides_eval_ψ_grad_ψ() << '\n'
       << "              get_box_C: " << problem.provides_get_box_C() << '\n'
       << "              get_box_D: " << problem.provides_get_box_D() << '\n';
}

/// @}

} // namespace alpaqa