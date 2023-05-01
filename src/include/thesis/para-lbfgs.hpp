#include<alpaqa/accelerators/lbfgs.hpp>
#include<alpaqa/implementation/accelerators/lbfgs.tpp>

#include<Kokkos_Core.hpp>


namespace alpaqa {

template <Config Conf = DefaultConfig>
class ParaLBFGS {
  public:
    USING_ALPAQA_CONFIG(Conf);

    using Params = LBFGSParams<config_t>;

    /// The sign of the vectors @f$ p @f$ passed to the @ref update method.
    enum class Sign {
        Positive, ///< @f$ p \sim \nabla \psi(x) @f$
        Negative, ///< @f$ p \sim -\nabla \psi(x) @f$
    };

    ParaLBFGS() = default;
    ParaLBFGS(Params params) : params(params) {}
    ParaLBFGS(Params params, length_t n) : params(params) { resize(n); }

    /// Check if the new vectors s and y allow for a valid BFGS update that
    /// preserves the positive definiteness of the Hessian approximation.
    static bool update_valid(const Params &params, real_t yᵀs, real_t sᵀs,
                             real_t pᵀp);

    /// Update the inverse Hessian approximation using the new vectors
    /// sₖ = xₙₑₓₜ - xₖ and yₖ = pₙₑₓₜ - pₖ.
    bool update_sy(crvec s, crvec y, real_t pₙₑₓₜᵀpₙₑₓₜ, bool forced = false);
    /// @see @ref update_sy
    bool update_sy_impl(const auto &s, const auto &y, real_t pₙₑₓₜᵀpₙₑₓₜ,
                        bool forced = false);

    /// Update the inverse Hessian approximation using the new vectors xₙₑₓₜ
    /// and pₙₑₓₜ.
    bool update(crvec xₖ, crvec xₙₑₓₜ, crvec pₖ, crvec pₙₑₓₜ,
                Sign sign = Sign::Positive, bool forced = false);

    /// Apply the inverse Hessian approximation to the given vector q.
    /// Initial inverse Hessian approximation is set to @f$ H_0 = \gamma I @f$.
    /// If @p γ is negative, @f$ H_0 = \frac{s^\top y}{y^\top y} I @f$.
    bool apply(rvec q, real_t γ = -1) const;
    bool apply(rvec q, real_t γ = -1, index_t nthrds) const {
        if (idx == 0 && not full)
            return false;

        // If the step size is negative, compute it as sᵀy/yᵀy
        if (params.stepsize == LBFGSStepSize::BasedOnCurvature || γ < 0) {
            auto new_idx = pred(idx);
            real_t yᵀy   = y(new_idx).squaredNorm();
            γ            = 1 / (ρ(new_idx) * yᵀy);
        }

        foreach_rev([&](index_t i) {
            α(i) = ρ(i) * s(i).dot(q);
            q -= α(i) * y(i);
        });

        // r ← H_0 q
        q *= γ;

        foreach_fwd([&](index_t i) {
            real_t β = ρ(i) * y(i).dot(q);
            q -= (β - α(i)) * s(i);
        });

        return true;
    }

    /// Apply the inverse Hessian approximation to the given vector q, applying
    /// only the columns and rows of the Hessian in the index set J.
    bool apply_masked(rvec q, real_t γ, crindexvec J) const;
    /// @copydoc apply_masked(rvec, real_t, crindexvec) const
    bool apply_masked(rvec q, real_t γ, const std::vector<index_t> &J) const;
    /// @copydoc apply_masked(rvec, real_t, crindexvec) const
    bool apply_masked_impl(rvec q, real_t γ, const auto &J) const;

    /// Throw away the approximation and all previous vectors s and y.
    void reset();
    /// Re-allocate storage for a problem with a different size. Causes
    /// a @ref reset.
    void resize(length_t n);

    /// Scale the stored y vectors by the given factor.
    void scale_y(real_t factor);

    /// Get a string identifier for this accelerator.
    std::string get_name() const {
        return "LBFGS<" + std::string(config_t::get_name()) + '>';
    }
    /// Get the parameters.
    const Params &get_params() const { return params; }

    /// Get the size of the s and y vectors in the buffer.
    length_t n() const { return sto.n(); }
    /// Get the number of previous vectors s and y stored in the buffer.
    length_t history() const { return sto.history(); }
    /// Get the next index in the circular buffer of previous s and y vectors.
    index_t succ(index_t i) const { return i + 1 < history() ? i + 1 : 0; }
    /// Get the previous index in the circular buffer of s and y vectors.
    index_t pred(index_t i) const { return i > 0 ? i - 1 : history() - 1; }
    /// Get the number of previous s and y vectors currently stored in the
    /// buffer.
    length_t current_history() const { return full ? history() : idx; }

    auto s(index_t i) { return sto.s(i); }
    auto s(index_t i) const { return sto.s(i); }
    auto y(index_t i) { return sto.y(i); }
    auto y(index_t i) const { return sto.y(i); }
    real_t &ρ(index_t i) { return sto.ρ(i); }
    real_t &ρ(index_t i) const { return sto.ρ(i); }
    real_t &α(index_t i) { return sto.α(i); }
    real_t &α(index_t i) const { return sto.α(i); }

    /// Iterate over the indices in the history buffer, oldest first.
    template <class F>
    void foreach_fwd(const F &fun) const {
        if (full)
            for (index_t i = idx; i < history(); ++i)
                fun(i);
        if (idx)
            for (index_t i = 0; i < idx; ++i)
                fun(i);
    }    

    /// Iterate over the indices in the history buffer, newest first.
    template <class F>
    void foreach_rev(const F &fun) const {
        if (idx)
            for (index_t i = idx; i-- > 0;)
                fun(i);
        if (full)
            for (index_t i = history(); i-- > idx;)
                fun(i);
    }

    real_t dot_product [&](index_t i, index_t N){
        if (i == N) {
            α(i) = ρ(i) * s(i).segment((nx+nu)*N,nu).dot(q.segment((nx+nu)*N,nu));
        } else {
            α(i) = ρ(i) * s(i).segment((nx+nu)*i,nu+nx).dot(q.segment((nx+nu)*i,nu+nx));
        }
    } 

    void sum_reduction [&](index_t i, index_t N){
        if (i == N) {
            α(i) = ρ(i) * s(i).segment((nx+nu)*N,nu).dot(q.segment((nx+nu)*N,nu));
        } else {
            α(i) = ρ(i) * s(i).segment((nx+nu)*i,nu+nx).dot(q.segment((nx+nu)*i,nu+nx));
        }
    }

  private:
    LBFGSStorage<config_t> sto;
    index_t idx = 0;
    bool full   = false;
    Params params;
};


}