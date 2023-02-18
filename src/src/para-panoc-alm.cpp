#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <Thesis/para-panoc.hpp>
#include <Thesis/para-alm.hpp>
#include <Thesis/para-panoc-alm.hpp>

namespace alpaqa {

// clang-format off
ALPAQA_EXPORT_TEMPLATE(class, ParaPANOCSolver, LBFGSDirection<DefaultConfig>);
ALPAQA_EXPORT_TEMPLATE(class, ParaPANOCSolver, LBFGSDirection<EigenConfigf>);
ALPAQA_EXPORT_TEMPLATE(class, ParaPANOCSolver, LBFGSDirection<EigenConfigd>);
ALPAQA_EXPORT_TEMPLATE(class, ParaPANOCSolver, LBFGSDirection<EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, ParaPANOCSolver, LBFGSDirection<EigenConfigq>);
#endif

ALPAQA_EXPORT_TEMPLATE(class, ParaALMSolver, ParaPANOCSolver<LBFGSDirection<DefaultConfig>>);
ALPAQA_EXPORT_TEMPLATE(class, ParaALMSolver, ParaPANOCSolver<LBFGSDirection<EigenConfigf>>);
ALPAQA_EXPORT_TEMPLATE(class, ParaALMSolver, ParaPANOCSolver<LBFGSDirection<EigenConfigd>>);
ALPAQA_EXPORT_TEMPLATE(class, ParaALMSolver, ParaPANOCSolver<LBFGSDirection<EigenConfigl>>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, ParaALMSolver, ParaPANOCSolver<LBFGSDirection<EigenConfigq>>);
#endif
// clang-format on

} // namespace alpaqa