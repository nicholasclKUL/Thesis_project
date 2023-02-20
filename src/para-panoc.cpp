#include </home/nicholas/Thesis_project/src/src/para-panoc.tpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCParams, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCParams, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCParams, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCParams, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCStats, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCStats, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCStats, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCStats, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCStats, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCProgressInfo, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCProgressInfo, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCProgressInfo, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCProgressInfo, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, ParaPANOCProgressInfo, EigenConfigq);
#endif

} // namespace alpaqa