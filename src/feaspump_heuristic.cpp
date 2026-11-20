/**
 * @file feaspump_heuristic.cpp
 * @brief Factory registration for the FeasibilityPump adapter
 *
 * See docs/spec-heuristic-interface.md section 7.1.
 *
 * FeasibilityPumpHeuristic is header-only forwarding; this translation unit exists to register it
 * with HeuristicsFactory. Registration uses the same static-recorder idiom as the transformers
 * (src/transformers.cpp) and var rankers (src/ranking.cpp), which works because libkp is
 * whole-archive linked into kp and results.
 */

#include <iostream>

#include "kernelpump/feaspump_heuristic.h"

using namespace dominiqs;

// auto registration

class HEUR_FACTORY_RECORDER
{
public:
	HEUR_FACTORY_RECORDER()
	{
#ifndef SILENT_EXEC
		std::cout << "Registering FeasibilityHeuristics...";
#endif //< SILENT_EXEC
		HeuristicsFactory::getInstance().registerClass<FeasibilityPumpHeuristic>("feaspump");
#ifndef SILENT_EXEC
		std::cout << "done" << std::endl;
#endif //< SILENT_EXEC
	}
};

HEUR_FACTORY_RECORDER my_heur_factory_recorder;
