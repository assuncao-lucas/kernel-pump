/**
 * @file feasibility_heuristic.cpp
 * @brief Free helpers for the FeasibilityHeuristic interface
 *
 * See docs/spec-heuristic-interface.md. These live apart from any single implementation because
 * both drivers need them: main.cpp resolves `method=<name>` for a standalone run, and KernelPump
 * resolves `kp.heuristic=<name>` for the per-bucket heuristic.
 */

#include <algorithm>
#include <cctype>
#include <iterator>
#include <string>
#include <vector>

#include "kernelpump/feasibility_heuristic.h"

namespace
{
	std::string upcase(const std::string &s)
	{
		std::string out = s;
		std::transform(out.begin(), out.end(), out.begin(),
					   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
		return out;
	}
} // namespace

namespace dominiqs
{

	std::string resolveHeuristicName(const std::string &name)
	{
		if (name.empty())
			return std::string();

		std::vector<std::string> ids;
		HeuristicsFactory::getInstance().getIDs(std::back_inserter(ids));

		const std::string wanted = upcase(name);
		for (const auto &id : ids)
			if (upcase(id) == wanted)
				return id;
		return std::string();
	}

	std::string registeredHeuristicNames()
	{
		std::vector<std::string> ids;
		HeuristicsFactory::getInstance().getIDs(std::back_inserter(ids));
		std::sort(ids.begin(), ids.end());

		std::string out;
		for (const auto &id : ids)
			out += (out.empty() ? "" : ", ") + id;
		return out;
	}

} // namespace dominiqs
