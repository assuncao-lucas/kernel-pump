/**
 * @file consolelog.h
 * @brief Console logging functions
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * 2023
 */

#ifndef CONSOLELOG_H
#define CONSOLELOG_H

#include <fmt/format.h>
#include <fmt/color.h>
#include <unistd.h>

// macro type savers
#define LOG_ITEM(name, value) consoleLog("{} = {}", name, value)
#define LOG_CONFIG(what) LOG_ITEM("fp." #what, what)

template <typename... Args>
void consoleLog(Args &&...args)
{
#ifndef SILENT_EXEC
	fmt::print(std::forward<Args>(args)...);
	fmt::print("\n");
#endif //< SILENT_EXEC
}

template <typename... Args>
void consoleInfo(Args &&...args)
{
#ifndef SILENT_EXEC
	if (isatty(STDOUT_FILENO))
		fmt::print(fmt::fg(fmt::color::green), std::forward<Args>(args)...);
	else
		fmt::print(std::forward<Args>(args)...);
	fmt::print("\n");
#endif //< SILENT_EXEC
}

template <typename... Args>
void consoleWarn(Args &&...args)
{
#ifndef SILENT_EXEC
	if (isatty(STDOUT_FILENO))
		fmt::print(fmt::fg(fmt::color::yellow), std::forward<Args>(args)...);
	else
		fmt::print(std::forward<Args>(args)...);
	fmt::print("\n");
#endif //< SILENT_EXEC
}

template <typename... Args>
void consoleError(Args &&...args)
{
#ifndef SILENT_EXEC
	if (isatty(STDOUT_FILENO))
		fmt::print(fmt::fg(fmt::color::red), std::forward<Args>(args)...);
	else
		fmt::print(std::forward<Args>(args)...);
	fmt::print("\n");
#endif //< SILENT_EXEC
}

#ifdef DEBUG_LOG

enum class DebugLevel
{
	None = 0,
	Minimal = 1,
	Normal = 2,
	Verbose = 3,
	VeryVerbose = 4,
	All = 5
};

#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL DebugLevel::Normal
#endif //< DEBUG_LEVEL

template <typename... Args>
void consoleDebug(DebugLevel level, Args &&...args)
{
	if (level > DEBUG_LEVEL)
		return;
	consoleLog(std::forward<Args>(args)...);
}

#else //< DEBUG_LOG

#define consoleDebug(level, ...) \
	do                           \
	{                            \
	} while (false)

#endif //<DEBUG_LOG

#endif /* CONSOLELOG_H */
