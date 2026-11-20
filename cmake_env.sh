#!/bin/bash
# Source this before invoking cmake:  . ./cmake_env.sh
#
# Why this exists
# ---------------
# Some vendor packages install their private dependency libraries into a directory that they then
# add to the *global* loader path via /etc/ld.so.conf.d/. When one of those bundled libraries is
# older than the system's, it shadows the system copy for every process on the machine. If the
# bundled library in turn needs a soname the distribution has since dropped, any program linked
# against it stops starting at all.
#
# Seen on this project's dev machine (Ubuntu 26.04): Dell's dca-enabler ships
# /etc/ld.so.conf.d/dcaenabler.conf pointing at /usr/lib/dcaenabler/, which contains an Ubuntu-22
# libcurl.so.4 built against OpenLDAP 2.5. Ubuntu 26.04 ships OpenLDAP 2.6 (libldap.so.2), so that
# libcurl can no longer resolve libldap-2.5.so.0 -- and cmake, which links libcurl, dies with
#
#   cmake: error while loading shared libraries: libldap-2.5.so.0: cannot open shared object file
#
# (/usr/bin/curl is hit by the same thing.) The system's own libcurl is fine; it is simply being
# shadowed.
#
# This is a machine configuration bug, not a project bug, and the real fix needs root:
#
#   sudo mv /etc/ld.so.conf.d/dcaenabler.conf /etc/ld.so.conf.d/dcaenabler.conf.disabled
#   sudo ldconfig
#
# Until that is done, the workaround below lets the build scripts run. It is conditional: on a
# healthy machine it changes nothing, because cmake already starts.

if ! command -v cmake >/dev/null 2>&1; then
	echo "cmake_env.sh: cmake not found on PATH" >&2
elif ! cmake --version >/dev/null 2>&1; then
	_kp_syslib="/usr/lib/$(uname -m)-linux-gnu"
	if [ -d "$_kp_syslib" ]; then
		export LD_LIBRARY_PATH="$_kp_syslib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
		if cmake --version >/dev/null 2>&1; then
			echo "cmake_env.sh: cmake could not start; prepended $_kp_syslib to LD_LIBRARY_PATH."
			echo "cmake_env.sh: a library on the global loader path is shadowing a system library."
			echo "cmake_env.sh: see the comments in cmake_env.sh for the permanent fix."
		else
			echo "cmake_env.sh: cmake still cannot start. Missing libraries:" >&2
			ldd "$(command -v cmake)" 2>/dev/null | grep 'not found' >&2
		fi
	fi
	unset _kp_syslib
fi
