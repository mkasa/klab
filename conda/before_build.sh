# This script must be SOURCED, e.g.
#
#     . conda/before_build.sh
#
# It prepares the shell for running conda-build: it drops the compiler flags
# the surrounding environment may have set, activates the build environment
# and exports the flags conda-build expects.

# Make sure that this script is sourced.
# See https://stackoverflow.com/questions/2683279/how-to-detect-if-a-script-is-being-sourced for details.
is_sourced() {
    if [ -n "$ZSH_VERSION" ]; then
        case $ZSH_EVAL_CONTEXT in *:file:*) return 0;; esac
    else  # Add additional POSIX-compatible shell names here, if needed.
        # In a LOGIN shell $0 is '-bash' (not 'bash'), so a leading '-' has to
        # be stripped before matching.
        _bb_shell_name=${0##*/}
        _bb_shell_name=${_bb_shell_name#-}
        case $_bb_shell_name in dash|bash|ksh|sh|zsh) return 0;; esac
    fi
    return 1
}

if ! is_sourced ; then
    echo "ERROR: This script must be sourced (e.g. '. conda/before_build.sh')." >&2
    # Not sourced, so this is a real subprocess and exiting is safe.
    exit 2
fi

# From here on the script is known to be sourced: use 'return', never 'exit',
# so that a failure does not kill the user's interactive shell.
if [ -z "$CONDA_EXE" ]; then
    echo "ERROR: This script must be executed from a conda environment." >&2
    return 2
fi

# 'module' is a shell function provided by environment-modules, which is not
# installed everywhere.
if command -v module >/dev/null 2>&1; then
    module unload gcc
fi

unset CFLAGS
unset CPPFLAGS
unset CXXFLAGS
unset LDFLAGS
unset FCFLAGS
unset JAVA_HOME
unset FC
export CONDA_BUILD=1

# The name of the environment used for building. Override with
#     CONDA_BUILD_ENV=<name> . conda/before_build.sh
CONDA_BUILD_ENV=${CONDA_BUILD_ENV:-build}
if conda env list 2>/dev/null | awk '{print $1}' | grep -qx -- "$CONDA_BUILD_ENV"; then
    conda activate "$CONDA_BUILD_ENV"
else
    echo "WARNING: conda environment '$CONDA_BUILD_ENV' does not exist;" >&2
    echo "         staying in the currently activated environment." >&2
fi

# The flags below used to be written with a literal, UNEXPANDED '${PREFIX}',
# which left '-isystem /include', '-L/lib' and '-Wl,-rpath,/lib' behind, i.e.
# the HOST system directories were prepended to the include/link path and
# baked into RPATH. Resolve the prefix explicitly instead.
_bb_prefix=${PREFIX:-$CONDA_PREFIX}
if [ -z "$_bb_prefix" ]; then
    echo "ERROR: neither PREFIX nor CONDA_PREFIX is set; cannot set the build flags." >&2
    unset _bb_shell_name
    return 2
fi

# '-march=nocona -mtune=haswell' is x86-64 only and fails outright on
# aarch64 / Apple Silicon, so it is only used where it is valid.
_bb_arch_flags=""
case "$(uname -m)" in
    x86_64|amd64) _bb_arch_flags="-march=nocona -mtune=haswell" ;;
esac

# The two -fdebug-prefix-map entries conda-build normally emits map the source
# directory and the prefix onto stable paths, so fall back sensibly when the
# corresponding variables are not set.
_bb_src=${SRC_DIR:-$_bb_prefix}
_bb_pkg="${PKG_NAME:-unknown}-${PKG_VERSION:-0}"

export CPPFLAGS="-DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem ${_bb_prefix}/include"
export CFLAGS="${_bb_arch_flags} -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem ${_bb_prefix}/include -fdebug-prefix-map=${_bb_src}=/usr/local/src/conda/${_bb_pkg} -fdebug-prefix-map=${_bb_prefix}=/usr/local/src/conda-prefix"
# NOTE: most of these -Wl options are GNU ld specific and are rejected by the
# macOS linker; this script targets Linux builds.
export LDFLAGS="-Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,${_bb_prefix}/lib -Wl,-rpath-link,${_bb_prefix}/lib -L${_bb_prefix}/lib"

unset _bb_shell_name
unset _bb_prefix
unset _bb_arch_flags
unset _bb_src
unset _bb_pkg
