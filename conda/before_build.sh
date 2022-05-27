
if [[ ! -n CONDA_EXE ]]; then
    echo "ERROR: This script must be executed from a conda environment."
    exit 2
fi

# Make sure that this script is sourced.
# See https://stackoverflow.com/questions/2683279/how-to-detect-if-a-script-is-being-sourced for details.
is_sourced() {
    if [ -n "$ZSH_VERSION" ]; then
        case $ZSH_EVAL_CONTEXT in *:file:*) return 0;; esac
    else  # Add additional POSIX-compatible shell names here, if needed.
        case ${0##*/} in dash|bash|ksh|sh) return 0;; esac
    fi
    return 1
}

if is_sourced ; then
    module unload gcc
    unset CFLAGS
    unset CPPFLAGS
    unset CXXFLAGS
    unset LDFLAGS
    unset FCFLAGS
    unset JAVA_HOME
    unset FC
    export CONDA_BUILD=1
    conda activate build
    unset FC
    export CPPFLAGS="-DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /include"
    export CFLAGS="-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /include -fdebug-prefix-map==/usr/local/src/conda/- -fdebug-prefix-map==/usr/local/src/conda-prefix"
    export LDFLAGS="-Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/lib -Wl,-rpath-link,/lib -L/lib"
else
    echo "ERROR: This script must be sourced."
fi
