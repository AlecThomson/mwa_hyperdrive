#!/bin/bash

# Make a hyperdrive binary portable by copying the shared libraries it needs
# (except the core system libraries that every Linux has) next to it, and
# pointing its rpath at them. This is the same idea as Python's auditwheel.
#
# Usage: bundle_libs.sh <path/to/hyperdrive> <output dir>
# The output dir gets bin/hyperdrive and lib/*.so*.

set -euo pipefail

BINARY="$1"
OUT="$2"

mkdir -p "$OUT/bin" "$OUT/lib"
cp "$BINARY" "$OUT/bin/hyperdrive"

# Libraries that must come from the host system (glibc and friends, and the
# C++/GCC runtimes, which are backwards compatible).
EXCLUDE='^(linux-vdso|ld-linux|libc\.so|libm\.so|libmvec\.so|libpthread\.so|libdl\.so|librt\.so|libutil\.so|libresolv\.so|libstdc\+\+\.so|libgcc_s\.so)'

# Copy dependencies, then their dependencies, until nothing new is found.
copy_deps() {
    ldd "$1" | awk '/=> \// {print $1, $3}' | while read -r name path; do
        if [[ "$name" =~ $EXCLUDE ]]; then
            continue
        fi
        if [ ! -e "$OUT/lib/$name" ]; then
            echo "Bundling $name ($path)"
            cp -L "$path" "$OUT/lib/$name"
            chmod u+w "$OUT/lib/$name"
            copy_deps "$OUT/lib/$name"
        fi
    done
}
copy_deps "$OUT/bin/hyperdrive"

patchelf --set-rpath '$ORIGIN/../lib' "$OUT/bin/hyperdrive"
for lib in "$OUT"/lib/*.so*; do
    patchelf --set-rpath '$ORIGIN' "$lib"
done

# Check that everything resolves to the bundled copies or the excluded
# system libraries.
if ldd "$OUT/bin/hyperdrive" | grep -q "not found"; then
    ldd "$OUT/bin/hyperdrive"
    echo "Some libraries weren't found" >&2
    exit 1
fi
ldd "$OUT/bin/hyperdrive"
