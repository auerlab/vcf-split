#!/bin/sh -e

if [ -e .git ]; then
    v=$(git describe --tags | cut -d - -f 1-2)
    commit=${v#*-}
    version=${v%-*}
elif [ -n "$PORTVERSION" ]; then
    version=$PORTVERSION
elif [ -n "$PKGVERSION" ]; then
    version=$PKGVERSION
else
    version=Unknown
fi
echo $version

