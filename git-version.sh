#!/bin/sh -e

v=$(git describe --tags | cut -d - -f 1-2)
commit=${v#*-}
version=${v%-*}
echo $version-$(($commit + 1))

