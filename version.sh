#!/bin/sh -e

if [ -e .git ]; then
    version=$(git describe --tags | cut -d - -f 1-2 | tr - .)
elif [ -n "$VERSION" ]; then
    version=$VERSION
else
    version="Unknown-version"
fi
echo $version

