#! /bin/sh
set -e
platforms='aarch64 x86_64'
root=$(realpath $(dirname $0)/..)
for platform in $platforms; do
  echo platform=$platform
  (
    cd $root/docker
    docker build . -t iterative-solver-$platform --platform linux/$platform
  )
  mkdir -p build-$platform
  docker run --rm -v $root:$root -v $root/build-$platform:$root/build -w $root \
    --platform linux/$platform iterative-solver-$platform \
    python -m build
done
