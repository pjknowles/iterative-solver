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
  docker run --rm -v $PWD:$PWD -v $PWD/build-$platform:$PWD/build -w $PWD \
    --platform linux/$platform iterative-solver-$platform \
    pip install dist/*${platform}*whl
done
