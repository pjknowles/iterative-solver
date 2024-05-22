#! /bin/sh
set -e
platforms='aarch64 x86_64'
root=$(realpath $(dirname $0)/..)
for platform in $platforms; do
  echo platform=$platform
  (
    cd $root/python
    docker build . -t iterative-solver-$platform --platform linux/$platform
  )
  docker run --rm -v $root:$root -w $root \
    --platform linux/$platform iterative-solver-$platform \
    bash --login -c " pip install -v -e python; python -c 'import iterative_solver; print(iterative_solver.__version__)';  python python/test_optimize.py   "
done
