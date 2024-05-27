#! /bin/sh
set -e
platforms="$@"
if [ -z "$platforms" ]; then
  platforms='aarch64 x86_64'
fi
root=$(
  cd -- $(dirname "$0")/.. >/dev/null 2>&1
  pwd -P
)
echo root=$root
for platform in $platforms; do
  echo platform=$platform
  docker build $root/python -t iterative-solver-$platform --platform linux/$platform
  docker run --rm -v $root:$root -w $root \
    --platform linux/$platform iterative-solver-$platform \
    bash --login -c " bash python/build.sh; cd /tmp; python -c 'import iterative_solver; print(iterative_solver.__version__)';  python $root/python/test/test_optimize.py   "
  ls -lR $CONDA_PREFIX/lib/python*/site-packages/iterative_solver
done
