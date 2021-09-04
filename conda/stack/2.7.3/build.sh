

mkdir -p ${PREFIX}/bin
cp stack ${PREFIX}/bin/

mkdir -p ${PREFIX}/share/haskell-stack
cp -r doc LICENSE *.md ${PREFIX}/share/haskell-stack/
