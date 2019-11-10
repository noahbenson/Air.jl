#! /bin/bash

AIRPATH="`dirname \"$0\"`"
DIR0="$PWD"
cd "$AIRPATH"
AIRPATH="$PWD"
cd "$DIR0"

julia -q -p auto <<EOF
using Pkg
let airpath = "${AIRPATH}", ps = Pkg.PackageSpec(path=airpath)
    Pkg.develop(ps)
    Pkg.test("Air")
end
exit()
EOF



