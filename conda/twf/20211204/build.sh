
export GOPATH=$PWD
go install github.com/wvanlint/twf/cmd/twf\@latest
mkdir -p ${PREFIX}/bin
cp ${GOPATH}/bin/twf ${PREFIX}/bin/twf

