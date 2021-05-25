#! /bin/sh

./@compile-it
cd main
./@compile-it
cd ../tools
./@compile-it
cd ../

rm ../bin/swift_readpl4glisser
cp main/swift_readpl4glisser ../bin/.

exit

