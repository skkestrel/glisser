#! /bin/sh

./@compile-it
cd main
./@compile-it
cd ../tools
./@compile-it
cd ../

exit

