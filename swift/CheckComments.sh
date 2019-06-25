#! /bin/bash

grep -r -l "	1" * 2> /dev/null | grep "\.f$"

exit
