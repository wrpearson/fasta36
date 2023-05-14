#!/bin/csh

## should start from fasta_dir/src

set bin = ../bin
if (! -d ../bin ) mkdir $bin
# if (! -d ../bin/ppc) mkdir $bin/ppc
# if (! -d ../bin/i386) mkdir $bin/i386
if (! -d ../bin/x86_64) mkdir $bin/x86_64
if (! -d ../bin/arm64) mkdir $bin/arm64

## cd ../src
# rm *.o
# make -f ../make/Makefile.os_x_ppc all
# make -f ../make/Makefile.os_x_ppc uinstall

# rm *.o
# make -f ../make/Makefile.os_x86 all
# make -f ../make/Makefile.os_x86 uinstall

rm *.o
make -f ../make/Makefile.os_x86_64 all
make -f ../make/Makefile.os_x86_64 uinstall

rm *.o
make -f ../make/Makefile.os_x_arm64 all
make -f ../make/Makefile.os_x_arm64 uinstall

rm *.o

pushd ../bin
foreach n ( x86_64/* arm64/*)
set f=$n:t
#lipo -create ppc/$f i386/$f x86_64/$f -output $f
lipo -create x86_64/$f arm64/$f -output $f
echo "Universal $f built"
end

# rm -rf x86_64 arm64
popd

echo "Done!"
