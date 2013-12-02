#!/bin/csh

set bin = ../bin
if (! -d ../bin ) mkdir $bin
#if (! -d ../bin/ppc) mkdir $bin/ppc
if (! -d ../bin/i386) mkdir $bin/i386
if (! -d ../bin/x86_64) mkdir $bin/x86_64

# cd ../src
# rm *.o
# make -f ../make/Makefile.os_x all
# make -f ../make/Makefile.os_x uinstall

rm *.o
make -f ../make/Makefile.os_x86 all
make -f ../make/Makefile.os_x86 uinstall

rm *.o
make -f ../make/Makefile.os_x86_64 all
make -f ../make/Makefile.os_x86_64 uinstall
rm *.o
cd ../bin
foreach n ( i386/* )
set f=$n:t
#lipo -create ppc/$f i386/$f x86_64/$f -output $f
lipo -create i386/$f x86_64/$f -output $f
echo "Universal $f built"
end
#rm -rf ppc/ i386/ x86_64/
#rm -rf i386/ x86_64/
echo "Done!"
