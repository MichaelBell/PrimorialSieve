#!/bin/bash

set -e

ESDK=${EPIPHANY_HOME}
ELIBS=${ESDK}/tools/host/lib
EINCS=${ESDK}/tools/host/include
ELDF=link/fast.ldf

# Create the binaries directory
mkdir -p bin/

CROSS_PREFIX=
case $(uname -p) in
	arm*)
		# Use native arm compiler (no cross prefix)
		CROSS_PREFIX=
		;;
	   *)
		# Use cross compiler
		CROSS_PREFIX="arm-linux-gnueabihf-"
		;;
esac

# Build HOST side application
${CROSS_PREFIX}gcc -O3 -std=gnu99 src/primorial.c -o bin/primorial.elf -I ${EINCS} -L ${ELIBS} -le-hal  -lm -lpthread #-le-loader

# Build DEVICE side program
e-gcc -O3 -g -T ${ELDF} -std=c99 src/e_primorial.c -o bin/e_primorial.elf -le-lib -lm -ffast-math -mfp-mode=int
e-gcc -O3 -S -T ${ELDF} -std=c99 src/e_primorial.c -o bin/e_primorial.s -le-lib -lm -ffast-math -mfp-mode=int

# Convert ebinary to SREC file
e-objcopy --srec-forceS3 --output-target srec bin/e_primorial.elf bin/e_primorial.srec


