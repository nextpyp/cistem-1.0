#!/bin/bash

# assume libcore.a is in the same path as this script

# create folders for libraries if they don't exist
folders=("libcore" "libstdc" "libwx_baseu-3.0" "libwx_baseu_net-3" "libfftw" "libtiff") 
for folder in ${folders[*]}
do
    if [ ! -d $folder ]
    then 
        mkdir $folder
    fi

    nfiles=$(ls $folder | wc -l)
    if [ $nfiles -gt 0 ]
    then
        rm $folder/*
    fi

done

if [ ! -f ./libcore.a ]
then
    echo "cisTEM-1.0 core library not found! Please compile cisTEM-1.0 first"
    exit
fi


# remove previously built object file or statically linked library if they exist
test -f programs/refine3d_cspt/refine3d_cspt-refine3d_cspt.o && rm programs/refine3d_cspt/refine3d_cspt-refine3d_cspt.o
test -f programs/refine3d/cistem.a && rm programs/refine3d/cistem.a

# copy statically linked libraries and convert them into object files 
cd libcore && cp ../libcore.a . && ar x libcore.a && cd ..
cd libstdc && cp /usr/lib/gcc/x86_64-linux-gnu/7/libstdc++.a . && ar x libstdc++.a && cd ..
cd libwx_baseu-3.0 && cp /usr/local/lib/libwx_baseu-3.0.a . && ar x libwx_baseu-3.0.a && cd ..
cd libwx_baseu_net-3 && cp /usr/local/lib/libwx_baseu_net-3.0.a . && ar x libwx_baseu_net-3.0.a && cd ..
cd libtiff && cp /usr/local/lib/libtiff.a . && ar x libtiff.a && cd ..
cd libfftw && cp /usr/lib/x86_64-linux-gnu/libfftw3f.a . && ar x libfftw3f.a && cd ..

# resolve name conflicts for libtiff 
symbols=$(ls .symbols)
for symb in $symbols
do 
    ofile=$(basename $symb .symbols)
    objcopy --redefine-syms=.symbols/${symb} libtiff/${ofile} libtiff/${ofile}.new
    rm libtiff/${ofile} 
    mv libtiff/${ofile}.new libtiff/${ofile}
done

if [ ! -d programs/refine3d_cspt/.deps ]
then 
    mkdir programs/refine3d_cspt/.deps
fi

# compile refine3d source into object file
icpc -DPACKAGE_NAME=\"cisTEM\" -DPACKAGE_TARNAME=\"cistem\" -DPACKAGE_VERSION=\"1.0.0-beta\" -DPACKAGE_STRING=\"cisTEM\ 1.0.0-beta\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE_URL=\"\" -DPACKAGE=\"cistem\" -DVERSION=\"1.0.0-beta\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_DLFCN_H=1 -DLT_OBJDIR=\".libs/\" -Dwx_is_available=1 -I.  -I/usr/local/lib/wx/include/gtk2-unicode-3.0 -I/usr/local/include/wx-3.0 -D_FILE_OFFSET_BITS=64 -DWXUSINGDLL -D__WXGTK__ -DwxUSE_GUI=0 -O3 -no-prec-div -no-prec-sqrt -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DMKL -mkl=sequential -I/usr/local/lib/wx/include/gtk2-unicode-3.0 -I/usr/local/include/wx-3.0 -D_FILE_OFFSET_BITS=64 -DWXUSINGDLL -D__WXGTK__ -DwxUSE_GUI=0 -no-multibyte-chars -fabi-version=11 -wd10237 -O3 -no-prec-div -no-prec-sqrt -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DMKL -MT programs/refine3d_cspt/refine3d_cspt-refine3d_cspt.o -MD -MP -MF programs/refine3d_cspt/.deps/refine3d_cspt-refine3d_cspt.Tpo -c -o programs/refine3d_cspt/refine3d_cspt-refine3d_cspt.o `test -f 'programs/refine3d_cspt/refine3d_cspt.cpp' || echo './'`programs/refine3d_cspt/refine3d_cspt.cpp


# merge all the necessary object files into one statically linked library 
if [ -f programs/refine3d_cspt/refine3d_cspt-refine3d_cspt.o ]
then
    ar rcs cspt_lib.a programs/refine3d_cspt/refine3d_cspt-refine3d_cspt.o libcore/*.o libstdc/*.o libwx_baseu-3.0/*.o libwx_baseu_net-3/*.o libtiff/*.o libfftw/*.o 
    echo "cisTEM-1.0 library (cspt_lib.a) built successfully!"
else
    echo "refine_cspt was not compiled. Please retry."
fi
