# cistem-1.0

This is a fork from version 1.0 of the [cisTEM](https://github.com/timothygrant80/cisTEM) repository.

Main additions:
- New `reconstruct3d_stats` version with support for constrained classification (`cistem1_based` branch)
- New `reconstruct3d_stats_DW` version with support for dose weighting (`main` branch)
- New `refine3d_no_prior` version that skips the use of priors during refinement (`cistem1_based` branch)

## Building from source

1. First, launch cisTEM's singularity container and enable the Intel compiler, for example:
<pre>
singularity shell -B /opt/apps cistem_dev_env_latest.sif
source  /opt/apps/rhel8/intel-2020/compilers_and_libraries/linux/bin/compilervars.sh intel64
</pre>

2. Build wxWidget which is the prerequisite library to build cisTEM software. First, download the source for [wxWidget 3.0.4](https://github.com/wxWidgets/wxWidgets/releases/tag/v3.0.4):
<pre>
wget https://github.com/wxWidgets/wxWidgets/releases/download/v3.0.4/wxWidgets-3.0.4.zip
unzip wxWidgets-3.0.4.zip -d wxWidgets-3

cd wxWidgets-3

./configure CXX=icpc CC=icc --disable-precomp-headers --prefix=$(pwd) --with-libnotify=no --disable-shared --without-gtkprint --with-libjpeg=builtin --with-libpng=builtin --with-libtiff=builtin --with-zlib=builtin --with-expat=builtin --disable-compat28 --without-liblzma --without-libjbig CFLAGS="-no-multibyte-chars" CXXFLAGS="-no-multibyte-chars"

make -j4
make install
</pre>

3. Build cisTEM using the wxWidgets library you just built (specify the path using the option `--with-xw-config`):
<pre>
git clone git@github.com:nextpyp/cistem-1.0.git
cd cisTEM-1.0

./configure CFLAGS="-no-multibyte-chars" CXXFLAGS="-no-multibyte-chars" --enable-staticmode --with-wx-config=/path_to_your_wxwidgets/wx-config --prefix=$(pwd)

make -j4
make install
</pre>
