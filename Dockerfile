# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang
ARG APPNAME=triqs_dft_tools

RUN apt-get install -y meson ninja-build python3-setuptools gfortran

# triqs_dft_tools imports triqs_dftkit at runtime (its converters re-export
# triqs_dftkit.{wien2k,vasp,...}). Build and install it into the TRIQS prefix
# before building this app so the import resolves.
ARG DFTKIT_BRANCH=unstable
RUN git clone --depth 1 -b $DFTKIT_BRANCH https://github.com/TRIQS/dftkit $SRC/dftkit \
 && mkdir -p $BUILD/dftkit && cd $BUILD/dftkit \
 && cmake $SRC/dftkit -DTRIQS_ROOT=${INSTALL} -DBuild_Tests=OFF \
 && (make -j4 || make -j1 VERBOSE=1) \
 && make install

COPY --chown=build . $SRC/$APPNAME
RUN mkdir $BUILD/$APPNAME && chown build $BUILD/$APPNAME

ARG BUILD_ID
ARG CMAKE_ARGS
USER build
WORKDIR $BUILD/$APPNAME
RUN cmake $SRC/$APPNAME -DTRIQS_ROOT=${INSTALL} $CMAKE_ARGS && make -j4 || make -j1 VERBOSE=1
USER root
RUN make install
