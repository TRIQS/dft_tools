# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang
ARG APPNAME=triqs_dft_tools

RUN apt-get install -y meson ninja-build python3-setuptools

COPY --chown=build . $SRC/$APPNAME
RUN mkdir $BUILD/$APPNAME && chown build $BUILD/$APPNAME

ARG BUILD_ID
ARG CMAKE_ARGS
USER build
WORKDIR $BUILD/$APPNAME
RUN cmake $SRC/$APPNAME -DTRIQS_ROOT=${INSTALL} $CMAKE_ARGS && make -j4 || make -j1 VERBOSE=1
USER root
RUN make install
