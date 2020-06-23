# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang
ARG APPNAME=triqs_dft_tools

COPY requirements.txt /src/$APPNAME/requirements.txt
RUN pip3 install -r /src/$APPNAME/requirements.txt

COPY --chown=build . $SRC/$APPNAME
WORKDIR $BUILD/$APPNAME
RUN chown build .
USER build
ARG BUILD_DOC=0
ARG BUILD_ID
RUN cmake $SRC/$APPNAME -DTRIQS_ROOT=${INSTALL} -DBuild_Documentation=${BUILD_DOC} -DBuild_Deps=Always && make -j2 || make -j1 VERBOSE=1
USER root
RUN make install
