# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:master-ubuntu-clang

ARG APPNAME
COPY . $SRC/$APPNAME
WORKDIR $BUILD/$APPNAME
RUN chown build .
USER build
ARG BUILD_DOC=0
RUN cmake $SRC/$APPNAME -DTRIQS_ROOT=${INSTALL} -DBuild_Documentation=${BUILD_DOC} && make -j2 && make test
USER root
RUN make install
