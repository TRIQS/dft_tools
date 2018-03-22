# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:master-ubuntu-clang

COPY . ${SRC}/app4triqs
WORKDIR ${BUILD}/app4triqs
RUN chown build .
USER build
RUN cmake ${SRC}/app4triqs -DTRIQS_ROOT=${INSTALL} && make -j2 && make test
USER root
RUN make install
