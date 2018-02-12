# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:master-ubuntu-clang

COPY . ${SRC}/dft_tools
WORKDIR ${BUILD}/dft_tools
RUN chown build .
USER build
RUN cmake ${SRC}/dft_tools -DTRIQS_ROOT=${INSTALL} && make -j2 && make test
USER root
RUN make install
