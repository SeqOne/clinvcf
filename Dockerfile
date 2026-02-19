FROM debian:bookworm-slim
ARG HTSLIB_VERSION=1.10.2

RUN apt-get update && apt-get install -y --no-install-recommends \
  ca-certificates \
  build-essential \
  curl \
  make \
  libbz2-dev \
  libdeflate-dev \
  liblzma-dev \
  zlib1g-dev \
  bzip2 \
  libcurl4-openssl-dev \
  procps \
  && rm -rf /var/lib/apt/lists/*



# INSTALL HTSLIB
RUN curl -fsSL -o htslib-${HTSLIB_VERSION}.tar.bz2 https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
  tar xf htslib-${HTSLIB_VERSION}.tar.bz2 && \
  cd htslib-${HTSLIB_VERSION} && \
  ./configure --with-libdeflate && \
  make && \
  make install && \
  make clean


ENV LD_LIBRARY_PATH=/usr/bin/htslib-$HTSLIB_VERSION

COPY clinvcf /usr/bin/

ENTRYPOINT []
