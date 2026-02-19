FROM debian:bookworm-slim AS builder
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
  && rm -rf /var/lib/apt/lists/*

# Build htslib
RUN curl -fsSL -o htslib-${HTSLIB_VERSION}.tar.bz2 https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
  tar xf htslib-${HTSLIB_VERSION}.tar.bz2 && \
  cd htslib-${HTSLIB_VERSION} && \
  ./configure -with-libdeflate && \
  make && \
  make install

# --- Runtime stage ---
FROM debian:bookworm-slim

# Runtime dependencies only
RUN apt-get update && apt-get install -y --no-install-recommends \
  libpcre3 \
  libbz2-1.0 \
  liblzma5 \
  zlib1g \
  libcurl4 \
  libdeflate0 \
  procps \
  && rm -rf /var/lib/apt/lists/*

# Copy htslib shared libraries into /usr/lib (always in the dynamic linker search path)
COPY --from=builder /usr/local/lib/libhts* /usr/lib/
RUN ldconfig

COPY clinvcf /usr/bin/

ENTRYPOINT []
