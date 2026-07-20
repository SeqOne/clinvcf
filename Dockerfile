# Base images are pulled through the GitLab dependency proxy in CI.
# CI_DEPENDENCY_PROXY_GROUP_IMAGE_PREFIX is injected by CI; locally it defaults
# to Docker Hub so `docker build .` works without extra flags.
ARG CI_DEPENDENCY_PROXY_GROUP_IMAGE_PREFIX=docker.io/library

# --- Build stage: toolchain + compile clinvcf ---
FROM ${CI_DEPENDENCY_PROXY_GROUP_IMAGE_PREFIX}/debian:bookworm-slim AS builder
ARG NIM_VERSION=1.6.14
ARG HTSLIB_VERSION=1.10.2

RUN apt-get update && apt-get install -y --no-install-recommends \
  ca-certificates \
  build-essential \
  curl \
  make \
  git \
  bzip2 \
  libbz2-dev \
  libdeflate-dev \
  liblzma-dev \
  zlib1g-dev \
  libcurl4-openssl-dev \
  libpcre3-dev \
  && rm -rf /var/lib/apt/lists/*

# Build htslib
RUN curl -fsSL -o htslib-${HTSLIB_VERSION}.tar.bz2 https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
  tar xf htslib-${HTSLIB_VERSION}.tar.bz2 && \
  cd htslib-${HTSLIB_VERSION} && \
  ./configure --with-libdeflate && \
  make && \
  make install && \
  make clean && \
  ldconfig

# Install Nim (version-pinned) via choosenim
ENV CHOOSENIM_CHOOSE_VERSION=${NIM_VERSION}
RUN curl https://nim-lang.org/choosenim/init.sh -sSf | sh -s -- -y
ENV PATH="/root/.nimble/bin:$PATH"

# Compile clinvcf
WORKDIR /src
COPY . .
RUN nimble build -d:release -y

# --- Runtime stage ---
FROM ${CI_DEPENDENCY_PROXY_GROUP_IMAGE_PREFIX}/debian:bookworm-slim

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

COPY --from=builder /src/clinvcf /usr/bin/

ENTRYPOINT []