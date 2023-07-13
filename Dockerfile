FROM debian:stable-slim as builder
ARG HTSLIB_VERSION=1.10

RUN apt-get update && apt-get install --no-install-recommends -y \
    libpcre3 libpcre3-dev \
    make \
    wget \
    libbz2-dev \
    bzip2 \
    ca-certificates \
    liblzma-dev \
    zlib1g-dev  libcurl4-gnutls-dev gcc

# Install HTSLIB
RUN cd /usr/bin && \
    wget https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2 && \
    tar -vxjf htslib-$HTSLIB_VERSION.tar.bz2 && \
    cd htslib-$HTSLIB_VERSION && \
    make

# FROM debian:stable-slim
# ARG HTSLIB_VERSION=1.10

# COPY --from=builder /usr/bin/htslib-$HTSLIB_VERSION /usr/bin/htslib-$HTSLIB_VERSION

ENV LD_LIBRARY_PATH=/usr/bin/htslib-$HTSLIB_VERSION


ARG PRIVATE_CI_ACCESS_TOKEN
ARG VERSION

COPY clinvcf /usr/bin/

ENTRYPOINT ["clinvcf"]
