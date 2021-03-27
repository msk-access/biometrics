################## BASE IMAGE ######################

FROM python:3.6-slim

################## ARGUMENTS/Environments ##########

ARG BUILD_DATE
ARG BUILD_VERSION
ARG LICENSE="Apache-2.0"
ARG BIOMETRICS_VERSION
ARG VCS_REF
################## METADATA ########################
LABEL org.opencontainers.image.vendor="MSKCC"
LABEL org.opencontainers.image.authors="Charlie Murphy (murphyc4@mskcc.org)"

LABEL org.opencontainers.image.created=${BUILD_DATE} \
    org.opencontainers.image.version=${BUILD_VERSION} \
    org.opencontainers.image.licenses=${LICENSE} \
    org.opencontainers.image.version.biometrics=${BIOMETRICS_VERSION} \
    org.opencontainers.image.source.biometrics="https://pypi.org/project/biometrics/" \
    org.opencontainers.image.vcs-url="https://github.com/msk-access/biometrics.git" \
    org.opencontainers.image.vcs-ref=${VCS_REF}

LABEL org.opencontainers.image.description="This container uses python3.6 as the base image to build \
    biometrics version ${BIOMETRICS_VERSION}"

################## INSTALL ##########################

WORKDIR /app
ADD . /app

RUN apt-get update \
  && apt-get install gcc g++ zlib1g-dev -y \
  && pip install -r requirements.txt \
  && pip install .
