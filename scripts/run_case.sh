#!/usr/bin/env bash
set -euo pipefail

CONFIG=${1:-configs/default.ini}
OPENLB_DIR=${OPENLB_DIR:-}
CASE_NAME=droplet_particle_jump3d
BUILD_DIR=build/${CASE_NAME}

if [[ -z "${OPENLB_DIR}" ]]; then
  echo "OPENLB_DIR must point to an OpenLB source tree" >&2
  exit 2
fi
if [[ ! -f "${OPENLB_DIR}/src/olb3D.h" && ! -f "${OPENLB_DIR}/olb3D.h" ]]; then
  echo "Could not find OpenLB headers under OPENLB_DIR=${OPENLB_DIR}" >&2
  exit 2
fi

mkdir -p "${BUILD_DIR}"
CXX=${CXX:-g++}
CXXFLAGS=${CXXFLAGS:-"-O3 -std=c++14"}
INCLUDES="-I${OPENLB_DIR}/src -I${OPENLB_DIR}"
LIBS=${OPENLB_LIBS:-""}

${CXX} ${CXXFLAGS} ${INCLUDES} src/${CASE_NAME}.cpp ${LIBS} -o "${BUILD_DIR}/${CASE_NAME}"
"${BUILD_DIR}/${CASE_NAME}" "${CONFIG}"
