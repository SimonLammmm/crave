#!/usr/bin/env bash
# Build, save and publish the CRAVE image.
#
# amd64 only, unlike Exorcise. Exorcise needs a multi-architecture manifest
# because BLAT and twoBitToFa are architecture-specific binaries and its users
# run it on Apple silicon as well as on servers. CRAVE is an R Shiny app served
# over HTTP: it is deployed to a server, not run on a laptop, and none of its
# dependencies are architecture-sensitive. One image keeps the build to minutes
# rather than the hours QEMU emulation would cost.
#
# If that ever stops being true, the multi-architecture recipe is in Exorcise's
# docker/build-and-push.sh.
#
# Run from the repository root:
#
#   ./build-and-push.sh build     # build and tag locally
#   ./build-and-push.sh push      # push :VERSION and :latest
#   ./build-and-push.sh release   # build then push
#   ./build-and-push.sh save      # write a tarball next to the source
#   ./build-and-push.sh inspect   # show what is published

set -euo pipefail

REPO="${CRAVE_REPO:-simonlammmm/crave}"
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PLATFORM="linux/amd64"

# The version lives in the app, not here, so the image tag cannot disagree with
# what the running app reports on its home page.
CONSTANTS="${ROOT}/shiny-server/R/01_constants.R"
VERSION="$(sed -n 's/^CRAVE_VERSION <- "\(.*\)"$/\1/p' "${CONSTANTS}" | head -n 1)"
if [[ -z "${VERSION}" ]]; then
  echo "ERROR: could not read CRAVE_VERSION from ${CONSTANTS}" >&2
  exit 1
fi

IMAGE_DIR="${CRAVE_IMAGE_DIR:-${ROOT}/images}"

usage() {
  cat <<EOF
Usage: ./build-and-push.sh COMMAND

Commands:
  build     Build ${REPO}:${VERSION} and ${REPO}:latest for ${PLATFORM}
  push      Push both tags
  release   build, then push
  save      Save the built image to ${IMAGE_DIR}/crave-${VERSION}_amd64.tar.gz
  inspect   Show the published manifest
  version   Print the version and exit

Version ${VERSION}, read from shiny-server/R/01_constants.R.
Override the repository with CRAVE_REPO.
EOF
}

# --load puts the result in the local image store, which --platform alone does
# not do. Possible here only because this is a single-platform build.
cmd_build() {
  cd "${ROOT}"
  echo "Building ${REPO}:${VERSION} for ${PLATFORM}..."
  docker buildx build \
    --platform "${PLATFORM}" \
    -t "${REPO}:${VERSION}" \
    -t "${REPO}:latest" \
    --load \
    .
  echo "Built ${REPO}:${VERSION}"
}

cmd_push() {
  for tag in "${VERSION}" latest; do
    if ! docker image inspect "${REPO}:${tag}" >/dev/null 2>&1; then
      echo "ERROR: ${REPO}:${tag} is not built. Run './build-and-push.sh build' first." >&2
      exit 1
    fi
  done
  docker push "${REPO}:${VERSION}"
  docker push "${REPO}:latest"
  echo "Pushed ${REPO}:${VERSION} and ${REPO}:latest"
}

cmd_save() {
  mkdir -p "${IMAGE_DIR}"
  local out="${IMAGE_DIR}/crave-${VERSION}_amd64.tar.gz"
  if ! docker image inspect "${REPO}:${VERSION}" >/dev/null 2>&1; then
    echo "ERROR: ${REPO}:${VERSION} is not built. Run './build-and-push.sh build' first." >&2
    exit 1
  fi
  docker save "${REPO}:${VERSION}" | gzip > "${out}"
  echo "Saved ${out}"
}

cmd_inspect() {
  docker buildx imagetools inspect "${REPO}:${VERSION}"
}

case "${1:-}" in
  build)   cmd_build ;;
  push)    cmd_push ;;
  release) cmd_build; cmd_push; cmd_inspect ;;
  save)    cmd_save ;;
  inspect) cmd_inspect ;;
  version) echo "${VERSION}" ;;
  -h | --help | help | "") usage ;;
  *) echo "Unknown command: $1" >&2; echo >&2; usage >&2; exit 2 ;;
esac
