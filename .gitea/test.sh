#!/usr/bin/env bash
# test one image reference. Usage: .gitea/test.sh <image-ref>
set -euo pipefail

ref="${1:?usage: $0 <image-ref>}"
name="test-$$"
# curl runs in its own container inside the app container's network
curl_image="${CURL_IMAGE:-curlimages/curl:latest}"

docker pull -q "$curl_image" >/dev/null

trap 'docker rm -f "$name" >/dev/null 2>&1 || true' EXIT

[ "$(docker inspect -f '{{.Config.User}}' "$ref")" = "app" ] \
  || { echo "image does not default to the app user" >&2; exit 1; }

in_app_net() { docker run --rm --network "container:$name" "$curl_image" "$@"; }

# only checks for a 2xx. Grep the body for something your app renders
# to also catch error pages served with a 200.
wait_for() {
  for _ in $(seq 1 30); do
    in_app_net -fsS -o /dev/null "$1" 2>/dev/null && return 0
    # in_app_net -fsS "$1" 2>/dev/null | grep -qF "My Page Title" && return 0
    sleep 2
  done
  docker logs "$name" || true
  echo "app did not answer at $1" >&2
  return 1
}

docker run -d --name "$name" "$ref" >/dev/null
wait_for http://127.0.0.1:8080/
docker rm -f "$name" >/dev/null

# APP_BASE_PATH should move where the app answers, not break it
docker run -d --name "$name" -e APP_BASE_PATH=/sub "$ref" >/dev/null
wait_for http://127.0.0.1:8080/sub/
docker rm -f "$name" >/dev/null

echo "test passed"
