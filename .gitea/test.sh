#!/usr/bin/env bash
# test one image reference. Usage: .gitea/test.sh <image-ref>
set -euo pipefail

ref="${1:?usage: $0 <image-ref>}"
name="test-$$"
# the image ships without data: a small database is built from the repo's test VCFs
# into a volume (not a bind mount, as the runner's paths may not be visible to the daemon)
volume="pathstr-test-db-$$"
builder="pathstr-test-build-$$"
# curl runs in its own container inside the app container's network
curl_image="${CURL_IMAGE:-curlimages/curl:latest}"

docker pull -q "$curl_image" >/dev/null

cleanup() {
  docker rm -f "$name" "$builder" >/dev/null 2>&1 || true
  docker volume rm -f "$volume" >/dev/null 2>&1 || true
}
trap cleanup EXIT

[ "$(docker inspect -f '{{.Config.User}}' "$ref")" = "app" ] \
  || { echo "image does not default to the app user" >&2; exit 1; }

in_app_net() { docker run --rm --network "container:$name" "$curl_image" "$@"; }

# build the test database (needs network access to fetch the STRchive catalog).
# runs as root only because the fresh volume is root-owned; the app itself runs as "app"
docker volume create "$volume" >/dev/null
docker create --name "$builder" --user 0 -v "$volume:/data" -w /data/build "$ref" \
  python /app/app.py --vcf test/genotypes_test.fofn --sample_info data/pathSTR_samples.tsv \
  --save_db /data/pathSTR_database --store_only >/dev/null
docker cp test "$builder:/data/build/test"
docker cp data "$builder:/data/build/data"
docker start -a "$builder" >/dev/null 2>&1 \
  || { docker logs "$builder" | tail -50; echo "building the test database failed" >&2; exit 1; }
docker rm -f "$builder" >/dev/null

run_app() { docker run -d --name "$name" -v "$volume:/data:ro" "$@" "$ref" >/dev/null; }

# grep the page title too, so an error page served with a 200 does not pass
wait_for() {
  for _ in $(seq 1 30); do
    in_app_net -fsS "$1" 2>/dev/null | grep -qF "<title>pathSTR</title>" && return 0
    sleep 2
  done
  docker logs "$name" || true
  echo "app did not answer at $1" >&2
  return 1
}

run_app
wait_for http://127.0.0.1:8080/
docker rm -f "$name" >/dev/null

# APP_BASE_PATH should move where the app answers, not break it
run_app -e APP_BASE_PATH=/sub
wait_for http://127.0.0.1:8080/sub/
# Dash serves its index page for any unknown path, so a 200 proves nothing here:
# check that the API and assets answer with their own content types under the prefix
content_type() { in_app_net -fsS -o /dev/null -w '%{content_type}' "$1"; }
case "$(content_type http://127.0.0.1:8080/sub/_dash-layout)" in
  application/json*) ;;
  *) echo "dash routes not served under /sub/" >&2; exit 1 ;;
esac
case "$(content_type http://127.0.0.1:8080/sub/assets/custom.css)" in
  text/css*) ;;
  *) echo "assets not served under /sub/" >&2; exit 1 ;;
esac
docker rm -f "$name" >/dev/null

echo "test passed"
