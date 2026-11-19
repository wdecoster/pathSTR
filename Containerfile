FROM python:3.14-slim AS builder

COPY --from=ghcr.io/astral-sh/uv:0.12.12 /uv /bin/uv

RUN apt-get update \
 && apt-get install -y --no-install-recommends g++ \
 && rm -rf /var/lib/apt/lists/*

# build against the base image's interpreter
# so the venv it links to is the one the runtime stage has
ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PYTHON=/usr/local/bin/python3 \
    UV_PYTHON_DOWNLOADS=never

WORKDIR /app

COPY pyproject.toml uv.lock ./
RUN uv sync --frozen --no-dev


FROM python:3.14-slim
LABEL org.opencontainers.image.source="https://github.com/wdecoster/pathSTR"

# the base image lags behind Debian's security updates
RUN apt-get update \
 && apt-get upgrade -y \
 && rm -rf /var/lib/apt/lists/*

# nothing installs packages at runtime so remove it
RUN rm -rf /usr/local/lib/python3.*/site-packages/pip* /usr/local/bin/pip*

ENV PYTHONUNBUFFERED=1 \
    PATH=/app/.venv/bin:$PATH

WORKDIR /app

# same path as the builder, so the venv's shebangs still resolve
COPY --from=builder /app/.venv /app/.venv
COPY app.py .
COPY pathSTR/ ./pathSTR/
COPY assets/ ./assets/

RUN useradd --create-home --uid 10001 app
USER app

EXPOSE 8080

CMD ["gunicorn", "--bind", "0.0.0.0:8080", "--workers", "2", "--access-logfile", "-", "app:create_wsgi_app()"]
