FROM python:3.12-slim-bookworm AS native
RUN apt-get update && apt-get install -y --no-install-recommends \
    git gcc make libc6-dev zlib1g-dev ca-certificates \
    && rm -rf /var/lib/apt/lists/*
ARG MINIMAP2_COMMIT=3c28777e7e2dcc90f825de1b9f17a89cca7d4452
RUN git clone https://github.com/lh3/minimap2.git /build/minimap2 \
    && cd /build/minimap2 && git checkout "$MINIMAP2_COMMIT" \
    && test "$(git rev-parse HEAD)" = "$MINIMAP2_COMMIT" \
    && make -j2

FROM python:3.12-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates libgomp1 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=native /build/minimap2/minimap2 /usr/local/bin/minimap2
WORKDIR /opt/genocache
COPY requirements-runtime.txt pyproject.toml ./
COPY genocache ./genocache
RUN pip install --no-cache-dir -r requirements-runtime.txt \
    && pip install --no-cache-dir --no-deps .
RUN useradd --uid 10001 --create-home worker && mkdir -p /work && chown worker:worker /work
USER worker
WORKDIR /work
ENV PYTHONUNBUFFERED=1 GENOCACHE_ROOT=/work/cache
CMD ["genocache", "--help"]
