# Use a slim Python base image
FROM python:3.10-slim-bookworm

LABEL org.opencontainers.image.source="https://github.com/Singular-Genomics/G4X-helpers"

# Add uv to the image
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

# Install build dependencies
RUN apt-get update -qq \
 && apt-get install -y --no-install-recommends \
      libexpat1 \
      libopenjp2-7 \
      gcc \
      g++ \
      libgdal-dev gdal-bin \
 && rm -rf /var/lib/apt/lists/*

# Install the project into `/app`
WORKDIR /app

# Enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1

# Copy from the cache instead of linking since it's a mounted volume
ENV UV_LINK_MODE=copy

RUN mkdir -p /usr/local/tmp && chmod a+rwx /usr/local/tmp
RUN mkdir -p /usr/local/tmp/fontconfig && chmod a+rwx /usr/local/tmp/fontconfig

ENV MPLCONFIGDIR=/usr/local/tmp
ENV NUMBA_CACHE_DIR=/usr/local/tmp
ENV XDG_CACHE_HOME=/usr/local/tmp

# Install the project's dependencies using the lockfile and settings
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project --no-dev

    # Then, add the rest of the project source code and install it
# Installing separately from its dependencies allows optimal layer caching
COPY . /app
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev

# Place executables in the environment at the front of the path
ENV PATH="/app/.venv/bin:$PATH"

# Reset the entrypoint, don't invoke `uv`
ENTRYPOINT []

CMD ["bash"]
