# --- Stage 1: Build ---
FROM ghcr.io/prefix-dev/pixi:latest AS build

# Set the working directory
WORKDIR /app

# Install dependencies
COPY pyproject.toml pixi.lock README.md LICENSE src ./
RUN pixi install --locked

# create the shell-hook bash script to activate the environment
RUN pixi shell-hook -e default -s bash > /shell-hook
RUN echo "#!/bin/bash" > /app/entrypoint.sh
RUN cat /shell-hook >> /app/entrypoint.sh
# extend the shell-hook script to run the command passed to the container
RUN echo 'exec "$@"' >> /app/entrypoint.sh


# --- Stage 2: Runtime ---
FROM debian:bookworm-slim AS runtime

WORKDIR /app
# only copy the default environment into container
# please note that the "prefix" (path) needs to stay
# the same as in the build container
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default
COPY --from=build --chmod=0755 /app/entrypoint.sh /app/entrypoint.sh

# Extend PYTHONPATH appropriately
ENV PYTHONPATH="/app/src"
# copy in the application code
COPY ./src /app/src

# Set entrypoint.
ENTRYPOINT ["/app/entrypoint.sh"]
