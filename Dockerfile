FROM nginxinc/nginx-unprivileged

ENV RUSTUP_TOOLCHAIN="nightly"

RUN set -eux; \
    \
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

WORKDIR /app

ENV PATH $PATH:/app

RUN cargo build --release; \
    rm -rf src/; \
    cp /app/target/release/hgb /app/hgb

ENTRYPOINT ["/app/hgb"]
