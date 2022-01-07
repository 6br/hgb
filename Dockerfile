FROM nginxinc/nginx-unprivileged

ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH \
    RUSTFLAGS="-C target-feature=-avx2" \
    RUSTUP_TOOLCHAIN="nightly"

RUN set -eux; \
    \
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -y; \
    chmod -R a+w $RUSTUP_HOME $CARGO_HOME;

WORKDIR /app

ENV PATH $PATH:/app

RUN cargo build --release; \
    rm -rf src/; \
    cp /app/target/release/hgb /app/hgb

ENTRYPOINT ["/app/hgb"]
