FROM rust:1.56.0

WORKDIR /app

COPY . .

ENV RUSTFLAGS -C target-feature=+crt-static

RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/app/target \
    cargo build --release && cp /app/target/release/hgb /app/hgb

ENTRYPOINT ["/app/hgb"]