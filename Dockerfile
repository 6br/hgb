FROM rust:1.56.0

WORKDIR /app

COPY . .

RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/app/target \
    cargo build --release

COPY /app/target/release/hgb /app/hgb

ENTRYPOINT ["/app/hgb"]