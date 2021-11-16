FROM rustlang/rust:nightly

WORKDIR /app

COPY . .

RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/app/target \
    cargo build --release && cp /app/target/release/hgb /app/hgb

ENTRYPOINT ["/app/hgb"]