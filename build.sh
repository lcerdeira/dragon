#!/bin/bash
cargo build --release 2>&1
echo "EXIT_CODE=$?"
