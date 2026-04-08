#!/usr/bin/env bash
# AWS Index Builder Setup for Dragon
# ===================================
# Provisions an EC2 instance to build a pre-built Dragon index
# over the full GTDB prokaryotic genome database (2.34M genomes).
#
# Recommended instance: r6i.2xlarge (8 vCPU, 64 GB RAM, ~$0.504/hr)
# With --low-memory mode: r6i.xlarge (4 vCPU, 32 GB RAM, ~$0.252/hr)
#
# Estimated costs:
#   Index build: ~2-4 hours on r6i.2xlarge = $1-2
#   Storage: ~100 GB EBS gp3 for genomes + index = $8/month
#
# Usage:
#   # 1. Launch an EC2 instance (Amazon Linux 2023 or Ubuntu 22.04)
#   # 2. SSH in and run:
#   bash setup_index_builder.sh
#
#   # 3. After setup completes, build the index:
#   bash build_gtdb_index.sh
#
#   # 4. Download the index to your local machine:
#   scp -r ec2-user@<IP>:~/dragon_index/ .

set -euo pipefail

echo "====================================="
echo "Dragon AWS Index Builder Setup"
echo "====================================="

# Detect OS
if [ -f /etc/os-release ]; then
    . /etc/os-release
    OS=$ID
else
    OS="unknown"
fi

echo "Detected OS: $OS"

# Install system dependencies
echo ""
echo "--- Installing system dependencies ---"
if [[ "$OS" == "amzn" || "$OS" == "amazon" ]]; then
    sudo yum update -y
    sudo yum install -y gcc gcc-c++ make cmake git curl wget tar gzip
elif [[ "$OS" == "ubuntu" || "$OS" == "debian" ]]; then
    sudo apt-get update
    sudo apt-get install -y build-essential cmake git curl wget
else
    echo "Unsupported OS: $OS. Install gcc, cmake, git, curl manually."
fi

# Install Rust
echo ""
echo "--- Installing Rust ---"
if ! command -v rustc &>/dev/null; then
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source "$HOME/.cargo/env"
fi
rustc --version
cargo --version

# Clone and build Dragon
echo ""
echo "--- Building Dragon ---"
if [ ! -d "$HOME/Dragon" ]; then
    git clone https://github.com/yourusername/Dragon.git "$HOME/Dragon"
fi
cd "$HOME/Dragon"
cargo build --release
echo "Dragon binary: $(realpath target/release/dragon)"

# Create workspace
mkdir -p "$HOME/dragon_workspace"
mkdir -p "$HOME/dragon_index"

echo ""
echo "====================================="
echo "Setup complete!"
echo ""
echo "Next steps:"
echo "  1. Run: bash ~/Dragon/aws/build_gtdb_index.sh"
echo "  2. Or build manually:"
echo "     ~/Dragon/target/release/dragon download --database gtdb --output ~/dragon_workspace/genomes"
echo "     ~/Dragon/target/release/dragon index -i ~/dragon_workspace/genomes -o ~/dragon_index -k 15 --low-memory --max-ram 56"
echo "====================================="
