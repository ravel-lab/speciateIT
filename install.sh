#!/usr/bin/env bash

# Requirements
# curl jq unzip bash

set -e  # stop if error

CWD=$(pwd)
DB_PATH="$CWD/vSpeciateDB_models"
ARTICLE_ID="25254229"

mkdir -p "$DB_PATH"
cd "$DB_PATH"

echo "Fetching file list from Figshare..."

# Get files list and download each
curl -s "https://api.figshare.com/v2/articles/$ARTICLE_ID" \
| jq -r '.files[] | "\(.name) \(.download_url)"' \
| while read -r name url; do
    echo "Downloading $name ..."
    curl -sS -L "$url" -o "$name"

    # unzip only if zip file
    if [[ "$name" == *.zip ]]; then
        echo "Unzipping $name ..."
        unzip -q -o "$name"
        rm $name
    fi
done

# Remove __MACOSX dir
rm -rf "$DB_PATH/__MACOSX"

echo "Download + unzip complete."

# Detect OS
OS=""
uname_out=$(uname)
if [ "$uname_out" == "Darwin" ]; then
    OS="macosx"

elif [ "$uname_out" == "Linux" ]; then
    OS="linux"

else
    echo "Unsupported OS: $uname_out"
    exit 1
fi
echo "OS: $OS"

# Add classify to PATH
CLASSIFY_PATH="$CWD/bin/$OS"
export PATH="$CLASSIFY_PATH:$PATH"

RC_FILE="$HOME/.bashrc"
if [ -n "$ZSH_VERSION" ]; then
    RC_FILE="$HOME/.zshrc"
fi

if ! grep -Fxq "export PATH=\"$CLASSIFY_PATH:\$PATH\"" "$RC_FILE"; then
    echo "" >> "$RC_FILE"
    echo "# Added by vSpeciate installer" >> "$RC_FILE"
    echo "export PATH=\"$CLASSIFY_PATH:\$PATH\"" >> "$RC_FILE"
    echo "$CLASSIFY_PATH added to $RC_FILE"
fi

echo "Installation complete. You can now run a test with:"
echo "  classify -d vSpeciateDB_models/vSpeciateIT_V3V4 -i test.fasta -o test"
