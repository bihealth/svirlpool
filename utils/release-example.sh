#!/bin/bash
# This is an example script showing how to create a manual release.
# DO NOT RUN THIS SCRIPT DIRECTLY - use it as a reference only.
# Follow RELEASE_CHECKLIST.md for the actual release process.

set -e

# Configuration - UPDATE THESE VALUES
NEW_VERSION="0.2.0"  # Change to your desired version
RELEASE_DATE=$(date +%Y-%m-%d)

echo "==============================================="
echo "Manual Release Script for Svirlpool"
echo "==============================================="
echo "This script demonstrates the release process."
echo "NEW VERSION: ${NEW_VERSION}"
echo "DATE: ${RELEASE_DATE}"
echo "==============================================="
echo ""
echo "WARNING: This is an example/reference script."
echo "Please follow RELEASE_CHECKLIST.md manually instead."
echo ""
read -p "Press Ctrl+C to exit, or Enter to see the steps..." 

echo ""
echo "Step 1: Update version in pyproject.toml"
echo "  sed -i 's/^version = .*/version = \"${NEW_VERSION}\"/' pyproject.toml"

echo ""
echo "Step 2: Update version in src/svirlpool/version.py"
echo "  sed -i 's/__version__ = .*/__version__ = \"${NEW_VERSION}\"/' src/svirlpool/version.py"

echo ""
echo "Step 3: Update version in .release-please-manifest.json"
echo "  sed -i 's/\"\\.\": .*/\".\": \"${NEW_VERSION}\"/' .release-please-manifest.json"

echo ""
echo "Step 4: Update CHANGELOG.md"
echo "  # Edit CHANGELOG.md manually to add new release section"
echo "  # Example format:"
echo "  ## [${NEW_VERSION}](https://github.com/bihealth/svirlpool/compare/v0.1.0...v${NEW_VERSION}) (${RELEASE_DATE})"
echo "  ### Features"
echo "  * List your features here"
echo "  ### Bug Fixes"
echo "  * List your bug fixes here"

echo ""
echo "Step 5: Commit changes"
echo "  git add pyproject.toml src/svirlpool/version.py .release-please-manifest.json CHANGELOG.md"
echo "  git commit -m \"chore: release v${NEW_VERSION}\""

echo ""
echo "Step 6: Push to main"
echo "  git push origin main"

echo ""
echo "Step 7: Create and push tag"
echo "  git tag -a v${NEW_VERSION} -m \"Release v${NEW_VERSION}\""
echo "  git push origin v${NEW_VERSION}"

echo ""
echo "Step 8: Create GitHub release"
echo "  Option A (using gh CLI):"
echo "    gh release create v${NEW_VERSION} --title \"v${NEW_VERSION}\" --notes-file CHANGELOG.md --target main"
echo ""
echo "  Option B (using web UI):"
echo "    1. Go to: https://github.com/bihealth/svirlpool/releases/new"
echo "    2. Select tag: v${NEW_VERSION}"
echo "    3. Set title: v${NEW_VERSION}"
echo "    4. Copy changelog to description"
echo "    5. Click 'Publish release'"

echo ""
echo "Step 9: Verify release"
echo "  # Check GitHub release page"
echo "  # Verify GitHub Actions workflow runs"
echo "  # Test Docker image:"
echo "  docker pull ghcr.io/bihealth/svirlpool:v${NEW_VERSION}"
echo "  docker run --rm ghcr.io/bihealth/svirlpool:v${NEW_VERSION} svirlpool --version"

echo ""
echo "==============================================="
echo "Release process completed!"
echo "==============================================="
echo ""
echo "Remember: This was just a demonstration."
echo "Follow RELEASE_CHECKLIST.md for actual releases."
