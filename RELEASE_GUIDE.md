# Release Guide for Svirlpool

This guide explains how to create a manual release for the svirlpool project.

## Overview

The svirlpool project uses [Release Please](https://github.com/googleapis/release-please) for automated releases. However, if you need to create a manual release, follow the steps below.

## Current State

- **Current Version**: 0.1.0
- **Last Release**: v0.0.0 (published January 7, 2026)
- **Release Automation**: GitHub Actions workflow `.github/workflows/release-please.yml`
- **Version Files**:
  - `pyproject.toml` (line 9)
  - `src/svirlpool/version.py` (line 1)
  - `.release-please-manifest.json`

## Prerequisites

Before creating a release, ensure:

1. All changes are merged to the `main` branch
2. All tests pass
3. The changelog is up to date
4. You have write access to the repository

## Step-by-Step Manual Release Process

### Step 1: Update Version Numbers

You need to update the version in three files. Choose your new version number (e.g., `0.2.0` for a minor release, `0.1.1` for a patch, or `1.0.0` for a major release).

**1.1. Update `pyproject.toml`:**
```bash
# Edit line 9 in pyproject.toml
version = "0.2.0"  # Change to your new version
```

**1.2. Update `src/svirlpool/version.py`:**
```bash
# Edit line 1 in src/svirlpool/version.py
__version__ = "0.2.0"  # Change to your new version
```

**1.3. Update `.release-please-manifest.json`:**
```bash
# Edit the version in .release-please-manifest.json
{
  ".": "0.2.0"  # Change to your new version
}
```

### Step 2: Update the Changelog

Edit `CHANGELOG.md` to add a new section for your release:

```markdown
## [0.2.0](https://github.com/bihealth/svirlpool/compare/v0.1.0...v0.2.0) (2026-02-03)

### Features

* List your new features here

### Bug Fixes

* List your bug fixes here

### Documentation

* List documentation changes here
```

### Step 3: Commit and Push Changes

```bash
git add pyproject.toml src/svirlpool/version.py .release-please-manifest.json CHANGELOG.md
git commit -m "chore: release v0.2.0"
git push origin main
```

### Step 4: Create a Git Tag

```bash
# Create an annotated tag
git tag -a v0.2.0 -m "Release v0.2.0"

# Push the tag to GitHub
git push origin v0.2.0
```

### Step 5: Create a GitHub Release

You have two options for creating the GitHub release:

#### Option A: Using GitHub CLI (gh)

```bash
gh release create v0.2.0 \
  --title "v0.2.0" \
  --notes-file CHANGELOG.md \
  --target main
```

#### Option B: Using GitHub Web Interface

1. Go to https://github.com/bihealth/svirlpool/releases/new
2. Select the tag: `v0.2.0`
3. Set the release title: `v0.2.0`
4. Copy the relevant section from `CHANGELOG.md` into the description
5. Click "Publish release"

### Step 6: Verify the Release

After creating the release:

1. **Check the GitHub Release page**: https://github.com/bihealth/svirlpool/releases
2. **Verify Docker Image Build**: The `release-please.yml` workflow should trigger automatically and build the Docker image
   - Go to: https://github.com/bihealth/svirlpool/actions
   - Check that the "release-please" workflow runs successfully
   - Verify the Docker image is available at: `ghcr.io/bihealth/svirlpool:v0.2.0`

3. **Test the Docker image**:
```bash
docker pull ghcr.io/bihealth/svirlpool:v0.2.0
docker run --rm ghcr.io/bihealth/svirlpool:v0.2.0 svirlpool --version
```

## Semantic Versioning Guide

Svirlpool follows [Semantic Versioning](https://semver.org/):

- **MAJOR version** (X.0.0): Incompatible API changes
- **MINOR version** (0.X.0): New functionality in a backward-compatible manner
- **PATCH version** (0.0.X): Backward-compatible bug fixes

Current version: `0.1.0` (pre-1.0 development phase)

### During Pre-1.0 Development:
- Breaking changes can bump the MINOR version
- New features bump the MINOR version
- Bug fixes bump the PATCH version

## Automated Release Process (Normal Flow)

By default, releases are automated via Release Please:

1. Push commits to `main` with conventional commit messages:
   - `feat:` for new features (bumps minor version)
   - `fix:` for bug fixes (bumps patch version)
   - `feat!:` or `BREAKING CHANGE:` for breaking changes (bumps major version)

2. Release Please creates/updates a release PR automatically

3. When you merge the release PR:
   - Version numbers are updated
   - Changelog is generated
   - Git tag is created
   - GitHub release is published
   - Docker image is built and pushed

## Troubleshooting

### Issue: Tag already exists

```bash
# Delete the local tag
git tag -d v0.2.0

# Delete the remote tag (careful!)
git push origin :refs/tags/v0.2.0

# Recreate and push
git tag -a v0.2.0 -m "Release v0.2.0"
git push origin v0.2.0
```

### Issue: Docker build fails

Check the GitHub Actions logs:
1. Go to: https://github.com/bihealth/svirlpool/actions
2. Find the failed workflow run
3. Review the logs for errors
4. Fix the issues and re-run the workflow if possible

### Issue: Version mismatch

Ensure all three version files are updated consistently:
- `pyproject.toml`
- `src/svirlpool/version.py`
- `.release-please-manifest.json`

## Additional Resources

- [Release Please Documentation](https://github.com/googleapis/release-please)
- [Conventional Commits](https://www.conventionalcommits.org/)
- [Semantic Versioning](https://semver.org/)
- [GitHub Releases Documentation](https://docs.github.com/en/repositories/releasing-projects-on-github)
- [Example Release Script](utils/release-example.sh) - Reference only, do not run directly

## Questions or Issues?

If you encounter any issues with the release process, please:
1. Check the GitHub Actions logs
2. Review this guide
3. Contact the maintainers via GitHub issues
