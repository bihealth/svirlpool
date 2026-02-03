# Quick Release Checklist

Use this checklist when creating a manual release for svirlpool.

## Pre-Release Checklist

- [ ] All changes merged to `main`
- [ ] All tests passing
- [ ] Determine version number (current: 0.1.0)
  - Patch (0.1.X): Bug fixes only
  - Minor (0.X.0): New features, backward-compatible
  - Major (X.0.0): Breaking changes

## Release Steps

### 1. Update Version Files

- [ ] Update `pyproject.toml` line 9: `version = "X.Y.Z"`
- [ ] Update `src/svirlpool/version.py` line 1: `__version__ = "X.Y.Z"`
- [ ] Update `.release-please-manifest.json`: `".": "X.Y.Z"`

### 2. Update Changelog

- [ ] Edit `CHANGELOG.md` with new release section
- [ ] Include features, bug fixes, and documentation changes

### 3. Commit and Tag

```bash
git add pyproject.toml src/svirlpool/version.py .release-please-manifest.json CHANGELOG.md
git commit -m "chore: release vX.Y.Z"
git push origin main
git tag -a vX.Y.Z -m "Release vX.Y.Z"
git push origin vX.Y.Z
```

### 4. Create GitHub Release

Option A - Using GitHub CLI:
```bash
gh release create vX.Y.Z --title "vX.Y.Z" --notes-file CHANGELOG.md --target main
```

Option B - Web UI:
- [ ] Go to: https://github.com/bihealth/svirlpool/releases/new
- [ ] Select tag: vX.Y.Z
- [ ] Title: vX.Y.Z
- [ ] Copy changelog content to description
- [ ] Click "Publish release"

### 5. Verify Release

- [ ] Check release page: https://github.com/bihealth/svirlpool/releases
- [ ] Verify GitHub Actions workflow runs successfully
- [ ] Confirm Docker image built: `ghcr.io/bihealth/svirlpool:vX.Y.Z`
- [ ] Test Docker image:
  ```bash
  docker pull ghcr.io/bihealth/svirlpool:vX.Y.Z
  docker run --rm ghcr.io/bihealth/svirlpool:vX.Y.Z svirlpool --version
  ```

## Notes

- See `RELEASE_GUIDE.md` for detailed instructions
- Normal releases should use the automated Release Please workflow
- Only create manual releases when necessary
