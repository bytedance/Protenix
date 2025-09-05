# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Test cases for installation and compatibility issues by Shad Nygren (@ShadNygren)
  - `test_installation.py` - Tests for DeepSpeed/Pydantic compatibility (issue #182)
  - `test_triton_compatibility.py` - Tests for Triton GPU compatibility (issue #185)
  - `test_esm_loading.py` - Tests for ESM weights loading with PyTorch 2.6+ (issue #176)
- CHANGELOG.md to track project changes (Apache 2.0 compliance)

### Fixed
- Issue #182: DeepSpeed and Pydantic v2.x compatibility by upgrading to DeepSpeed 0.17.5
- Issue #185: Triton kernel compatibility with consumer GPUs (RTX 3090/4090)
  - Added PyTorch fallback when Triton is unavailable or unsupported
  - Enables Protenix to run on consumer hardware without code changes
  - Modified `protenix/model/tri_attention/__init__.py` to provide transparent fallback

### Contributors
- Shad Nygren, Virtual Hipster Corporation (@ShadNygren)

## [v0.6.1] - 2025-08-20

### Fixed
- ESM models loading in PyTorch 2.6 and newer versions (commit dbbee14)

## [v0.6.0] - Previous releases

See commit history for earlier changes.