"""Tests for TemplateHitProcessor._fetch_or_read_cif with remote fetching."""

import os
import tempfile

import pytest
import requests

from protenix.data.template.template_utils import TemplateHitProcessor


@pytest.fixture
def tmp_mmcif_dir(tmp_path):
    """Create a temporary mmcif directory."""
    d = tmp_path / "mmcif"
    d.mkdir()
    return str(d)


class TestFetchOrReadCifLocal:
    """Tests for local file reading (no network)."""

    def test_reads_existing_local_file(self, tmp_mmcif_dir):
        cif_path = os.path.join(tmp_mmcif_dir, "1abc.cif")
        with open(cif_path, "w") as f:
            f.write("LOCAL_CIF_CONTENT")

        proc = TemplateHitProcessor(mmcif_dir=tmp_mmcif_dir, fetch_remote=False)
        result = proc._fetch_or_read_cif("1abc")
        assert result == "LOCAL_CIF_CONTENT"

    def test_local_preferred_over_remote(self, tmp_mmcif_dir):
        """Even with fetch_remote=True, local file should be used if it exists."""
        cif_path = os.path.join(tmp_mmcif_dir, "1abc.cif")
        with open(cif_path, "w") as f:
            f.write("LOCAL_CIF_CONTENT")

        proc = TemplateHitProcessor(mmcif_dir=tmp_mmcif_dir, fetch_remote=True)
        result = proc._fetch_or_read_cif("1abc")
        assert result == "LOCAL_CIF_CONTENT"

    def test_raises_when_missing_and_fetch_remote_false(self, tmp_mmcif_dir):
        proc = TemplateHitProcessor(mmcif_dir=tmp_mmcif_dir, fetch_remote=False)
        with pytest.raises(FileNotFoundError, match="CIF not found"):
            proc._fetch_or_read_cif("9zzz")


class TestFetchOrReadCifRemote:
    """Tests that actually hit the PDBe API (requires network)."""

    @pytest.fixture(autouse=True)
    def _check_network(self):
        """Skip if PDBe is unreachable."""
        try:
            requests.head("https://www.ebi.ac.uk/pdbe/", timeout=5)
        except requests.ConnectionError:
            pytest.skip("PDBe unreachable, skipping remote tests")

    def test_fetches_from_pdbe_and_caches(self, tmp_mmcif_dir):
        proc = TemplateHitProcessor(mmcif_dir=tmp_mmcif_dir, fetch_remote=True)

        result = proc._fetch_or_read_cif("1a2b")

        assert len(result) > 1000, "Downloaded CIF should be a substantial file"
        assert "data_" in result, "CIF content should contain a data_ block"

        cached_path = os.path.join(tmp_mmcif_dir, "1a2b.cif")
        assert os.path.exists(cached_path), "CIF should be cached locally"
        with open(cached_path) as f:
            assert f.read() == result, "Cached content should match returned content"

    def test_second_call_uses_cache(self, tmp_mmcif_dir):
        proc = TemplateHitProcessor(mmcif_dir=tmp_mmcif_dir, fetch_remote=True)

        result1 = proc._fetch_or_read_cif("1a2b")
        result2 = proc._fetch_or_read_cif("1a2b")
        assert result1 == result2

    def test_creates_mmcif_dir_if_missing(self, tmp_path):
        new_dir = str(tmp_path / "nonexistent" / "mmcif")
        assert not os.path.exists(new_dir)

        proc = TemplateHitProcessor(mmcif_dir=new_dir, fetch_remote=True)
        result = proc._fetch_or_read_cif("1a2b")

        assert os.path.isdir(new_dir)
        assert len(result) > 1000

    def test_invalid_pdb_id_raises(self, tmp_mmcif_dir):
        proc = TemplateHitProcessor(mmcif_dir=tmp_mmcif_dir, fetch_remote=True)
        with pytest.raises(requests.HTTPError):
            proc._fetch_or_read_cif("0000")


class TestTemplateHitProcessorInit:
    """Tests that fetch_remote plumbing works correctly."""

    def test_default_fetch_remote_is_false(self, tmp_mmcif_dir):
        proc = TemplateHitProcessor(mmcif_dir=tmp_mmcif_dir)
        assert proc._fetch_remote is False

    def test_fetch_remote_set_true(self, tmp_mmcif_dir):
        proc = TemplateHitProcessor(mmcif_dir=tmp_mmcif_dir, fetch_remote=True)
        assert proc._fetch_remote is True
