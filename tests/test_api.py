import json
from pathlib import Path
from unittest.mock import patch

import pytest
import requests

from pgx_starfinder.pgx_starfinder_app.pgx_api import pharmvar_api
from pgx_starfinder.pgx_starfinder_app.vv_api import vv_api


MOCK_DIR = Path(__file__).parent / "mock_data"

API_pgx = "pgx_starfinder.pgx_starfinder_app.pgx_api"
API_vv = "pgx_starfinder.pgx_starfinder_app.vv_api"

MOCKS = {
    "cpic_allele": json.loads((MOCK_DIR / "cyp2c19_2_cpic.json").read_text()),
    "pharmvar_allele": json.loads((MOCK_DIR / "cyp2c19_2_pharmvar.json").read_text()),
    "clinpgx_gene": json.loads((MOCK_DIR / "clinpgx_annotation" / "clinpgx_summary.json").read_text()),
    "clinpgx_guideline": json.loads((MOCK_DIR / "clinpgx_annotation" / "clinpgx_guideline_example.json").read_text()),
    "vv_correct_38": json.loads((MOCK_DIR / "vv_variable" / "vv_variant_grch38.json").read_text()),
    "vv_correct_37": json.loads((MOCK_DIR / "vv_variable" / "vv_variant_grch37.json").read_text()),
    "vv_incorrect": json.loads((MOCK_DIR / "vv_variable" / "vv_variant_grch37_wrong.json").read_text()),
    "vv_NG": json.loads((MOCK_DIR / "vv_variable" / "vv_NG.json").read_text()),
    "vv_NG_outdated": json.loads((MOCK_DIR / "vv_variable" / "vv_NG_outdated.json").read_text()),
    "vv_LRG": json.loads((MOCK_DIR / "vv_variable" / "vv_LRG.json").read_text()),
    "vv_LRG_outdated": json.loads((MOCK_DIR / "vv_variable" / "vv_LRG_outdated.json").read_text()),
}

# Fake responses

class FakeResponse:
    def __init__(self, payload=None, status_code=200):
        self.payload = payload
        self.status_code = status_code
        # VariantValidator login reads the token from headers and the session from cookies.
        self.headers = {"Authorization": "Bearer fake-id-token"}
        self.cookies = {"session": "fake-session"}

    def json(self):
        return self.payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.exceptions.HTTPError(
                f"{self.status_code} error",
                response=self,
            )


# Each VariantValidator mock, found by the variant it was made for (its "submitted_variant").
VV_MOCKS = {}
for name in ("vv_correct_38", "vv_correct_37", "vv_incorrect", "vv_NG", "vv_NG_outdated", "vv_LRG", "vv_LRG_outdated"):
    for record in MOCKS[name].values():
        if isinstance(record, dict) and "submitted_variant" in record:
            VV_MOCKS[record["submitted_variant"]] = MOCKS[name]


def fake_vv_lookup(data):
    """Returns the VariantValidator mock for each variant in the request.

    The real request sends the variants as JSON in data=..., e.g.
    [{"variant_description": ["NM_000088.3:c.589G>T"], ...}].
    A variant with no mock gets the "could not validate" mock (vv_incorrect).
    """
    variants = json.loads(data)[0]["variant_description"]
    return [VV_MOCKS.get(variant, MOCKS["vv_incorrect"]) for variant in variants]


def fake_get(url, headers=None, timeout=None, **kwargs):
    # VariantValidator: login (GET), then the variant lookup (POST, patch requests.post with fake_get too).
    if "shaipup.com/market/shaip/variantvalidator" in url:
        return FakeResponse(fake_vv_lookup(kwargs["data"]))

    matches = {
        "cpicpgx.org/v1/allele?": MOCKS["cpic_allele"],
        "pharmvar.org/api-service/alleles/": MOCKS["pharmvar_allele"],
        "clinpgx.org/v1/data/summaryAnnotation": MOCKS["clinpgx_gene"],
        "clinpgx.org/v1/data/guidelineAnnotation/": MOCKS["clinpgx_guideline"],
        "shaipup.com/market/token": {"refreshToken": "fake-refresh-token"},
    }

    for key, value in matches.items():
        if key in url:
            return FakeResponse(value)

    return FakeResponse(None, status_code=404)

# Fixtures: switch the fake APIs on for a test

@pytest.fixture
def mock_pgx_apis():
    """Replaces PharmVar, CPIC and ClinPGx requests with the mock data."""
    with patch(f"{API_pgx}.requests.get", side_effect=fake_get), \
         patch("time.sleep"): # skip pauses between ClinPGx requests
        yield


@pytest.fixture
def mock_vv_api():
    """Replaces the VariantValidator login (GET) and lookup (POST) with the mock data."""
    with patch(f"{API_vv}.requests.get", side_effect=fake_get), \
         patch(f"{API_vv}.requests.post", side_effect=fake_get):
        yield


def submitted_variant(mock_name):
    """Returns the variant a VariantValidator mock was made for, so tests always match it."""
    for record in MOCKS[mock_name].values():
        if isinstance(record, dict) and "submitted_variant" in record:
            return record["submitted_variant"]


# Tests for pharmvar_api()

def test_pharmvar_returns_allele(mock_pgx_apis):
    """
    Tests if a known star allele returns a valid result from the mock data.
    """
    result = pharmvar_api("CYP2C19*2")

    assert result["Valid"] is True
    assert result["Allele"] == "CYP2C19*2"
    assert result["Variants"]


def test_pharmvar_not_found():
    """
    Tests if a 404 from PharmVar comes back as Valid=False with a reason.
    """
    def fake_404(url, **kwargs):
        return FakeResponse(None, status_code=404)

    with patch(f"{API_pgx}.requests.get", side_effect=fake_404):
        result = pharmvar_api("CYP2C19*999")

    assert result["Valid"] is False
    assert "not found" in result["Reason"].lower()


# Tests for vv_api()

def test_vv_outdated_lrg_gives_note(mock_vv_api):
    """
    Tests if an LRG input mapping to an outdated transcript is standardised with a note.
    """
    variant = submitted_variant("vv_LRG_outdated")

    result = vv_api([variant])[variant]

    assert result["Valid"] is True
    assert result["GRCh38_HGVS"].startswith("NC_")
    assert "older than MANE Select" in result["Reason"]


def test_vv_invalid_variant(mock_vv_api):
    """
    Tests if a variant VariantValidator cannot validate returns Valid=False.
    """
    variant = submitted_variant("vv_incorrect")

    result = vv_api([variant])[variant]

    assert result["Valid"] is False
    assert result["Reason"]


# TODO: CPIC, ClinPGx and data_unification tests.
# Use mock_pgx_apis for success cases, and patch(...) with a failing
# fake (as in test_pharmvar_not_found) for error cases.