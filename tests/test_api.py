import json
from pathlib import Path
import requests

MOCK_DIR = Path(__file__).parent / "fixtures"

MOCKS = {
    "cpic_allele": json.loads((MOCK_DIR / "cyp2c19_2_cpic.json").read_text()),
    "cpic_pair": json.loads((MOCK_DIR / "cyp2c19_cpic_drug.json").read_text()),
    "cpic_gene": json.loads((MOCK_DIR / "cyp2c19_cpic.json").read_text()),
    "pharmvar_allele": json.loads((MOCK_DIR / "cyp2c19_2_pharmvar.json").read_text()),
    "pharmvar_gene": json.loads((MOCK_DIR / "cyp2c19_pharmvar.json").read_text()),
}


class FakeResponse:
    def __init__(self, payload=None, status_code=200):
        self.payload = payload
        self.status_code = status_code

    def json(self):
        return self.payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.exceptions.HTTPError(
                f"{self.status_code} error",
                response=self,
            )


def fake_get(url, headers=None, timeout=None, **kwargs):
    matches = {
        "cpicpgx.org/allele?": MOCKS["cpic_allele"],
        "cpicpgx.org/pair?": MOCKS["cpic_pair"],
        "pharmvar.org/alleles/": MOCKS["pharmvar_allele"],
        "pharmvar.org/genes/": MOCKS["pharmvar_gene"],
    }

    for key, value in matches.items():
        if key in url:
            return FakeResponse(value)

    return FakeResponse(None, status_code=404)