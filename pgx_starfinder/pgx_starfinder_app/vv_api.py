import os
import json
import base64
from pathlib import Path
import logging
from dotenv import load_dotenv
import requests
from pgx_starfinder.pgx_starfinder_app.errors import handle_http_error

logger = logging.getLogger(__name__)

# VariantValidator (via TerraCipher SHAIP)
vv_token_url = "https://shaipup.com/market/token"
vv_refresh_url = "https://shaipup.com/market/refreshtoken"
vv_base_url = "https://shaipup.com/market/shaip"
vv_owner = "VariantValidator"
vv_shaip = "rest-variantvalidator2"
load_dotenv(Path(__file__).resolve().parents[2] / ".env")
vv_email = os.environ.get("SHAIP_EMAIL")
vv_password = os.environ.get("SHAIP_PASSWORD")
vv_variant_endpoint = "variantvalidator"

def get_initial_tokens(vv_email, vv_password):
    """Gets the login tokens needed to talk to VariantValidator."""
    # Think of this as logging in. We need a session token so the API knows
    # who we are and allows the next request to go through.
    # Codes were taken from the API documentation.
    credentials = f"{vv_email}:{vv_password}"
    encoded_credentials = base64.b64encode(credentials.encode()).decode()
    headers = {"Authorization": f"Basic {encoded_credentials}"}

    try:
        response = requests.get(vv_token_url, headers=headers, timeout=30)
        response.raise_for_status() # Raises an exception for bad responses
    except requests.exceptions.HTTPError as error:
        return handle_http_error(error, api_name="VariantValidator", context="initial token retrieval")

    id_token = response.headers.get('Authorization', '').replace('Bearer ', '')
    refresh_token = response.json().get('refreshToken')
    session_cookie = response.cookies.get('session')

    return id_token, refresh_token, session_cookie


def refresh_id_token(refresh_token, session_cookie):
    """Uses the Refresh Token to get a new ID Token."""
    # Codes were taken from the API documentation.
    headers = {"Authorization": f"Bearer {refresh_token}"}
    cookies = {"session": session_cookie}

    response = requests.get(vv_refresh_url, headers=headers, cookies=cookies, timeout=30)
    response.raise_for_status()

    new_id_token = response.json().get('userIDToken')
    return new_id_token


def call_algorithm_api(id_token, algorithm_url, payload, session_cookie):
    # Codes were taken from the API documentation.
    """Calls the target algorithm API with the ID Token."""
    headers = {
        "Authorization": f"Bearer {id_token}",
        "Content-Type": "application/json",
    }
    cookies = {"session": session_cookie}

    # VariantValidator can take well over 30 s for a batch, so allow longer here.
    response = requests.post(algorithm_url, headers=headers, cookies=cookies, data=json.dumps(payload), timeout=120)
    response.raise_for_status()
    return response.json()


def standardise_variant(variant, records):
    """
    Picks the GRCh38 version of a variant and the latest MANE transcript.
    """
    result = {"Valid": False, "Reason": None}
    found_record = None
    mane_record = None

    for record in records:
        if not record.get("gene_symbol"):
            continue

        if found_record is None:
            found_record = record

        #annotations contains "mane_select": true or false. If the key is 
        # missing, .get() returns None, which counts as false.
        #It is tested that the out of date MANE Select transcript is not marked as 
        # "mane_select": true, so a separate found_record is needed.

        annotations = record.get("annotations") or {}
        if annotations.get("mane_select"):
            mane_record = record
            break

    if not found_record:
        logger.info("No valid result for %s on GRCh38", variant)
        result["Reason"] = f"VariantValidator could not validate {variant}, check input"
        logger.warning("Could not validate %s", variant)
        return result


    all_loci = found_record.get("primary_assembly_loci") or {}
    grch38_locus = all_loci.get("grch38") or {}
    grch38_hgvs = grch38_locus.get("hgvs_genomic_description")

    grch37_locus = all_loci.get("grch37") or {}
    grch37_hgvs = grch37_locus.get("hgvs_genomic_description")

    gene_id = found_record.get("gene_ids") or {}
    hgnc_id = gene_id.get("hgnc_id") or None

    protein_info = found_record.get("hgvs_predicted_protein_consequence") or {}
    protein_consequences = protein_info.get("tlr")


    notes=[]
    input_reference = variant.split(":")[0]
    grch37_reference = (grch37_hgvs or "").split(":")[0]
    grch38_reference = (grch38_hgvs or "").split(":")[0]

    if input_reference == grch37_reference and input_reference != grch38_reference:
        notes.append(f"Input was on GRCh37; converted to GRCh38 ({grch38_hgvs}).")

    mane_transcript = None

    if mane_record:
        mane_hgvs = mane_record.get("hgvs_transcript_variant") or ""
        mane_transcript = mane_hgvs.split(":")[0]

    #outdated MANE Select could also return None mane_record, but the API will include 
    #a warning in the "validation_warnings" list. If we find a warning, we can extract 
    # the latest MANE Select transcript from it.
    else:
        for warning in found_record.get("validation_warnings") or []:
            if warning.startswith("TranscriptVersionWarning"):
                mane_transcript = warning.rsplit("(", 1)[-1].rstrip(")").strip()
                mane_idx = (found_record.get("hgvs_transcript_variant") or "").split(":")[-1]
                mane_hgvs = f"{mane_transcript}:{mane_idx}"
                break

    if not mane_transcript or not grch38_hgvs:
        result["Reason"] = f"Could not get a GRCh38 position or MANE Select transcript for {variant}, check input"
        logger.warning("Incomplete VariantValidator result for %s", variant)
        return result

    used_hgvs = found_record.get("hgvs_transcript_variant") or ""
    used_transcript = used_hgvs.split(":")[0]

    #Split each transcript into name and version: "NM_000088.3" -> "NM_000088", "3"
    used_name, _, used_version = used_transcript.partition(".")
    mane_name, _, mane_version = mane_transcript.partition(".")

    #Outdated = same transcript as MANE Select, but a lower version number.
    if used_name == mane_name and used_version.isdigit() and mane_version.isdigit():
        if int(used_version) < int(mane_version):
            if used_transcript in variant:
                notes.append(f"Input transcript {used_transcript} is older than MANE Select {mane_transcript}.")
            else:
                notes.append(f"Input {variant} maps to {used_transcript}, which is older than MANE Select {mane_transcript}.")

    result["Valid"] = True
    result["Reason"] = " ".join(notes) if notes else None
    result["GRCh38_HGVS"] = grch38_hgvs
    result["MANE_Select_Transcript"] = mane_hgvs
    result["HGNC_ID"] = hgnc_id
    result["Protein_Consequences"] = protein_consequences

    logger.info("Standardised %s -> %s | %s", variant, grch38_hgvs, mane_transcript)
    return result


def vv_api(variant_hgvs):
    """
    Queries VariantValidator for one HGVS description or a list of them.

    Input: any HGVS description (NC_, NG_, NM_, ENST or LRG), or a list of them.
    All variants are sent in one request.

    Returns: for each variant, GRCh38_HGVS (the variant on the GRCh38 chromosome),
    MANE_Select_Transcript (the current MANE Select transcript) and HGNC_ID.
    A single string input gives one result dictionary; a list gives a dictionary
    of results keyed by each input variant.
    """
    # If the input is missing or not text, stop early.
    # This avoids sending useless requests and makes the error clearer.
    if not isinstance(variant_hgvs, (str, list)) or not variant_hgvs:
        logger.warning("VariantValidator lookup skipped: input is not a string or list (%r)", variant_hgvs)
        return {"Valid": False, "Reason": "Please enter an HGVS variant description."}

    single = isinstance(variant_hgvs, str)
    if single:
        variant_hgvs = [variant_hgvs]
    variants = [str(v).strip() for v in variant_hgvs]

    results = {}
    to_query = []

    for variant in variants:
        if not variant:
            logger.warning("VariantValidator lookup skipped: empty variant hgvs")
            results[variant] = {"Valid": False, "Reason": "Please enter an HGVS variant description."}
        elif ":" not in variant:
            logger.warning("VariantValidator lookup skipped: no ':' in %r", variant)
            results[variant] = {"Valid": False, "Reason": f"Not an HGVS description: {variant}, expected e.g. NM_000088.3:c.589G>T"}
        else:
            to_query.append(variant)

    # If any step below fails, every variant still waiting gets this error.
    error_result = None
    variant_results = {}

    if to_query:
        try:
            tokens = get_initial_tokens(vv_email, vv_password)
            # If the login fails, we get a dictionary with Valid=False and a Reason.
            if not isinstance(tokens, tuple) or len(tokens) != 3:
                error_result = {"Valid": False, "Reason": "VariantValidator login returned an unexpected result, contact IT."}
                logger.error("Unexpected get_initial_tokens result: %r", type(tokens))
            else:
                id_token, refresh_token, session_cookie = tokens
                algorithm_url = f"{vv_base_url}/{vv_variant_endpoint}?owner={vv_owner}&shaip={vv_shaip}"

                payload = [{
                    "variant_description": to_query,
                    "genome_build": "GRCh38",
                    "select_transcripts": ["mane_select"],
                }]

                api_response = call_algorithm_api(id_token, algorithm_url, payload, session_cookie)
                api_response = api_response if isinstance(api_response, list) else [api_response]

                for response in api_response:
                    #skip the "flag" and "metadata" keys, and any non-dict values.
                    for key, record in response.items():
                        if key in ("flag", "metadata") or not isinstance(record, dict):
                            continue

                        # Records from a batch come back together; submitted_variant says which input each belongs to.
                        submitted = record.get("submitted_variant")
                        if submitted not in variant_results:
                            variant_results[submitted] = []
                        variant_results[submitted].append(record)

        except requests.exceptions.HTTPError as error:
            error_result = handle_http_error(error, api_name="VariantValidator", context=f"variant lookup for {to_query}")

        except requests.RequestException as error:
            logger.error("VariantValidator request failed for %s: %s", to_query, error)
            error_result = {"Valid": False, "Reason": f"Network failure during VariantValidator consult: {error}, please try again"}

        except Exception as error:
            logger.exception("Reading VariantValidator response failed for %s", to_query)
            error_result = {"Valid": False, "Reason": f"Failed to extract keys from VariantValidator variant: {error}"}
    # Now we have either a valid variant_results dictionary or an error_result to apply to all variants.
    for variant in to_query:
        if error_result:
            results[variant] = dict(error_result)
            continue
        try:
            results[variant] = standardise_variant(variant, variant_results.get(variant, []))
        except Exception as error:
            results[variant] = {"Valid": False, "Reason": f"Failed to extract keys from VariantValidator variant: {error}"}
            logger.exception("Parsing VariantValidator response failed for %s", variant)

    # One input gives one result; a list gives results keyed by each input, in input order.
    if single:
        return results[variants[0]]
    return {variant: results[variant] for variant in variants}


if __name__ == "__main__":
    # Run from the repository root: python -m pgx_starfinder.pgx_starfinder_app.vv_api
    logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(name)s | %(message)s")
    print(vv_api("NM_000088.3:c.589G>T"))
    print(vv_api(["NM_000088.3:c.589G>T", "NC_000017.10:g.41243432del"]))
    print(vv_api("CYP2C19*2"))
