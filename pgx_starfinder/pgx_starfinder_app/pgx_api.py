import base64
import html
import json
import logging
import os
import re
import time
from pathlib import Path

import pandas as pd
import requests
from dotenv import load_dotenv

# This file talks to external pharmacogenomics databases.
# In simple words: it checks a star allele, gathers useful facts,
# and returns a clean result. If something fails, it explains the problem
# instead of crashing or giving a confusing blank answer.

logger = logging.getLogger(__name__)

# API keys and logins are read from the .env file at the repository root.
load_dotenv(Path(__file__).resolve().parents[2] / ".env")

# PharmVar base URL and API key.
pm_url = "https://www.pharmvar.org/api-service"
pm_key = os.environ.get("PHARMVAR_API_KEY")

# CPIC base URL.
cpic_url = "https://api.cpicpgx.org/v1"

# ClinPGx base URL.
clinpgx_url = "https://api.clinpgx.org/v1"

# PubMed link prefix for citations.
pubmed_url = "https://pubmed.ncbi.nlm.nih.gov/"

# VariantValidator (via TerraCipher SHAIP)
vv_token_url = "https://shaipup.com/market/token"
vv_refresh_url = "https://shaipup.com/market/refreshtoken"
vv_base_url = "https://shaipup.com/market/shaip"
vv_owner = "VariantValidator"
vv_shaip = "rest-variantvalidator2"
vv_email = os.environ.get("SHAIP_EMAIL")
vv_password = os.environ.get("SHAIP_PASSWORD")
vv_g2t_endpoint = "gene2transcripts-v2"
vv_variant_endpoint = "variantvalidator"


def handle_http_error(error, api_name="API", context="consult"):
    """
    Takes an HTTPError object and returns a standardised result dictionary.
    """
    response = getattr(error, "response", None)
    status_code = response.status_code if response is not None else None

    # A few common HTTP problems need a friendlier explanation.

    if status_code in (401, 403):
        reason = f"{api_name} API unauthorised, contact IT."
        logger.error("%s returned %s during %s: API key missing or invalid",
                     api_name, status_code, context)
    elif status_code == 404:
        reason = f"Endpoint or item not found in {api_name}, check input."
        logger.warning("%s returned 404 during %s: %s", api_name, context, error)
    else:
        reason = f"HTTP Error {status_code} during {context}: {error}, please try again."
        logger.error("%s returned HTTP %s during %s: %s",
                     api_name, status_code, context, error)

    return {"Valid": False, "Reason": reason}


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
    result["MANE_Select_Transcript"] = mane_transcript
    result["HGNC_ID"] = hgnc_id

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


def pharmvar_api(star_allele):
    """
    Queries PharmVar to check star allele information. Possible inputs include
    a specific star allele or a PharmVar allele ID (output: highly specific variant, validated HGVS and
    MANE Select transcript).
    Never fails silently: a failed call comes back as Valid=False with a Reason.
    """
    # If the input is missing or not a text string, we stop early.
    if pd.isna(star_allele) or not isinstance(star_allele, str):
        logger.warning("PharmVar lookup skipped: input is not a string (%r)", star_allele)
        return None

    result = {"Valid": False, "Reason": "Please enter a star allele."}

    # Clean the input so we compare names consistently.
    star_allele = str(star_allele).strip().upper()

    is_star_allele = "*" in star_allele
    is_PV = star_allele.startswith("PV") and len(star_allele) > 2 and star_allele[2:].isdigit()

    # Star alleles and PharmVar IDs both use the alleles endpoint; gene symbols are not supported.
    if is_star_allele or is_PV:
        url = f"{pm_url}/alleles/{star_allele}"
    else:
        result["Valid"] = False
        result["Reason"] = f"Please enter a star allele (e.g., CYP2C19*2) or PharmVar allele ID (e.g., PV12345) rather than a gene symbol ({star_allele})."
        logger.error("PharmVar gene-level lookup not implemented for %s", star_allele)
        return result

    headers = {
        "accept": "application/json",
        "api-key": pm_key,
    }

    logger.debug("Querying PharmVar: %s", url)

    try:
        # Send the request to PharmVar.
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()
        data = response.json()

    except requests.exceptions.HTTPError as error:
        return handle_http_error(error, api_name="PharmVar", context="allele consult")

    except requests.RequestException as error:
        result["Reason"] = f"Network failure during PharmVar consult: {error}, please try again"
        logger.error("PharmVar request failed for %s: %s", star_allele, error)
        return result

    try:
        # PharmVar may return no matching record.
        if not data:
            result["Reason"] = f"No records found in PharmVar for: {star_allele}, check input"
            logger.warning("No PharmVar records for %s", star_allele)
            return result

        result["Valid"] = True
        result["Reason"] = None

        # For one allele, the API usually returns a list with one main record, take the main record and discard the rest.
        data = data[0] if isinstance(data, list) else data
        result["Allele"] = data.get("alleleName")
        result["Type"] = data.get("alleleType")
        result["Haplotype_HGVS"] = data.get("hgvs")
        result["Variants"] = []

        variant_list = []
        for variant in data.get("variants") or []:
            ref_collections = variant.get("referenceCollections") or []
            # Only include variants that have a GRCh38 reference collection.
            if "GRCh38" in ref_collections:
                var_info = {
                    "Variant_HGVS": variant.get("hgvs"),
                    "rsID": variant.get("rsId"),
                    "Impact": variant.get("impact"),
                    "GnomAD_Frequency": None,
                }
                variant_list.append(var_info["Variant_HGVS"])

                # Some variant records include frequency in GnomAD. Keep the first one found.
                for freq in variant.get("variantFrequency") or []:
                    if freq.get("source") == "GnomAD":
                        var_info["GnomAD_Frequency"] = freq.get("frequency")
                        break

                result["Variants"].append(var_info)

        # Check the GRCh38 variants with VariantValidator (one request for all of them)
        # to add the gene's HGNC ID and each variant's MANE Select transcript.
        result["HGNC_ID"] = None
        if variant_list:
            vv_results = vv_api(variant_list)
            for var_info in result["Variants"]:
                vv_result = vv_results.get(var_info["Variant_HGVS"]) or {}
                var_info["MANE_Select"] = vv_result.get("MANE_Select_Transcript")
                if result["HGNC_ID"] is None:
                    result["HGNC_ID"] = vv_result.get("HGNC_ID")
                if not vv_result.get("Valid"):
                    logger.warning("No VariantValidator result for %s: %s",
                                    var_info["Variant_HGVS"], vv_result.get("Reason"))

        return result

    except Exception as error:
        result["Valid"] = False
        result["Reason"] = f"Failed to extract keys from PharmVar: {error}"
        logger.exception("Parsing PharmVar response failed for %s", star_allele)
        return result


def cpic_api(star_allele):
    """
    Queries CPIC to check clinical information of the star allele and its academic summary
    and publications.
    Possible inputs include a specific star allele (output: Academic finding summary and
    list of publication links).
    Never fails silently: a failed call comes back as Valid=False with a Reason.
    """
    if pd.isna(star_allele) or not isinstance(star_allele, str):
        logger.warning("CPIC lookup skipped: input is not a string (%r)", star_allele)
        return None

    result = {"Valid": False, "Reason": "Please enter a star allele."}

    # Clean the user input so CYP2C19*2 and cyp2c19*2 behave the same way.
    star_allele = str(star_allele).strip().upper()

    headers = {"accept": "application/json"}

    # Example: CYP2C19*2 should be split into gene = CYP2C19 and allele = *2.
    gene = star_allele.split("*")[0]
    allele_idx = f'*{star_allele.split("*")[1]}'
    url_allele = f"{cpic_url}/allele?name=eq.{allele_idx}&genesymbol=eq.{gene}"

    logger.debug("Querying CPIC allele: %s", url_allele)

    try:
        # First fetch the allele-specific record.
        allele_response = requests.get(url_allele, headers=headers, timeout=30)
        allele_response.raise_for_status()
        allele_info = allele_response.json()

    except requests.exceptions.HTTPError as error:
        return handle_http_error(error, api_name="CPIC", context="allele consult")

    except requests.RequestException as error:
        result["Reason"] = f"Network failure during CPIC allele consult: {error}, please try again"
        logger.error("CPIC allele request failed for %s: %s", star_allele, error)
        return result

    try:
        if not allele_info:
            result["Reason"] = f"No CPIC records found for star allele: {star_allele}, check input"
            logger.warning("No CPIC allele record for %s", star_allele)
            return result

        # Take the first matching result and keep the most useful fields.
        record = allele_info[0]
        result["Allele"] = star_allele
        result["Valid"] = True
        result["Reason"] = None
        result["Academic_Findings"] = record.get("findings", None)
        cite_list = record.get("citations") or []
        result["Citation_Link"] = [f"{pubmed_url}{pmid}" for pmid in cite_list]
        return result

    except Exception as error:
        result["Valid"] = False
        result["Reason"] = f"Data parsing for CPIC star allele failed. Reason: {error}"
        logger.exception("Parsing CPIC allele response failed for %s", star_allele)
        return result


def clinpgx_api(star_allele):
    """
    Queries ClinPGx for the drugs linked to a star allele.
    Possible inputs include a specific star allele (output: one entry per drug, with its
    annotations - related disease, level of evidence, allele effect and type of evidence -
    and the details of any related prescribing guidelines).
    Never fails silently: a failed call comes back as Valid=False with a Reason.
    """
    if pd.isna(star_allele) or not isinstance(star_allele, str):
        logger.warning("ClinPGx lookup skipped: input is not a string (%r)", star_allele)
        return None

    result = {"Valid": False, "Reason": "Please enter a star allele."}

    # Clean the user input so CYP2C19*2 and cyp2c19*2 behave the same way.
    star_allele = str(star_allele).strip().upper()

    gene, _, number = star_allele.partition("*")
    allele = f"*{number}"

    try:
        response = requests.get(
            f"{clinpgx_url}/data/summaryAnnotation",
            params={"location.genes.symbol": gene},
            headers={"accept": "application/json"},
            timeout=30,
        )
        response.raise_for_status()
        clinpgx_data = response.json().get("data") or []

    except requests.exceptions.HTTPError as error:
        return handle_http_error(error, api_name="ClinPGx", context="allele consult")

    except requests.RequestException as error:
        result["Reason"] = f"Network failure during ClinPGx allele consult: {error}, please try again"
        logger.error("ClinPGx allele request failed for %s: %s", star_allele, error)
        return result

    try:
        if not clinpgx_data:
            result["Reason"] = f"No ClinPGx records found for star allele: {star_allele}, check input"
            logger.warning("No ClinPGx allele record for %s", star_allele)
            return result

        drugs = []
        for annotation in clinpgx_data:
            for allele_phenotype in annotation.get("allelePhenotypes") or []:
                if allele_phenotype.get("allele") == allele:
                    drugs.append({
                        "Drug": ", ".join(c.get("name") for c in annotation.get("relatedChemicals") or []),
                        "Related_Disease": ", ".join(c.get("name") for c in annotation.get("relatedDiseases") or []),
                        "Level": (annotation.get("levelOfEvidence") or {}).get("term"),
                        "Allele_Effect": allele_phenotype.get("phenotype"),
                        "Type_of_Evidence": ", ".join(annotation.get("types") or []),
                        "Guideline_IDs": [g.get("id") for g in annotation.get("relatedGuidelines") or []],
                    })

        if not drugs:
            return {"Valid": False, "Reason": f"No ClinPGx annotations found for {star_allele}"}

        #Merge entries with the same drug name: one entry per drug, its annotations inside.
        merged = {}
        for drug in drugs:
            name = drug["Drug"]
            if name not in merged:
                merged[name] = {"Drug": name, "Annotations": [], "Guideline_IDs": []}

            merged[name]["Annotations"].append({
                "Related_Disease": drug["Related_Disease"],
                "Level": drug["Level"],
                "Allele_Effect": drug["Allele_Effect"],
                "Type_of_Evidence": drug["Type_of_Evidence"],
            })

            #Add each guideline ID only once per drug.
            for guideline_id in drug["Guideline_IDs"]:
                if guideline_id not in merged[name]["Guideline_IDs"]:
                    merged[name]["Guideline_IDs"].append(guideline_id)

        drugs = list(merged.values())

        #Each guideline is fetched once, even if several drugs share it.
        guideline_cache = {}

        for drug in drugs:
            drug["Guideline_Details"] = []
            #Fetch each guideline's details from ClinPGx, but only once per guideline ID.
            for guideline_id in drug["Guideline_IDs"]:
                if guideline_id not in guideline_cache:
                    guideline_link = f"{clinpgx_url}/data/guidelineAnnotation/{guideline_id}"
                    try:
                        # ClinPGx has a rate limit of 5 requests per second, so we wait 1 s between requests.
                        time.sleep(1)
                        guideline_response = requests.get(
                            guideline_link,
                            headers={"accept": "application/json"},
                            timeout=30,
                        )
                        # If we hit the rate limit, wait for the time specified in the Retry-After header (default 5 s) and retry.
                        if guideline_response.status_code == 429:
                            wait = int(guideline_response.headers.get("Retry-After", 5))
                            logger.info("ClinPGx rate limit hit; waiting %s s before retrying %s", wait, guideline_id)
                            time.sleep(wait)
                            guideline_response = requests.get(guideline_link, headers={"accept": "application/json"}, timeout=30)

                        guideline_response.raise_for_status()
                        guideline_data = guideline_response.json().get("data") or {}

                    except (requests.exceptions.RequestException, ValueError) as error:
                        logger.warning("Failed to fetch guideline data for %s: %s", guideline_link, error)
                        guideline_cache[guideline_id] = None
                        continue

                    #Turn the HTML summary into plain text:
                    summary_html = (guideline_data.get("summaryMarkdown") or {}).get("html") or ""
                    summary = re.sub(r"<[^>]+>", " ", summary_html)
                    summary = html.unescape(summary)
                    summary = re.sub(r"\s+", " ", summary).strip()

                    guideline_cache[guideline_id] = {
                        "Title": guideline_data.get("name"),
                        "Source": guideline_data.get("source"),
                        "Summary": summary or None,
                        "Citation_Link": ", ".join(
                            c.get("_sameAs") for c in guideline_data.get("literature") or [] if c.get("_sameAs")
                        ),
                    }

                #Outside the "if": runs for new and already-fetched guidelines alike.
                if guideline_cache.get(guideline_id):
                    drug["Guideline_Details"].append(guideline_cache[guideline_id])

        result["Valid"] = True
        result["Reason"] = None
        result["Allele"] = star_allele
        result["Drugs"] = drugs
        logger.info("ClinPGx: %s found for %d drug(s)", star_allele, len(drugs))
        return result

    except Exception as error:
        result["Valid"] = False
        result["Reason"] = f"Data parsing for drug info failed. Reason: {error}"
        logger.exception("Parsing ClinPGx drug-pair response failed for %s", gene)
        return result


def data_unification(star_allele):
    """
    Combine PharmVar, VariantValidator, CPIC, and ClinPGx results into one easy-to-read output.

    This is the final step: gather both the genetic information and the clinical
    drug information, then merge them into a neat summary. If one database fails,
    we still keep the data from the other database when possible.
    """
    if pd.isna(star_allele) or not isinstance(star_allele, str):
        logger.warning("Search skipped: input is not a string (%r)", star_allele)
        return None

    # PharmVar first: it also accepts PharmVar IDs (e.g. PV12345), and the allele
    # name it returns (e.g. CYP2C19*2) is what CPIC and ClinPGx need.
    pharmvar_data = pharmvar_api(star_allele)
    allele_name = pharmvar_data.get("Allele") or star_allele.strip().upper()

    # CPIC and ClinPGx only work with star allele names (GENE*number), so a gene
    # symbol or an unknown PharmVar ID is not sent to them.
    if "*" in allele_name:
        cpic_data = cpic_api(allele_name)
        clinpgx_data = clinpgx_api(allele_name)
    else:
        cpic_data = {"Valid": False, "Reason": f"Not a star allele: {allele_name}, e.g. CYP2C19*2"}
        clinpgx_data = cpic_data

    # Each source's error message has its own field, so data fields only ever
    # hold data (or None / an empty list) and a template can loop over them safely.
    reasons = {
        "PharmVar_Missing_Reason": None if pharmvar_data["Valid"] else pharmvar_data["Reason"],
        "CPIC_Missing_Reason": None if cpic_data["Valid"] else cpic_data["Reason"],
        "ClinPGx_Missing_Reason": None if clinpgx_data["Valid"] else clinpgx_data["Reason"],
    }

    if not (pharmvar_data["Valid"] or cpic_data["Valid"] or clinpgx_data["Valid"]):
        # If every source failed, return only the reasons.
        logger.warning("No result for %s: %s", star_allele, reasons)
        return {"Valid": False, **reasons}

    failed = [name for name, data in (("PharmVar", pharmvar_data), ("CPIC", cpic_data), ("ClinPGx", clinpgx_data))
              if not data["Valid"]]
    if failed:
        logger.warning("Partial result for %s: %s failed", star_allele, ", ".join(failed))

    unified_result = {
        "Valid": True,
        "Allele_Data": {
            "Allele": allele_name,
            "Type": pharmvar_data.get("Type"),
            "Haplotype_HGVS": pharmvar_data.get("Haplotype_HGVS"),
            "HGNC_ID": pharmvar_data.get("HGNC_ID"),
            "Variants": pharmvar_data.get("Variants") or [],
            "Academic_Findings": cpic_data.get("Academic_Findings"),
            "Citation_Link": cpic_data.get("Citation_Link") or [],
        },
        # Drug information from ClinPGx: one entry per drug, with its annotations and guidelines.
        "Drug_Data": clinpgx_data.get("Drugs") or [],
        **reasons,
    }
    return unified_result


if __name__ == "__main__":
    # Set up basic console logging to see the messages while testing.
    logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(name)s | %(message)s")

    # Simple test lines for quick checking during development.

    #print(data_unification("CYP2C19*2"))
    #print(data_unification("CYP2C19*9"))
    #print(data_unification("CYP2C19*9"))
    print(data_unification("CYP2C19*2"))
