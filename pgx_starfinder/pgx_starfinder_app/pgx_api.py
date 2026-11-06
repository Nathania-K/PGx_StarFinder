import html
import logging
import os
import re
import time
from pathlib import Path

import pandas as pd
import requests
from dotenv import load_dotenv

from pgx_starfinder.pgx_starfinder_app.errors import handle_http_error
from pgx_starfinder.pgx_starfinder_app.vv_api import vv_api

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
        result["PharmVar_ID"] = data.get("pvId")
        result["Variants"] = []

        variant_list = []
        for variant in data.get("variants") or []:
            ref_collections = variant.get("referenceCollections") or []
            # Only include variants that have a GRCh38 reference collection.
            if "GRCh38" in ref_collections:
                var_info = {
                    "Variant_HGVS": variant.get("hgvs") if variant.get("hgvs") else None,
                    "rsID": variant.get("rsId") if variant.get("rsId") else None,
                    "rsID_Link": f"https://www.ncbi.nlm.nih.gov/snp/{variant.get('rsId')}" if variant.get("rsId") else None,
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
        result["Allele_Description"] = None
        allele_description = []
        if variant_list:
            vv_results = vv_api(variant_list)
            for var_info in result["Variants"]:
                vv_result = vv_results.get(var_info["Variant_HGVS"]) or {}
                var_info["Protein_Consequences"] = vv_result.get("Protein_Consequences")
                var_info["MANE_Select_Transcript"] = vv_result.get("MANE_Select_Transcript")
                allele_description.append(var_info["MANE_Select_Transcript"])
                if result["HGNC_ID"] is None:
                    result["HGNC_ID"] = vv_result.get("HGNC_ID")
                if not vv_result.get("Valid"):
                    logger.warning("No VariantValidator result for %s: %s",
                                    var_info["Variant_HGVS"], vv_result.get("Reason"))
            #merge different "NM_000088.4:c.589G>T" into one comprehensive cell content.
            base = None
            transcripts = []
            for allele in allele_description:
                if not allele:
                    continue
                if base is None:
                    base = allele.rsplit(".", 1)[0]
                transcript = allele.rsplit(".", 1)[-1]
                transcripts.append(transcript)
            description_merge = ";".join(transcripts)
            merge = f"{base}.[{description_merge}]"
            result["Allele_Description"] = merge if base else None



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
            return {"Valid": False, "Reason": f"No ClinPGx annotations found for {star_allele}, check input"}

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

        # Display evidence levels from highest to lowest, with 1A first, putting anything unexpected at the end.
        evidence_order = {"1A": 0, "1B": 1, "2A": 2, "2B": 3, "3": 4, "4": 5}

        drugs = list(merged.values())
        for drug in drugs:
            drug["Annotations"].sort(
                key=lambda annotation: evidence_order.get(
                    str(annotation.get("Level") or "").strip().upper(),
                    len(evidence_order),
                )
            )
        drugs.sort(
            key=lambda drug: min(
                (   
                    evidence_order.get(
                        str(annotation.get("Level") or "").strip().upper(),
                        len(evidence_order),
                    )
                    for annotation in drug["Annotations"]
                ),
                default=len(evidence_order),
            )
        )

        #Each guideline is fetched once, even if several drugs share it.
        guideline_cache = {}

        for drug in drugs:
            # make the drug display in order of lowest level of evidence to highest level of evidence
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
            "Type": pharmvar_data.get("Type") or pharmvar_data.get("Reason"),
            "PharmVar_ID": pharmvar_data.get("PharmVar_ID") or pharmvar_data.get("Reason"),
            "Allele_Description": pharmvar_data.get("Allele_Description") or pharmvar_data.get("Reason"),
            "HGNC_ID": pharmvar_data.get("HGNC_ID") or pharmvar_data.get("Reason"),
            "Variants": pharmvar_data.get("Variants") or [],
            "Academic_Findings": cpic_data.get("Academic_Findings") or cpic_data.get("Reason"),
            "Citation_Link": cpic_data.get("Citation_Link") or [],
        },
        # Drug information from ClinPGx: one entry per drug, with its annotations and guidelines.
        "Drug_Data": clinpgx_data.get("Drugs") or clinpgx_data.get("Reason"),
        **reasons,
    }
    return unified_result


if __name__ == "__main__":
    # Run from the repository root: python -m pgx_starfinder.pgx_starfinder_app.pgx_api
    logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(name)s | %(message)s")
    print(pharmvar_api("CYP2C19*2"))
    print(cpic_api("CYP2C19*2"))
    print(clinpgx_api("CYP2C19*2"))
    print(data_unification("CYP2C19*2"))
    print(data_unification("CYP2C19*9"))
