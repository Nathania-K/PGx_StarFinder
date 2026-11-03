
import os
import logging
import requests
from dotenv import load_dotenv
import pandas as pd
#from pgx_starfinder.pgx_starfinder_app.errors import handle_http_error

# This file connects to external pharmacogenomics databases.
# The goal is to check gene/allele information and return a clean, simple result.
# If something goes wrong, the functions return a friendly error message instead of failing silently.

logger = logging.getLogger(__name__)

# PharmVar base URL and API key.
pm_url = "https://www.pharmvar.org/api-service"
load_dotenv()
pm_key = os.environ.get("PHARMVAR_API_KEY")

# CPIC base URL.
cpic_url = "https://api.cpicpgx.org/v1/"

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

def pharmvar_api(star_allele):
    """
    Queries PharmVar to check star allele information. Possible inputs include
    a specific star allele (output: highly specific variant, validation, HGVS,
    and population frequency data) or gene symbol (output: a list of recorded
    star alleles).
    Never fails silently: a failed call comes back as Valid=False with a Reason.
    """
    # If the input is missing or not a text string, we stop early.
    if pd.isna(star_allele) or not isinstance(star_allele, str):
        logger.warning("PharmVar lookup skipped: input is not a string (%r)", star_allele)
        return None

    result = {"Valid": False, "Reason":  "Please enter a gene or star allele."}

    # Clean the input so we compare names consistently.
    star_allele = str(star_allele).strip().upper()

    is_star_allele = "*" in star_allele

    # PharmVar uses different endpoints for a whole gene versus a specific allele.
    if is_star_allele:
        url = f"{pm_url}/alleles/{star_allele}"
    else:
        url = f"{pm_url}/genes/{star_allele}"

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

        if is_star_allele:
            # For one allele, the API usually returns a list with one main record, take the main record and discard the rest.
            data = data[0] if isinstance(data, list) else data
            result["Allele"] = data.get("alleleName")
            result["Type"] = data.get("alleleType")
            result["Haplotype_HGVS"] = data.get("hgvs")
            result["Variants"] = []

            for variant in data.get("variants") or []:
                ref_collections = variant.get("referenceCollections") or []
                if "GRCh38" in ref_collections:
                    var_info = {
                        "Variant_HGVS": variant.get("hgvs"),
                        "rsID": variant.get("rsId"),
                        "Impact": variant.get("impact"),
                        "GnomAD_Frequency": None,
                    }

                # Some variant records include frequency in GnomAD. Keep the first one found.
                    for freq in variant.get("variantFrequency") or []:
                        if freq.get("source") == "GnomAD":
                            var_info["GnomAD_Frequency"] = freq.get("frequency")
                            break

                    result["Variants"].append(var_info)

        else:
            # For a gene name, summarise how many alleles are documented.
            result["Gene"] = data.get("geneSymbol")
            alleles = data.get("alleles", [])
            result["Total_Alleles_Documented"] = len(alleles)

            core_alleles = []
            for allele in alleles:
                if allele.get("alleleType") == "Core":
                    core_alleles.append(allele.get("alleleName"))

            result["Core_Alleles"] = core_alleles

        return result

    except Exception as error:
        result["Valid"] = False
        result["Reason"] = f"Failed to extract keys from PharmVar: {error}"
        logger.exception("Parsing PharmVar response failed for %s", star_allele)
        return result


def cpic_api(star_allele):
    """
    Queries CPIC to check clinical information of the star allele and its related gene.
    Possible inputs include a specific star allele (output: Academic finding summary,
    list of publication link, drug pair with its originated gene, clinical actionability
    and FDA label) or gene symbol (output: drug pair, clinical actionability
    and FDA label).
    Never fails silently: a failed call comes back as Valid=False with a Reason.
    """
    if pd.isna(star_allele) or not isinstance(star_allele, str):
        logger.warning("CPIC lookup skipped: input is not a string (%r)", star_allele)
        return None

    result = {"Valid": False, "Reason": "Please enter a gene or star allele."}

    # Clean the user input so CYP2C19*2 and cyp2c19*2 behave the same way.
    star_allele = str(star_allele).strip().upper()
    is_star_allele = "*" in star_allele

    headers = {"accept": "application/json"}

    if is_star_allele:
        # Example: CYP2C19*2 should be split into gene = CYP2C19 and allele = *2.
        gene = star_allele.split("*")[0]
        allele_idx = f'*{star_allele.split("*")[1]}'
        url_allele = f"https://api.cpicpgx.org/v1/allele?name=eq.{allele_idx}&genesymbol=eq.{gene}"
        url_gene = f"https://api.cpicpgx.org/v1/pair?genesymbol=eq.{gene}"

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
            result["Academic_Findings"] = record.get("findings", None)
            cite_list = record.get("citations") or []
            pubmed_str = "https://pubmed.ncbi.nlm.nih.gov/"
            result["Citation_Link"] = [f"{pubmed_str}{pmid}" for pmid in cite_list]

        except Exception as error:
            result["Valid"] = False
            result["Reason"] = f"Data parsing for CPIC star allele failed. Reason: {error}"
            logger.exception("Parsing CPIC allele response failed for %s", star_allele)
            return result

    else:
        gene = star_allele
        url_gene = f"https://api.cpicpgx.org/v1/pair?genesymbol=eq.{gene}"

    # DEBUG: the gene-level drug query is made for both allele and gene searches.
    logger.debug("Querying CPIC drug pairs: %s", url_gene)

    # For both gene-only and allele-specific queries, we also fetch drug information by gene.
    try:
        response = requests.get(url_gene, headers=headers, timeout=30)
        response.raise_for_status()
        gene_dict = response.json()

    except requests.exceptions.HTTPError as error:
        return handle_http_error(error, api_name="CPIC", context="gene consult")

    except requests.RequestException as error:
        result["Valid"] = False
        result["Reason"] = f"Network failure during CPIC gene consult: {error}, please try again"
        logger.error("CPIC drug-pair request failed for %s: %s", gene, error)
        return result

    try:
        if not gene_dict:
            result["Valid"] = False
            result["Reason"] = f"No drug records found in CPIC for gene: {gene}, check input"
            logger.warning("No CPIC drug pairs for gene %s", gene)
            return result

        # Store the gene-level drug info in a simple list.
        result["Valid"] = True
        result["Gene"] = gene
        result["Drug INFO"] = []

        for entry in gene_dict:
            drug_info = {
                "Drug ID": entry.get("drugid"),
                "Clinical Actionability": entry.get("cpiclevel"),
                "FDA label": entry.get("pgxtesting"),
            }
            result["Drug INFO"].append(drug_info)

        logger.info("CPIC: %s found with %d drug pair(s)", star_allele, len(result["Drug INFO"]))

        return result

    except Exception as error:
        result["Valid"] = False
        result["Reason"] = f"Data parsing for drug info failed. Reason: {error}"
        logger.exception("Parsing CPIC drug-pair response failed for %s", gene)
        return result


def data_unification(star_allele):
    """
    Combine PharmVar and CPIC results into one easy-to-read output.

    This is the final step: gather both the genetic information and the clinical
    drug information, then merge them into a neat summary. If one database fails,
    we still keep the data from the other database when possible.
    """
    if pd.isna(star_allele) or not isinstance(star_allele, str):
        logger.warning("Search skipped: input is not a string (%r)", star_allele)
        return None

    # Ask both databases separately.
    pharmvar_data = pharmvar_api(star_allele)
    cpic_data = cpic_api(star_allele)

    is_star_allele = "*" in star_allele
    gene_symbol = star_allele.split("*")[0] if is_star_allele else star_allele

    if pharmvar_data["Valid"] is True or cpic_data["Valid"] is True:
        if not pharmvar_data["Valid"]:
            logger.warning("Partial result for %s: PharmVar failed (%s)", star_allele, pharmvar_data["Reason"])
        elif not cpic_data["Valid"]:
            logger.warning("Partial result for %s: CPIC failed (%s)", star_allele, cpic_data["Reason"])

        # Data fields only ever hold data (or None / an empty list when missing),
        # and each database's error message has its own field. A template can
        # then loop over lists safely and show the reasons separately; mixing
        # them made a template loop over an error message one letter at a time.
        unified_result = {
            "Target_Gene": gene_symbol,
            "Gene_Data": None,
            "Allele_Data": None,
            "Gene_Drug_Table": [],
            "PharmVar_Missing_Reason": None if pharmvar_data["Valid"] else pharmvar_data["Reason"],
            "CPIC_Missing_Reason": None if cpic_data["Valid"] else cpic_data["Reason"],
        }

        # Format drug information into a readable list.
        drug_raw_data = cpic_data.get("Drug INFO")
        if isinstance(drug_raw_data, list):
            expanded_drugs = []

            for drug in drug_raw_data:
                drug_name = drug.get('Drug ID', 'Unknown')
                cpic_level = drug.get('Clinical Actionability', 'N/A')
                fda_label = drug.get('FDA label', 'None')

                formatted_string = f"{drug_name} (CPIC Level: {cpic_level} | FDA: {fda_label})"
                expanded_drugs.append(formatted_string)

            unified_result["Gene_Drug_Table"] = expanded_drugs

        if is_star_allele:
            # If the user gave a specific allele, keep allele-level details together.
            unified_result["Allele_Data"] = {
                "Allele": pharmvar_data.get("Allele"),
                "Type": pharmvar_data.get("Type"),
                "Haplotype_HGVS": pharmvar_data.get("Haplotype_HGVS"),
                "Variants": pharmvar_data.get("Variants") or [],
                "Academic_Findings": cpic_data.get("Academic_Findings"),
                "Citation_Link": cpic_data.get("Citation_Link") or [],
            }

        else:
            # For gene-only searches, keep the gene summary together.
            # Sort the core alleles by their numeric suffix for better readability
            core_alleles = pharmvar_data.get("Core_Alleles") or []
            if core_alleles:
                core_alleles.sort(key=lambda x: int(''.join(c for c in x.split('*')[-1] if c.isdigit()) or 0))
            unified_result["Gene_Data"] = {
                "Total_Alleles_Documented": pharmvar_data.get("Total_Alleles_Documented"),
                "Core_Alleles": core_alleles,
            }
    else:
        # If both databases failed, return a combined error message.
        unified_result = {
            "Valid": False,
            "PharmVar_Missing_Reason": pharmvar_data.get("Reason"),
            "CPIC_Missing_Reason": cpic_data.get("Reason")
        }
        logger.warning("No result for %s: %s", star_allele, unified_result["Reason"])
    return unified_result

if __name__ == "__main__":
    # Set up basic console logging to see the messages while testing.
    logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(name)s | %(message)s")

    # Simple test lines for quick checking during development.
    print(data_unification("CYP2C19*2"))
    print(data_unification("CYP2C19"))
