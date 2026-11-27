"""
Django views for displaying the search form and rendering unified PharmVar and CPIC results.
"""
import logging

from django.http import HttpRequest, HttpResponse
from django.shortcuts import render
from datetime import timedelta
from django.utils import timezone

from .pgx_api import data_unification
from .models import star_allele_search, star_allele_result
from pprint import pformat

# Sets up logging for view module.
logger = logging.getLogger(__name__)

logger.info("Initialising Django view layer")

HOME_TEMPLATE = "pgx_starfinder/index.html"
RESULTS_TEMPLATE = "pgx_starfinder/results.html"
LINKS_TEMPLATE = "pgx_starfinder/links.html"

# Defines and render application homepage.
def home(request: HttpRequest) -> HttpResponse:
    """
    Renders the application's homepage.

    Args:
        request (HttpRequest): Incoming HTTP request for homepage.

    Returns:
        HttpResponse: Rendered homepage.
    """
    logger.debug("Rendering homepage.")

    return render(request, HOME_TEMPLATE)


# Defines and renders the star allele search and results.
def search(request: HttpRequest) -> HttpResponse:
    """
    Processes star allele search requests. Model developed using code written by Codex AI

    Args:
        request (HttpRequest): Incoming HTTP request contiaining submitted form data

    Returns:
        HttpReponse:  Rendered HTML response containing either search results or an
        explainatory error message.
    """

    # Ensures only POST requests allows searches to be performed.
    if request.method != "POST":
        logger.debug("Ignoring non-POST request.")

        return render(request, HOME_TEMPLATE)

    # prepares parameter input by cleaning whitespaces to allow for star allele search.
    search_term = request.POST.get("star_allele", "").strip()

    # Redisplays the homepage with an error when no search term was submitted.
    if not search_term:
        return render(
            request,
            HOME_TEMPLATE,
            {
                "error_message": "Please enter a star allele.",
                "search_term": search_term,
            },
            status=400,
        )

    logger.info("searching for star allele: %s", search_term)

    # Enables old searches to be invalid when API structure is changed
    # Change this when data_unification's output changes
    CURRENT_SCHEMA_VERSION ="v2"
    # Stored searches can be reused for brief period, after which new search is made
    RESULT_LIFETIME = timedelta(days=7)

    result_instance = None
    normalised_term = search_term.upper()
    now = timezone.now()

    # Search database for valid result
    stored_result = (
        star_allele_result.objects.filter(
            normalised_term=normalised_term,
            schema_version=CURRENT_SCHEMA_VERSION,
            reusable=True,
            expires_at__gt=now,
        )
        .order_by("-fetched_at")
        .first()
    )

    # Either use valid returned database result or run API retrieval
    if stored_result:
        result = stored_result.search_return
        logger.info("database record found for %s", normalised_term)
        result_instance = stored_result
        served_from_database = True
    else:
        logger.info("No existing record, executing API retrieval")
        result = data_unification(search_term)
        served_from_database = False

    # When result is returned from API, check if it is valid to store
    if not served_from_database:
        # Permanent "not found" reasons (e.g. "check input", "Not a star allele") are
        # fine to remember, but temporary failures (network, HTTP, parsing, auth) are
        # not, so the next search tries the APIs again.
        reason_options = ("try again", "contact it", "failed")
        has_temporary_failure = False

        if result:
            for key in (
                "PharmVar_Missing_Reason", 
                "CPIC_Missing_Reason", 
                "ClinPGx_Missing_Reason",
            ):

                reason = (result.get(key) or "").lower()
            
                for option in reason_options:
                    if option in reason:
                        has_temporary_failure = True
                        break

                if has_temporary_failure:
                    logger.info("API return has failures, result return not stored in database")
                    break

        # Create result record if valid to store
        if result and not has_temporary_failure:
            result_instance = star_allele_result.objects.create(
                normalised_term=normalised_term,
                schema_version=CURRENT_SCHEMA_VERSION,
                search_return=result,
                reusable=True,
                expires_at=(now + RESULT_LIFETIME),
            )
            logger.info("search return for %s added to database", normalised_term)
               
    # Create search record and note where served from (db vs API)
    star_allele_search.objects.create(
        submitted_term=search_term,
        normalised_term=normalised_term,
        result=result_instance,
        served_from_database=served_from_database,
    )
    logger.info("new search record stored for %s", search_term)

    # Formats the nested dictionary over multiple lines for a more readable HTML display.
    result_display = pformat(
        result,
        sort_dicts=False,
        width=1
    )

    # Renders formatted result on the results page.
    return render(
        request,
        RESULTS_TEMPLATE,
        {
            "search_term": search_term,
            "result": result,
            "result_display":result_display,
        },
    )

def links(request: HttpRequest) -> HttpResponse:
    """
    Renders the links page containing useful links.

    Args:
        request (HttpRequest): Incoming HTTP request for Useful Links page.

    Returns:
        HttpResponse: Rendered Useful Links page.
    """
    
    logger.debug("Rendering useful links page.")    
    
    return render(request, LINKS_TEMPLATE)
