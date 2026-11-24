"""
Django views for displaying the search form and rendering unified PharmVar and CPIC results.
"""
import logging

from django.http import HttpRequest, HttpResponse
from django.shortcuts import render
from django.core.cache import cache

from .pgx_api import data_unification
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
    Processes star allele search requests. 

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

    #Add cache function, can resume even if the website restarts
    CACHE_VERSION = "v2"            #change this when data_unification's output changes
    CACHE_SECONDS = 60 * 60 * 24 * 7   #keep results for one week

    # Same key for "cyp2c19*2", " CYP2C19*2 " and so on.
    cache_key = f"pgx:{CACHE_VERSION}:{search_term.strip().upper()}"

    result = cache.get(cache_key)

    if result is None:
        # Combines the PharmVar, CPIC and ClinPGx results for the submitted allele.
        result = data_unification(search_term)

        # Permanent "not found" reasons (e.g. "check input", "Not a star allele") are
        # fine to remember, but temporary failures (network, HTTP, parsing, auth) are
        # not, so the next search tries the APIs again.
        temp_marker = ("try again", "contact it", "failed")
        has_temporary_failure = False

        if result:
            for key in ("PharmVar_Missing_Reason", "CPIC_Missing_Reason", "ClinPGx_Missing_Reason"):
                reason = result.get(key) or ""
                reason = reason.lower()

                for marker in temp_marker:
                    if marker in reason:
                        has_temporary_failure = True
                        break

                if has_temporary_failure:
                    break

        if result and not has_temporary_failure:
            cache.set(cache_key, result, timeout=CACHE_SECONDS)

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
