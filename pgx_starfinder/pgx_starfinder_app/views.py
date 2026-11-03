"""
Django views for displaying the search form and rendering unified PharmVar and CPIC results.
"""
import logging

from django.http import HttpRequest, HttpResponse
from django.shortcuts import render

from .pgx_api import data_unification

from pprint import pformat

# Sets up logging for view module.
logger = logging.getLogger(__name__)

logger.info("Initialising Django view layer")

HOME_TEMPLATE = "pgx_starfinder/index.html"
RESULTS_TEMPLATE = "pgx_starfinder/results.html"

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

    # Combines the PharmVar and CPIC results for the submitted allele.
    result = data_unification(search_term)

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
            "result_display":result_display,
        },
    )
