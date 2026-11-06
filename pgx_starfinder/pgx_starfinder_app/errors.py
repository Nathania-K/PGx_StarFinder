"""
Shared error handling for the external API modules.

The module:
Converts HTTP errors into a standard result dictionary.

Used by pgx_api.py and vv_api.py so every
failed call returns the same Valid/Reason structure.
"""
import logging

logger = logging.getLogger(__name__)


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
