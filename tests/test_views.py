"""Tests for the Django views and HTML pages.

Covers the home and search views in
pgx_starfinder/pgx_starfinder_app/views.py and the templates in
pgx_starfinder/pgx_starfinder_app/templates/pgx_starfinder/.
"""

from django.urls import reverse


def test_home_page_returns_200(client):
    """views.home returns the home page with status 200."""
    response = client.get(reverse("home"))

    assert response.status_code == 200


def test_home_page_uses_index_template(client):
    """views.home renders the pgx_starfinder/index.html template."""
    response = client.get(reverse("home"))

    template_names = [template.name for template in response.templates]

    assert "pgx_starfinder/index.html" in template_names


def test_home_page_has_search_form(client):
    """The home page (views.home) has the search box used by views.search.
    The search box is an HTML <input> element with name="star_allele".
    """
    response = client.get(reverse("home"))

    assert 'name="star_allele"' in response.content.decode()


def test_empty_search_returns_400(client, monkeypatch):
    """views.search returns 400 for an empty search and skips the API.

    Tests the empty-input branch of search() in
    pgx_starfinder/pgx_starfinder_app/views.py. data_unification is
    replaced so no real API call can happen.
    """
    # This list records every call to the pretend API function below.
    calls = []
    # monkeypatch swaps the real API function for a stand-in during this test.
    monkeypatch.setattr(
        # The place to swap: the name data_unification as views.py sees it.
        "pgx_starfinder.pgx_starfinder_app.views.data_unification",
        # The stand-in: "lambda" makes oneline function. "term" is
        # its input (the search text). The stand-in writes its input
        # into the calls list. At the end the test checks the list is empty,
        # which proves it was never called.
        lambda term: calls.append(term),
    )
    # submit the search form with only spaces in the box.
    # "star_allele" is the input's name= in index.html, read in views.py
    # with request.POST.get(...).
    # The view strips the spaces,# so the search text ends up empty.
    response = client.post(reverse("search"), {"star_allele": "   "})

    assert response.status_code == 400
    # The message text is set in views.py, in search(), as "error_message".
    assert response.context["error_message"] == "Please enter a star allele."
    # The pretend API function was never called, so the list is empty.
    assert calls == []
