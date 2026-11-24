"""Provides root URL config for PGx_StarFinder."""

from django.urls import include, path

urlpatterns = [
    path("", include("pgx_starfinder.pgx_starfinder_app.urls"))
]