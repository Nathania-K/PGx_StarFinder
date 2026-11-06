"""
Urls to determine which display view to handle requests.
Checks URL patterns in order (first matching used).
"""

from django.urls import path
from pgx_starfinder.pgx_starfinder_app import views

urlpatterns = [

    # Homepage: 
    # url: '/' as Django does not use leading '/'
    # Displays the search form. 
    path("", views.home, name="home"),

    # Star allele search endpoint: 
    # url: '/search/'
    # Processes POST requests from the home form. 
    path("search/", views.search, name="search"),

    # Link page:
    # url: '/links/'
    # Displays the links page which may be useful to CS.
    path("links/", views.links, name="links"),

]