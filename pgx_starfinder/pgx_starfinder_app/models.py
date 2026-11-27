from django.db import models

class star_allele_result(models.Model):
    """Storage of search return json. Model developed using code written by Codex AI"""

    normalised_term = models.CharField(max_length=100, db_index=True)
    schema_version = models.CharField(max_length=20)
    search_return = models.JSONField(null=True)

    # False for temporary network/auth/parsing failures
    reusable = models.BooleanField(default=True)

    fetched_at = models.DateTimeField(auto_now_add=True)
    expires_at = models.DateTimeField(db_index=True)

    class Meta:
        indexes = [
            models.Index(
                fields=[
                    "normalised_term",
                    "schema_version",
                    "reusable",
                    "expires_at"
                ]
            )
        ]
        ordering = ["-fetched_at"]
    
    def __str__(self):
        """return human-readable description of the cached result"""
        return f"{self.normalised_term} ({self.schema_version})"


class star_allele_search(models.Model):
    """Records searches made through app. Model developed using code written by Codex AI"""

    submitted_term = models.CharField(max_length=100)
    normalised_term = models.CharField(max_length=100, db_index=True)
    searched_at = models.DateTimeField(auto_now_add=True)

    result = models.ForeignKey(
        star_allele_result,
        on_delete=models.SET_NULL,
        null=True,
        related_name="searches",
    )

    # Shows whether an API call was avoided
    served_from_database = models.BooleanField(default=False)

    def __str__(self):
        """return human-readable description of the cached result"""
        return f"{self.submitted_term} ({self.searched_at})"