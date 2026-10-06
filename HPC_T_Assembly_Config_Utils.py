def busco_config_for_lineage(template, lineage):
    """Return a request-specific BUSCO command without mutating its template."""
    return template.replace("{buscolineage}", lineage)
