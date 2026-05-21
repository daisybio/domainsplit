process EXTRACT_UNIQUE_DOMAINS {
    tag { "${domainsplit_db.simpleName}" }
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db

    output:
    path "pfam_ids.txt", emit: pfam_ids
    path "versions.yml", emit: versions

    script:
    """
    sqlite3 "${domainsplit_db}" <<'SQL' > pfam_ids.txt
    .mode list
    .headers off
    SELECT DISTINCT d.pfam_id
    FROM domain d
    WHERE d.id IN (
        SELECT domain_id_a FROM domain_domain_interaction
        UNION
        SELECT domain_id_b FROM domain_domain_interaction
    )
    AND d.pfam_id IS NOT NULL
    AND d.pfam_id != '';
    SQL

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqlite3: \$(sqlite3 --version | awk '{print \$1}')
    END_VERSIONS
    """
}
