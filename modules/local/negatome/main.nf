process DOWNLOAD_NEGATOME {
    tag "negatome"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    val url

    output:
    path "combined_pfam.txt", emit: negatome
    path "versions.yml",      emit: versions

    script:
    """
    #!/usr/bin/env python3
    import shutil
    import sys
    import ssl
    import os
    import urllib.request
    import urllib.error

    URL = "${url}"
    OUT = "combined_pfam.txt"

    def download(url, ctx):
        req = urllib.request.Request(url, headers={"User-Agent": "domainsplit-pipeline"})
        with urllib.request.urlopen(req, context=ctx, timeout=120) as resp, open(OUT, "wb") as fh:
            while True:
                chunk = resp.read(64 * 1024)
                if not chunk:
                    break
                fh.write(chunk)

    if URL.startswith(("http://", "https://", "ftp://", "file://")):
        try:
            download(URL, ssl.create_default_context())
        except (urllib.error.URLError, ssl.SSLError) as exc:
            print(
                f"WARNING: Negatome HTTPS cert validation failed for {URL} ({exc!r}); "
                "retrying with unverified SSL context.",
                file=sys.stderr,
                flush=True,
            )
            download(URL, ssl._create_unverified_context())
    elif os.path.exists(URL):
        shutil.copy(URL, OUT)
    else:
        raise SystemExit(f"url_negatome '{URL}' is neither a supported URL nor a local file")

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
    """
}
