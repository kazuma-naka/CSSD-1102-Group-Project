def cid_to_smiles(cid: int) -> str:
    """
    Given the cid, return the canonical smiles string

    >>> cid_to_smiles(962)
    'O'
    """
    # Build the PubChem PUG-REST URL to retrieve CanonicalSMILES
    # Example:
    # https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/962/property/CanonicalSMILES/txt
    url = generate_url(
        # input part: identify the compound by CID
        f"compound/cid/{cid}",
        "property/CanonicalSMILES",      # operation: request the CanonicalSMILES property
        "txt",                           # output format: plain text
    )

    # Send the request to PubChem and get the raw response text
    smiles_text = make_query(url)

    # Remove any leading/trailing whitespace and newline characters
    smiles = smiles_text.strip()

    # Return the cleaned Canonical SMILES string
    return smiles
